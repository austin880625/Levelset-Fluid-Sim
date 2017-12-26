/*
 *  solver.cpp
 *  smoke
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "utility.h"

static char subcell = 0;
static char solver_mode = 0;

// Clamped Fetch
static FLOAT x_ref( FLOAT **A, FLOAT **x, int fi, int fj, int i, int j, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	if( A[i][j] < 0.0 ) return x[i][j];
	return subcell ? A[i][j]/fmin(1.0e-6,A[fi][fj])*x[fi][fj] : 0.0;
}

// Ans = Ax
static void compute_Ax( FLOAT **A, FLOAT **x, FLOAT **ans, int n ) {
	FLOAT h2 = 1.0/(n*n);
	FOR_EVERY_CELL(n) {
		if( A[i][j] < 0.0 ) {
			ans[i][j] = (4.0*x[i][j]-x_ref(A,x,i,j,i+1,j,n)-x_ref(A,x,i,j,i-1,j,n)-x_ref(A,x,i,j,i,j+1,n)-x_ref(A,x,i,j,i,j-1,n))/h2;
		} else {
			ans[i][j] = 0.0;
		}
	} END_FOR
}

// ans = x^T * x
static FLOAT product( FLOAT **A, FLOAT **x, FLOAT **y, int n ) {
	FLOAT ans = 0.0;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) ans += x[i][j]*y[i][j];
		}
	}
	return ans;
}

// x = 0
static void clear( FLOAT **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = 0.0;
		}
	}
}

static void flip( FLOAT **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = -x[i][j];
		}
	}
}

// x <= y
static void copy( FLOAT **x, FLOAT **y, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = y[i][j];
		}
	}
}
				 
// Ans = x + a*y
static void op( FLOAT **A, FLOAT **x, FLOAT **y, FLOAT **ans, FLOAT a, int n ) {
	static FLOAT **tmp = alloc2D<FLOAT>(n);
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) tmp[i][j] = x[i][j]+a*y[i][j];
		}
	}
	copy(ans,tmp,n);
}

// r = b - Ax
static void residual( FLOAT **A, FLOAT **x, FLOAT **b, FLOAT **r, int n ) {
	compute_Ax(A,x,r,n);
	op( A, b, r, r, -1.0, n );
}

static inline FLOAT square( FLOAT a ) {
	return a*a;
}

static FLOAT A_ref( FLOAT **A, int i, int j, int qi, int qj, int n ) {
	if( i<0 || i>n-1 || j<0 || j>n-1 || A[i][j]>0.0 ) return 0.0;
	if( qi<0 || qi>n-1 || qj<0 || qj>n-1 || A[qi][qj]>0.0 ) return 0.0;
	return -1.0;
}

static FLOAT A_diag( FLOAT **A, int i, int j, int n ) {
	FLOAT diag = 4.0;
	if( A[i][j] > 0.0 ) return diag;
	int q[][2] = { {i-1,j}, {i+1,j}, {i,j-1}, {i,j+1} };
	for( int m=0; m<4; m++ ) {
		int qi = q[m][0];
		int qj = q[m][1];
		if( qi<0 || qi>n-1 || qj<0 || qj>n-1 ) diag -= 1.0;
		else if( A[qi][qj] > 0.0 && subcell ) {
			diag -= A[qi][qj]/fmin(1.0e-6,A[i][j]);
		}
	}
	return diag;
}

static FLOAT P_ref( FLOAT **P, int i, int j, int n ) {
	if( i<0 || i>n-1 || j<0 || j>n-1 ) return 0.0;
	return P[i][j];
}

static void buildPreconditioner( FLOAT **P, FLOAT **A, int n ) {
	clear(P,n);
	FLOAT t = solver_mode == 2 ? 0.97 : 0.0;
	FLOAT a = 0.25;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) {
				FLOAT left = A_ref(A,i-1,j,i,j,n)*P_ref(P,i-1,j,n);
				FLOAT bottom = A_ref(A,i,j-1,i,j,n)*P_ref(P,i,j-1,n);
				FLOAT mleft = A_ref(A,i-1,j,i,j,n)*A_ref(A,i,j-1,i,j,n)*square(P_ref(P,i-1,j,n));
				FLOAT mbottom = A_ref(A,i,j-1,i,j,n)*A_ref(A,i-1,j,i,j,n)*square(P_ref(P,i,j-1,n));
				
				FLOAT diag = A_diag( A, i, j, n );
				FLOAT e = diag - square(left) - square(bottom) - t*( mleft + mbottom );
				if( e < a*diag ) e = diag;
				P[i][j] = 1.0/sqrtf(e);
			}
		}
	}
}

static void applyPreconditioner( FLOAT **z, FLOAT **r, FLOAT **P, FLOAT **A, int n ) {
	if( solver_mode == 0 ) {
		copy(z,r,n);
		return;
	}
	
	static FLOAT **q = alloc2D<FLOAT>(n);
	clear(q,n);
	
	// Lq = r
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( A[i][j] < 0.0 ) {
				FLOAT left = A_ref(A,i-1,j,i,j,n)*P_ref(P,i-1,j,n)*P_ref(q,i-1,j,n);
				FLOAT bottom = A_ref(A,i,j-1,i,j,n)*P_ref(P,i,j-1,n)*P_ref(q,i,j-1,n);
				
				FLOAT t = r[i][j] - left - bottom;
				q[i][j] = t*P[i][j];
			}
		}
	}
	
	// L^T z = q
	for( int i=n-1; i>=0; i-- ) {
		for( int j=n-1; j>=0; j-- ) {
			if( A[i][j] < 0.0 ) {
				FLOAT right = A_ref(A,i,j,i+1,j,n)*P_ref(P,i,j,n)*P_ref(z,i+1,j,n);
				FLOAT top = A_ref(A,i,j,i,j+1,n)*P_ref(P,i,j,n)*P_ref(z,i,j+1,n);
				
				FLOAT t = q[i][j] - right - top;
				z[i][j] = t*P[i][j];
			}
		}
	}
}

static void conjGrad( FLOAT **A, FLOAT **P, FLOAT **x, FLOAT **b, int n ) {
	// Pre-allocate Memory
	static FLOAT **r = alloc2D<FLOAT>(n);
	static FLOAT **z = alloc2D<FLOAT>(n);
	static FLOAT **s = alloc2D<FLOAT>(n);
	
	clear(x,n);									// p = 0
	copy(r,b,n);								// r = b
	applyPreconditioner(z,r,P,A,n);				// Apply Conditioner z = f(r)
	copy(s,z,n);								// s = z
	
	FLOAT a = product( A, z, r, n );			// a = z . r
	for( int k=0; k<n*n; k++ ) {
		compute_Ax( A, s, z, n );				// z = applyA(s)
		FLOAT alpha = a/product( A, z, s, n );	// alpha = a/(z . s)
		op( A, x, s, x, alpha, n );				// p = p + alpha*s
		op( A, r, z, r, -alpha, n );			// r = r - alpha*z;
		FLOAT error2 = product( A, r, r, n );	// error2 = r . r
		if( error2/(n*n) < 1.0e-6 ) break;
		applyPreconditioner(z,r,P,A,n);			// Apply Conditioner z = f(r)
		FLOAT a2 = product( A, z, r, n );		// a2 = z . r
		FLOAT beta = a2/a;
		op( A, z, s, s, beta, n );				// s = z + beta*s
		a = a2;
	}
}

FLOAT solver::solve( FLOAT **A, FLOAT **x, FLOAT **b, int n, char subcell_aware, char solver_type ) {
	static FLOAT **r = alloc2D<FLOAT>(n);
	static FLOAT **P = alloc2D<FLOAT>(n);
	clear(r,n);
	
	// Save Mode
	subcell = subcell_aware;
	solver_mode = solver_type;
	
	// Flip Divergence
	flip(b,n);
	
	// Build Modified Incomplete Cholesky Precondioner Matrix
	if( solver_mode >= 1 ) buildPreconditioner(P,A,n);
	
	// Conjugate Gradient Method
	conjGrad(A,P,x,b,n);

	residual(A,x,b,r,n);
	FLOAT res = sqrt(product( A, r, r, n ))/(n*n);
	// printf( "Residual = %e\n", res );
	return res;
}