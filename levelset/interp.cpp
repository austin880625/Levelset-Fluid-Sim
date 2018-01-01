/*
 *  interp.cpp
 */

#include "interp.h"

static int method = 0;
/*
static FLOAT clamp( FLOAT v, FLOAT min, FLOAT max ) {
	if( v < min ) return min;
	if( v > max ) return max;
	return v;
}

static int iclamp( int v, int min, int max ) {
	if( v < min ) return min;
	if( v > max ) return max;
	return v;
}

static FLOAT spline_cubic(const FLOAT a[4], FLOAT x) {
	if( a[0] == a[1] && a[0] == a[2] && a[0] == a[3] ) return a[0];
	int i, j;
	FLOAT alpha[4], l[4], mu[4], z[4];
	FLOAT b[4], c[4], d[4];
	for(i = 1; i < 3; i++)
		alpha[i] = 3. * (a[i+1] - a[i]) - 3. * (a[i] - a[i-1]);
	l[0] = 1.;
	mu[0] = 0.;
	z[0] = 0.;
	for(i = 1; i < 3; i++){
		l[i] = 4. - mu[i-1];
		mu[i] = 1. / l[i];
		z[i] = (alpha[i] - z[i-1]) / l[i];
	}
	l[3] = 1.;
	z[3] = 0.;
	c[3] = 0.;
	for(j = 2; 0 <= j; j--){
		c[j] = z[j] - mu[j] * c[j+1];
		b[j] = a[j+1] - a[j] - (c[j+1] + 2. * c[j]) / 3.;
		d[j] = (c[j+1] - c[j]) / 3.;
	}
	return clamp( a[1] + b[1] * x + c[1] * x * x + d[1] * x * x * x, fmin(a[1],a[2]), fmax(a[1],a[2]) );
}

FLOAT interp::spline( FLOAT **d, FLOAT x, FLOAT y, int w, int h ) {
	FLOAT f[16];
	FLOAT xn[4];
	
	if (x<0.0) x=0.0; if (x>=w) x=w-0.01;
	if (y<0.0) y=0.0; if (y>=h) y=h-0.01;
	
	for( int j=0; j<4; j++ ) {
		for( int i=0; i<4; i++ ) {
			int hd = (int)x - 1 + i;
			int v = (int)y - 1 + j;
			f[4*j+i] = d[iclamp(hd,0,w-1)][iclamp(v,0,h-1)];
		}
	}
	
	for( int j=0; j<4; j++ ) xn[j] = spline_cubic( &f[4*j], x - (int)x );	
	return spline_cubic( xn, y - (int)y );
}
*/
FLOAT interp::linear ( FLOAT ***d, FLOAT x, FLOAT y, FLOAT z, int w, int h, int dep ) {
	x = fmax(0.0,fmin(w,x));
	y = fmax(0.0,fmin(h,y));
	z = fmax(0.0,fmin(dep,z));
	int i = min(x,w-2);
	int j = min(y,h-2);
	int k = min(z,dep-2);
	
	//return ( ((i+1-x)*d[i][j][k]+(x-i)*d[i+1][j][k])*(j+1-y) + ((i+1-x)*d[i][j+1][k]+(x-i)*d[i+1][j+1][k])*(y-j) )*(k+1-z) + ( ((i+1-x)*d[i][j][k+1]+(x-i)*d[i+1][j][k+1])*(j+1-y) + ((i+1-x)*d[i][j+1][k+1]+(x-i)*d[i+1][j+1][k+1])*(y-j) )*(z-k);
	return ( ((k+1-z)*d[i][j][k]+(z-k)*d[i][j][k+1])*(j+1-y) + ((k+1-z)*d[i][j+1][k]+(z-k)*d[i][j+1][k+1])*(y-j) )*(i+1-x) + ( ((k+1-z)*d[i+1][j][k]+(z-k)*d[i+1][j][k+1])*(j+1-y) + ((k+1-z)*d[i+1][j+1][k]+(z-k)*d[i+1][j+1][k+1])*(y-j) )*(x-i);
}

FLOAT interp::interp ( FLOAT ***d, FLOAT x, FLOAT y, FLOAT z, int w, int h, int dep ) {
	FLOAT r = 0.0;
	switch (method) {
		case 0:
			r = linear(d,x,y,z,w,h,dep);
			break;
		case 1:
			//r = spline(d,x,y,w,h);
			break;
	} 
	return r;
}

void interp::setInterpMethod( int num ) {
	method = num;
}
