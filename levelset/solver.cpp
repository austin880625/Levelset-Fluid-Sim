/*
 *  solver.cpp
 *  smoke
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include "solver.h"
#include "utility.h"

static char subcell = 0;
static char solver_mode = 0;
namespace {
	int n;
	cl_command_queue Q;
	cl_context context;
	cl_program program;
	cl_mem A_cl, p_cl, r_cl, z_cl, s_cl;
	cl_mem tmp_prod;
	cl_kernel kern_compute_Ax, kern_op, kern_product, kern_copy, kern_clear;
}
static void print_cl_mem(const char * s, cl_mem* tmp, int offset, int len)
{
	return ;
	static float* buf = new float[len];
	clEnqueueReadBuffer(Q, *tmp, CL_TRUE, offset*sizeof(float), len*sizeof(float), buf, 0, NULL, NULL);	
	clFinish(Q);
	printf("print_cl_mem %s : ", s);
	for(int i = 0; i < len; i++)
		printf("%f ", buf[i]);
	puts("");
}

// Clamped Fetch
static FLOAT x_ref( FLOAT ***A, FLOAT ***x, int fi, int fj, int fk, int i, int j, int k, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	k = min(max(0,k),n-1);
	if( A[i][j][k] < 0.0 ) return x[i][j][k];
	return subcell ? A[i][j][k]/fmin(1.0e-6,A[fi][fj][fk])*x[fi][fj][fk] : 0.0;
}

// Ans = Ax
static void compute_Ax( FLOAT ***A, FLOAT ***x, FLOAT ***ans, int n ) {
	FLOAT h2 = 1.0/(n*n);
#pragma omp parallel for
	FOR_EVERY_CELL(n) {
		if( A[i][j][k] < 0.0 ) {
			ans[i][j][k] = (6.0*x[i][j][k]-x_ref(A,x,i,j,k,i+1,j,k,n)-x_ref(A,x,i,j,k,i-1,j,k,n)-x_ref(A,x,i,j,k,i,j+1,k,n)-x_ref(A,x,i,j,k,i,j-1,k,n)-x_ref(A,x,i,j,k,i,j,k-1,n)-x_ref(A,x,i,j,k,i,j,k+1,n))/h2;
		} else {
			ans[i][j][k] = 0.0;
		}
	} END_FOR
}

// ans = x^T * x
static FLOAT product( cl_mem *A, cl_mem *x, cl_mem *y, size_t ng ) {
	FLOAT ans = 0.0;
	size_t dotGlobalws = 16*((n*n*n)/16+1);
	size_t dotLocalws = 16;
	cl_int error;
	static float *partial = (float*)malloc((ng)*sizeof(float));
	//printf("partial: %d\n",partial);
#pragma omp parallel for
	for(size_t i=0;i<ng;i++){
		partial[i]=0;
	}
	clSetKernelArg(kern_product, 0, sizeof(cl_mem), (void*)A);
	clSetKernelArg(kern_product, 1, sizeof(cl_mem), (void*)x);
	clSetKernelArg(kern_product, 2, sizeof(cl_mem), (void*)y);
	clSetKernelArg(kern_product, 3, sizeof(cl_mem), (void*)&tmp_prod);
	clSetKernelArg(kern_product, 4, sizeof(int), (void*)&n);
	error = clSetKernelArg(kern_product, 5, sizeof(float)*16, NULL);
	error = clEnqueueNDRangeKernel(Q, kern_product, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
	printf("%d\n",error);
	clFinish(Q);
	error=clEnqueueReadBuffer(Q, tmp_prod, CL_TRUE, 0, sizeof(float)*ng, partial, 0,NULL,NULL);
	clFinish(Q);
	
	printf("%d\n",error);
#pragma omp parallel for reduction(+:ans)
	for( size_t i=0; i<ng; i++ ) {
		if(partial[i]==partial[i]){
			ans+=partial[i];
			//printf("ipartial: %f\n",partial[i]);
		}
	}
	//printf("prod %f\n",ans);
	return ans;
}
static FLOAT product2( cl_mem *A, cl_mem *x, cl_mem *y, size_t ng ){
	static FLOAT *xx = (float*)malloc(n*n*n*sizeof(FLOAT));
	static FLOAT *yy = (float*)malloc(n*n*n*sizeof(FLOAT));
	static FLOAT *AA = (float*)malloc(n*n*n*sizeof(FLOAT));
	clEnqueueReadBuffer(Q, *A, CL_TRUE, 0, sizeof(float)*n*n*n, AA, 0,NULL,NULL);
	clEnqueueReadBuffer(Q, *x, CL_TRUE, 0, sizeof(float)*n*n*n, xx, 0,NULL,NULL);
	clEnqueueReadBuffer(Q, *y, CL_TRUE, 0, sizeof(float)*n*n*n, yy, 0,NULL,NULL);
	clFinish(Q);
	FLOAT ans = 0.0;
	for(int i=0;i<n*n*n;i++){
		if(AA[i]<0.0)ans+=xx[i]*yy[i];
	}
	return ans;
}
// x = 0
static void clear( cl_mem *x, int n ) {
	clSetKernelArg(kern_clear, 0, sizeof(cl_mem), (void*)x);
	clSetKernelArg(kern_clear, 1, sizeof(int), (void*)&n);
	size_t dotGlobalws = 16*((n*n*n)/16+1);
	size_t dotLocalws = 16;
	int err = clEnqueueNDRangeKernel(Q, kern_clear, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
	//printf("clear err: %d\n", err);
	print_cl_mem("after clear", x, 0, n*n*n);
	clFinish(Q);

}

static void flip( FLOAT ***x, int n ) {
#pragma omp parallel for
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			for(int k=0; k<n; k++){
				x[i][j][k] = -x[i][j][k];
			}
		}
	}
}

// x <= y
static void copy( FLOAT ***x, FLOAT ***y, int n ) {
#pragma omp parallel for
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			for(int k=0; k<n; k++){
				x[i][j][k] = y[i][j][k];
			}
		}
	}
}
				 
// Ans = x + a*y
static void op( FLOAT ***A, FLOAT ***x, FLOAT ***y, FLOAT ***ans, FLOAT a, int n ) {
	static FLOAT ***tmp = alloc3D<FLOAT>(n);
#pragma omp parallel for
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			for(int k=0; k<n; k++){
				if( A[i][j][k] < 0.0 ) tmp[i][j][k] = x[i][j][k]+a*y[i][j][k];
				else tmp[i][j][k] = 0.0;
			}
		}
	}
	copy(ans,tmp,n);
}

// r = b - Ax
static void residual( FLOAT ***A, FLOAT ***x, FLOAT ***b, FLOAT ***r, int n ) {
	compute_Ax(A,x,r,n);
	op( A, b, r, r, -1.0, n );
}

static inline FLOAT square( FLOAT a ) {
	return a*a;
}
static void applyPreconditioner( FLOAT ***z, FLOAT ***r, FLOAT ***P, FLOAT ***A, int n ) {
	if( solver_mode == 0 ) {
		copy(z,r,n);
		return;
	}
}

static void conjGrad( FLOAT ***A, FLOAT ***x, FLOAT ***b, int n ) {
	// Pre-allocate Memory
	//printf("%d",solver_mode);
	//static FLOAT ***r = alloc3D<FLOAT>(n);
	//static FLOAT ***z = alloc3D<FLOAT>(n);
	//static FLOAT ***s = alloc3D<FLOAT>(n);
	//printf("starting CG\n");
	cl_int error = CL_SUCCESS;
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			//printf("%d %d %f\n",n,sizeof(A),*(A[i][j]+n-1));
			error=clEnqueueWriteBuffer(Q, A_cl, CL_TRUE, (i*n*n+j*n)*sizeof(float), n*sizeof(float), (void*)A[i][j], 0,NULL,NULL);
			//printf("A_cl finished with code %d\n",error);
			error=clEnqueueWriteBuffer(Q, r_cl, CL_TRUE, (i*n*n+j*n)*sizeof(float), n*sizeof(float), b[i][j], 0,NULL,NULL);
			//printf("r_cl finished with code %d\n",error);
			error=clEnqueueWriteBuffer(Q, z_cl, CL_TRUE, (i*n*n+j*n)*sizeof(float), n*sizeof(float), b[i][j], 0,NULL,NULL);
			//printf("z_cl finished with code %d\n",error);
			error=clEnqueueWriteBuffer(Q, s_cl, CL_TRUE, (i*n*n+j*n)*sizeof(float), n*sizeof(float), b[i][j], 0,NULL,NULL);
			clFinish(Q);
			//printf("s_cl finished with code %d\n",error);
		}
	}
	size_t num_g = 16*(n/16+1);
	size_t globalws[3] = {num_g,num_g,num_g};
	size_t localws[3] = {16,16,16};
	size_t dotGlobalws = 16*((n*n*n)/16+1);
	size_t dotLocalws = 16;
	
	size_t dotnum_g = (n*n*n/16+1);
	//puts("starting iteration");
	clear(&p_cl,n);
	print_cl_mem("p_before", &p_cl, 0, n*n*n);
	FLOAT a = product( &A_cl, &z_cl, &r_cl, dotnum_g );			// a = z . r
	for( int k=0; k<50; k++ ) {
		//compute_Ax( A, s, z, n );				// z = applyA(s)
		
		error=clSetKernelArg(kern_compute_Ax, 0, sizeof(cl_mem), (void*)&A_cl);
		clSetKernelArg(kern_compute_Ax, 1, sizeof(cl_mem), (void*)&s_cl);
		clSetKernelArg(kern_compute_Ax, 2, sizeof(cl_mem), (void*)&z_cl);
		clSetKernelArg(kern_compute_Ax, 3, sizeof(int), (void*)&n);
		printf("%d\n",error);
		
		error=clEnqueueNDRangeKernel(Q, kern_compute_Ax, 3, NULL, (size_t*)globalws, (size_t*)localws, 0,NULL,NULL);
		printf("%d\n", error);
		clFinish(Q);
		print_cl_mem("z_cl", &z_cl, 0, n*n*n);
		print_cl_mem("p_cl", &p_cl, 0, n*n*n);
		print_cl_mem("s_cl", &s_cl, 0, n*n*n);
		FLOAT alpha = a/product( &A_cl, &z_cl, &s_cl, dotnum_g );	// alpha = a/(z . s)
		//printf("kern_compute_Ax %d %f %f\n",error, a, alpha);
		printf("%d\n",error);
		clSetKernelArg(kern_op, 0, sizeof(cl_mem), (void*)&A_cl);
		clSetKernelArg(kern_op, 1, sizeof(cl_mem), (void*)&p_cl);
		clSetKernelArg(kern_op, 2, sizeof(cl_mem), (void*)&s_cl);
		clSetKernelArg(kern_op, 3, sizeof(cl_mem), (void*)&p_cl);
		clSetKernelArg(kern_op, 4, sizeof(float), (void*)&alpha);
		clSetKernelArg(kern_op, 5, sizeof(int), (void*)&n);
		error = clEnqueueNDRangeKernel(Q, kern_op, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
		printf("%d\n",error);
		clFinish(Q);
		print_cl_mem("p_cl", &p_cl, 0, n*n*n);
		//printf("op %d\n",error);
		//op( A, x, s, x, alpha, n );				// p = p + alpha*s
		
		alpha=-alpha;
		print_cl_mem("r_cl", &r_cl, 0, n*n*n);
		print_cl_mem("z_cl", &z_cl, 0, n*n*n);
		clSetKernelArg(kern_op, 0, sizeof(cl_mem), (void*)&A_cl);
		clSetKernelArg(kern_op, 1, sizeof(cl_mem), (void*)&r_cl);
		clSetKernelArg(kern_op, 2, sizeof(cl_mem), (void*)&z_cl);
		clSetKernelArg(kern_op, 3, sizeof(cl_mem), (void*)&r_cl);
		clSetKernelArg(kern_op, 4, sizeof(float), (void*)&alpha);
		clSetKernelArg(kern_op, 5, sizeof(int), (void*)&n);
		error=clEnqueueNDRangeKernel(Q, kern_op, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
		printf("%d\n",error);
		clFinish(Q);
		print_cl_mem("r_cl", &r_cl, 0, n*n*n);
		//printf("op %d\n",error);
		//op( A, r, z, r, -alpha, n );			// r = r - alpha*z;
		FLOAT error2 = product( &A_cl, &r_cl, &r_cl, dotnum_g );	// error2 = r . r
		//printf("%f\n", error2);
		if( error2/(n*n*n) < 1.0e-6 ) break;
		
		clSetKernelArg(kern_copy, 0, sizeof(cl_mem), (void*)&z_cl);
		clSetKernelArg(kern_copy, 1, sizeof(cl_mem), (void*)&r_cl);
		clSetKernelArg(kern_copy, 2, sizeof(int), (void*)&n);
		error = clEnqueueNDRangeKernel(Q, kern_copy, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
		printf("%d\n",error);
		clFinish(Q);
		//printf("copy %d\n",error);
		//applyPreconditioner(z,r,P,A,n);			// Apply Conditioner z = f(r)
		
		FLOAT a2 = product( &A_cl, &z_cl, &r_cl, dotnum_g );		// a2 = z . r
		FLOAT beta = a2/a;
		clSetKernelArg(kern_op, 0, sizeof(cl_mem), (void*)&A_cl);
		clSetKernelArg(kern_op, 1, sizeof(cl_mem), (void*)&z_cl);
		clSetKernelArg(kern_op, 2, sizeof(cl_mem), (void*)&s_cl);
		clSetKernelArg(kern_op, 3, sizeof(cl_mem), (void*)&s_cl);
		clSetKernelArg(kern_op, 4, sizeof(float), (void*)&beta);
		clSetKernelArg(kern_op, 5, sizeof(int), (void*)&n);
		error = clEnqueueNDRangeKernel(Q, kern_op, 1, NULL, (size_t*)&dotGlobalws, (size_t*)&dotLocalws, 0,NULL,NULL);
		printf("%d\n",error);
		clFinish(Q);
		//op( A, z, s, s, beta, n );				// s = z + beta*s
		a = a2;
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			clEnqueueReadBuffer(Q, p_cl, CL_TRUE, (i*n*n+j*n)*sizeof(float), sizeof(float)*n, x[i][j], 0,NULL,NULL);
			clFinish(Q);
		}
	}
}

FLOAT solver::solve( FLOAT ***A, FLOAT ***x, FLOAT ***b, char subcell_aware, char solver_type ) {
	//static FLOAT ***r = alloc3D<FLOAT>(n);
	//static FLOAT ***P = alloc3D<FLOAT>(n);
	//clear(r,n);
	
	// Save Mode
	subcell = subcell_aware;
	solver_mode = solver_type;
	
	// Flip Divergence
	flip(b,n);
	
	// Build Modified Incomplete Cholesky Precondioner Matrix
	//if( solver_mode >= 1 ) buildPreconditioner(P,A,n);
	
	// Conjugate Gradient Method
	conjGrad(A,x,b,n);
	//residual(A,x,b,r,n);
	//FLOAT res = sqrt(product( A, r, r, n ))/(n*n*n);
	//printf( "Residual = %e\n", res );
	return 0.0;
}

void loadKernel(cl_device_id dev){
	static const int MAX_SOURCE_SIZE = 4048576;
	FILE *kernel_fp = fopen("levelset/kernels/solver.cl", "r");
	printf("File opend as %d\n");
	char * kernelSource = (char *)malloc(MAX_SOURCE_SIZE*sizeof(char));
	unsigned int sourceSize = fread( kernelSource, 1, MAX_SOURCE_SIZE, kernel_fp );
	fclose(kernel_fp);
	kernelSource[sourceSize++] = '\0';

	cl_int error = CL_SUCCESS;
	program = clCreateProgramWithSource(context, 1, (const char **)&kernelSource, NULL, &error);
	printf("program creating error code: %d\n",error);
	error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	printf("program building error code: %d\n",error);
	if(error < 0) {
		size_t log_size;
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 
			0, NULL, &log_size);
		char *program_log = (char*) malloc(log_size + 1);
		program_log[log_size] = '\0';
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 
			log_size + 1, program_log, NULL);
		printf("%s\n", program_log);
		free(program_log);
		exit(1);
	}

	kern_compute_Ax = clCreateKernel(program, "compute_Ax", &error);
	kern_op = clCreateKernel(program, "op", &error);
	kern_product = clCreateKernel(program, "product", &error);
	kern_copy = clCreateKernel(program, "copy", &error);
	kern_clear = clCreateKernel(program, "clear", &error);
	printf("%d\n",error);
}
void solver::setCL(int _n){
	n=_n;
	cl_uint platformIdCount = 0;
	clGetPlatformIDs(0, NULL, &platformIdCount);

	std::vector<cl_platform_id> platformIds(platformIdCount);
	clGetPlatformIDs(platformIdCount, platformIds.data(), NULL);

	cl_uint deviceIdCount = 0;
	clGetDeviceIDs(platformIds[0], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceIdCount);
	
	std::vector<cl_device_id> deviceIds(deviceIdCount);
	clGetDeviceIDs(platformIds[0], CL_DEVICE_TYPE_ALL, deviceIdCount, deviceIds.data(), NULL);

	const cl_context_properties contextProperties [] = {CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(platformIds[0]), 0, 0};

	cl_int error = CL_SUCCESS;
	context = clCreateContext(contextProperties, deviceIdCount, deviceIds.data(), NULL, NULL, &error);
	error = CL_SUCCESS;
	printf("num of dev: %u\n",deviceIds.size());
	Q = clCreateCommandQueue(context, deviceIds[1], 0, &error);
	
	printf(">>> %d\n",error);

	loadKernel(deviceIds[1]);
	printf("Kernel loaded\n");
	int n3=n*n*n;
	A_cl = clCreateBuffer(context, CL_MEM_READ_ONLY, n3*sizeof(float), NULL, &error);
	p_cl = clCreateBuffer(context, CL_MEM_WRITE_ONLY, n3*sizeof(float), NULL, &error);
	r_cl = clCreateBuffer(context, CL_MEM_READ_WRITE, n3*sizeof(float), NULL, &error);
	z_cl = clCreateBuffer(context, CL_MEM_READ_WRITE, n3*sizeof(float), NULL, &error);
	s_cl = clCreateBuffer(context, CL_MEM_READ_WRITE, n3*sizeof(float), NULL, &error);
	tmp_prod = clCreateBuffer(context, CL_MEM_WRITE_ONLY, (n3/16+1)*sizeof(float), NULL, &error);
}
