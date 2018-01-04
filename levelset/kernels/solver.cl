float x_ref( global float *A, global float *x, int fi, int fj, int fk, int i, int j, int k, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	k = min(max(0,k),n-1);
	if( A[i*n*n+j*n+k] < 0.0 ) return x[i*n*n+j*n+k];
	return A[i*n*n+j*n+k]/min((float)1.0e-6,A[fi*n*n+fj*n+fk])*x[fi*n*n+fj*n+fk];
}
__kernel void compute_Ax( __global float *A, __global float *x, __global float *ans, int n ) {
	int i=get_global_id(0);
	int j=get_global_id(1);
	int k=get_global_id(2);
	if(i>=n||j>=n||k>=n)return ;
	float h2 = 1.0/(n*n);
	if( A[i*n*n+j*n+k] < 0.0 ) {
		ans[i*n*n+j*n+k] = (6.0*x[i*n*n+j*n+k]-x_ref(A,x,i,j,k,i+1,j,k,n)-x_ref(A,x,i,j,k,i-1,j,k,n)-x_ref(A,x,i,j,k,i,j+1,k,n)-x_ref(A,x,i,j,k,i,j-1,k,n)-x_ref(A,x,i,j,k,i,j,k-1,n)-x_ref(A,x,i,j,k,i,j,k+1,n))/h2;
	} else {
		ans[i*n*n+j*n+k] = 0.0;
	}
}
__kernel void flip( __global float *x, int n ) {
	int i=get_global_id(0);
	if(i>=n*n*n)return;
	x[i] = -x[i];
}
__kernel void clear( __global float *x, int n ) {
	int i=get_global_id(0);
	if(i>=n*n*n)return;
	x[i] = 0;
}
__kernel void copy( __global float *x, __global float *y, int n ) {
	int i=get_global_id(0);
	if(i>=n*n*n)return;
	x[i] = y[i];
}
__kernel void op( __global float *A, __global float *x, __global float *y, __global float *ans, float a, int n ) {
	int i=get_global_id(0);
	if(i>=n*n*n)return ;
	if( A[i] < 0.0 ) ans[i] = x[i]+a*y[i];
	else ans[i] = 0.0;
}
__kernel void product(__global float *A, __global float *x, __global float *y, __global float *res, int n, __local float *partial){
	int gi=get_global_id(0);
	int li=get_local_id(0);
	int gs=get_local_size(0);
	partial[li]=0;

	if(gi>=n*n*n){ return ;}
	if(A[gi]<0.0)partial[li]=x[gi]*y[gi];
	barrier(CLK_LOCAL_MEM_FENCE);

	for(int i=gs/2; i>0; i>>=1){
		if(li<i){
			partial[li]+=partial[li+i];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if(li==0){
		res[get_group_id(0)]=partial[0];
	}
}


