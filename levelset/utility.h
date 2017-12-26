/*
 *  utility.h
 *  smoke
 *
 */

#include "common.h"

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW(N)	for( int xn=0; xn<(N+1)*(N); xn++ ) { int i=xn%(N+1); int j=xn/(N+1);
#define FOR_EVERY_Y_FLOW(N)	for( int yn=0; yn<(N+1)*(N); yn++ ) { int i=yn%(N); int j=yn/(N);
#define FOR_EVERY_CELL(N)	for( int ci=0; ci<(N)*(N); ci++ ) { int i=ci%(N); int j=ci/(N);
#define END_FOR }

#include <omp.h>
#ifdef _OPENMP
#define OPENMP_FOR		_Pragma("omp parallel for" )
#define OPENMP_SECTION  _Pragma("omp section" )
#define OPENMP_BEGIN	_Pragma("omp parallel" ) {
#define OPENMP_END		}
#define OPENMP_FOR_P	_Pragma("omp for" )
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif

template <class T> T ** alloc2D( int n ) {
	T **ptr = new T *[n+1];
	for( int i=0; i<n; i++ ) {
		ptr[i] = new T[n+1];
	}
	ptr[n] = 0;
	return ptr;
}

template <class T> void free2D( T **ptr ) {
	for( int i=0; ptr[i]; i++ ) delete [] ptr[i];
	delete [] ptr;
}

static inline FLOAT hypot2( FLOAT a, FLOAT b ) {
    return a*a+b*b;
}

static inline FLOAT length2( FLOAT *p0, FLOAT *p1 ) {
    return hypot2(p0[0]-p1[0],p0[1]-p1[1]);
}

