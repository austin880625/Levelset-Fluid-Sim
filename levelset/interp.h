/*
 *  interp.h
 */

#include "common.h"
namespace interp {
	FLOAT spline( FLOAT **d, FLOAT x, FLOAT y, int w, int h );
	FLOAT linear ( FLOAT **d, FLOAT x, FLOAT y, int w, int h );
	FLOAT interp ( FLOAT **d, FLOAT x, FLOAT y, int w, int h );
	///
	void setInterpMethod( int num );
}