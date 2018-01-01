/*
 *  interp.h
 */

#include "common.h"
namespace interp {
	//FLOAT spline( FLOAT **d, FLOAT x, FLOAT y, int w, int h );
	FLOAT linear ( FLOAT ***d, FLOAT x, FLOAT y, FLOAT z, int w, int h, int dep );
	FLOAT interp ( FLOAT ***d, FLOAT x, FLOAT y, FLOAT z, int w, int h, int dep );
	///
	void setInterpMethod( int num );
}
