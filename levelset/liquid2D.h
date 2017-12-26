/*
 *  liquid2D.h
 */

#include "common.h"

namespace liquid2D {
	void init( int n );
	void display(int width, int height);
	void keyDown( char key );
	void mouse( FLOAT x, FLOAT y, int state );
	void motion( FLOAT x, FLOAT y, FLOAT dx, FLOAT dy );
}
