/*
 *  levelset2D.h
 */

#include "common.h"

namespace levelset2D {
	void init( int n );
	void buildLevelset( bool (*func)(FLOAT x, FLOAT y), FLOAT tolerance );
	void advect( void (*func)( FLOAT x, FLOAT y, FLOAT &u, FLOAT &v, FLOAT &dt ) );
	void redistance( FLOAT tolerance );
	void extrapolate( FLOAT **q, char **region );
	void display( bool cell_centered );
	void keyDown( char key );
	///
	FLOAT getLevelSet( int i, int j );
	void getLevelSet( FLOAT **dists );
	void setLevelSet( int i, int j, FLOAT d );
	FLOAT getVolume();
	///
	void setVisibility( bool show_grid, bool show_dist, bool show_region );
}
