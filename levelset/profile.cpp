/*
 *  profile.cpp
 */

#include "profile.h"
#include "levelset2D.h"

namespace {
	int gn;
	FLOAT maxdist = 0.15;
}

static bool sphere( FLOAT x, FLOAT y ) {
	return hypot(x-0.5,y-0.75) < 0.15;
}

void profile::init( int n ) {
	gn = n;
	levelset2D::init(n);
	
	// Build Simple LevelSet
	levelset2D::buildLevelset(sphere,maxdist);
}

static void flow( FLOAT x, FLOAT y, FLOAT &u, FLOAT &v, FLOAT &dt ) {
	u = -(y-0.5);
	v = x-0.5;
	dt = 0.01;
}

void profile::display() {
	// Advect
	levelset2D::advect(flow);
	levelset2D::redistance(maxdist);	
	levelset2D::display(true);
}

void profile::keyDown( unsigned char key ) {
	if( key == 'r' ) {
		levelset2D::redistance(maxdist);
	}
}