/*
 *  levelset2D.cpp
 */

#include "levelset2D.h"
#include <stdio.h>
#include <vector>
#include <queue>
#include "interp.h"
#include "utility.h"
#include <math.h>

using namespace std;

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <sys/time.h>
#elif defined(WIN32)
#include "glut.h"
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <sys/time.h>
#endif

namespace {
	typedef struct {
		char known;
		char estimated;
		FLOAT dist;
		FLOAT pos[2];
		int gp[2];
	} grid;
	grid ** grids = NULL;
	FLOAT ** grad[2] = { NULL, NULL };
	FLOAT **test_q = NULL;
	static int gn;
	struct gridComparison{
		bool operator () ( grid * left, grid * right){
			return fabs(left->dist) > fabs(right->dist); 
		}
	};
	typedef priority_queue<grid *,vector<grid*>,gridComparison> grid_queue;
	///
	char show_grid = 1;
	char show_dist = 1;
	char show_region = 1;
}

void levelset2D::init( int n ) {
	// Allocate LevelSet Grids
	grids = alloc2D<grid>(n);
	grad[0] = alloc2D<FLOAT>(n);
	grad[1] = alloc2D<FLOAT>(n);
	test_q = alloc2D<FLOAT>(n);
	gn = n;
	FOR_EVERY_CELL(gn) {
		grids[i][j].dist = 1.0;
		grids[i][j].known = false;
		grids[i][j].gp[0] = i;
		grids[i][j].gp[1] = j;
		grad[0][i][j] = grad[1][i][j] = 0.0;
	} END_FOR;
}

// Helper Function Used In FastMarch()
static bool update_distance( int i, int j, FLOAT pos[2], FLOAT &dist ) {
	FLOAT x = i/(FLOAT)(gn-1);
	FLOAT y = j/(FLOAT)(gn-1);
	dist = 1.0;
	bool updated = false;
	
	// For Neigboring Grids
	int query[][2] = { {i+1,j},{i-1,j},{i,j+1},{i,j-1},{i-1,j-1},{i-1,j+1},{i+1,j-1},{i+1,j+1} };
	for( int q=0; q<8; q++ ) {
		int qi = query[q][0];
		int qj = query[q][1];
		if( qi>=0 && qi<gn && qj>=0 && qj<gn ) {
			FLOAT sgn = grids[qi][qj].dist > 0 ? 1.0 : -1.0;
			// If The Neighborhood Is A Known Grid
			if( grids[qi][qj].known ) {
				// Compute The Distance From The Point Of Its Neighbors
				FLOAT d = hypot(grids[qi][qj].pos[0]-x,grids[qi][qj].pos[1]-y);
				// If The Distance Is Closer Than Holding One, Just Replace
				if( d < fabs(dist) ) {
					dist = sgn*d;
					for( int n=0; n<2; n++ ) pos[n] = grids[qi][qj].pos[n];
					updated = true;
				}
			}
		}
	}
	return updated;
}

static void fastMarch( FLOAT tolerance ) {
	// Unknowns
	grid_queue unknowns;
	
	// Compute First Estimate of Distances
	FOR_EVERY_CELL(gn) {
		FLOAT pos[2];
		FLOAT dist;
		if( ! grids[i][j].known ) {
			grids[i][j].estimated = false;
			if( update_distance(i,j,pos,dist)) {
				grids[i][j].dist = dist;
				if( fabs(dist) < tolerance ) {
					grids[i][j].pos[0] = pos[0];
					grids[i][j].pos[1] = pos[1];
					grids[i][j].estimated = true;
					unknowns.push(&grids[i][j]);
				}
			} else {
				grids[i][j].dist = 1.0;
			}
		}
	} END_FOR;
	
	// While Unknown Grid Exists
	while( ! unknowns.empty() ) {
		// Pop Out Top
		grid *mingrid = unknowns.top();
		int i = mingrid->gp[0];
		int j = mingrid->gp[1];
		unknowns.pop();
		grids[i][j].estimated = false;
		mingrid->known = true;
		
		// Dilate...
		int query[][2] = { {i+1,j},{i-1,j},{i,j+1},{i,j-1},{i-1,j-1},{i-1,j+1},{i+1,j-1},{i+1,j+1} };
		for( int q=0; q<8; q++ ) {
			int qi = query[q][0];
			int qj = query[q][1];
			if( qi>=0 && qi<gn && qj>=0 && qj<gn ) {
				if( ! grids[qi][qj].estimated && ! grids[qi][qj].known ) {
					FLOAT pos[2];
					FLOAT dist;
					if( update_distance(qi,qj,pos,dist)) {
						if( fabs(dist) < tolerance ) {
							grids[qi][qj].dist = dist;
							grids[qi][qj].pos[0] = pos[0];
							grids[qi][qj].pos[1] = pos[1];
							grids[qi][qj].estimated = true;
							unknowns.push(&grids[qi][qj]);
						}
					}
				}
			}
		}
	}

#if 1
	// Fill Inner Region With Scanline Order
	for( int dir=0; dir<2; dir++ ) {
		OPENMP_FOR for( int j=0; j<gn; j++ ) {
			FLOAT sgn = 1.0;
			
			// Forward Search
			for( int i=0; i<gn; i++ ) {
				int idx[][2] = { {i,j}, {j,i} };
				
				if( grids[idx[dir][0]][idx[dir][1]].known ) {
					sgn = grids[idx[dir][0]][idx[dir][1]].dist < 0.0 ? -1.0 : 1.0;
				} else if( sgn < 0.0 ) {
					grids[idx[dir][0]][idx[dir][1]].dist = -1.0;
					grids[idx[dir][0]][idx[dir][1]].known = true;
				}
			}
			
			// Backward Search
			for( int i=gn-1; i>=0; i-- ) {
				int idx[][2] = { {i,j}, {j,i} };
				
				if( grids[idx[dir][0]][idx[dir][1]].known ) {
					sgn = grids[idx[dir][0]][idx[dir][1]].dist < 0.0 ? -1.0 : 1.0;
				} else if( sgn < 0.0 ) {
					grids[idx[dir][0]][idx[dir][1]].dist = -1.0;
					grids[idx[dir][0]][idx[dir][1]].known = true;
				}
			}
		}
	}
#endif
}

inline int clamp( int a ) {
	return max(min(a,gn-1),0);
}

static void computeDistGradient() {
	FOR_EVERY_CELL(gn) {
		grad[0][i][j] = grad[1][i][j] = 0.0;
		if( grids[i][j].known &&
		   grids[clamp(i+1)][j].known && grids[clamp(i-1)][j].known &&
		   grids[i][clamp(j+1)].known && grids[i][clamp(j-1)].known ) {
			grad[0][i][j] = (grids[clamp(i+1)][j].dist-grids[clamp(i-1)][j].dist)*(gn-1);
			grad[1][i][j] = (grids[i][clamp(j+1)].dist-grids[i][clamp(j-1)].dist)*(gn-1);
		}
		FLOAT d = hypot(grad[0][i][j],grad[1][i][j]);
		if( d > 0.0 ) {
			grad[0][i][j] /= d;
			grad[1][i][j] /= d;
		}
	} END_FOR;
}

inline FLOAT sd( int i, int j ) {
	return (i>=0 && i<gn && j>=0 && j<gn) ? grids[i][j].dist : 1.0;
}

void levelset2D::extrapolate( FLOAT **q, char **region ) {	
	// Unknowns
	grid_queue unknowns;
	
	// Computed Field
	static char **computed = alloc2D<char>(gn);
	
	// Push
	FOR_EVERY_CELL(gn) {
		if( ! region[i][j] && grids[i][j].known ) {
			unknowns.push(&grids[i][j]);
		}
		computed[i][j] = region[i][j];
	} END_FOR;

	// While Unknowns Exists
	while( ! unknowns.empty() ) {
		// Pop Out Top
		grid *mingrid = unknowns.top();
		int i = mingrid->gp[0];
		int j = mingrid->gp[1];
		
		int sum = 0;
		FLOAT sumq = 0.0;
		int query[][2] = { {i+1,j},{i-1,j},{i,j+1},{i,j-1},{i-1,j-1},{i-1,j+1},{i+1,j-1},{i+1,j+1} };
		for( int nq=0; nq<8; nq++ ) {
			int qi = query[nq][0];
			int qj = query[nq][1];
			if( qi>=0 && qi<gn && qj>=0 && qj<gn ) {
				if( computed[qi][qj] ) {
					sumq += q[qi][qj];
					sum ++;
				}
			}
		}
		if( sum ) {
			q[i][j] = sumq / sum;
		}
		unknowns.pop();
		computed[i][j] = 1;
	}
}

static void intersect( FLOAT x1, FLOAT y1, FLOAT x2, FLOAT y2, FLOAT &x, FLOAT &y  ) {
	FLOAT d = hypot2(x2-x1,y2-y1);
	FLOAT u = ((x-x1)*(x2-x1) + (y-y1)*(y2-y1))/d;
	u = fmin(1.0,fmax(0.0,u));
	x = x1 + u*(x2-x1);
	y = y1 + u*(y2-y1);
}

void levelset2D::redistance( FLOAT tolerance ) {
	FLOAT w = 1.0/(gn-1);
	
	// Make Everything Known First
	FOR_EVERY_CELL(gn) {
		if( ! grids[i][j].known ) grids[i][j].dist = 1.0;
		grids[i][j].known = true;
	} END_FOR;
	
	FOR_EVERY_CELL(gn) {
		bool farCell = true;
		int query[][2] = { {i+1,j},{i+1,j+1},{i,j+1},{i-1,j+1},{i-1,j},{i-1,j-1},{i,j-1},{i+1,j-1} };
		int qnum = 8;
		FLOAT pos[qnum][2];
		bool fnd[qnum];
		for( int q=0; q<qnum; q++ ) {
			int qi = query[q][0];
			int qj = query[q][1];
			fnd[q] = false;
			if( qi>=0 && qi<gn && qj>=0 && qj<gn ) {
				if( grids[qi][qj].known && grids[qi][qj].dist * grids[i][j].dist < 0.0 ) {
					farCell = false;
					
					// Calculate New Cross Position
					FLOAT y0 = grids[i][j].dist;
					FLOAT y1 = grids[qi][qj].dist;
					FLOAT a = y0/(y0-y1);
					FLOAT p0[2] = { w*i, w*j };
					FLOAT p1[2] = { w*qi, w*qj };
					pos[q][0] = (1.0-a)*p0[0]+a*p1[0];
					pos[q][1] = (1.0-a)*p0[1]+a*p1[1];
					fnd[q] = true;
				}
			}
		}
		if( farCell ) {
			grids[i][j].known = false;
		} else {
			grids[i][j].known = true;
			FLOAT mind = 9999.0;
			for ( int q=0; q<qnum; q++ ) {
				for( int rp=1; rp<3; rp++ ) {
					if( fnd[(q+rp)%qnum] && fnd[q] ) {
						FLOAT x = i*w;
						FLOAT y = j*w;
						intersect( pos[(q+rp)%qnum][0], pos[(q+rp)%qnum][1], pos[q][0], pos[q][1], x, y );
						FLOAT d = hypot(x-i*w,y-j*w);
						if( d < mind ) {
							mind = d;
							grids[i][j].pos[0] = x;
							grids[i][j].pos[1] = y;
							grids[i][j].dist = (grids[i][j].dist > 0.0 ? 1.0 : -1.0)*fmax(d,1.0e-8);
						}
					}
				}
			}
			
			for ( int q=0; q<qnum; q++ ) {
				if( !fnd[q] && fnd[(q+1)%qnum] && !fnd[(q+2)%qnum] ) {
					FLOAT x = pos[(q+1)%qnum][0];
					FLOAT y = pos[(q+1)%qnum][1];
					FLOAT d = hypot(x-i*w,y-j*w);
					if( d < mind ) {
						mind = d;
						grids[i][j].pos[0] = x;
						grids[i][j].pos[1] = y;
						grids[i][j].dist = (grids[i][j].dist > 0.0 ? 1.0 : -1.0)*fmax(d,1.0e-8);
					}
				}
			}
		}
	} END_FOR;
	
	FOR_EVERY_CELL(gn) {
		if( ! grids[i][j].known ) grids[i][j].dist = 1.0;
		
	} END_FOR;
	fastMarch(tolerance);
	
#if 0 // Just For A Debug
	static char **region = alloc2D<char>(gn);
	FOR_EVERY_CELL(gn) {
		region[i][j] = grids[i][j].dist <= 0.0;
		test_q[i][j] = region[i][j] ? 0.5 : 0.0;
	} END_FOR;
	extrapolate( test_q, region );
#endif
}

void levelset2D::buildLevelset( bool (*func)(FLOAT x, FLOAT y), FLOAT tolerance ) {
	// Initialize Distances
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		FLOAT x = i/(FLOAT)(gn-1);
		FLOAT y = j/(FLOAT)(gn-1);
		FLOAT sx = x;
		FLOAT sy = y;
		FLOAT hit;
		hit = func(sx,sy);
		grids[i][j].known = true;
		grids[i][j].dist = hit ? -1.0 : 1.0;
		grids[i][j].pos[0] = 0.0;
		grids[i][j].pos[1] = 0.0;
	} END_FOR;
	
	// Redistance...
	redistance(tolerance);
}

void levelset2D::advect( void (*func)( FLOAT x, FLOAT y, FLOAT &u, FLOAT &v, FLOAT &dt ) ) {	
	static FLOAT ** source_dists = alloc2D<FLOAT>(gn);
	static FLOAT ** swap_dists = alloc2D<FLOAT>(gn);
	
	// Copy Dists
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		source_dists[i][j] = grids[i][j].dist;
	} END_FOR;
	
	// Advect
	FLOAT w = 1.0/(gn-1);
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		FLOAT x, y;
		FLOAT u, v, dt;
		x = i*w;
		y = j*w;
		func( x, y, u, v, dt );
		// Semi-Lagragian
		x -= dt*u;
		y -= dt*v;
		x = fmin(1.0,fmax(0.0,x));
		y = fmin(1.0,fmax(0.0,y));
		// Interpolate Dists
		if( grids[(int)((gn-1)*x)][(int)((gn-1)*y)].known ) {
			swap_dists[i][j] = y<1.0 ? interp::interp( source_dists, (gn-1)*x, (gn-1)*y, gn, gn ) : 1.0;
		} else {
			swap_dists[i][j] = source_dists[i][j];
		}
	} END_FOR;
	
	// Swap
	FOR_EVERY_CELL(gn) {
		grids[i][j].dist = swap_dists[i][j];
	} END_FOR;
}

void levelset2D::keyDown( char key ) {
}

static void calcMarchingPoints( int i, int j, FLOAT p[8][2], int &pnum ) {
	pnum = 0;
	FLOAT w = 1.0/(gn-1);
	int quads[][2] = { {i, j}, {i+1, j}, {i+1, j+1}, {i, j+1} };
	for( int n=0; n<4; n++ ) {
		if( grids[quads[n][0]][quads[n][1]].known ) {
			// Inside Liquid
			if( grids[quads[n][0]][quads[n][1]].dist < 0.0 ) {
				p[pnum][0] = w*quads[n][0];
				p[pnum][1] = w*quads[n][1];
				pnum ++;
			}
			// If Line Crossed
			if( grids[quads[n][0]][quads[n][1]].dist * grids[quads[(n+1)%4][0]][quads[(n+1)%4][1]].dist < 0 ) {
				// Calculate Cross Position
				FLOAT y0 = grids[quads[n][0]][quads[n][1]].dist;
				FLOAT y1 = grids[quads[(n+1)%4][0]][quads[(n+1)%4][1]].dist;
				FLOAT a = y0/(y0-y1);
				FLOAT p0[2] = { w*quads[n][0], w*quads[n][1] };
				FLOAT p1[2] = { w*quads[(n+1)%4][0], w*quads[(n+1)%4][1] };
				p[pnum][0] = (1.0-a)*p0[0]+a*p1[0];
				p[pnum][1] = (1.0-a)*p0[1]+a*p1[1];
				pnum ++;
			}
		}
	}	
}

static void drawMarchingCube() {
	// Paint Distance Field
	glColor4f(0.5,0.6,1.0,0.5);
	FOR_EVERY_CELL(gn-1) {
		FLOAT p[8][2];
		int pnum;
		glBegin(GL_TRIANGLE_FAN);
		calcMarchingPoints( i, j, p, pnum );
		for( int m=0; m<pnum; m++ ) {
			glVertex2f(p[m][0],p[m][1]);
		}
		glEnd();
	} END_FOR;
}

void levelset2D::display( bool cell_centered ) {	
	FLOAT w = 1.0/(gn-1);

	// Cell Centered
	glPushMatrix();
	if( cell_centered ) {
		FLOAT s = (gn-1)/(FLOAT)gn;
		glTranslated( w/2, w/2, 0.0 );
		glScaled(s,s,s); 
	}
	
	if( show_dist ) {
		// Paint Distance Field
		FOR_EVERY_CELL(gn-1) {
			int quads[] = { i, j, i+1, j, i+1, j+1, i, j+1 };
			glBegin(GL_QUADS);
			for( int q=0; q<4; q++ ) {
				int fi = quads[2*q+0];
				int fj = quads[2*q+1];
				FLOAT alpha = grids[fi][fj].known ? 1.0 : 0.0;
				FLOAT d = grids[fi][fj].dist;
				FLOAT s = 5.0;
				glColor4f( s*d*(d>0), 0.0, s*(-d)*(d<0), alpha );
				glVertex2d( fi*w, fj*w );
			}
			glEnd();
		} END_FOR;
	}
	
	if( show_region ) drawMarchingCube();
	
#if 0
	// Paint Test Quantity
	FOR_EVERY_CELL(gn-1) {
		int quads[] = { i, j, i+1, j, i+1, j+1, i, j+1 };
		glBegin(GL_QUADS);
		for( int q=0; q<4; q++ ) {
			int fi = quads[2*q+0];
			int fj = quads[2*q+1];
			FLOAT d = grids[fi][fj].dist;
			FLOAT s = 1.0/max_d;
			FLOAT q = test_q[fi][fj];
			glColor4f( 1.0, 1.0, 1.0, q );
			glVertex2d( fi*w, fj*w );
		}
		glEnd();
	} END_FOR;
#endif
	
	if( show_grid ) {
		glPointSize(2);
		// Plot Grid Points
		FOR_EVERY_CELL(gn) {
			grids[i][j].dist < 0 ? glColor4f(0.0,0.0,1.0,1.0) : glColor4f(1.0,0.0,0.0,1.0);
			if( ! grids[i][j].known ) glColor4f(0.5,0.5,0.5,0.5);
			glBegin(GL_POINTS);
			glVertex2f(i*w,j*w);
			glEnd();
		} END_FOR;
		glPointSize(1);
	}
	
	// Draw Surface Reference Line
#if 0
	FOR_EVERY_CELL(gn) {
		if( grids[i][j].known ) {
			FLOAT alpha = 0.5;
			if( grids[i][j].dist > 0.0 ) glColor4f(1.0,1.0,0.0,alpha);
			else glColor4f(0.0,1.0,1.0,alpha);
			
			FLOAT *pos = grids[i][j].pos;
			
			glBegin(GL_LINES);
			glVertex2f(i*w,j*w);
			glVertex2f(pos[0],pos[1]);
			glEnd();
		}
	} END_FOR;
#endif
	
#if 0
	// Draw Gradient Of Distance Field
	glColor4f(1.0,1.0,1.0,0.5);
	FOR_EVERY_CELL(gn) {
		glBegin(GL_LINES);
		glVertex2f(i*w,j*w);
		glVertex2f(i*w+grad[0][i][j]*w,j*w+grad[1][i][j]*w);
		glEnd();
	} END_FOR;
#endif
	
#if 0
	glColor4f(1.0,1.0,1.0,1.0);
	FLOAT p[3][2] = { {0.1,0.1}, {0.4,0.8}, {0.7,0.3} };
	glBegin(GL_LINE_LOOP);
	for( int i=0; i<3; i++ ) {
		glVertex2fv(p[i]);
	}
	glEnd();
	
	FLOAT x = p[2][0];
	FLOAT y = p[2][1];
	intersect( p[0][0], p[0][1], p[1][0], p[1][1], x, y  );
	glBegin(GL_LINES);
	glVertex2f( x, y );
	glVertex2fv(p[2]);
	glEnd();
#endif
	
	glPopMatrix();
}

FLOAT levelset2D::getLevelSet( int i, int j ) {
	if( i<0 || i>gn-1 || j<0 || j>gn-1 ) return 1.0;
	return grids[i][j].dist;
}

void levelset2D::getLevelSet( FLOAT **dists ) {
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		dists[i][j] = grids[i][j].known ? grids[i][j].dist : 1.0;
	} END_FOR;
}

FLOAT levelset2D::getVolume() {
	FLOAT volume = 0.0;
	FOR_EVERY_CELL(gn-1) {
		FLOAT p[8][2];
		int pnum;
		calcMarchingPoints( i, j, p, pnum );
		for( int m=0; m<pnum; m++ ) {
			volume += p[m][0]*p[(m+1)%pnum][1]-p[m][1]*p[(m+1)%pnum][0];
		}
	} END_FOR;
	
	return volume*0.5;
}

void levelset2D::setLevelSet( int i, int j, FLOAT d ) {
	grids[i][j].dist = d;
	grids[i][j].known = true;
}

void levelset2D::setVisibility( bool grid, bool dist, bool region ) {
	show_grid = grid;
	show_dist = dist;
	show_region = region;
}
