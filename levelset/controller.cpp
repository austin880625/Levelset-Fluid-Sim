/*
 *  controller.cpp
 *
 */

#include "controller.h"
#include "liquid2D.h"
#include "profile.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <sys/time.h>
#elif defined(WIN32)
#include "glut.h"
#include <windows.h>
#else
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <sys/time.h>
#endif

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define TEST	0

void controller::init( int gn ) {
	// Turn On Blending
	//glEnable(GL_BLEND);
	//glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	
	// Initialize Liquid
	liquid2D::init(gn);
}
/*
void controller::reshape( int w, int h ) {
	double margin = 0.02;
	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(-margin,1.0+margin,-margin,1.0+margin,-1.0,1.0);
}
*/
void controller::display() {
	liquid2D::display();
}

void controller::keyDown( unsigned char key ) {
	liquid2D::keyDown(key);
}

void controller::mouse( double x, double y, int state ) {
#if TEST
#else
	liquid2D::mouse(x,y,state);
#endif	
}

void controller::motion( double x, double y, double dx, double dy ) {
#if TEST
#else
	liquid2D::motion(x,y,dx,dy);
#endif	
}












