/*
 *  controller.h
 *
 */
#include <GLFW/glfw3.h>
namespace controller {
	void init( int gn );
	void reshape( int w, int h );
	void display(GLFWwindow *window);
	void mouse( double x, double y, int state );
	void motion( double x, double y, double dx, double dy );
	void keyDown( unsigned char key );
}
