#include "controller.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
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

int win_x = 512;
int win_y = 512;
int prev_x, prev_y;
int mstat = 1;
GLFWwindow *pWindow;

static void display(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glfwPollEvents();
	controller::display(pWindow);
	glfwSwapBuffers(pWindow);
}

static void init( int gn ) {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	controller::init(gn);
}

static void reshape(GLFWwindow *window, int w, int h) {
	win_x = w;
	win_y = h;
	controller::reshape(w,h);
}

static void idle ( void ) {
	//glutPostRedisplay ();
}

static void keyboard( GLFWwindow *window, int key_int, int scancode, int action, int mods ) {
	if(action!=GLFW_PRESS)return ;
	unsigned char key = (unsigned char)key_int;
	//printf("%c\n",key);
	if( key == '\e' ) exit(0);
	controller::keyDown(key);
	//glutPostRedisplay ();
}

static void mouse ( GLFWwindow *window, int button, int state, int mods ) {
	//prev_x = x;
	//prev_y = y;
	mstat = state;
	controller::mouse( prev_x/(GLdouble)win_x, 1.0 - prev_y/(GLdouble)win_y,  state );
}

static void motion ( GLFWwindow *window, double x, double y ) {
	if( mstat == 1 ) {
		controller::motion( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y,
						 (x-prev_x)/(GLdouble)win_x, -(y-prev_y)/(GLdouble)win_y );
	}
	prev_x = x;
	prev_y = y;
}


int main (int argc, char * argv[]) {
	
	int grid_size = 64;
	if( argc == 2  ) {
		sscanf( argv[1], "%d", &grid_size );
	}
	
	//glutInit(&argc, argv);
	if(!glfwInit()){
		fprintf(stderr, "Error initializing GLFW\n");
		exit(1);
	}
#if _OPENMP
	printf( "Number of threads: %d\n", omp_get_num_procs() );
#endif
	/*
	glutInitDisplayMode(GLUT_RGBA | GL_DOUBLE);
	glutInitWindowPosition ( 100, 100 );
	glutInitWindowSize ( win_x, win_y );
	glutCreateWindow(argv[0]);
	*/
	pWindow = glfwCreateWindow(win_x, win_y, argv[0], NULL, NULL);
	if(!pWindow){
		fprintf(stderr, "failed to open GLFW window\n");
		getchar();
		glfwTerminate();
		exit(1);
	}
	glfwMakeContextCurrent(pWindow);

	//glutIdleFunc(idle);
	glfwSetKeyCallback(pWindow, keyboard);
	glfwSetMouseButtonCallback(pWindow, mouse);
	glfwSetCursorPosCallback (pWindow, motion);
	glfwSetWindowSizeCallback(pWindow, reshape);
	init(grid_size);
	reshape(pWindow, win_x, win_y);
	while(!glfwWindowShouldClose(pWindow)){
		display();
	}
	return 0;
}
