#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GLM_ENABLE_EXPERIMENTAL

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include "common/shader.hpp"
#include "levelset/utility.h"
#include "levelset/controller.h"

GLFWwindow *window;
int screen_width = 1024;
int screen_height = 600;
int gn = 32;

using namespace glm;

int main()
{
	if(glfwInit()==0)
	{
		fprintf(stderr, "GLFW initialization failed. exit.");
		getchar();
		return -1;
	}
	
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	window = glfwCreateWindow(screen_width, screen_height, "Trivago", NULL, NULL);

	glfwMakeContextCurrent(window);
	glewExperimental = true; // Needed for core profile
	if(glewInit() != GLEW_OK){
		fprintf(stderr, "GLEW initialization failed. exit.");
		getchar();
		glfwTerminate();
	}
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	GLuint programID = LoadShaders( "SimpleVertexShader.vertexshader", "SimpleFragmentShader.fragmentshader" );

	/*
	int maxn=129;
	sp = (GLfloat***)malloc((maxn+1)*sizeof(GLfloat**));sp[maxn]=NULL;
	for(int i=0; i<maxn; i++){
		sp[i] = (GLfloat**)malloc((maxn+1)*sizeof(GLfloat*));sp[i][maxn]=NULL;
		for(int j=0; j<maxn; j++){
			sp[i][j] = (GLfloat*)malloc(maxn*sizeof(GLfloat));
			for(int k=0; k<maxn; k++){
				GLfloat d=14.0f/(maxn-1); 
				GLfloat x=-7.0f+i*d, y=-7.0f+j*d, z=-7.0f+k*d;
				//GLfloat x=-1.0f+i*d, y=-1.0f+j*d, z=-1.0f+k*d;
				//sp[i][j][k] = x*x+y*y+z*z-1;
				//printf("%f %f %f\n",x,y,z);
				float r=sqrt(x*x+y*y+z*z);
				float phi = 0.204*r*exp(-r/2)*(0.488*z/r);
				sp[i][j][k] = phi*phi>0.001 ? -1 :1;
				//printf("%f\n",sp[i][j][k]);
			}
		}
	}
	*/
	controller::init(gn);

	glClearColor(0, 0, 0, 0);
	
	mat4 projection = perspective(radians(45.0f), (float)screen_width/(float)screen_height, 0.1f, 100.0f);	
	
	GLuint mvpMatrixID = glGetUniformLocation(programID, "mvp");
	GLuint mMatrixID = glGetUniformLocation(programID, "m");
	GLuint vMatrixID = glGetUniformLocation(programID, "v");
	GLuint lightVecID = glGetUniformLocation(programID, "lightPos");
	int angle=90;
	do{
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		
		//angle = (angle+2)%360;
		mat4 view = lookAt(vec3(0,2,4), vec3(0,0,0), vec3(0,1,0));
		mat4 model = mat4(1.0f);
		mat4 mvp = projection*view*model;

		glUniformMatrix4fv(mvpMatrixID, 1, GL_FALSE, &mvp[0][0]);
		glUniformMatrix4fv(mMatrixID, 1, GL_FALSE, &model[0][0]);
		glUniformMatrix4fv(vMatrixID, 1, GL_FALSE, &view[0][0]);
		vec3 light(5,5,5);
		glUniform3f(lightVecID, light.x, light.y, light.z);
		/*
		num_tri = (num_tri+1)%5+1;
		for(int i=1; i<num_tri; i++){
			for(int j=0; j<6; j++){
				g_vertex_buffer_data[9*i+j] = g_vertex_buffer_data[9*i-6+j];
			}
			g_vertex_buffer_data[9*i+6] = g_vertex_buffer_data[9*i]+2.0f;
			g_vertex_buffer_data[9*i+7] = g_vertex_buffer_data[9*i+1];
			g_vertex_buffer_data[9*i+8] = 0.0f;
		}
		*/
		glUseProgram(programID);
		
		controller::display();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}while( glfwGetKey(window, GLFW_KEY_ESCAPE)!=GLFW_PRESS&&glfwWindowShouldClose(window)==0);

	glfwTerminate();
	return 0;
}
