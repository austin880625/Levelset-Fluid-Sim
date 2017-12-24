#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include <stdio.h>

namespace CubeMarching{
	void init();
	void genVertices(GLfloat ***grid, GLfloat length, int w, int h, int d);
	void draw();
}
