CXX = g++
LIBS = -lGL -lGLU -lGLEW -lX11 -lXxf86vm -lXrandr -lpthread -ldl -lXi  -lXinerama -lXcursor
OPT = --std=c++11
GLFW = -lglfw
GLFW3 = -lglfw3

all: gg2
gg1: CubeMarching.cpp common/shader.cpp tri4.cpp
	$(CXX) -o bin/fluid-sim CubeMarching.cpp common/shader.cpp tri4.cpp $(LIBS) $(GLFW) $(OPT)
gg2: CubeMarching.cpp common/shader.cpp tri4.cpp
	$(CXX) -o bin/fluid-sim CubeMarching.cpp common/shader.cpp tri4.cpp $(LIBS) $(GLFW3) $(OPT)
