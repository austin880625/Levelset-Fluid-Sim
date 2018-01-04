CXX = g++
LIBS = -lGL -lGLU -lGLEW -lX11 -lXxf86vm -lXrandr -lpthread -ldl -lXi  -lXinerama -lXcursor -fopenmp -lOpenCL
OPT = --std=c++11
GLFW = -lglfw
GLFW3 = -lglfw3

all: gg2

interp: levelset/interp.cpp
	$(CXX) -c -o interp.o levelset/interp.cpp

liquid: levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp
	make interp
	$(CXX) -c interp.o levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp $(LIBS) $(OPT)

gg1: common/shader.cpp tri4.cpp
	make liquid
	$(CXX) -o bin/fluid-sim *.o common/shader.cpp levelset/controller.cpp tri4.cpp $(LIBS) $(GLFW) $(OPT)
gg2: common/shader.cpp tri4.cpp
	make liquid
	$(CXX) -o bin/fluid-sim *.o common/shader.cpp levelset/controller.cpp tri4.cpp $(LIBS) $(GLFW3) $(OPT)
