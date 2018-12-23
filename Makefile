CXX = g++
CXX_MAC = g++-7
LIBS = -lGL -lGLU -lGLEW -lX11 -lXxf86vm -lXrandr -lpthread -ldl -lXi  -lXinerama -lXcursor
LIBS_MAC = -framework OpenGL -lGLEW -lpthread -ldl
OPT = --std=c++11
GLFW = -lglfw
GLFW3 = -lglfw3

all: gg2

interp: levelset/interp.cpp
	$(CXX) -c -o interp.o levelset/interp.cpp

liquid: levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp
	make interp
	$(CXX) -c levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp levelset/common.cpp $(LIBS) $(OPT)

gg1: common/shader.cpp tri4.cpp
	make liquid
	mkdir -p bin
	$(CXX) -o bin/fluid-sim *.o common/shader.cpp levelset/controller.cpp tri4.cpp $(LIBS) $(GLFW) $(OPT)
gg2: common/shader.cpp tri4.cpp
	make liquid
	mkdir -p bin
	$(CXX) -o bin/fluid-sim *.o common/shader.cpp levelset/controller.cpp tri4.cpp $(LIBS) $(GLFW3) $(OPT)

liquid_mac: levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp
	make interp
	$(CXX_MAC) -c levelset/solver.cpp levelset/CubeMarching.cpp levelset/levelset2D.cpp levelset/liquid2D.cpp levelset/common.cpp $(LIBS_MAC) $(OPT)
mac: common/shader.cpp tri4.cpp
	make liquid_mac
	mkdir -p bin
	$(CXX_MAC) -o bin/fluid-sim *.o common/shader.cpp levelset/controller.cpp tri4.cpp $(LIBS_MAC) $(GLFW) $(OPT)
