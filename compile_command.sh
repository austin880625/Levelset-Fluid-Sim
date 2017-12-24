#!/bin/sh
g++ common/shader.cpp CubeMarching.cpp tri4.cpp -lGL -lGLU -lglfw -lGLEW -lX11 -lXxf86vm -lXrandr -lpthread -ldl -lXi  -lXinerama -lXcursor --std=c++11
