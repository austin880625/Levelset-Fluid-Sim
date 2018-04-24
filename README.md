# Levelset Fluid Simulation
A fluid simulation program that uses levelset and cube-marching algorithm, and openMP, openCL for parallelization.

Modified from https://code.google.com/archive/p/levelset2d/

Versions using different parallel libraries can be found in corresponding branches.

## Usage
There are two targets in the makefile: gg1 and gg2, depends on the name of your OpenGL3 library in the system.
Some of feature in the original project is removed due to lack of time to modifying

## Screenshots
![Alt text](/screenshot.png "screenshot")

## About architecture of the code
The code is messy since we started this project as we just learned OpenGL and OpenCL. And hence we do some clarification(may not be useful, XD)

The `tri4.cpp` is the actual `main.cpp` file, there is only two shaders, and the `main.cpp` in the `levelset` directory is not used.
