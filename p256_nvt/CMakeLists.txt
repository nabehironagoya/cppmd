cmake_minimum_required(VERSION 2.8.12)
project(cppmd CXX)

add_executable(cppmd
               mdrun.cpp 
			   string.cpp 
			   gro_util.cpp 
			   math.cpp)

set(CMAKE_CXX_FLAGS "-Wall -fopenmp -std=c++11")

option(OPENMP "use openmp" OFF)
if(OPENMP)
		add_definitions(-DOPENMP)
endif()

