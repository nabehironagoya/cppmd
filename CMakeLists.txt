cmake_minimum_required(VERSION 2.8.12)
project(cppmd CXX)

add_executable(cppmd mdrun.cpp string.cpp read_gro.cpp math.cpp)
