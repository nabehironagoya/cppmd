cmake_minimum_required(VERSION 2.8.12)
project(cppmd CXX)

add_executable(cppmd mdrun.cpp string.cpp gro_util.cpp math.cpp)
