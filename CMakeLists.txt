cmake_minimum_required(VERSION 3.5)
project(clothoid)

set(CMAKE_CXX_STANDARD 23)
#set(CMAKE_CXX_FLAGS "-fconcepts-diagnostics-depth=2")

find_package(Eigen3 REQUIRED)

add_executable(clothoid main.cpp
        clothoid.h
        gauss-legendre.h)

target_include_directories(clothoid PUBLIC "${EIGEN3_INCLUDE_DIRS}")