cmake_minimum_required(VERSION 3.5)
project(clothoid)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)

add_executable(clothoid main.cpp
        include/clothoid.h
        include/gauss-legendre.h
        include/vec2.h
        include/distance.h
        include/shape.h
        include/common.h
)

target_include_directories(clothoid PUBLIC "./include" "${EIGEN3_INCLUDE_DIRS}")