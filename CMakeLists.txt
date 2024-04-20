cmake_minimum_required(VERSION 3.24)
project(23_Home_2)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
        Matrix.hpp
        main.cpp)
add_executable(23_Home_2 ${SOURCE_FILES})

include_directories(.)