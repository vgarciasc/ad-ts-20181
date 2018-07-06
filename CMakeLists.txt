cmake_minimum_required(VERSION 3.10)
project(AD)

set(CMAKE_CXX_STANDARD 17)

add_executable(AD main.cpp Event.cpp Event.h Stats.cpp Stats.h)