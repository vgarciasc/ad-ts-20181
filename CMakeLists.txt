cmake_minimum_required(VERSION 3.0)
project(AD)

set(CMAKE_CXX_STANDARD 17)

add_executable(AD main.cpp Event.cpp Event.h Packet.cpp Packet.h)