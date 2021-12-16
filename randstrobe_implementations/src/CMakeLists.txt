cmake_minimum_required(VERSION 3.19)
project(src)

set(CMAKE_CXX_STANDARD 14)

add_executable(src
               main.cpp
               index.cpp
        index.hpp
        robin_hood.h
        xxhash.c
        xxhash.h)


include_directories(src)

