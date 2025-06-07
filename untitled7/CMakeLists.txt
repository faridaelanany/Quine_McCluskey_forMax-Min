cmake_minimum_required(VERSION 3.29)
project(untitled7)

set(CMAKE_CXX_STANDARD 20)

add_executable(untitled7 cmake-build-debug/main.cpp
        cmake-build-debug/qm.cpp
        cmake-build-debug/qm.h
        qm-test.cpp
)




