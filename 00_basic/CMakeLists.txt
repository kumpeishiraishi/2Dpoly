cmake_minimum_required(VERSION 3.13)

project(my_program CXX)

add_executable(a.out main.cpp)

target_compile_options(a.out PRIVATE -O3)
target_compile_features(a.out PRIVATE cxx_std_17)
