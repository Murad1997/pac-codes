cmake_minimum_required(VERSION 3.5)

project(llrscl_decoder LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
file(GLOB hdr src/*.cpp)
#set(Armadillo_USE_STATIC_LIBS OFF)
#find_package(Armadillo REQUIRED COMPONENTS )
#include_directories(${Armadillo_INCLUDE_DIR})

add_executable(llrscl ${hdr})

#target_compile_options(astar PRIVATE -larmadillo)
target_include_directories(llrscl PRIVATE include)
target_link_libraries(llrscl -larmadillo)

