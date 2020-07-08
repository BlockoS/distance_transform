cmake_minimum_required (VERSION 3.12)

project(distance C)

find_package(Doxygen)

file(WRITE ${CMAKE_BINARY_DIR}/distance.c "#define DISTANCE_TRANSFORM_IMPLEMENTATION\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/distance.h\"")
add_library(distance STATIC ${CMAKE_BINARY_DIR}/distance.c)
target_compile_options(distance PUBLIC -Wall -Wshadow -Wextra)

add_subdirectory(test)
