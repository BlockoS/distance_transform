cmake_minimum_required (VERSION 3.20)

project(distance C)

file(WRITE ${CMAKE_BINARY_DIR}/distance.c "#define DISTANCE_TRANSFORM_IMPLEMENTATION\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/distance.h\"")
add_library(distance STATIC ${CMAKE_BINARY_DIR}/distance.c distance.h)
target_include_directories(distance PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_compile_options(distance PRIVATE $<IF:$<C_COMPILER_ID:MSVC>,/W4,-Wall -Wshadow -Wextra>)
set_property(TARGET distance PROPERTY C_STANDARD 11)
add_subdirectory(test)
