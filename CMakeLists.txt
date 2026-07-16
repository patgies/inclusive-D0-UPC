
# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(rcbkdipole CXX C)

# Optimization on (for debug flags, change to Debug)
# Release or Debug
set(CMAKE_BUILD_TYPE Release)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
execute_process(COMMAND ${LHAPDF_CONFIG} --incdir
        OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LHAPDF_CONFIG} --libdir
        OUTPUT_VARIABLE LHAPDF_LIBRARY_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${LHAPDF_INCLUDE_DIR})

add_subdirectory(src)

set_target_properties(amplitude PROPERTIES
         LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

set_target_properties(dipole PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

