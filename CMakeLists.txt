
# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(rcbkdipole CXX C)

# std::make_unique (gamma_aa.cpp, amplitudelib.cpp) needs C++14; set this
# explicitly instead of relying on the compiler's default, since that
# default varies (e.g. Puhti's gcc/11.3.0 module defaults to C++17, but
# older system compilers may default to pre-C++11).
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization on (for debug flags, change to Debug)
# Release or Debug
set(CMAKE_BUILD_TYPE Release)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_subdirectory(src)

set_target_properties(amplitude PROPERTIES
         LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

set_target_properties(dipole PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

