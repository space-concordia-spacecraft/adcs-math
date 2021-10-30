# CMAKE minimum version
cmake_minimum_required(VERSION 3.17)

# Setup project
project(ADCS_MATH)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)

set(ADCS_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

set(ADCS_INCLUDE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")


add_executable(ADCS_MATH "${ADCS_SRC}")

target_include_directories(ADCS_MATH PRIVATE "${ADCS_INCLUDE}")