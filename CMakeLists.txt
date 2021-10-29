# CMAKE minimum version
cmake_minimum_required(VERSION 3.17)

# Setup project
project("ADCS Math")

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(POLICY CMP0091)
    # https://stackoverflow.com/a/56490614
    cmake_policy(SET CMP0091 NEW)
endif()

# Clump source files from src
file(GLOB ADCS_SRC "src/*.cpp" "src/**/*.cpp")

# Specify include directories
include_directories(include)

file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/assets)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/assets)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/assets DESTINATION ${CMAKE_CURRENT_BINARY_DIR})