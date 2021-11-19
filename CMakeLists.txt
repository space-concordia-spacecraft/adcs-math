# CMAKE minimum version
cmake_minimum_required(VERSION 3.17)

# Setup project
project(ADCS_MATH)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 14)

set(ADCS_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/AstroLib.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/DTM_12.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/EopSpw.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/MathTimeLib.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/MSIS_Vers.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/sgp4/SGP4.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/sgp4/SGP4TJK.cpp")

set(ADCS_INCLUDE
        "${CMAKE_CURRENT_SOURCE_DIR}/extern/kronos-osal/types"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/utils"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/sgp4")


add_executable(ADCS_MATH "${ADCS_SRC}")

target_include_directories(ADCS_MATH PRIVATE "${ADCS_INCLUDE}")