cmake_minimum_required(VERSION 3.5)

# set the project name and version
project(viterbi3 CXX)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

set(CMAKE_C_COMPILER /opt/intel/oneapi/compiler/latest/bin/icx)
set(CMAKE_CXX_COMPILER /opt/intel/oneapi/compiler/latest/bin/icpx)

# Set the output directory for the release build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ../bin/release)

# Set the output directory for the debug build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ../bin/debug)

message(STATUS "Building: ${CMAKE_BUILD_TYPE}")

# system
set(CMAKE_SYSTEM_NAME Linux)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)

# c++
set(CMAKE_CXX_FLAGS_DEBUG "-qipp=common -axCORE-AVX2 -xCORE-AVX2 -D_DEBUG -debug parallel -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ipo -qipp=common -axCORE-AVX2 -xCORE-AVX2")

# Linker
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-qipp=common -debug parallel -D_DEBUG -O0 -g")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-O3 -qipp=common")

# #source files
FILE(GLOB  Sources ${PROJECT_SOURCE_DIR}/src/*.c??)

add_custom_target(update_version
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# add the executable
add_executable(${PROJECT_NAME} ${Sources})


SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
# libs
target_include_directories(${PROJECT_NAME} PUBLIC src /opt/intel/oneapi/ipp/latest/include/ipp)
target_link_libraries(${PROJECT_NAME} uhd pthread boost_system boost_program_options)
