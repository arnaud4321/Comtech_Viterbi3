cmake_minimum_required(VERSION 3.5)

# Compilateur Intel oneAPI optionnel (à définir avant project())
set(INTEL_COMPILER_ROOT "/opt/intel/oneapi/compiler/latest" CACHE PATH "Intel oneAPI compiler installation directory")
if(EXISTS "${INTEL_COMPILER_ROOT}/bin/icpx")
	set(CMAKE_C_COMPILER "${INTEL_COMPILER_ROOT}/bin/icx")
	set(CMAKE_CXX_COMPILER "${INTEL_COMPILER_ROOT}/bin/icpx")
else()
	message(STATUS "Intel compiler not found at ${INTEL_COMPILER_ROOT}, using default compiler")
endif()

# set the project name and version
project(viterbi3 CXX)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} - ${CMAKE_CXX_COMPILER}")

# Pour clangd / "Go to Implementation" (Ctrl+F12)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF()

# Set the output directory for the release build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ../bin/release)

# Set the output directory for the debug build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ../bin/debug)

message(STATUS "Building: ${CMAKE_BUILD_TYPE}")

# system
set(CMAKE_SYSTEM_NAME Linux)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)

# c++ (flags selon le compilateur détecté : Intel icpx ou GCC)
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
	set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG -O0 -g -xCORE-AVX2 -mrdrnd")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xCORE-AVX2 -ipo -mrdrnd")
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-ipo")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG -O0 -g -mavx2 -mfma -mrdrnd")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2 -mfma -mrdrnd")
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE "")
endif()
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "")

# #source files
FILE(GLOB  Sources ${PROJECT_SOURCE_DIR}/src/*.c??)

add_custom_target(update_version
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# add the executable
add_executable(${PROJECT_NAME} ${Sources})


SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
# libs
set(IPP_ROOT "/opt/intel/oneapi/ipp/latest" CACHE PATH "Intel IPP installation directory")
target_include_directories(${PROJECT_NAME} PUBLIC src ${IPP_ROOT}/include/ipp)
target_link_directories(${PROJECT_NAME} PRIVATE ${IPP_ROOT}/lib/intel64)
target_link_libraries(${PROJECT_NAME} uhd pthread boost_system boost_program_options ipps ippcore)

# Lien symbolique pour que clangd trouve compile_commands.json à la racine (Ctrl+F12)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/compile_commands.json ${CMAKE_SOURCE_DIR}/compile_commands.json
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
