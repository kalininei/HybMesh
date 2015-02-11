cmake_minimum_required(VERSION 2.8)

set (PROJECT HYGRID)

option(USE_CLANG "build application with clang" ON)

if(USE_CLANG)
	set (CMAKE_C_COMPILER   "/usr/bin/clang")
	set (CMAKE_CXX_COMPILER "/usr/bin/clang++")
	set (CMAKE_AR      "/usr/bin/llvm-ar")
	set (CMAKE_LINKER  "/usr/bin/llvm-ld")
	set (CMAKE_NM      "/usr/bin/llvm-nm")
	set (CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
	set (CMAKE_RANLIB  "/usr/bin/llvm-ranlib")	
endif(USE_CLANG)

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Installation prefix")   

project(${PROJECT})

if(USE_CLANG)
	set (CMAKE_C_FLAGS                "-Wall -std=c99")
	set (CMAKE_C_FLAGS_DEBUG          "-g")
	set (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
	set (CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
	set (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
	
	set (CMAKE_CXX_FLAGS                "-Wall -std=c++11")
	set (CMAKE_CXX_FLAGS_DEBUG          "-O0 -g -gdwarf-2")
	set (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
	set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
else(USE_CLANG)
	set(CMAKE_CXX_FLAGS "-Wall -std=c++11")
endif(USE_CLANG)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

find_package(VTK REQUIRED)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#headers from include folder
set (CommonInclude ${CMAKE_SOURCE_DIR}/include)

#included programs
#shared libraries
set(CROSSGRID_SONAME crossgrid)
#binaries
#2d grid builder
set(CROSSGRID_EXNAME  crossgrid_test) 

add_subdirectory(crossgrid)
add_subdirectory(crossgrid/tests)



