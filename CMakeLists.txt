# ======================== Project structure:
# [ROOT]
#  | - [src]
#  |    | - CMakeLists.txt
#  |    | - $MainGui.py    //main python executable. Should return version on --version argument.
#  |    | - $MainConf.py   //configuration module
#  |    | - [$PyPack]   //python package
#  |    |     | - ...   //python source files
#  |    | - [libs]      //c libraries
#  |          | - CMakeLists.txt
#  |          | - [...] //directories with cpp libraries. Each contains CMakeLists.txt.   
#  | - [BuildScripts]
#  |     | - RunPyinstaller.cmake.in   //script for building windows binary from ${MainGui}
#  |     | - make_installer.nsi.in     //script for building windows installer
#  |     | - setup.py.in               //python setup script
#  |     | - setup.cfg.in              //python setup script options
#  | - [misc]
#  |     | - $MainGui.ico  //icon file         
#  | - CMakeLists.txt (This. Root cmake file)

# ========================= Installation strategy
# Unix:
# -- $MainGui -> ${CMAKE_INSTALL_PREFIX}/bin
# -- $MainConf and $PyPack will be installed as standard python packages
# -- cpp libs -> ${CMAKE_INSTALL_PREFIX}/lib/${APPNAME}/*.so
# Windows
# -- $MainGui -> ${CMAKE_INSTALL_PREFIX}/py/$MainGui.py
# -- cpp libs -> ${CMAKE_INSTALL_PREFIX}/lib/*.dll, *.dll.a
# -- PyInstaller output binaries, MinGW runtime -> ${CMAKE_INSTALL_PREFIX}/bin
# -- NSIS installer -> ${CMAKE_INSTALL_PREFIX} (Run "make pack" after "make install")

cmake_minimum_required(VERSION 3.1)

project(HybMesh)

# ========================== Project defines
# name of resulting installed executable
set(APPNAME HybMesh)
# names of files and directories in the file tree
set(MainGui HybMesh)
set(MainConf HybMesh_conf)
set(PyPack HybMeshGui)
# description
set(APPDESCRIPTION "composite grid generator")

# ========================== Python check
find_program(PYTHON "python")
if (NOT PYTHON)
	message(FATAL_ERROR "Python was not found")
endif()

# ============================ System checks
set(CXX_BITNESS "32")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(CXX_BITNESS "64")
endif()
execute_process(COMMAND ${PYTHON}
	-c "import sys; print 64 if sys.maxsize > 2**32 else 32"
	RESULT_VARIABLE _PRESULT
	OUTPUT_VARIABLE PYTHON_BITNESS)
if (NOT _PRESULT EQUAL 0)
	message(FATAL_ERROR "Cannot detect python bitness")
endif()
if (NOT CXX_BITNESS EQUAL PYTHON_BITNESS)
	message(FATAL_ERROR "It seems python bitness ${PYTHON_BITNESS} differs from compiler bitness ${CXX_BITNESS}")
endif()

# ========================== Windows installer options
if (WIN32)
	find_program(PYINSTALLER "pyinstaller" PATHS ${PYINST_HINT_PATH})
	find_program(NSIS "makensis" PATHS ${NSIS_HINT_PATH})
	if (NOT PYINSTALLER)
		message(STATUS "PyInstaller was not found")
	endif()
	if (NOT NSIS)
		message(STATUS "NSIS was not found")
	endif()
	if (PYINSTALLER)
		option(USE_PYINSTALLER "Build windows binaries on make install" ON)
		if (NSIS)
			option(USE_NSIS "Build windows installer on make pack" ON)
		endif()
	endif()
endif(WIN32)

# ============================ Version
# provides ${APPVERSION} string
execute_process(COMMAND
	${PYTHON} ${CMAKE_SOURCE_DIR}/src/setup.py --version
	RESULT_VARIABLE _PRESULT
	OUTPUT_VARIABLE APPVERSION)
if (NOT _PRESULT EQUAL 0)
	message(FATAL_ERROR "setup.py doesn't process --version argument")
endif()
string(STRIP ${APPVERSION} APPVERSION)
string(LENGTH "${APPVERSION}" APPVERSION_LENGTH)
if(APPVERSION_LENGTH EQUAL 0 OR APPVERSION_LENGTH GREATER 8)
	message(FATAL_ERROR "${MainGui}.py --version returns invalid version string")
endif()

# =========================== linking libraries flags
#libGmsh: global definition of BUILD_GMSH, GMSH_TARGET, GMSH_INCLUDE
find_library(FIND_LIBGMSH Gmsh)
if (FIND_LIBGMSH)
	option(USE_INTERNAL_LIBGMSH "use internal build of libgmsh" OFF)
	message(STATUS "LibGmsh found at ${FIND_LIBGMSH}")
else()
	message(STATUS "LibGmsh not found")
endif()
if (NOT FIND_LIBGMSH OR USE_INTERNAL_LIBGMSH)
	set(GMSH_TARGET gmsh_shared)
	set(BUILD_GMSH ON)
	set(_GPATH "src/libs/external/gmsh")
	set(GMSH_INCLUDE
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Common"
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Geo"
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Graphics"
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Numeric"
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Plugin"
		"${CMAKE_SOURCE_DIR}/${_GPATH}/Solver"
		"${CMAKE_BINARY_DIR}/${_GPATH}/Common"
		)
else()
	message(STATUS "Using system libGmsh")
	set(BUILD_GMSH OFF)
	set(GMSH_TARGET ${FIND_LIBGMSH})
	find_path(GMSH_INCLUDE NAMES Gmsh.h GModel.h PATH_SUFFIXES gmsh)
endif()
#libpolyclipping: global definition of BUILD_CLIPPER, CLIPPER_TARGET, CLIPPER_INCLUDE
find_library(FIND_LIBCLIPPER polyclipping)
if (FIND_LIBCLIPPER)
	option(USE_INTERNAL_LIBCLIPPER "use internal build of libpolyclipping" OFF)
	message(STATUS "clipper found at ${FIND_LIBGMSH}")
else()
	message(STATUS "clipper not found")
endif()
if (NOT FIND_LIBCLIPPER OR USE_INTERNAL_LIBCLIPPER)
	set(CLIPPER_TARGET polyclipping)
	set(BUILD_CLIPPER on)
	set(CLIPPER_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/external/clipper")
else()
	message(STATUS "Using system libpolyclipping")
	set(BUILD_CLIPPER off)
	set(CLIPPER_TARGET ${FIND_LIBCLIPPER})
	find_path(CLIPPER_INCLUDE NAMES clipper.hpp)
endif()
#libbgeom2d: global definition of BGEOM2D_INCLUDE, BGEOM2D_TARGET
set(BGEOM2D_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/bgeom2d")
set(BGEOM2D_TARGET bgeom2d)
#libhybmesh_contours2d: global definition of HYBMESH_CONTOURS2D_INCLUDE, HYBMESH_CONTOURS2D_TARGET
set(HYBMESH_CONTOURS2D_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hybmesh_contours2d")
set(HYBMESH_CONTOURS2D_TARGET hybmesh_contours2d)
#libcrossgird: global definition of CROSSGRID_INCLUDE, CROSSGRID_TARGET
set(CROSSGRID_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/crossgrid")
set(CROSSGRID_TARGET crossgrid)
#libhmblay: global definition of HMBLAY_INCLUDE, HMBLAY_TARGET
set(HMBLAY_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmblay")
set(HMBLAY_TARGET hmblay)
#libhmblay: global definition of HMFEM_INCLUDE, HMFEM_TARGET
set(HMFEM_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmfem")
set(HMFEM_TARGET hmfem)
#libhmcport: global definition of HMCPORT_INCLUDE, HMCPORT_TARGET
set(HMCPORT_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmcport")
set(HMCPORT_TARGET hmcport)


# =========================== configuration file 
if (WIN32)
	set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")
	set(HYBMESH_HOMEDIR "$ENV{PROGRAMDATA}\\HybMesh")
else()
	set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/${APPNAME}")
	set(HYBMESH_HOMEDIR "~/.HybMesh")
endif()
configure_file(${CMAKE_SOURCE_DIR}/BuildScripts/config_hybmesh.in
	config_hybmesh)
if (UNIX)
	install(FILES "${CMAKE_BINARY_DIR}/config_hybmesh"
		DESTINATION /etc/HybMesh)
else()
	install(FILES "${CMAKE_BINARY_DIR}/config_hybmesh"
		DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
	install(FILES "${CMAKE_BINARY_DIR}/config_hybmesh"
		DESTINATION "${CMAKE_INSTALL_PREFIX}/py")
endif()


# ============================ build source
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)
