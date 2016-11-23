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

project(HybMesh C CXX Fortran)

# ========================== Project defines
# name of resulting installed executable
set(APPNAME hybmesh)
# names of files and directories in the file tree
set(MainGui hybmesh-gui)
set(MainBin hybmesh)
set(PyPack hybmeshpack)
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
set(APPVERSION 0.4.5)

# =========================== linking libraries flags
#libGmsh: global definition of BUILD_GMSH, GMSH_TARGET, GMSH_INCLUDE
find_library(FIND_LIBGMSH NAMES Gmsh)
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
set(CLIPPER_TARGET polyclipping)
set(BUILD_CLIPPER on)
set(CLIPPER_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/external/clipper")

#libgpc: global definition of GPC_TARGET, GPC_INCLUDE
set(GPC_TARGET gpc)
set(GPC_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/external/gpc")

#suitesparse: defines
# SUITESPARSE_FOUND: TRUE iff SuiteSparse and all dependencies have been found.
# SUITESPARSE_INCLUDE_DIRS: Include directories for all SuiteSparse components.
# SUITESPARSE_LIBRARIES: Libraries for all SuiteSparse component libraries and dependencies.
# ------------- Detect SuiteSparse libraries:
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
if (WIN32)
	set(SUITESPARSE_WINDIR ${CMAKE_SOURCE_DIR}/src/libs/external/winlib/)
	set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_WINDIR}include)
	set(SSTMP "")
	list(APPEND SSTMP spqr.2.0.5 cholmod.3.0.9 ccolamd.2.9.4 camd.2.4.4 umfpack.5.7.4 
	                          colamd.2.9.4 amd.2.4.4 suitesparseconfig.4.5.1 blas
							  lapack metis)
	foreach(T ${SSTMP})
		list(APPEND SUITESPARSE_LIBRARIES "${SUITESPARSE_WINDIR}lib${T}.dll")
	endforeach()
	set(SUITESPARSE_FOUND TRUE)
else(WIN32)
	find_package(SuiteSparse REQUIRED)
endif()
#libxml2
if (WIN32)
	set(LibXML2_WINDIR ${CMAKE_SOURCE_DIR}/src/libs/external/winlib/)
	set(LibXML2_INCLUDE_DIRS ${LibXML2_WINDIR}include)
	set(LibXML2_LIBRARIES "${LibXML2_WINDIR}libxml2-2.dll;${LibXML2_WINDIR}libiconv-2.dll")
	set(LibXML2_FOUND TRUE)
else(WIN32)
	find_package(LibXML2 REQUIRED)
endif()

#libhmproject: HMPROJECT_INCLUDE, HMPROJECT_TARGET
set(HMPROJECT_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmproject")
set(HMPROJECT_TARGET hmproject)
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
#libhmblay: global definition of HMGRID3D_INCLUDE, HMGRID3D_TARGET
set(HMGRID3D_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmgrid3d")
set(HMGRID3D_TARGET hmgrid3d)
#libhmmath: global definition of HMMATH_INCLUDE, HMMATH_TARGET
set(HMMATH_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmmath")
set(HMMATH_TARGET hmmath)
#libhmnumeric: global definition of HMNUMERIC_INCLUDE, HMNUMERIC_TARGET
set(HMNUMERIC_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmnumeric")
set(HMNUMERIC_TARGET hmnumeric)
#libhmmapping: global definition of HMMAPPING_INCLUDE, HMMAPPING_TARGET
set(HMMAPPING_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmmapping")
set(HMMAPPING_TARGET hmmapping)
#libscpack: SCPACK_TARGET
set(SCPACK_TARGET scpack)
#libscpack: DSCPACK_TARGET
set(DSCPACK_TARGET dscpack)
#libhmcport: global definition of HMCPORT_INCLUDE, HMCPORT_TARGET
set(HMCPORT_INCLUDE "${CMAKE_SOURCE_DIR}/src/libs/hmcport")
set(HMCPORT_TARGET hmcport)


# ============================ build source
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)
