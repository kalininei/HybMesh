# Gmsh - Copyright (C) 1997-2015 C. Geuzaine, J.-F. Remacle
#
# See the LICENSE.txt file for license information. Please report all
# bugs and problems to the public mailing list <gmsh@geuz.org>.

#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#project(gmsh CXX C)

set(ENABLE_BLAS_LAPACK ON)  #this should be on always
set(ENABLE_BUILD_LIB OFF)
set(ENABLE_BUILD_SHARED ON)

set(GMSH_MAJOR_VERSION 2)
set(GMSH_MINOR_VERSION 9)
set(GMSH_PATCH_VERSION 1)
set(GMSH_EXTRA_VERSION "" CACHE STRING "Gmsh extra version string")

set(GMSH_VERSION "${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}")
set(GMSH_VERSION "${GMSH_VERSION}.${GMSH_PATCH_VERSION}${GMSH_EXTRA_VERSION}")
set(GMSH_SHORT_LICENSE "GNU General Public License")

set(GMSH_API
  ${CMAKE_CURRENT_BINARY_DIR}/Common/GmshConfig.h 
  ${CMAKE_CURRENT_BINARY_DIR}/Common/GmshVersion.h
  Common/Gmsh.h Common/Context.h Common/GmshDefines.h Common/GmshMessage.h
    Common/VertexArray.h Common/Octree.h Common/OctreeInternals.h 
    Common/OS.h Common/StringUtils.h Common/OpenFile.h
    Common/onelab.h Common/GmshSocket.h Common/onelabUtils.h
  Numeric/Numeric.h Numeric/GaussIntegration.h Numeric/polynomialBasis.h
    Numeric/JacobianBasis.h Numeric/MetricBasis.h Numeric/bezierBasis.h Numeric/fullMatrix.h
    Numeric/FuncSpaceData.h
    Numeric/simpleFunction.h Numeric/cartesian.h Numeric/ElementType.h
    Numeric/BasisFactory.h
  Geo/GModel.h Geo/GEntity.h Geo/GPoint.h Geo/GVertex.h Geo/GEdge.h 
    Geo/GFace.h Geo/GRegion.h Geo/GEdgeLoop.h Geo/GEdgeCompound.h 
    Geo/GFaceCompound.h Geo/GRegionCompound.h Geo/GRbf.h Geo/MVertex.h
    Geo/MVertexBoundaryLayerData.h Geo/MEdge.h 
    Geo/MFace.h Geo/MElement.h Geo/MElementOctree.h Geo/MPoint.h Geo/MLine.h 
    Geo/MTriangle.h Geo/MQuadrangle.h Geo/MTetrahedron.h Geo/MHexahedron.h 
    Geo/MPrism.h Geo/MPyramid.h Geo/MElementCut.h Geo/MElementOctree.h 
    Geo/discreteVertex.h Geo/discreteEdge.h Geo/discreteFace.h Geo/discreteRegion.h 
    Geo/SPoint2.h Geo/SPoint3.h Geo/SVector3.h Geo/STensor3.h Geo/SBoundingBox3d.h 
    Geo/Pair.h Geo/Range.h Geo/SOrientedBoundingBox.h 
    Geo/CellComplex.h Geo/ChainComplex.h Geo/Cell.h Geo/Homology.h Geo/Chain.h 
    Geo/GenericVertex.h Geo/GenericEdge.h Geo/GenericFace.h Geo/GenericRegion.h
    Geo/partitionEdge.h Geo/CGNSOptions.h Geo/gmshLevelset.h Geo/boundaryLayersData.h
  Mesh/meshGEdge.h Mesh/meshGFace.h Mesh/meshGFaceOptimize.h 
    Mesh/meshGFaceElliptic.h Mesh/meshPartition.h Mesh/meshGFaceDelaunayInsertion.h 
    Mesh/simple3D.h Mesh/meshPartitionOptions.h Mesh/directions3D.h Mesh/yamakawa.h
    Mesh/Voronoi3D.h Mesh/Levy3D.h Mesh/periodical.h Mesh/meshMetric.h
  Numeric/mathEvaluator.h
  Solver/dofManager.h Solver/femTerm.h Solver/laplaceTerm.h Solver/elasticityTerm.h
    Solver/crossConfTerm.h Solver/orthogonalTerm.h
    Solver/linearSystem.h Solver/linearSystemGMM.h Solver/linearSystemCSR.h 
    Solver/linearSystemFull.h Solver/elasticitySolver.h Solver/sparsityPattern.h 
    Solver/groupOfElements.h Solver/linearSystemPETSc.h Solver/linearSystemMUMPS.h
    Solver/thermicSolver.h
  Post/PView.h Post/PViewData.h Plugin/PluginManager.h Post/OctreePost.h 
  Post/PViewDataList.h Post/PViewDataGModel.h Post/PViewOptions.h Post/ColorTable.h
   Numeric/nodalBasis.h
  Graphics/drawContext.h
  contrib/kbipack/gmp_normal_form.h contrib/kbipack/gmp_matrix.h 
    contrib/kbipack/gmp_blas.h contrib/kbipack/mpz.h
  contrib/DiscreteIntegration/Integration3D.h
  contrib/HighOrderMeshOptimizer/OptHOM.h contrib/HighOrderMeshOptimizer/OptHomMesh.h
  contrib/HighOrderMeshOptimizer/OptHomRun.h contrib/HighOrderMeshOptimizer/ParamCoord.h
  contrib/HighOrderMeshOptimizer/OptHomFastCurving.h contrib/HighOrderMeshOptimizer/SuperEl.h
  contrib/HighOrderMeshOptimizer/OptHomIntegralBoundaryDist.h
  contrib/HighOrderMeshOptimizer/CADDistances.h
  contrib/HighOrderMeshOptimizer/OptHomObjContribScaledJac.h
  contrib/HighOrderMeshOptimizer/OptHomObjContribMetricMin.h
  contrib/HighOrderMeshOptimizer/OptHomObjContribCADDist.h
  contrib/MeshOptimizer/MeshOptPatch.h contrib/MeshOptimizer/MeshOpt.h
  contrib/MeshOptimizer/MeshOptCommon.h contrib/MeshOptimizer/MeshOptimizer.h
  contrib/MeshOptimizer/MeshOptObjContribFunc.h contrib/MeshOptimizer/MeshOptObjContrib.h
  contrib/MeshOptimizer/MeshOptObjContribScaledNodeDispSq.h
  contrib/MeshOptimizer/MeshOptObjectiveFunction.h contrib/MeshOptimizer/MeshOptVertexCoord.h 
  contrib/MeshQualityOptimizer/MeshQualityObjContribIdealJac.h
  contrib/MeshQualityOptimizer/MeshQualityObjContribInvCond.h
  contrib/MeshQualityOptimizer/MeshQualityOptimizer.h
  contrib/MathEx/mathex.h)

get_property(IAMCHILD DIRECTORY  PROPERTY PARENT_DIRECTORY)
if(IAMCHILD)
  set(GMSH_API ${GMSH_API} PARENT_SCOPE)
endif(IAMCHILD)

execute_process(COMMAND date "+%Y%m%d" OUTPUT_VARIABLE DATE 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT DATE)
  set(DATE "unknown")
endif(NOT DATE)
set(GMSH_DATE "${DATE}")

if(NOT HOSTNAME)
  set(HOSTNAME "unknown")
endif(NOT HOSTNAME)
set(GMSH_HOST "${HOSTNAME}")

if(NOT PACKAGER)
  set(PACKAGER "unknown")
endif(NOT PACKAGER)
set(GMSH_PACKAGER "HybMeshBuilder")

set(GMSH_OS "${CMAKE_SYSTEM_NAME}")

include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckCXXCompilerFlag)

macro(append_gmsh_src DIRNAME FILES)
  foreach(FILE ${FILES})
    list(APPEND LIST ${DIRNAME}/${FILE})
  endforeach(FILE)
  set(GMSH_SRC ${GMSH_SRC};${LIST} PARENT_SCOPE)
  set(GMSH_DIRS ${GMSH_DIRS};${DIRNAME} PARENT_SCOPE)
endmacro(append_gmsh_src)

macro(find_all_libraries VARNAME LISTNAME PATH SUFFIX)
  set(${VARNAME})
  list(LENGTH ${LISTNAME} NUM_LIST)
  foreach(LIB ${${LISTNAME}})
    find_library(FOUND_LIB ${LIB} PATHS ${PATH} PATH_SUFFIXES ${SUFFIX})
    if(FOUND_LIB)
      list(APPEND ${VARNAME} ${FOUND_LIB})
    endif(FOUND_LIB)
    unset(FOUND_LIB CACHE)
  endforeach(LIB)
  list(LENGTH ${VARNAME} NUM_FOUND_LIBRARIES)
  if(NUM_FOUND_LIBRARIES LESS NUM_LIST)
    set(${VARNAME})
  endif(NUM_FOUND_LIBRARIES LESS NUM_LIST)
endmacro(find_all_libraries)

macro(set_config_option VARNAME STRING)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "Found " ${STRING})
endmacro(set_config_option)

# check if the machine is 64 bits (this is more reliable than using
# CMAKE_SIZEOF_VOID_P, which does not seem to work e.g. on some Suse
# machines)
check_type_size("void*" SIZEOF_VOID_P)
if(SIZEOF_VOID_P EQUAL 8)
  set_config_option(HAVE_64BIT_SIZE_T "64Bit")
endif(SIZEOF_VOID_P EQUAL 8)

if(ENABLE_BLAS_LAPACK)
  if(BLAS_LAPACK_LIBRARIES)
    # use libs as specified in the BLAS_LAPACK_LIBRARIES variable
    set_config_option(HAVE_BLAS "Blas(Custom)")
    set_config_option(HAVE_LAPACK "Lapack(Custom)")
    set(LAPACK_LIBRARIES ${BLAS_LAPACK_LIBRARIES})
  else(BLAS_LAPACK_LIBRARIES)
    if (NOT LAPACK_LIBRARIES)
      # try with generic names
      set(GENERIC_LIBS_REQUIRED lapack blas)
      find_all_libraries(LAPACK_LIBRARIES GENERIC_LIBS_REQUIRED ${GMSH_LAPACK_PATH} "")
      if(LAPACK_LIBRARIES)
        set_config_option(HAVE_BLAS "Blas(Generic)")
        set_config_option(HAVE_LAPACK "Lapack(Generic)")
      endif()
    endif()
    #if(NOT LAPACK_LIBRARIES)
      ## try OpenBLAS
      #set(OPENBLAS_LIBS_REQUIRED openblas gfortran-3 quadmath-0)
      #find_all_libraries(LAPACK_LIBRARIES OPENBLAS_LIBS_REQUIRED "" "")
      #if(LAPACK_LIBRARIES)
        #set_config_option(HAVE_BLAS "Blas(OpenBlas)")
        #set_config_option(HAVE_LAPACK "Lapack(OpenBlas)")
      #endif(LAPACK_LIBRARIES)  
    #endif()
    if(NOT LAPACK_LIBRARIES)
      message(FATAL_ERROR "Could not find Blas or Lapack")
    endif()
  endif(BLAS_LAPACK_LIBRARIES)
endif(ENABLE_BLAS_LAPACK)

set(LINK_LIBRARIES ${EXTERNAL_LIBRARIES} ${LAPACK_LIBRARIES})

# check existance of headers
check_include_file(dlfcn.h DLFCN_H)
if(DLFCN_H)
  set_config_option(HAVE_DLOPEN "Dlopen")
  list(APPEND EXTERNAL_LIBRARIES ${CMAKE_DL_LIBS})
endif(DLFCN_H)
check_include_file(linux/joystick.h LINUX_JOYSTICK_H)
if(LINUX_JOYSTICK_H)
  set_config_option(HAVE_LINUX_JOYSTICK "LinuxJoystick")
endif(LINUX_JOYSTICK_H)


add_subdirectory(Common)
add_subdirectory(Numeric)
add_subdirectory(Geo)
add_subdirectory(Mesh)
set_config_option(HAVE_MESH "Mesh")


# disable compile optimization on some known problematic files
check_cxx_compiler_flag("-O0" NOOPT)
if(NOOPT)
  file(GLOB_RECURSE NOOPT_SRC Numeric/robustPredicates.cpp Mesh/BDS.cpp
                              Parser/Gmsh.tab.cpp contrib/Tetgen*/predicates.cxx)
  foreach(FILE ${NOOPT_SRC})
    get_source_file_property(PROP ${FILE} COMPILE_FLAGS)
    if(PROP)
      set_source_files_properties(${FILE} PROPERTIES COMPILE_FLAGS "${PROP} -O0")
    else(PROP)
      set_source_files_properties(${FILE} PROPERTIES COMPILE_FLAGS "-O0")
    endif(PROP)
  endforeach(FILE)
endif(NOOPT)

list(SORT CONFIG_OPTIONS)
set(GMSH_CONFIG_OPTIONS "")
foreach(OPT ${CONFIG_OPTIONS})
  set(GMSH_CONFIG_OPTIONS "${GMSH_CONFIG_OPTIONS} ${OPT}")
endforeach(OPT)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Common/GmshConfig.h.in 
               ${CMAKE_CURRENT_BINARY_DIR}/Common/GmshConfig.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Common/GmshVersion.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/Common/GmshVersion.h)

# process cmake environment variables so we can append them to the -I include
# commands. This is not recommended (we should only use the cache variables) but
# it is very convenient: otherwise we have to remember providing the
# -D... options to cmake for each new build.
set(ENV_CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
set(ENV_CMAKE_INCLUDE_PATH $ENV{CMAKE_INCLUDE_PATH})
if(UNIX)
  if(ENV_CMAKE_PREFIX_PATH)
    string(REPLACE ":" ";" ENV_CMAKE_PREFIX_PATH ${ENV_CMAKE_PREFIX_PATH})
  endif(ENV_CMAKE_PREFIX_PATH)
  if(ENV_CMAKE_INCLUDE_PATH)
    string(REPLACE ":" ";" ENV_CMAKE_INCLUDE_PATH ${ENV_CMAKE_INCLUDE_PATH})
  endif(ENV_CMAKE_INCLUDE_PATH)
endif(UNIX)
list(APPEND EXTERNAL_INCLUDES ${CMAKE_INCLUDE_PATH} ${ENV_CMAKE_INCLUDE_PATH})
list(APPEND EXTERNAL_INCLUDES ${CMAKE_PREFIX_PATH} ${ENV_CMAKE_PREFIX_PATH})
foreach(DIR ${CMAKE_PREFIX_PATH} ${ENV_CMAKE_PREFIX_PATH})
  list(APPEND EXTERNAL_INCLUDES ${DIR}/include)
endforeach(DIR)

if(EXTERNAL_INCLUDES)
  list(REMOVE_DUPLICATES EXTERNAL_INCLUDES)
endif(EXTERNAL_INCLUDES)

# we could specify include dirs more selectively, but this is simpler
include_directories(Common Fltk Geo Graphics Mesh Solver Numeric Parser
  Plugin Post Qt ${EXTERNAL_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR}/Common)

# required for mingw64
if (WIN32)
  LIST(APPEND LINK_LIBRARIES wsock32 psapi winmm)
endif()

# set this for external codes that might include this CMakeList file
set(GMSH_EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDES} CACHE 
    STRING "External include directories" FORCE)
set(GMSH_EXTERNAL_LIBRARIES ${LINK_LIBRARIES} CACHE 
    STRING "External libraries" FORCE)

# =================================== add tetgen
set(HAVE_TETGEN "Tetgen1.5")
include_directories(${TETGEN_INCLUDE})
LIST(APPEND LINK_LIBRARIES ${TETGEN_TARGET})

# =================================== static library target
if(ENABLE_BUILD_LIB)
  add_library(gmsh_static STATIC ${GMSH_SRC})
  set_target_properties(gmsh_static PROPERTIES OUTPUT_NAME GmshLocal)
  target_link_libraries(gmsh_static ${LINK_LIBRARIES})
endif(ENABLE_BUILD_LIB)

# =================================== shared library target
if(ENABLE_BUILD_SHARED)
  add_library(${GMSH_TARGET} SHARED ${GMSH_SRC})
  set_target_properties(${GMSH_TARGET} PROPERTIES OUTPUT_NAME GmshLocal)
  set_target_properties(${GMSH_TARGET} PROPERTIES 
    VERSION ${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}.${GMSH_PATCH_VERSION}
    SOVERSION ${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION})
  if(HAVE_LAPACK AND LAPACK_FLAGS)
    set_target_properties(${GMSH_TARGET} PROPERTIES LINK_FLAGS ${LAPACK_FLAGS})
  endif(HAVE_LAPACK AND LAPACK_FLAGS)
  target_link_libraries(${GMSH_TARGET} ${LINK_LIBRARIES})
  target_compile_definitions(${GMSH_TARGET} PRIVATE HAVE_TETGEN=${HAVE_TETGEN} TETLIBRARY)
  install(TARGETS ${GMSH_TARGET} 
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
  )
  if (WIN32)
    install(FILES ${LAPACK_LIBRARIES} DESTINATION ${LIB_INSTALL_DIR})
  endif()
endif()

message(STATUS "")
message(STATUS "Gmsh ${GMSH_VERSION} has been configured for ${GMSH_OS}")
message(STATUS "")
message(STATUS " * Build options:" ${GMSH_CONFIG_OPTIONS})
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS "")
