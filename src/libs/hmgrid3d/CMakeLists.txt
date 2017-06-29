set (HEADERS
	hmgrid3d.hpp
	buildgrid3d.hpp
	revolve_grid3d.hpp
	tetrahedral.hpp
	tetramesh_preproc.hpp
	merge3d.hpp
	pyramid_layer.hpp
)

set (SOURCES
	hmgrid3d.cpp
	buildgrid3d.cpp
	revolve_grid3d.cpp
	tetrahedral.cpp
	tetramesh_preproc.cpp
	merge3d.cpp
	pyramid_layer.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HMGRID3D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMGRID3D_TARGET} ${GMSH_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HMMATH_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HYBMESH_SURFACES3D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HMGRID2D_TARGET})

include_directories(${GMSH_INCLUDE})
include_directories(${HMPROJECT_INCLUDE})
include_directories(${HMMATH_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${HYBMESH_SURFACES3D_INCLUDE})
include_directories(${HMGRID2D_INCLUDE})

install(TARGETS ${HMGRID3D_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)
add_subdirectory(tests)
