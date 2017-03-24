#library file
set (HEADERS
	primitives2d.hpp
	contour.hpp
	contour_tree.hpp
	contabs2d.hpp
	finder2d.hpp
	treverter2d.hpp
	buildcont.hpp
	clipper_core.hpp
	clipdomain.hpp
	debug2d.hpp
	gpc_core.hpp
	assemble2d.hpp
	modcont.hpp
	partcont.hpp
	coarsencont.hpp
	export2d_hm.hpp
	export2d_vtk.hpp
	export2d_tecplot.hpp
	export2d_gmsh.hpp
	export2d_fluent.hpp
	import2d_hm.hpp
)

set (SOURCES
	primitives2d.cpp
	contour.cpp
	contour_tree.cpp
	contabs2d.cpp
	finder2d.cpp
	treverter2d.cpp
	buildcont.cpp
	clipper_core.cpp
	clipdomain.cpp
	debug2d.cpp
	gpc_core.cpp
	assemble2d.cpp
	modcont.cpp
	partcont.cpp
	coarsencont.cpp
	export2d_hm.cpp
	export2d_vtk.cpp
	export2d_tecplot.cpp
	export2d_gmsh.cpp
	export2d_fluent.cpp
	import2d_hm.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HYBMESH_CONTOURS2D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${GPC_TARGET})
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${CLIPPER_TARGET})
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${HMMATH_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${CLIPPER_INCLUDE})
include_directories(${GPC_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${HYBMESH_CONTOURS2D_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

add_subdirectory(tests)
