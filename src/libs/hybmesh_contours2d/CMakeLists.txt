#library file
set (HEADERS
	primitives2d.hpp
	contour.hpp
	tree.hpp
	contabs2d.hpp
	finder2d.hpp
	treverter2d.hpp
	constructor.hpp
	clipper_core.hpp
	contclipping.hpp
	debug2d.hpp
	gpc_core.hpp
	cont_assembler.hpp
	algos.hpp
	cont_partition.hpp
	cont_repart.hpp
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
	tree.cpp
	contabs2d.cpp
	finder2d.cpp
	treverter2d.cpp
	constructor.cpp
	clipper_core.cpp
	contclipping.cpp
	debug2d.cpp
	gpc_core.cpp
	cont_assembler.cpp
	algos.cpp
	cont_partition.cpp
	cont_repart.cpp
	export2d_hm.cpp
	export2d_vtk.cpp
	export2d_tecplot.cpp
	export2d_gmsh.cpp
	export2d_fluent.cpp
	import2d_hm.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

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

#C++11
set_property(TARGET ${HYBMESH_CONTOURS2D_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HYBMESH_CONTOURS2D_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
