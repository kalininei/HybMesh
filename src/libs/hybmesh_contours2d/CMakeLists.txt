#library file
set (HEADERS
	c_hybmesh_contours2d.h
	hybmesh_contours2d.hpp
	collections.hpp
	containers.hpp
	edges.hpp
	tree.hpp
	contour.hpp
	constructor.hpp
	clipper_core.hpp
	contclipping.hpp
	debug_cont2d.hpp
)

set (SOURCES
	hybmesh_contours2d.cpp
	constructor.cpp
	collections.cpp
	contour.cpp
	clipper_core.cpp
	tree.cpp
	contclipping.cpp
	debug_cont2d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HYBMESH_CONTOURS2D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${CLIPPER_TARGET})
target_link_libraries(${HYBMESH_CONTOURS2D_TARGET} ${BGEOM2D_TARGET})

include_directories(${CommonInclude})
include_directories(${CLIPPER_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})

install(TARGETS ${HYBMESH_CONTOURS2D_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${HYBMESH_CONTOURS2D_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HYBMESH_CONTOURS2D_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
