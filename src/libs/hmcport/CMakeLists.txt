#library file
set (HEADERS
	hmcport.h
)

set (SOURCES
	hmcport.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMCPORT_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMCPORT_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMBLAY_TARGET})

include_directories(${CommonInclude})
include_directories(${GMSH_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMBLAY_INCLUDE})

install(TARGETS ${HMCPORT_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${HMCPORT_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMCPORT_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)
