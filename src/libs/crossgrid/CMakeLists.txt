set(CROSSGRID_SONAME crossgrid)

#library file
set (HEADERS
	grid.h
	bgeom.h
	fileproc.h
	trigrid.h
	buffergrid.h
	wireframegrid.h
	contours.h
	intrusion.h
)

set (SOURCES
	crossgrid.cpp
	grid.cpp
	bgeom.cpp
	fileproc.cpp
	trigrid.cpp
	buffergrid.cpp
	wireframegrid.cpp
	contours.cpp
	intrusion.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${CROSSGRID_SONAME} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${CROSSGRID_SONAME} ${GMSH_TARGET})

include_directories(${CommonInclude})
include_directories(${GMSH_INCLUDE})

install(TARGETS ${CROSSGRID_SONAME}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${CROSSGRID_SONAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CROSSGRID_SONAME} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
