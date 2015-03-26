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
target_link_libraries(${CROSSGRID_SONAME} ${GMSH_LIBRARY})

include_directories(${CommonInclude})

install(TARGETS ${CROSSGRID_SONAME}
	LIBRARY DESTINATION ".")

