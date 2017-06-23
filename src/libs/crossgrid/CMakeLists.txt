#library file
set (HEADERS
	unite_grids.hpp
	buffergrid.hpp
	buildgrid.hpp
	trigrid.hpp
	pebi.hpp
	healgrid.hpp
	modgrid.hpp
	infogrid.hpp
	wireframegrid.hpp
	snap_grid2cont.hpp
	inscribe_grid.hpp
	debug_grid2d.hpp
)

set (SOURCES
	unite_grids.cpp
	buffergrid.cpp
	buildgrid.cpp
	trigrid.cpp
	pebi.cpp
	healgrid.cpp
	modgrid.cpp
	infogrid.cpp
	wireframegrid.cpp
	snap_grid2cont.cpp
	inscribe_grid.cpp
	debug_grid2d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${CROSSGRID_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${CROSSGRID_TARGET} ${GMSH_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${HMMATH_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${GMSH_INCLUDE})
include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${CROSSGRID_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

add_subdirectory(tests)
