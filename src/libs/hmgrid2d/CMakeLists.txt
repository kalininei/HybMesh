#library file
set (HEADERS
	buildgrid.hpp
	trigrid.hpp
	pebi.hpp
	healgrid.hpp
	modgrid.hpp
	infogrid.hpp
	snap_grid2cont.hpp
)

set (SOURCES
	buildgrid.cpp
	trigrid.cpp
	pebi.cpp
	healgrid.cpp
	modgrid.cpp
	infogrid.cpp
	snap_grid2cont.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HMGRID2D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMGRID2D_TARGET} ${GMSH_TARGET})
target_link_libraries(${HMGRID2D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMGRID2D_TARGET} ${HMMATH_TARGET})
target_link_libraries(${HMGRID2D_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMGRID2D_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${GMSH_INCLUDE})
include_directories(${HMPROJECT_INCLUDE})
include_directories(${HMMATH_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})

install(TARGETS ${GRID2D_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)
