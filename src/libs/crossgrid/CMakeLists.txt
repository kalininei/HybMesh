#library file
set (HEADERS
	crossgrid.h
	grid.h
	trigrid.h
	buffergrid.h
	wireframegrid.h
	contours.h
	intrusion.h
	procgrid.h
	debug_grid2d.h
	vtk_export_grid2d.hpp
	fluent_export_grid2d.hpp
	tecplot_export_grid2d.hpp
)

set (SOURCES
	crossgrid.cpp
	grid.cpp
	trigrid.cpp
	buffergrid.cpp
	wireframegrid.cpp
	contours.cpp
	intrusion.cpp
	procgrid.cpp
	debug_grid2d.cpp
	vtk_export_grid2d.cpp
	fluent_export_grid2d.cpp
	tecplot_export_grid2d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${CROSSGRID_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${CROSSGRID_TARGET} ${GMSH_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${CROSSGRID_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${GMSH_INCLUDE})
include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})

install(TARGETS ${CROSSGRID_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

#C++11
set_property(TARGET ${CROSSGRID_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CROSSGRID_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
