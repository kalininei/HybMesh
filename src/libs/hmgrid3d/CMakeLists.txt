set (HEADERS
	hmgrid3d.hpp
	primitives_grid3d.hpp
	construct_grid3d.hpp
	revolve_grid3d.hpp
	debug_grid3d.hpp
	surface_grid3d.hpp
	serialize_grid3d.hpp
	tetrahedral.hpp
	unite_grid3d.hpp
	fluent_export_grid3d.hpp
	vtk_export_grid3d.hpp
	tecplot_export_grid3d.hpp
)

set (SOURCES
	hmgrid3d.cpp
	primitives_grid3d.cpp
	construct_grid3d.cpp
	revolve_grid3d.cpp
	debug_grid3d.cpp
	surface_grid3d.cpp
	serialize_grid3d.cpp
	tetrahedral.cpp
	unite_grid3d.cpp
	fluent_export_grid3d.cpp
	vtk_export_grid3d.cpp
	tecplot_export_grid3d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMGRID3D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMGRID3D_TARGET} ${GMSH_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${CROSSGRID_TARGET})

include_directories(${HMGRID3D_TARGET} ${GMSH_INCLUDE})
include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})

install(TARGETS ${HMGRID3D_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

#C++11
set_property(TARGET ${HMGRID3D_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMGRID3D_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
