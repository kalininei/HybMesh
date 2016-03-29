set (HEADERS
	hmgrid3d.hpp
	primitives_grid3d.hpp
	construct_grid3d.hpp
	export_grid3d.hpp
	debug_grid3d.hpp
)

set (SOURCES
	hmgrid3d.cpp
	primitives_grid3d.cpp
	construct_grid3d.cpp
	export_grid3d.cpp
	debug_grid3d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMGRID3D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMGRID3D_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMGRID3D_TARGET} ${CROSSGRID_TARGET})

include_directories(${CommonInclude})
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
