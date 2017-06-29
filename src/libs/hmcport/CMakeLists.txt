#library file
set (HEADERS
	hmcport.h
	cport_hmxml.h
	cport_grid2d.h
	cport_grid3d.h
	cport_cont2d.h
	cport_surface3d.h
	c2cpp_helper.hpp
	tscaler.hpp
)

set (SOURCES
	hmcport.cpp
	cport_hmxml.cpp
	cport_grid2d.cpp
	cport_grid3d.cpp
	cport_cont2d.cpp
	cport_surface3d.cpp
	tscaler.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HMCPORT_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMCPORT_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMMATH_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMMAPPING_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMBLAY_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HYBMESH_SURFACES3D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMGRID3D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMGRID2D_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HMMATH_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMMAPPING_INCLUDE})
include_directories(${HMBLAY_INCLUDE})
include_directories(${HYBMESH_SURFACES3D_INCLUDE})
include_directories(${HMGRID3D_INCLUDE})
include_directories(${HMGRID2D_INCLUDE})


install(TARGETS ${HMCPORT_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)
