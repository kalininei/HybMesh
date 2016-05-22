#library file
set (HEADERS
	hmcport.h
	cport_grid2d.h
	cport_grid3d.h
	cport_cont2d.h
)

set (SOURCES
	hmcport.cpp
	cport_grid2d.cpp
	cport_grid3d.cpp
	cport_cont2d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMCPORT_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMCPORT_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMMAPPING_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMBLAY_TARGET})
target_link_libraries(${HMCPORT_TARGET} ${HMGRID3D_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMMAPPING_INCLUDE})
include_directories(${HMBLAY_INCLUDE})
include_directories(${HMGRID3D_INCLUDE})

install(TARGETS ${HMCPORT_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

#C++11
set_property(TARGET ${HMCPORT_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMCPORT_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)
