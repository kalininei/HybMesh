set(HYBMESH_SURFACES3D_EXNAME hybmesh_surfaces3d_test)

set (HEADERS
)
set (SOURCES
    hybmesh_surfaces3d_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HYBMESH_SURFACES3D_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HYBMESH_SURFACES3D_EXNAME} ${HYBMESH_SURFACES3D_TARGET})
target_link_libraries(${HYBMESH_SURFACES3D_EXNAME} ${HMGRID3D_TARGET})

include_directories(${CommonInclude})
include_directories(${HYBMESH_SURFACES3D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${HMGRID3D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMGRID2D_INCLUDE})
