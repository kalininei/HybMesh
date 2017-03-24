set(CONTOURS2D_EXNAME hybmesh_contours2d_test)

set (HEADERS
)
set (SOURCES
    hybmesh_contours2d_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${CONTOURS2D_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${CONTOURS2D_EXNAME} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${CommonInclude})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
