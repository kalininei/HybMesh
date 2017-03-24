set(HMGRID3D_EXNAME hmgrid3d_test)

set (HEADERS
)
set (SOURCES
    hmgrid3d_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HMGRID3D_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMGRID3D_EXNAME} ${HMGRID3D_TARGET})

include_directories(${CommonInclude})
include_directories(${HMGRID3D_INCLUDE})
