set(HMGRID3D_EXNAME hmgrid3d_test)

set (HEADERS
)
set (SOURCES
    hmgrid3d_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${HMGRID3D_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMGRID3D_EXNAME} ${HMGRID3D_TARGET})
target_link_libraries(${HMGRID3D_EXNAME} ${HMMAPPING_TARGET})

include_directories(${CommonInclude})
include_directories(${HMGRID3D_INCLUDE})
include_directories(${HMMAPPING_INCLUDE})

#C++11
set_property(TARGET ${HMGRID3D_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMGRID3D_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)

