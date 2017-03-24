set(CROSSGRID_EXNAME crossgrid_test)

set (HEADERS
)
set (SOURCES
    crossgrid_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${CROSSGRID_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${CROSSGRID_EXNAME} ${CROSSGRID_TARGET})
target_link_libraries(${CROSSGRID_EXNAME} ${HMCPORT_TARGET})
target_link_libraries(${CROSSGRID_EXNAME} ${HMPROJECT_TARGET})

include_directories(${CommonInclude})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMCPORT_INCLUDE})
