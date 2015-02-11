set (HEADERS
)
set (SOURCES
    crossgrid_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${CROSSGRID_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${CROSSGRID_EXNAME} ${CROSSGRID_SONAME})

include_directories(${CommonInclude})