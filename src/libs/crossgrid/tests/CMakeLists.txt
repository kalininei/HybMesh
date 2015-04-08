set(CROSSGRID_EXNAME crossgrid_test)

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

#C++11
set_property(TARGET ${CROSSGRID_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CROSSGRID_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)
