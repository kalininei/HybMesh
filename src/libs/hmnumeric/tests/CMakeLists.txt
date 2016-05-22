set(HMNUMERIC_EXNAME hmnumeric_test)

set (HEADERS
)
set (SOURCES
    hmnumeric_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${HMNUMERIC_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMNUMERIC_EXNAME} ${HMNUMERIC_TARGET})

include_directories(${CommonInclude})
include_directories(${HMNUMERIC_INCLUDE})

#C++11
set_property(TARGET ${HMNUMERIC_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMNUMERIC_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)

