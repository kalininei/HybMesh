set(HMNUMERIC_EXNAME hmnumeric_test)

set (HEADERS
)
set (SOURCES
    hmnumeric_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HMNUMERIC_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMNUMERIC_EXNAME} ${HMNUMERIC_TARGET})

include_directories(${CommonInclude})
include_directories(${HMNUMERIC_INCLUDE})
