set(HMBLAY_EXNAME hmblay_test)

set (HEADERS
)
set (SOURCES
    hmblay_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${HMBLAY_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMBLAY_EXNAME} ${HMBLAY_TARGET})

include_directories(${CommonInclude})
include_directories(${HMBLAY_INCLUDE})

#C++11
set_property(TARGET ${HMBLAY_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMBLAY_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)

