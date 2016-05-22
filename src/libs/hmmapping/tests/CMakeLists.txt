set(HMMAPPING_EXNAME hmmapping_test)

set (HEADERS
)
set (SOURCES
    hmmapping_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${HMMAPPING_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMMAPPING_EXNAME} ${HMMAPPING_TARGET})

include_directories(${CommonInclude})
include_directories(${HMMAPPING_INCLUDE})

#C++11
set_property(TARGET ${HMMAPPING_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMMAPPING_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)


