set(HMMAPPING_EXNAME hmmapping_test)

set (HEADERS
)
set (SOURCES
    hmmapping_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HMMAPPING_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMMAPPING_EXNAME} ${HMMAPPING_TARGET})

include_directories(${CommonInclude})
include_directories(${HMMAPPING_INCLUDE})

