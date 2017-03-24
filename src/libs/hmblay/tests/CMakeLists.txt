set(HMBLAY_EXNAME hmblay_test)

set (HEADERS
)
set (SOURCES
    hmblay_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HMBLAY_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMBLAY_EXNAME} ${HMBLAY_TARGET})

include_directories(${CommonInclude})
include_directories(${HMBLAY_INCLUDE})
