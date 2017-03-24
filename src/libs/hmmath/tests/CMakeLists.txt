set(HMMATH_EXNAME hmmath_test)

set (HEADERS
)
set (SOURCES
    hmmath_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

USE_CXX11()
add_executable (${HMMATH_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMMATH_EXNAME} ${HMMATH_TARGET})

include_directories(${CommonInclude})
include_directories(${HMMATH_INCLUDE})
