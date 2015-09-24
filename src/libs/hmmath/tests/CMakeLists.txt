set(HMMATH_EXNAME hmmath_test)

set (HEADERS
)
set (SOURCES
    hmmath_test.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (${HMMATH_EXNAME} ${HEADERS} ${SOURCES})

target_link_libraries(${HMMATH_EXNAME} ${HMMATH_TARGET})

include_directories(${CommonInclude})
include_directories(${HMMATH_INCLUDE})

#C++11
set_property(TARGET ${HMMATH_EXNAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMMATH_EXNAME} PROPERTY CXX_STANDARD_REQUIRED ON)

