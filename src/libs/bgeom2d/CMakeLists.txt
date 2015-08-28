#library file
set (HEADERS
	bgeom2d.h
)

set (SOURCES
	bgeom2d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${BGEOM2D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
include_directories(${CommonInclude})

install(TARGETS ${BGEOM2D_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${BGEOM2D_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BGEOM2D_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)
