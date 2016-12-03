#library file
set (HEADERS
	bgeom2d.h
	bgeom3d.h
)

set (SOURCES
	bgeom2d.cpp
	bgeom3d.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${BGEOM2D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${BGEOM2D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${BGEOM2D_TARGET} ${HMMATH_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${BGEOM2D_TARGET}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	)

#C++11
set_property(TARGET ${BGEOM2D_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BGEOM2D_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)
