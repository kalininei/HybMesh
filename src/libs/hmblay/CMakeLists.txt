#library file
set (HEADERS
	c_hmblay.h
	hmblay.hpp
	options.hpp
	bgrid.hpp
	connectors.hpp
	extpath.hpp
	canonic_bgrid.hpp
	bgrid_impose.hpp
)

set (SOURCES
	hmblay.cpp
	options.cpp
	bgrid.cpp
	connectors.cpp
	extpath.cpp
	canonic_bgrid.cpp
	bgrid_impose.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMBLAY_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMBLAY_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMBLAY_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMBLAY_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMBLAY_TARGET} ${HMMATH_TARGET})

include_directories(${CommonInclude})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${HMBLAY_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${HMBLAY_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMBLAY_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
