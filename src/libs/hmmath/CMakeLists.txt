#library file
set (HEADERS
	hmconformal.hpp
	hmcompute.hpp
	scpack_port.hpp
	hmfem.hpp
	femgrid43.hpp
	femmat.hpp
	femassembly.hpp
)

set (SOURCES
	hmconformal.cpp
	hmcompute.cpp
	scpack_port.cpp
	scpack_port.f
	hmfem.cpp
	femgrid43.cpp
	femmat.cpp
	femassembly.cpp
	$<TARGET_OBJECTS:${SCPACK_TARGET}>
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMMATH_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMMATH_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMMATH_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMMATH_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${CommonInclude})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})

install(TARGETS ${HMMATH_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${HMMATH_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMMATH_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
