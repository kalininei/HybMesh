#library file
set (HEADERS
	hmfem.hpp
	femgrid43.hpp
	femmat.hpp
	femassembly.hpp
)

set (SOURCES
	hmfem.cpp
	femgrid43.cpp
	femmat.cpp
	femassembly.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMFEM_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMFEM_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMFEM_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMFEM_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${CommonInclude})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})

install(TARGETS ${HMFEM_TARGET}
	DESTINATION ${LIB_INSTALL_DIR})

#C++11
set_property(TARGET ${HMFEM_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMFEM_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)
