#library file
set (HEADERS
	hmfem.hpp
	femgrid43.hpp
	femassembly.hpp
	hmfdm.hpp
)

set (SOURCES
	hmfem.cpp
	femgrid43.cpp
	femassembly.cpp
	hmfdm.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMNUMERIC_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMNUMERIC_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMNUMERIC_TARGET} ${HMMATH_TARGET})
target_link_libraries(${HMNUMERIC_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMNUMERIC_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMNUMERIC_TARGET} ${HYBMESH_CONTOURS2D_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${HMMATH_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})

install(TARGETS ${HMNUMERIC_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

#C++11
set_property(TARGET ${HMNUMERIC_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMNUMERIC_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
