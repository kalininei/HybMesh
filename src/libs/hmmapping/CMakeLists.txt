#library file
set (HEADERS
	hmconformal.hpp
	scpack_port.hpp
	dscpack_port.hpp
	confrect_fem.hpp
	hmmapping.hpp
	gridmap.hpp
	domapping.hpp
	mapped_contour.hpp
	rectangle_grid_builder.hpp
	circrect.hpp
)

set (SOURCES
	hmconformal.cpp
	scpack_port.cpp
	scpack_port.f
	dscpack_port.cpp
	dscpack_port.f
	confrect_fem.cpp
	$<TARGET_OBJECTS:${SCPACK_TARGET}>
	$<TARGET_OBJECTS:${DSCPACK_TARGET}>
	hmmapping.cpp
	gridmap.cpp
	domapping.cpp
	mapped_contour.cpp
	rectangle_grid_builder.cpp
	circrect.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HMMAPPING_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMMAPPING_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${HMNUMERIC_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${HMMATH_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMNUMERIC_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${HMMAPPING_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)
add_subdirectory(tests)
