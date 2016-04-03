#library file
set (HEADERS
	hmmapping.hpp
	gridmap.hpp
	mapped_contour.hpp
)

set (SOURCES
	hmmapping.cpp
	gridmap.cpp
	mapped_contour.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMMAPPING_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HMMAPPING_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${BGEOM2D_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${HYBMESH_CONTOURS2D_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${CROSSGRID_TARGET})
target_link_libraries(${HMMAPPING_TARGET} ${HMMATH_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})
include_directories(${HYBMESH_CONTOURS2D_INCLUDE})
include_directories(${CROSSGRID_INCLUDE})
include_directories(${HMMATH_INCLUDE})

install(TARGETS ${HMMAPPING_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

#C++11
set_property(TARGET ${HMMAPPING_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMMAPPING_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(tests)
