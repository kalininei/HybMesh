set (HEADERS
	primitives3d.hpp
	finder3d.hpp
	contabs3d.hpp
	debug3d.hpp
	surface.hpp
	surface_tree.hpp
	assemble3d.hpp
	serialize3d.hpp
	treverter3d.hpp
	
	export3d_fluent.hpp
	export3d_vtk.hpp
	export3d_gmsh.hpp
	export3d_tecplot.hpp
	export3d_hm.hpp
	import3d_hm.hpp
)

set (SOURCES
	primitives3d.cpp
	finder3d.cpp
	contabs3d.cpp
	debug3d.cpp
	surface.cpp
	surface_tree.cpp
	assemble3d.cpp
	serialize3d.cpp
	treverter3d.cpp

	export3d_fluent.cpp
	export3d_vtk.cpp
	export3d_gmsh.cpp
	export3d_tecplot.cpp
	export3d_hm.cpp
	import3d_hm.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HYBMESH_SURFACES3D_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
target_link_libraries(${HYBMESH_SURFACES3D_TARGET} ${HMPROJECT_TARGET})
target_link_libraries(${HYBMESH_SURFACES3D_TARGET} ${HMMATH_TARGET})
target_link_libraries(${HYBMESH_SURFACES3D_TARGET} ${BGEOM2D_TARGET})

include_directories(${HMPROJECT_INCLUDE})
include_directories(${HMMATH_INCLUDE})
include_directories(${BGEOM2D_INCLUDE})

install(TARGETS ${HYBMESH_SURFACES3D_TARGET}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)

add_subdirectory(tests)
