# ==================== find libgmsh or build it from internal source
if (BUILD_GMSH)
	add_subdirectory(gmsh)
endif()
add_subdirectory(clipper)
add_subdirectory(scpack)
add_subdirectory(dscpack)
add_subdirectory(gpc)
