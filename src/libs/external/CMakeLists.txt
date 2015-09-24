# ==================== find libgmsh or build it from internal source
if (BUILD_GMSH)
	add_subdirectory(gmsh)
endif()
if (BUILD_CLIPPER)
	add_subdirectory(clipper)
endif()
add_subdirectory(scpack)
