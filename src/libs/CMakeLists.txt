set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH}${LIB_INSTALL_DIR})

add_subdirectory(external)
add_subdirectory(hmproject)
add_subdirectory(bgeom2d)
add_subdirectory(hybmesh_contours2d)
add_subdirectory(crossgrid)
add_subdirectory(hmblay)
add_subdirectory(hmmath)
add_subdirectory(hmmapping)
add_subdirectory(hybmesh_surfaces3d)
add_subdirectory(hmgrid3d)
add_subdirectory(hmcport)
add_subdirectory(hmnumeric)


#copy mingw runtime libraries to binary directory
if (WIN32)
	find_library(RUNTIME_LIB1 "libstdc++-6.dll")
	find_library(RUNTIME_LIB2 "libwinpthread-1.dll")
	if (CXX_BITNESS STREQUAL "64")
		find_library(RUNTIME_LIB3 "libgcc_s_seh-1.dll")
	else()
		find_library(RUNTIME_LIB3 "libgcc_s_dw2-1.dll")
	endif()
	find_library(RUNTIME_LIB4 "libquadmath-0.dll")
	find_library(RUNTIME_LIB5 "libgfortran-3.dll")
	find_library(RUNTIME_LIB6 "libgomp-1.dll")
	
	set(MINGW_RUNTIME_LIBS ${RUNTIME_LIB1} ${RUNTIME_LIB2} ${RUNTIME_LIB3}
		${RUNTIME_LIB4} ${RUNTIME_LIB5} ${RUNTIME_LIB6} ${SUITESPARSE_LIBRARIES}
		${LIBXML2_LIBRARIES})
		
	# Copy runtime libs to assembly directory to be able to use hybmesh prior to installation.
	# This is needed to run bindparser.py. 
	add_custom_target(MINGW_RUNTIME_BIN_COPY ALL)
	foreach(_ARG ${MINGW_RUNTIME_LIBS})
		add_custom_command(TARGET MINGW_RUNTIME_BIN_COPY PRE_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${_ARG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
		)
	endforeach()
	
	# install
	install(FILES ${MINGW_RUNTIME_LIBS}
		DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
endif()
