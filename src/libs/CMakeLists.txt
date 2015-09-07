set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH}${LIB_INSTALL_DIR})

#headers from include folder
set (CommonInclude ${CMAKE_CURRENT_LIST_DIR}/include)

add_subdirectory(external)
add_subdirectory(bgeom2d)
add_subdirectory(hybmesh_contours2d)
add_subdirectory(crossgrid)
add_subdirectory(hmblay)
add_subdirectory(hmcport)


#copy mingw runtime libraries to binary directory
if (WIN32)
	find_library(RUNTIME_LIB1 "libstdc++-6.dll")
	find_library(RUNTIME_LIB2 "libwinpthread-1.dll")
	if (CXX_BITNESS STREQUAL "64")
		find_library(RUNTIME_LIB3 "libgcc_s_seh-1.dll")
	else()
		find_library(RUNTIME_LIB3 "libgcc_s_dw2-1.dll")
	endif()
	install(FILES ${RUNTIME_LIB1} ${RUNTIME_LIB2} ${RUNTIME_LIB3}
		DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
endif()
