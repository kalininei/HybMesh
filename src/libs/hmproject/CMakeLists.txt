#library file
set (HEADERS
	hmproject.h
	hmdebug.hpp
	hmtimer.hpp
	hmcallback.hpp
	hmtesting.hpp
	hmxmlreader.hpp
)

set (SOURCES
	hmproject.cpp
	hmdebug.cpp
	hmtimer.cpp
	hmcallback.cpp
	hmtesting.cpp
	hmxmlreader.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

USE_CXX11()
add_library(${HMPROJECT_TARGET} SHARED ${HEADERS} ${SOURCES})

target_link_libraries(${HMPROJECT_TARGET} ${LIBXML2_LIBRARIES})
target_link_libraries(${HMPROJECT_TARGET} ${GMSH_TARGET})

include_directories(${LIBXML2_INCLUDE_DIR})
include_directories(${GMSH_INCLUDE})

install(TARGETS ${HMPROJECT_TARGET}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	)
