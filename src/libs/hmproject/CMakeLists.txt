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

add_library(${HMPROJECT_TARGET} SHARED ${HEADERS} ${SOURCES})

target_link_libraries(${HMPROJECT_TARGET} ${LibXML2_LIBRARIES})
include_directories(${LibXML2_INCLUDE_DIRS})

install(TARGETS ${HMPROJECT_TARGET}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	)

#C++11
set_property(TARGET ${HMPROJECT_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMPROJECT_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

