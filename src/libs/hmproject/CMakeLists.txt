#library file
set (HEADERS
	hmproject.h
	hmdebug.hpp
	hmtimer.hpp
	hmcallback.hpp
	hmtesting.hpp
	piecewise.hpp
)

set (SOURCES
	hmproject.cpp
	hmdebug.cpp
	hmtimer.cpp
	hmcallback.cpp
	hmtesting.cpp
)

source_group ("Header Files" FILES ${HEADERS} ${HEADERS})
source_group ("Source Files" FILES ${SOURCES} ${SOURCES})

add_library(${HMPROJECT_TARGET} SHARED ${HEADERS} ${SOURCES})

# =========== linkage
install(TARGETS ${HMPROJECT_TARGET}
	LIBRARY DESTINATION ${LIB_INSTALL_DIR}
	RUNTIME DESTINATION ${LIB_INSTALL_DIR}
	)

#C++11
set_property(TARGET ${HMPROJECT_TARGET} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${HMPROJECT_TARGET} PROPERTY CXX_STANDARD_REQUIRED ON)

