set(BINDINGSOUT "${CMAKE_BINARY_DIR}/bindings")
set(BINDPARSER ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py)

if (INCLUDE_CPP_BINDINGS)
	set(CMDLINE cpp odir ${BINDINGSOUT}/cpp
		vers ${APPVERSION}
		exepath ${EXE_INSTALL_DIR}
		hmdir ${CMAKE_SOURCE_DIR}/src/py)
	add_custom_command(
		OUTPUT ${BINDINGSOUT}/cpp/Hybmesh.hpp
		COMMAND ${PYTHON} ${BINDPARSER} ${CMDLINE}
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cppbase.hpp
		        ${CMAKE_CURRENT_SOURCE_DIR}/commongen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/cppgen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/funcmask
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	add_custom_target(cpp_bindings_header ALL DEPENDS ${BINDINGSOUT}/cpp/Hybmesh.hpp)
	add_dependencies(cpp_bindings_header ${HMCPORT_TARGET})
	install(FILES ${BINDINGSOUT}/cpp/Hybmesh.hpp
	        DESTINATION ${SHARE_INSTALL_DIR}/cpp/)
endif()

if (INCLUDE_CS_BINDINGS)
	add_library(${CS_BINDINGS_TARGET} SHARED
			client.c core_hmconnection_mcspy.c)
	if (WIN32)
		set_target_properties(${CS_BINDINGS_TARGET} PROPERTIES PREFIX "")
	endif()
	set(CMDLINE cs odir ${BINDINGSOUT}/cs
		vers ${APPVERSION}
		exepath ${EXE_INSTALL_DIR}
		libpath ${SHARE_INSTALL_DIR}/cs
		hmdir ${CMAKE_SOURCE_DIR}/src/py)
	add_custom_command(
		OUTPUT ${BINDINGSOUT}/cs/Hybmesh.cs
		COMMAND ${PYTHON} ${BINDPARSER} ${CMDLINE} 
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/csbase.cs
		        ${CMAKE_CURRENT_SOURCE_DIR}/commongen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/csgen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/funcmask
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	add_custom_target(cs_bindings_header ALL DEPENDS ${BINDINGSOUT}/cs/Hybmesh.cs)
	add_dependencies(cs_bindings_header ${HMCPORT_TARGET})
	install(FILES ${BINDINGSOUT}/cs/Hybmesh.cs
	        DESTINATION ${SHARE_INSTALL_DIR}/cs)
	install(TARGETS ${CS_BINDINGS_TARGET}
	        LIBRARY DESTINATION ${SHARE_INSTALL_DIR}/cs
	        RUNTIME DESTINATION ${SHARE_INSTALL_DIR}/cs)
endif()

if (INCLUDE_JAVA_BINDINGS)
	add_library(${JAVA_BINDINGS_TARGET} SHARED
		client.c core_hmconnection_java.c
	)
	target_include_directories(${JAVA_BINDINGS_TARGET}
		PRIVATE ${JNI_INCLUDE_DIRS}
	)
	if (WIN32)
		set_target_properties(${JAVA_BINDINGS_TARGET} PROPERTIES PREFIX "")
	endif()
	set(CMDLINE java odir ${BINDINGSOUT}/java
		vers ${APPVERSION}
		exepath ${EXE_INSTALL_DIR}
		hmdir ${CMAKE_SOURCE_DIR}/src/py)
	add_custom_command(
		OUTPUT ${BINDINGSOUT}/java/Hybmesh.java
		COMMAND ${PYTHON} ${BINDPARSER} ${CMDLINE} 
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/javabase.java
		        ${CMAKE_CURRENT_SOURCE_DIR}/commongen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/javagen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/funcmask
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	add_custom_target(java_bindings_header ALL DEPENDS ${BINDINGSOUT}/java/Hybmesh.java)
	add_dependencies(java_bindings_header ${HMCPORT_TARGET})
	install(FILES ${BINDINGSOUT}/java/Hybmesh.java
	        DESTINATION ${SHARE_INSTALL_DIR}/java/)
	install(TARGETS ${JAVA_BINDINGS_TARGET}
	        LIBRARY DESTINATION ${SHARE_INSTALL_DIR}/java
	        RUNTIME DESTINATION ${SHARE_INSTALL_DIR}/java)
endif()


if (INCLUDE_PYTHON_BINDINGS)
	add_library(${PY_BINDINGS_TARGET} SHARED
		client.c core_hmconnection_mcspy.c
	)
	if (WIN32)
		set_target_properties(${PY_BINDINGS_TARGET} PROPERTIES PREFIX "")
	endif()
	set(CMDLINE py odir ${BINDINGSOUT}/py
		vers ${APPVERSION}
		exepath ${EXE_INSTALL_DIR}
		libpath ${SHARE_INSTALL_DIR}/py
		hmdir ${CMAKE_SOURCE_DIR}/src/py)
	add_custom_command(
		OUTPUT ${BINDINGSOUT}/py/Hybmesh.py
		COMMAND ${PYTHON} ${BINDPARSER} ${CMDLINE} 
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pybase.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/commongen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/pygen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/funcmask
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	add_custom_target(py_bindings_header ALL DEPENDS ${BINDINGSOUT}/py/Hybmesh.py)
	add_dependencies(py_bindings_header ${HMCPORT_TARGET})
	install(FILES ${BINDINGSOUT}/py/Hybmesh.py
	        DESTINATION ${SHARE_INSTALL_DIR}/py)
	install(TARGETS ${PY_BINDINGS_TARGET}
	        LIBRARY DESTINATION ${SHARE_INSTALL_DIR}/py
	        RUNTIME DESTINATION ${SHARE_INSTALL_DIR}/py)
endif()


if (INCLUDE_MATLAB_BINDINGS OR INCLUDE_OCTAVE_BINDINGS)
	set(CMDLINE m odir ${BINDINGSOUT}/m
		vers ${APPVERSION}
		exepath ${EXE_INSTALL_DIR}
		libpath ${SHARE_INSTALL_DIR}/m
		hmdir ${CMAKE_SOURCE_DIR}/src/py)
	add_custom_command(
		OUTPUT ${BINDINGSOUT}/m/Hybmesh.m
		COMMAND ${PYTHON} ${BINDPARSER} ${CMDLINE} 
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mbase.m
		        ${CMAKE_CURRENT_SOURCE_DIR}/commongen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/mgen.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/bindparser.py
		        ${CMAKE_CURRENT_SOURCE_DIR}/funcmask
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	add_custom_target(oct_bindings_header ALL DEPENDS ${BINDINGSOUT}/m/Hybmesh.m)
	add_dependencies(oct_bindings_header ${HMCPORT_TARGET})
	install(DIRECTORY ${BINDINGSOUT}/m
	        DESTINATION ${SHARE_INSTALL_DIR})
endif()

if (INCLUDE_MATLAB_BINDINGS)
	add_library(${MATLAB_BINDINGS_TARGET} SHARED
		client.c core_hmconnection_mcspy.c
	)
	install(TARGETS ${MATLAB_BINDINGS_TARGET}
	        LIBRARY DESTINATION ${SHARE_INSTALL_DIR}/m
	        RUNTIME DESTINATION ${SHARE_INSTALL_DIR}/m)
	install(FILES libcore_hmconnection_matlab.h
	        DESTINATION ${SHARE_INSTALL_DIR}/m)
endif()

if (INCLUDE_OCTAVE_BINDINGS)
	if (WIN32)
		set(OCTLIBRARY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${OCT_BINDINGS_TARGET}.oct")
		set(_F1 "${CMAKE_CURRENT_SOURCE_DIR}/core_hmconnection_oct.cpp")
		set(_F2 "${CMAKE_CURRENT_SOURCE_DIR}/client.c")
		# call visual basic script instead of usual octave.exe
		set(OCTCOMMAND "${OCTAVE_BIN_PATHS}/../octave.vbs")
		# do some obscure staff with paths slashing to make it work in mingw somehow
		string(REPLACE "\\" "/" OCTCOMMAND ${OCTCOMMAND})
		string(REPLACE "/" "\\" _F1 ${_F1})
		string(REPLACE "/" "\\" _F2 ${_F2})
		# add vbs invoker
		set(OCTCOMMAND cscript ${OCTCOMMAND})
		add_custom_command(
			OUTPUT ${OCTLIBRARY}
			COMMAND ${OCTCOMMAND}
			ARGS --no-gui --eval "\"mkoctfile('-s','${_F1}','${_F2}')\""
			DEPENDS ${_F1} ${_F2}
			WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
		)
		add_custom_target(${OCT_BINDINGS_TARGET} ALL DEPENDS ${OCTLIBRARY})
		install(FILES ${OCTLIBRARY}
			DESTINATION ${SHARE_INSTALL_DIR}/m)
	else()
		USE_CXX11()
		octave_add_oct(${OCT_BINDINGS_TARGET}
			SOURCES core_hmconnection_oct.cpp client.c
			EXTENSION oct
		)
		target_include_directories(${OCT_BINDINGS_TARGET}
			PRIVATE ${OCTAVE_INCLUDE_DIRS}
		)
		install(TARGETS ${OCT_BINDINGS_TARGET}
			LIBRARY DESTINATION ${SHARE_INSTALL_DIR}/m
			RUNTIME DESTINATION ${SHARE_INSTALL_DIR}/m)
	endif()
endif()
