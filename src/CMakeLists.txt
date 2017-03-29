# ================ Python distribution
set(PYROOTDIR "${CMAKE_SOURCE_DIR}/src/py")
set(SCRIPTSDIR "${CMAKE_SOURCE_DIR}/BuildScripts")
set(PYOUTDIR "${CMAKE_BINARY_DIR}/python-build")
# =========================== configuration file 
if (WIN32)
	set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")
	set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")
	set(EXE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
	set(HYBMESH_HOMEDIR "$ENV{PROGRAMDATA}\\hybmesh")
else()
	set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/${APPNAME}")
	set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${APPNAME}")
	set(EXE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
	set(HYBMESH_HOMEDIR "~/.hybmesh")
endif()

configure_file("${SCRIPTSDIR}/config_installed.py.in" "${PYOUTDIR}/config_installed.py" @ONLY)
configure_file("${SCRIPTSDIR}/setup.py.in" "${PYOUTDIR}/setup.py" @ONLY)
configure_file("${SCRIPTSDIR}/setup.cfg.in" "${PYOUTDIR}/setup.cfg" @ONLY)
configure_file("${SCRIPTSDIR}/check_python_vers.py.in" "${PYOUTDIR}/check_python_vers.py" @ONLY)

# ============= for assembling installer
set(DEPS "${PYROOTDIR}/${PyPack}/__init__.py")
add_custom_command(
	OUTPUT ${PYOUTDIR}/timestamp
	COMMAND ${PYTHON} "${PYOUTDIR}/setup.py" build
	COMMAND ${CMAKE_COMMAND} -E touch ${PYOUTDIR}/timestamp
	DEPENDS ${DEPS}
)
add_custom_target(pydistr ALL DEPENDS ${PYOUTDIR}/timestamp)

#installs python packages to system
install(CODE "execute_process(COMMAND ${PYTHON} ${PYOUTDIR}/setup.py install)")

if (WIN32)
	#install main python script to py directory
	install(PROGRAMS "${PYROOTDIR}/${MainBin}.py"
		DESTINATION py)
	
	# =========  run pyinstaller
	if (${USE_PYINSTALLER})
		set(PYINST_IN "${CMAKE_SOURCE_DIR}/BuildScripts/RunPyinstaller.cmake.in")
		if (CMAKE_BUILD_TYPE STREQUAL "Debug")
			set(PYINST_OPT "-y")
		else()
			#set(PYINST_OPT "-y -w --clean")
			set(PYINST_OPT "-y --clean")
		endif()
		set(PYINST_OUTDIR "${CMAKE_BINARY_DIR}/pyinstaller-build")
		set(PYINST_TARGET "${CMAKE_INSTALL_PREFIX}/py/${MainBin}.py")
		#set(PYINST_ICON "${CMAKE_SOURCE_DIR}/misc/${MainGui}.ico")
		configure_file(${PYINST_IN} 
			"${CMAKE_BINARY_DIR}/RunPyinstaller.cmake" @ONLY)
		install(SCRIPT "${CMAKE_BINARY_DIR}/RunPyinstaller.cmake")
		install(DIRECTORY "${PYINST_OUTDIR}/${MainBin}/"
			DESTINATION bin
		)
	endif()
	
	# ========= assemble installer
	if (${USE_NSIS})
		# configuration
		set(NSIS_INSTALLER_NAME "${APPNAME}-v${APPVERSION}-win64") 
		
		string(REPLACE "/" "\\" _TMP "${EXE_INSTALL_DIR}")
		set(NSIS_INSTALLER_BIN "/r \"${_TMP}\\*\"")

		string(REPLACE "/" "\\" _TMP "${LIB_INSTALL_DIR}")
		set(NSIS_INSTALLER_LIB "\"${_TMP}\\*.dll\"")

		string(REPLACE "/" "\\" _TMP "${SHARE_INSTALL_DIR}")
		set(NSIS_INSTALLER_SHARE "/r \"${_TMP}\\*\"")

		string(REPLACE "/" "\\" _TMP "${PYINST_ICON}")
		set(NSIS_INSTALLER_ICON "\"${_TMP}\"")

		string(REPLACE "/" "\\" _TMP ${CMAKE_SOURCE_DIR}/doc/build/html)
		set(NSIS_INSTALLER_HTML "/r \"${_TMP}\\*\"")

		string(REPLACE "/" "\\" _TMP ${CMAKE_BINARY_DIR}/build/lib/hybmeshpack)
		set(NSIS_HYBMESHPACK_FILES "/r \"${_TMP}\\*.py\"")

		string(REPLACE "/" "\\" _TMP ${PYOUTDIR})
		set(NSIS_HYBMESHPACK_SETUP "/r \"${_TMP}\\*.py\"")

		configure_file(
			"${CMAKE_SOURCE_DIR}/BuildScripts/make_installer.nsi.in" 
			"${CMAKE_BINARY_DIR}/make_installer.nsi" @ONLY)

		# run at make pack
		add_custom_target(pack
			COMMAND ${NSIS} make_installer.nsi
			COMMAND ${CMAKE_COMMAND} -E copy 
				"${NSIS_INSTALLER_NAME}.exe"
				"${CMAKE_INSTALL_PREFIX}"
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		)
	endif()
else()
	# console application
	install(PROGRAMS "${CMAKE_SOURCE_DIR}/src/py/${MainBin}.py" 
		DESTINATION ${EXE_INSTALL_DIR}
		RENAME ${APPNAME})
endif()


# ================ build libraries
add_subdirectory(libs)
add_subdirectory(bindings)
