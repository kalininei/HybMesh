# ================ Python distribution
set(PYROOTDIR "${CMAKE_SOURCE_DIR}/src")
set(PYOUTDIR "${CMAKE_BINARY_DIR}/python-build")
set(SETUP_PY_IN "${CMAKE_SOURCE_DIR}/BuildScripts/setup.py.in")
set(SETUP_PY "${PYOUTDIR}/setup.py")
set(SETUP_CFG_IN "${CMAKE_SOURCE_DIR}/BuildScripts/setup.cfg.in")
set(SETUP_CFG "${PYOUTDIR}/setup.cfg")
configure_file(${SETUP_PY_IN} ${SETUP_PY} @ONLY)
configure_file(${SETUP_CFG_IN} ${SETUP_CFG} @ONLY)

set(DEPS "${PYROOTDIR}/${PyPack}/__init__.py")
add_custom_command(
	OUTPUT ${PYOUTDIR}/timestamp
	COMMAND ${PYTHON} ${SETUP_PY} build
	COMMAND ${CMAKE_COMMAND} -E touch ${PYOUTDIR}/timestamp
	DEPENDS ${DEPS}
)
add_custom_target(pydistr ALL DEPENDS ${PYOUTDIR}/timestamp)

#installs python packages to system
install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")

if (WIN32)
	#install main python script to py directory
	install(PROGRAMS "${PYROOTDIR}/${MainGui}.py"
		DESTINATION py)
	
	# =========  run pyinstaller
	if (${USE_PYINSTALLER})
	    set(PYINST_IN "${CMAKE_SOURCE_DIR}/BuildScripts/RunPyinstaller.cmake.in")
		if (CMAKE_BUILD_TYPE STREQUAL "Debug")
			set(PYINST_OPT "-y")
		else()
			set(PYINST_OPT "-y -w --clean")
		endif()
		set(PYINST_OUTDIR "${CMAKE_BINARY_DIR}/pyinstaller-build")
		set(PYINST_TARGET "${CMAKE_INSTALL_PREFIX}/py/${MainGui}.py")
		set(PYINST_ICON "${CMAKE_SOURCE_DIR}/misc/${MainGui}.ico")
		configure_file(${PYINST_IN} 
			"${CMAKE_BINARY_DIR}/RunPyinstaller.cmake" @ONLY)

		install(SCRIPT "${CMAKE_BINARY_DIR}/RunPyinstaller.cmake")
		install(DIRECTORY "${PYINST_OUTDIR}/${MainGui}/"
			DESTINATION bin
		)
	endif()
	
	# ========= assemble installer
	if (${USE_NSIS})
		# configuration
		if(CXX_BITNESS EQUAL 32)
			set(NSIS_INSTALLER_BIT "32")
		else()
			set(NSIS_INSTALLER_BIT "64")
		endif()
		set(NSIS_INSTALLER_NAME "${APPNAME}-v${APPVERSION}-win${NSIS_INSTALLER_BIT}") 
		set(NSIS_INSTALLER_APPNAME "${APPNAME}")
		set(NSIS_INSTALLER_EXENAME "${MainGui}.exe")
		
		set(NSIS_INSTALLER_DESCRIPTION "${APPDESCRIPTION}")
		string(REPLACE "/" "\\" _TMP "${CMAKE_INSTALL_PREFIX}/bin")
		set(NSIS_INSTALLER_BIN "/r \"${_TMP}\\*\"")
		set(NSIS_INSTALLER_BINDEST "bin")
		string(REPLACE "/" "\\" _TMP "${CMAKE_INSTALL_PREFIX}/lib")
		set(NSIS_INSTALLER_LIB "\"${_TMP}\\*.dll\"")
		set(NSIS_INSTALLER_LIBDEST "lib")
		string(REPLACE "/" "\\" _TMP "${PYINST_ICON}")
		set(NSIS_INSTALLER_ICON "\"${_TMP}\"")
		
		configure_file("${CMAKE_SOURCE_DIR}/BuildScripts/make_installer.nsi.in" 
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
	install(PROGRAMS "${CMAKE_SOURCE_DIR}/src/${MainGui}.py" 
		DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
		RENAME ${APPNAME})
endif()


# ================ build libraries
add_subdirectory(libs)
