# This module defines the following variables:
#
# SUITESPARSE_FOUND: TRUE iff SuiteSparse and all dependencies have been found.
# SUITESPARSE_INCLUDE_DIRS: Include directories for all SuiteSparse components.
# SUITESPARSE_LIBRARIES: Libraries for all SuiteSparse component libraries and dependencies.
# SUITESPARSE_VERSION:
# SUITESPARSE_MAIN_VERSION:
# SUITESPARSE_SUB_VERSION:
# SUITESPARSE_SUBSUB_VERSION:
#
# The following variables control the behaviour of this module:
#
# SUITESPARSE_LOCAL_DISTR: if ON script will search suitesparse and metis only in directories defined below.
# Otherwise system paths will be examined.
#
# SUITESPARSE_LOCAL_INCDIR: List of directories in which to
# search for SuiteSparse includes if SUITESPARSE_LOCAL_DISTR is ON,
#
# SUITESPARSE_LOCAL_LIBDIR: List of directories in which to
# search for SuiteSparse (and metis) libraries if SUITESPARSE_LOCAL_DISTR is ON,
# 
#
# The following variables define the presence / includes & libraries for the
# SuiteSparse components searched for, the SUITESPARSE_XX variables are the
# union of the variables for all components.
#
# == Symmetric Approximate Minimum Degree (AMD)
# AMD_FOUND
# AMD_INCLUDE_DIR
# AMD_LIBRARY
#
# == Constrained Approximate Minimum Degree (CAMD)
# CAMD_FOUND
# CAMD_INCLUDE_DIR
# CAMD_LIBRARY
#
# == Column Approximate Minimum Degree (COLAMD)
# COLAMD_FOUND
# COLAMD_INCLUDE_DIR
# COLAMD_LIBRARY
#
# Constrained Column Approximate Minimum Degree (CCOLAMD)
# CCOLAMD_FOUND
# CCOLAMD_INCLUDE_DIR
# CCOLAMD_LIBRARY
#
# == Sparse Supernodal Cholesky Factorization and Update/Downdate (CHOLMOD)
# CHOLMOD_FOUND
# CHOLMOD_INCLUDE_DIR
# CHOLMOD_LIBRARY
#
# == UMFPACK
# UMFPACK_FOUND
# UMFPACK_INCLUDE_DIR
# UMFPACK_LIBRARY
#
# == Multifrontal Sparse QR (SuiteSparseQR)
# SUITESPARSEQR_FOUND
# SUITESPARSEQR_INCLUDE_DIR
# SUITESPARSEQR_LIBRARY
#
# == Common configuration for all but CSparse
# SUITESPARSE_CONFIG_FOUND
# SUITESPARSE_CONFIG_INCLUDE_DIR
# SUITESPARSE_CONFIG_LIBRARY
#
#
# == Serial Graph Partitioning and Fill-reducing Matrix Ordering (METIS)
# METIS_FOUND
# METIS_LIBRARY
#
# Called if we failed to find SuiteSparse or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.

MACRO(SUITESPARSE_REPORT_NOT_FOUND REASON_MSG)
	UNSET(SUITESPARSE_FOUND)
	UNSET(SUITESPARSE_INCLUDE_DIRS)
	UNSET(SUITESPARSE_LIBRARIES)
	UNSET(SUITESPARSE_VERSION)
	UNSET(SUITESPARSE_MAIN_VERSION)
	UNSET(SUITESPARSE_SUB_VERSION)
	UNSET(SUITESPARSE_SUBSUB_VERSION)

	# Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
	# use the camelcase library name, not uppercase.
	IF (SuiteSparse_FIND_QUIETLY)
		MESSAGE(STATUS "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
	ELSEIF (SuiteSparse_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
	ELSE()
		# Neither QUIETLY nor REQUIRED, use WARNING which emits a message
		# but continues configuration and allows generation.
		MESSAGE(WARNING "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
	ENDIF (SuiteSparse_FIND_QUIETLY)
ENDMACRO(SUITESPARSE_REPORT_NOT_FOUND)

MACRO (FIND_EXTERNAL_PACKAGE ALIAS WINLIBPATH)
	if(WIN32)
		SET(${ALIAS}_LIBRARIES ${WINLIBPATH} CACHE FILES "Custom ${ALIAS} libraries")
		if (EXISTS ${${ALIAS}_LIBRARIES})
			set(${ALIAS}_FOUND TRUE)
		else (EXISTS ${${ALIAS}_LIBRARIES})
			set(${ALIAS}_FOUND FALSE)
		endif(EXISTS ${${ALIAS}_LIBRARIES})
	else(WIN32)
		FIND_PACKAGE(${ALIAS} QUIET)
	endif(WIN32)
	
	IF (NOT ${ALIAS}_FOUND)
		SUITESPARSE_REPORT_NOT_FOUND("Did not find ${ALIAS} library (required for SuiteSparse).")
	ENDIF (NOT ${ALIAS}_FOUND)
ENDMACRO (FIND_EXTERNAL_PACKAGE ALIAS WINLIBPATH)

# BLAS
FIND_EXTERNAL_PACKAGE(BLAS ${SUITESPARSE_LOCAL_LIBDIR}/libblas.lib)
# LAPACK
FIND_EXTERNAL_PACKAGE(LAPACK ${SUITESPARSE_LOCAL_LIBDIR}/liblapack.lib)

# Specify search directories for include files and libraries (this is the union
# of the search directories for all OSs). Search user-specified hint
# directories first if supplied, and search user-installed locations first
# so that we prefer user installs to system installs where both exist.
SET(SUITESPARSE_CHECK_INCLUDE_DIRS)
SET(SUITESPARSE_CHECK_LIBRARY_DIRS)
if (SUITESPARSE_LOCAL_DISTR)
	SET(SUITESPARSE_CHECK_INCLUDE_DIRS ${SUITESPARSE_CHECK_INCLUDE_DIRS};${SUITESPARSE_LOCAL_INCDIR})
	SET(SUITESPARSE_CHECK_LIBRARY_DIRS ${SUITESPARSE_CHECK_LIBRARY_DIRS};${SUITESPARSE_LOCAL_LIBDIR})
	message(STATUS "Searching suitesparse headers   in ${SUITESPARSE_CHECK_INCLUDE_DIRS}")
	message(STATUS "Searching suitesparse libraries in ${SUITESPARSE_CHECK_LIBRARY_DIRS}")
	set(FINDOPT "NO_DEFAULT_PATH")
else()
	#add system paths here
	SET(SUITESPARSE_CHECK_INCLUDE_DIRS ${SUITESPARSE_CHECK_INCLUDE_DIRS};/usr/include/suitesparse)
	SET(SUITESPARSE_CHECK_LIBRARY_DIRS ${SUITESPARSE_CHECK_LIBRARY_DIRS};/usr/lib64/)
	message(STATUS "Searching suitesparse headers and libraries in system directories")
	set(FINDOPT "")
endif()

MACRO(SSMODULE_FIND_ROUTINE LIBPREF LIBNM LIBH)
	#unset to guarantee library path change without cache clear
	UNSET(${LIBPREF}_LIBRARY CACHE)
	UNSET(${LIBPREF}_INCLUDE_DIR CACHE)
	SET(${LIBPREF}_FOUND TRUE)

	FIND_LIBRARY(${LIBPREF}_LIBRARY NAMES ${LIBNM} PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS} ${FINDOPT})
	IF (EXISTS ${${LIBPREF}_LIBRARY})
		MESSAGE(STATUS "Found ${LIBPREF} library: ${${LIBPREF}_LIBRARY}")
	ELSE (EXISTS ${${LIBPREF}_LIBRARY})
		SUITESPARSE_REPORT_NOT_FOUND("Did not find ${LIBNM} library.")
		SET(${LIBPREF}_FOUND FALSE)
	ENDIF (EXISTS ${${LIBPREF}_LIBRARY})
	MARK_AS_ADVANCED(${LIBPREF}_LIBRARY)

	FIND_PATH(${LIBPREF}_INCLUDE_DIR NAMES ${LIBH} PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS} ${FINDOPT})
	IF (EXISTS ${${LIBPREF}_INCLUDE_DIR})
		MESSAGE(STATUS "Found ${LIBPREF} header in: ${${LIBPREF}_INCLUDE_DIR}")
	ELSE (EXISTS ${${LIBPREF}_INCLUDE_DIR})
		SUITESPARSE_REPORT_NOT_FOUND("Did not find ${LIBPREF} header.")
		SET(${LIBPREF}_FOUND FALSE)
	ENDIF (EXISTS ${${LIBPREF}_INCLUDE_DIR})
	MARK_AS_ADVANCED(${LIBPREF}_INCLUDE_DIR)
ENDMACRO()

#SuiteSparse libs
SSMODULE_FIND_ROUTINE(AMD amd amd.h)
SSMODULE_FIND_ROUTINE(CAMD camd camd.h)
SSMODULE_FIND_ROUTINE(COLAMD colamd colamd.h)
SSMODULE_FIND_ROUTINE(CCOLAMD ccolamd ccolamd.h)
SSMODULE_FIND_ROUTINE(CHOLMOD cholmod cholmod.h)
SSMODULE_FIND_ROUTINE(UMFPACK umfpack umfpack.h)
SSMODULE_FIND_ROUTINE(SUITESPARSEQR spqr SuiteSparseQR_C.h)
SSMODULE_FIND_ROUTINE(SUITESPARSE_CONFIG suitesparseconfig SuiteSparse_config.h)

#metis (optional)
UNSET(METIS_LIBRARY CACHE)
FIND_LIBRARY(METIS_LIBRARY NAMES metis PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS} ${FINDOPT})
IF (EXISTS ${METIS_LIBRARY})
	MESSAGE(STATUS "Found METIS library: ${METIS_LIBRARY}")
	SET(METIS_FOUND TRUE)
ELSE (EXISTS ${METIS_LIBRARY})
	MESSAGE(STATUS "Did not find METIS library (optional SuiteSparse dependency)")
	SET(METIS_FOUND FALSE)
ENDIF (EXISTS ${METIS_LIBRARY})
MARK_AS_ADVANCED(METIS_LIBRARY)

#librt (not in OSX)
IF (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
	FIND_LIBRARY(LIBRT_LIBRARY NAMES rt PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
	IF (LIBRT_LIBRARY)
		MESSAGE(STATUS "Adding librt: ${LIBRT_LIBRARY} to SuiteSparse_config libraries")
	ELSE (LIBRT_LIBRARY)
		MESSAGE(STATUS "Could not find librt. "
			"Assuming that SuiteSparse was compiled without timing.")
	ENDIF (LIBRT_LIBRARY)
	MARK_AS_ADVANCED(LIBRT_LIBRARY)
	LIST(APPEND SUITESPARSE_CONFIG_LIBRARY ${LIBRT_LIBRARY})
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)

# Extract the SuiteSparse version from the appropriate header
SET(SUITESPARSE_VERSION_FILE ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h)
FILE(READ ${SUITESPARSE_VERSION_FILE} SUITESPARSE_CONFIG_CONTENTS)
STRING(REGEX MATCH "#define SUITESPARSE_MAIN_VERSION [0-9]+"
	SUITESPARSE_MAIN_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
STRING(REGEX REPLACE "#define SUITESPARSE_MAIN_VERSION ([0-9]+)" "\\1"
	SUITESPARSE_MAIN_VERSION "${SUITESPARSE_MAIN_VERSION}")
STRING(REGEX MATCH "#define SUITESPARSE_SUB_VERSION [0-9]+"
	SUITESPARSE_SUB_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
STRING(REGEX REPLACE "#define SUITESPARSE_SUB_VERSION ([0-9]+)" "\\1"
	SUITESPARSE_SUB_VERSION "${SUITESPARSE_SUB_VERSION}")
STRING(REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION [0-9]+"
	SUITESPARSE_SUBSUB_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
STRING(REGEX REPLACE "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)" "\\1"
	SUITESPARSE_SUBSUB_VERSION "${SUITESPARSE_SUBSUB_VERSION}")
SET(SUITESPARSE_VERSION
	"${SUITESPARSE_MAIN_VERSION}.${SUITESPARSE_SUB_VERSION}.${SUITESPARSE_SUBSUB_VERSION}")

## Only mark SuiteSparse as found if all required dependencies have been found.
SET(SUITESPARSE_FOUND FALSE)
IF (AMD_FOUND AND CAMD_FOUND AND COLAMD_FOUND AND CCOLAMD_FOUND AND CHOLMOD_FOUND AND UMFPACK_FOUND
		AND SUITESPARSEQR_FOUND AND	SUITESPARSE_CONFIG_FOUND AND BLAS_FOUND AND LAPACK_FOUND)

	SET(SUITESPARSE_FOUND TRUE)

	LIST(APPEND SUITESPARSE_INCLUDE_DIRS
		${AMD_INCLUDE_DIR} ${CAMD_INCLUDE_DIR} ${UMFPACK_INCLUDE_DIR}
		${COLAMD_INCLUDE_DIR} ${CCOLAMD_INCLUDE_DIR}
		${CHOLMOD_INCLUDE_DIR} ${SUITESPARSEQR_INCLUDE_DIR} ${SUITESPARSE_CONFIG_INCLUDE_DIR})
	LIST(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)

	# Important: The ordering of these libraries is *NOT* arbitrary, as these
	# could potentially be static libraries their link ordering is important.
	LIST(APPEND SUITESPARSE_LIBRARIES
		${SUITESPARSEQR_LIBRARY} ${CHOLMOD_LIBRARY}
		${CCOLAMD_LIBRARY} ${CAMD_LIBRARY} ${UMFPACK_LIBRARY}
		${COLAMD_LIBRARY} ${AMD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY}
		${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
	IF(METIS_FOUND)
		LIST(APPEND SUITESPARSE_LIBRARIES ${METIS_LIBRARY})
	ENDIF(METIS_FOUND)

ELSE()
	SUITESPARSE_REPORT_NOT_FOUND("Failed to find some/all required components of SuiteSparse.")
ENDIF()

# Handle REQUIRED and QUIET arguments to FIND_PACKAGE
INCLUDE(FindPackageHandleStandardArgs)
# A change to CMake after release 2.8.10.2 means that
# FindPackageHandleStandardArgs() unsets <LibraryName>_FOUND without checking
# if it is one of the variables passed whose existence & validity is verified
# by FindPackageHandleStandardArgs() in conjunction with handling the REQUIRED
# and QUIET optional arguments, as such we use an intermediary variable.
SET(SUITESPARSE_FOUND_COPY ${SUITESPARSE_FOUND})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SuiteSparse REQUIRED_VARS SUITESPARSE_INCLUDE_DIRS 
	SUITESPARSE_LIBRARIES SUITESPARSE_FOUND_COPY VERSION_VAR SUITESPARSE_VERSION)
