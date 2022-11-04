# - Try to find arpack
# Once done, this will define
#
#  ARPACK_FOUND - system has htslib
#  ARPACK_INCLUDE_DIRS - the htslib include directories
#  ARPACK_LIBRARIES - link these to use htslib

set(ARPACK_SEARCH_DIRS
        ${ARPACK_SEARCH_DIRS}
	$ENV{CONDA_PREFIX}
        )

set(_arpack_ver_path "arpack-${ARPACK_FIND_VERSION}")
include(LibFindMacros)

# Dependencies
#libfind_package(ARPACK)

# Include dir
find_path(ARPACK_INCLUDE_DIR
        NAMES ${ARPACK_ADDITIONAL_HEADERS} arpack.h
        PATHS ${ARPACK_SEARCH_DIRS}
        PATH_SUFFIXES
        include include/arpack arpack/${_arpack_ver_path}/arpack
        )

# Finally the library itself
find_library(ARPACK_LIBRARY
        NAMES arpack libarpack.a
        PATHS ${ARPACK_SEARCH_DIRS}
        NO_DEFAULT_PATH
        PATH_SUFFIXES lib lib64 ${_arpack_ver_path}
        )

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(ARPACK_PROCESS_INCLUDES ARPACK_INCLUDE_DIR)
set(ARPACK_PROCESS_LIBS ARPACK_LIBRARY)
libfind_process(ARPACK)
message(STATUS "   ARPACK include dirs: ${ARPACK_INCLUDE_DIRS}")
message(STATUS "   ARPACK libraries: ${ARPACK_LIBRARIES}")
