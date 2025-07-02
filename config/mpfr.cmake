# You can force a path to mpfr.h using the environment variables MPFR, or
# MPFR_INCDIR and MPFR_LIBDIR
string(COMPARE NOTEQUAL "$ENV{MPFR}" "" HAS_MPFR_OVERRIDE)
if (HAS_MPFR_OVERRIDE)
    message(STATUS "Adding $ENV{MPFR} to the search path for Gnu MPFR")
    set(MPFR_INCDIR_HINTS "$ENV{MPFR}/include" ${MPFR_INCDIR_HINTS})
    set(MPFR_INCDIR_HINTS "$ENV{MPFR}"         ${MPFR_INCDIR_HINTS})
    set(MPFR_LIBDIR_HINTS "$ENV{MPFR}/lib"     ${MPFR_LIBDIR_HINTS})
    set(MPFR_LIBDIR_HINTS "$ENV{MPFR}/.libs"   ${MPFR_LIBDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPFR_INCDIR}" "" HAS_MPFR_INCDIR_OVERRIDE)
if (HAS_MPFR_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPFR_INCDIR} to the search path for Gnu MPFR")
    set(MPFR_INCDIR_HINTS "$ENV{MPFR_INCDIR}" ${MPFR_INCDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPFR_LIBDIR}" "" HAS_MPFR_LIBDIR_OVERRIDE)
if (HAS_MPFR_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPFR_LIBDIR} to the search path for Gnu MPFR")
    set(MPFR_LIBDIR_HINTS "$ENV{MPFR_LIBDIR}"     ${MPFR_LIBDIR_HINTS})
endif()

# First try overrides, really. We want cmake to shut up.
if (NOT MPFR_INCDIR)
    find_path   (MPFR_INCDIR mpfr.h PATHS ${MPFR_INCDIR_HINTS} DOC "Gnu MPFR headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif()
if (NOT MPFR_INCDIR)
    find_path   (MPFR_INCDIR mpfr.h HINTS ${MPFR_INCDIR_HINTS} DOC "Gnu MPFR headers"
        NO_DEFAULT_PATH
    )
endif()
if (NOT MPFR_INCDIR)
    find_path   (MPFR_INCDIR mpfr.h HINTS ${MPFR_INCDIR_HINTS} DOC "Gnu MPFR headers")
endif()

find_library(MPFR_LIB    mpfr   HINTS ${MPFR_LIBDIR_HINTS} DOC "Gnu MPFR library" NO_DEFAULT_PATH)
if(NOT MPFR_LIBDIR)
    find_library(MPFR_LIB    mpfr   HINTS ${MPFR_LIBDIR_HINTS} DOC "Gnu MPFR library")
endif()

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(MPFR_LIBDIR ${MPFR_LIB} PATH)
message(STATUS "MPFR_INCDIR=${MPFR_INCDIR}")
message(STATUS "MPFR_LIBDIR=${MPFR_LIBDIR}")
if(MPFR_INCDIR)
include_directories(${MPFR_INCDIR})
else()
    message(FATAL_ERROR "mpfr.h cannot be found. Please install Gnu MPFR, and specify its install prefix in local.sh")
endif()
if(MPFR_LIBDIR)
link_directories(${MPFR_LIBDIR})
else()
    message(FATAL_ERROR "MPFR library cannot be found. Please install Gnu MPFR, and specify its install prefix in local.sh")
endif()

set(WITH_MPFR 1 CACHE INTERNAL "MPFR is being used")
set(HAVE_MPFR 1 CACHE INTERNAL "MPFR is being used")
