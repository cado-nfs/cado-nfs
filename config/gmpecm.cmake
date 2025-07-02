
# You can force a path to ecm.h using the environment variables GMPECM, or
# GMPECM_INCDIR and GMPECM_LIBDIR

if (DEFINED ENV{GMPECM})
    message(STATUS "Adding $ENV{GMPECM} to the search path for GMP-ECM")
    set(GMPECM_INCDIR_HINTS ${GMPECM_INCDIR_HINTS} "$ENV{GMPECM}/include")
    set(GMPECM_INCDIR_HINTS ${GMPECM_INCDIR_HINTS} "$ENV{GMPECM}"        )
    set(GMPECM_LIBDIR_HINTS ${GMPECM_LIBDIR_HINTS} "$ENV{GMPECM}/lib"    )
    set(GMPECM_LIBDIR_HINTS ${GMPECM_LIBDIR_HINTS} "$ENV{GMPECM}/.libs"  )
endif()

if (DEFINED ENV{GMPECM_INCDIR})
    message(STATUS "Adding $ENV{GMPECM_INCDIR} to the search path for GMP-ECM")
    # prepend !
    set(GMPECM_INCDIR_HINTS "$ENV{GMPECM_INCDIR}" ${GMPECM_INCDIR_HINTS})
endif()

if (DEFINED ENV{GMPECM_LIBDIR})
    message(STATUS "Adding $ENV{GMPECM_LIBDIR} to the search path for GMP-ECM")
    # prepend !
    set(GMPECM_LIBDIR_HINTS "$ENV{GMPECM_LIBDIR}" ${GMPECM_LIBDIR_HINTS})
endif()

# Try in three passes, otherwise cmake gets in the way...
find_path   (GMPECM_INCDIR ecm.h HINTS ${GMPECM_INCDIR_HINTS} DOC "GMP-ECM headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
        )
if(NOT GMPECM_INCDIR)
    find_path   (GMPECM_INCDIR ecm.h HINTS ${GMPECM_INCDIR_HINTS} DOC "GMP-ECM headers"
        NO_DEFAULT_PATH
        )
endif()
if(NOT GMPECM_INCDIR)
    find_path   (GMPECM_INCDIR ecm.h HINTS ${GMPECM_INCDIR_HINTS} DOC "GMP-ECM headers")
endif()

find_library(GMPECM_LIB    ecm   HINTS ${GMPECM_LIBDIR_HINTS} DOC "GMP-ECM library"
    NO_DEFAULT_PATH
    )
if(NOT GMPECM_LIB)
find_library(GMPECM_LIB    ecm   HINTS ${GMPECM_LIBDIR_HINTS} DOC "GMP-ECM library")
endif()


# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "GMPECM_INCDIR=${GMPECM_INCDIR}")
message(STATUS "GMPECM_LIB=${GMPECM_LIB}")
string(COMPARE NOTEQUAL "${GMPECM_INCDIR}" GMPECM_INCDIR-NOTFOUND GMPECM_INCDIR_OK)
string(COMPARE NOTEQUAL "${GMPECM_LIB}" GMPECM_LIB-NOTFOUND GMPECM_LIBDIR_OK)

get_filename_component(GMPECM_LIBDIR ${GMPECM_LIB} PATH)

if(GMPECM_INCDIR_OK AND GMPECM_LIBDIR_OK)
include_directories(${GMPECM_INCDIR})
link_directories(${GMPECM_LIBDIR})
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-L${GMPECM_LIBDIR}")
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES ${GMPECM_INCDIR})
set(CMAKE_REQUIRED_LIBRARIES ecm ${gmp_libname} m ${pthread_libs})
CHECK_C_SOURCE_COMPILES("
    #include <ecm.h>
    #include <stdlib.h>
    #include <stdio.h>
    int main(void)
    {
      printf(\"ECM version: %s\", ecm_version());
      return EXIT_SUCCESS;
    }
    " HAVE_GMPECM)
endif()
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)
