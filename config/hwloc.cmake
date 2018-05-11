
# You can force a path to hwloc.h using the environment variables HWLOC, or
# HWLOC_INCDIR and HWLOC_LIBDIR (note that curl/curl.h is also searched)

if (DEFINED ENV{HWLOC})
    message(STATUS "Adding $ENV{HWLOC} to the search path for hwloc")
    set(HWLOC_INCDIR_HINTS ${HWLOC_INCDIR_HINTS} "$ENV{HWLOC}/include")
    set(HWLOC_INCDIR_HINTS ${HWLOC_INCDIR_HINTS} "$ENV{HWLOC}"        )
    set(HWLOC_LIBDIR_HINTS ${HWLOC_LIBDIR_HINTS} "$ENV{HWLOC}/lib"    )
    set(HWLOC_LIBDIR_HINTS ${HWLOC_LIBDIR_HINTS} "$ENV{HWLOC}/.libs"  )
endif()

if (DEFINED ENV{HWLOC_INCDIR})
    message(STATUS "Adding $ENV{HWLOC_INCDIR} to the search path for hwloc")
    # prepend !
    set(HWLOC_INCDIR_HINTS "$ENV{HWLOC_INCDIR}" ${HWLOC_INCDIR_HINTS})
endif()

if (DEFINED ENV{HWLOC_LIBDIR})
    message(STATUS "Adding $ENV{HWLOC_LIBDIR} to the search path for hwloc")
    # prepend !
    set(HWLOC_LIBDIR_HINTS "$ENV{HWLOC_LIBDIR}" ${HWLOC_LIBDIR_HINTS})
endif()

# Try in three passes, otherwise cmake gets in the way...
find_path   (HWLOC_INCDIR hwloc.h HINTS ${HWLOC_INCDIR_HINTS} DOC "hwloc headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
        )
if(NOT HWLOC_INCDIR)
find_path   (HWLOC_INCDIR hwloc.h HINTS ${HWLOC_INCDIR_HINTS} DOC "hwloc headers"
        NO_DEFAULT_PATH
        )
endif()
if(NOT HWLOC_INCDIR)
find_path   (HWLOC_INCDIR hwloc.h HINTS ${HWLOC_INCDIR_HINTS} DOC "hwloc headers")
endif()

find_library(HWLOC_LIB    hwloc   HINTS ${HWLOC_LIBDIR_HINTS} DOC "hwloc library"
    NO_DEFAULT_PATH
    )
if(NOT HWLOC_LIB)
find_library(HWLOC_LIB    hwloc   HINTS ${HWLOC_LIBDIR_HINTS} DOC "hwloc library")
endif()


# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "HWLOC_INCDIR=${HWLOC_INCDIR}")
message(STATUS "HWLOC_LIB=${HWLOC_LIB}")
string(COMPARE NOTEQUAL "${HWLOC_INCDIR}" HWLOC_INCDIR-NOTFOUND HWLOC_INCDIR_OK)
string(COMPARE NOTEQUAL "${HWLOC_LIB}" HWLOC_LIB-NOTFOUND HWLOC_LIBDIR_OK)

get_filename_component(HWLOC_LIBDIR ${HWLOC_LIB} PATH)

if(HWLOC_INCDIR_OK AND HWLOC_LIBDIR_OK)
include_directories(${HWLOC_INCDIR})
link_directories(${HWLOC_LIBDIR})
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-L${HWLOC_LIBDIR}")
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES ${HWLOC_INCDIR})
set(CMAKE_REQUIRED_LIBRARIES "-lhwloc")
CHECK_C_SOURCE_COMPILES("
    #include <hwloc.h>
    #ifdef HWLOC_API_VERSION
    #if HWLOC_API_VERSION < 0x00010400
    #error \"too old, never checked\"
    #endif
    #endif
    int main(void)
    {
      return hwloc_get_api_version();
    }
" HAVE_HWLOC)
endif()
