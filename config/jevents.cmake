
# You can force a path to jevents.h using the environment variables JEVENTS, or
# JEVENTS_INCDIR and JEVENTS_LIBDIR

if (DEFINED ENV{JEVENTS})
    message(STATUS "Adding $ENV{JEVENTS} to the search path for jevents")
    set(JEVENTS_INCDIR_HINTS ${JEVENTS_INCDIR_HINTS} "$ENV{JEVENTS}/include")
    set(JEVENTS_INCDIR_HINTS ${JEVENTS_INCDIR_HINTS} "$ENV{JEVENTS}"        )
    set(JEVENTS_LIBDIR_HINTS ${JEVENTS_LIBDIR_HINTS} "$ENV{JEVENTS}/lib"    )
    set(JEVENTS_LIBDIR_HINTS ${JEVENTS_LIBDIR_HINTS} "$ENV{JEVENTS}/lib64"    )
    set(JEVENTS_LIBDIR_HINTS ${JEVENTS_LIBDIR_HINTS} "$ENV{JEVENTS}/.libs"  )
endif()

if (DEFINED ENV{JEVENTS_INCDIR})
    message(STATUS "Adding $ENV{JEVENTS_INCDIR} to the search path for jevents")
    # prepend !
    set(JEVENTS_INCDIR_HINTS "$ENV{JEVENTS_INCDIR}" ${JEVENTS_INCDIR_HINTS})
endif()

if (DEFINED ENV{JEVENTS_LIBDIR})
    message(STATUS "Adding $ENV{JEVENTS_LIBDIR} to the search path for jevents")
    # prepend !
    set(JEVENTS_LIBDIR_HINTS "$ENV{JEVENTS_LIBDIR}" ${JEVENTS_LIBDIR_HINTS})
endif()

# Try in three passes, otherwise cmake gets in the way...
find_path   (JEVENTS_INCDIR jevents.h HINTS ${JEVENTS_INCDIR_HINTS} DOC "jevents headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
        )
if(NOT JEVENTS_INCDIR)
find_path   (JEVENTS_INCDIR jevents.h HINTS ${JEVENTS_INCDIR_HINTS} DOC "jevents headers"
        NO_DEFAULT_PATH
        )
endif()
if(NOT JEVENTS_INCDIR)
find_path   (JEVENTS_INCDIR jevents.h HINTS ${JEVENTS_INCDIR_HINTS} DOC "jevents headers")
endif()

find_library(JEVENTS_LIB    jevents   HINTS ${JEVENTS_LIBDIR_HINTS} DOC "jevents library"
    NO_DEFAULT_PATH
    )
if(NOT JEVENTS_LIB)
find_library(JEVENTS_LIB    jevents   HINTS ${JEVENTS_LIBDIR_HINTS} DOC "jevents library")
endif()


# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "JEVENTS_INCDIR=${JEVENTS_INCDIR}")
message(STATUS "JEVENTS_LIB=${JEVENTS_LIB}")
string(COMPARE NOTEQUAL "${JEVENTS_INCDIR}" JEVENTS_INCDIR-NOTFOUND JEVENTS_INCDIR_OK)
string(COMPARE NOTEQUAL "${JEVENTS_LIB}" JEVENTS_LIB-NOTFOUND JEVENTS_LIBDIR_OK)

get_filename_component(JEVENTS_LIBDIR ${JEVENTS_LIB} PATH)

if(JEVENTS_INCDIR_OK AND JEVENTS_LIBDIR_OK)
include_directories(${JEVENTS_INCDIR})
link_directories(${JEVENTS_LIBDIR})
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES "-ljevents")
set(CMAKE_REQUIRED_FLAGS "-L${JEVENTS_LIBDIR}")
set(CMAKE_REQUIRED_INCLUDES ${JEVENTS_INCDIR})
CHECK_C_SOURCE_COMPILES("
#include <rdpmc.h>

int main()
{
    struct rdpmc_ctx ctx;
    if (rdpmc_open(PERF_COUNT_HW_CPU_CYCLES, &ctx) < 0)
          return 1;
    rdpmc_close (&ctx);
    return 0;
}
" HAVE_JEVENTS)
endif()
