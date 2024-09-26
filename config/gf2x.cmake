
# You can force a path to gf2x.h using the environment variables GF2X, or
# GF2X_INCDIR and GF2X_LIBDIR
string(COMPARE NOTEQUAL "$ENV{GF2X}" "" HAS_GF2X_OVERRIDE)
if (HAS_GF2X_OVERRIDE)
    message(STATUS "Adding $ENV{GF2X} to the search path for gf2x")
    set(GF2X_INCDIR_HINTS "$ENV{GF2X}/include" ${GF2X_INCDIR_HINTS})
    set(GF2X_INCDIR_HINTS "$ENV{GF2X}"         ${GF2X_INCDIR_HINTS})
    set(GF2X_LIBDIR_HINTS "$ENV{GF2X}/lib"     ${GF2X_LIBDIR_HINTS})
    set(GF2X_LIBDIR_HINTS "$ENV{GF2X}/.libs"   ${GF2X_LIBDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{GF2X_INCDIR}" "" HAS_GF2X_INCDIR_OVERRIDE)
if (HAS_GF2X_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GF2X_INCDIR} to the search path for gf2x")
    set(GF2X_INCDIR_HINTS "$ENV{GF2X_INCDIR}" ${GF2X_INCDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{GF2X_LIBDIR}" "" HAS_GF2X_LIBDIR_OVERRIDE)
if (HAS_GF2X_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GF2X_LIBDIR} to the search path for gf2x")
    set(GF2X_LIBDIR_HINTS "$ENV{GF2X_LIBDIR}"     ${GF2X_LIBDIR_HINTS})
endif()

# First try overrides, really. We want cmake to shut up.
if (NOT GF2X_INCDIR)
    find_path   (GF2X_INCDIR gf2x.h PATHS ${GF2X_INCDIR_HINTS} DOC "gf2x headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif()
if (NOT GF2X_INCDIR)
    find_path   (GF2X_INCDIR gf2x.h HINTS ${GF2X_INCDIR_HINTS} DOC "gf2x headers"
        NO_DEFAULT_PATH
    )
endif()
if (NOT GF2X_INCDIR)
    find_path   (GF2X_INCDIR gf2x.h HINTS ${GF2X_INCDIR_HINTS} DOC "gf2x headers")
endif()

find_library(GF2X_LIB    gf2x   HINTS ${GF2X_LIBDIR_HINTS} DOC "gf2x library" NO_DEFAULT_PATH)
if(NOT GF2X_LIBDIR)
find_library(GF2X_LIB    gf2x   HINTS ${GF2X_LIBDIR_HINTS} DOC "gf2x library")
endif()

if(GF2X_INCDIR AND NOT EXISTS "${GF2X_INCDIR}/gf2x-cantor-field-impl.h")
    message(STATUS "An installed version of gf2x was found in ${GF2X_INCDIR}, but it lacks the the gf2x-cantor-impl.h file, which we need. Therefore we use our embedded version instead")
    set(GF2X_INCDIR GF2X_INCDIR-NOTFOUND)
    set(GF2X_LIB GF2X_LIB-NOTFOUND)
endif()

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(GF2X_LIBDIR ${GF2X_LIB} PATH)
message(STATUS "GF2X_INCDIR=${GF2X_INCDIR}")
message(STATUS "GF2X_LIBDIR=${GF2X_LIBDIR}")
if(GF2X_INCDIR)
    include_directories(${GF2X_INCDIR})
endif()
if(GF2X_LIBDIR)
    link_directories(${GF2X_LIBDIR})
endif()

string(COMPARE NOTEQUAL "${GF2X_INCDIR}" GF2X_INCDIR-NOTFOUND GF2X_INCDIR_OK)
string(COMPARE NOTEQUAL "${GF2X_LIB}" GF2X_LIB-NOTFOUND GF2X_LIBDIR_OK)

get_filename_component(GF2X_LIBDIR ${GF2X_LIB} PATH)

if(GF2X_INCDIR_OK AND GF2X_LIBDIR_OK)
    set(WITH_GF2X 1 CACHE INTERNAL "GF2X is being used")
    set(HAVE_GF2X 1 CACHE INTERNAL "GF2X is being used")
    if(HAVE_GF2X)
        message(STATUS "Using the gf2x library found on the system")
    else()
        message(FATAL_ERROR "A version of the gf2x library was found by cmake, but it apparently does not fit our criteria. This is a fatal error, as we cannot be sure that our embedded copy will work correctly in that case")
    endif()
elseif(GF2X_INCDIR_OK OR GF2X_LIBDIR_OK)
    message(FATAL_ERROR "A partly installed version of the gf2x library was found by cmake, but it lacks either the headers or the library. This is a fatal error, as we cannot be sure that our embedded copy will work correctly in that case")
else()
    message(STATUS "Using the embedded gf2x library")
endif()

set(gf2x_libname "gf2x")
