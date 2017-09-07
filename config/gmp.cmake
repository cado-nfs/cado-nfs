
# You can force a path to gmp.h using the environment variables GMP, or
# GMP_INCDIR and GMP_LIBDIR
string(COMPARE NOTEQUAL "$ENV{GMP}" "" HAS_GMP_OVERRIDE)
if (HAS_GMP_OVERRIDE)
    message(STATUS "Adding $ENV{GMP} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP}/include" ${GMP_INCDIR_HINTS})
    set(GMP_INCDIR_HINTS "$ENV{GMP}"         ${GMP_INCDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/lib"     ${GMP_LIBDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/.libs"   ${GMP_LIBDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{GMP_INCDIR}" "" HAS_GMP_INCDIR_OVERRIDE)
if (HAS_GMP_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_INCDIR} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP_INCDIR}" ${GMP_INCDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{GMP_LIBDIR}" "" HAS_GMP_LIBDIR_OVERRIDE)
if (HAS_GMP_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_LIBDIR} to the search path for Gnu MP")
    set(GMP_LIBDIR_HINTS "$ENV{GMP_LIBDIR}"     ${GMP_LIBDIR_HINTS})
endif()

# First try overrides, really. We want cmake to shut up.
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h PATHS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif()
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h HINTS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers"
        NO_DEFAULT_PATH
    )
endif()
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h HINTS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers")
endif()

find_library(GMP_LIB    gmp   HINTS ${GMP_LIBDIR_HINTS} DOC "Gnu MP library" NO_DEFAULT_PATH)
if(NOT GMP_LIBDIR)
find_library(GMP_LIB    gmp   HINTS ${GMP_LIBDIR_HINTS} DOC "Gnu MP library")
endif()

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(GMP_LIBDIR ${GMP_LIB} PATH)
message(STATUS "GMP_INCDIR=${GMP_INCDIR}")
message(STATUS "GMP_LIBDIR=${GMP_LIBDIR}")
if(GMP_INCDIR)
include_directories(${GMP_INCDIR})
else()
message(FATAL_ERROR "gmp.h cannot be found. Please install Gnu MP, and specify its install prefix in local.sh")
endif()
if(GMP_LIBDIR)
link_directories(${GMP_LIBDIR})
else()
message(FATAL_ERROR "GMP library cannot be found. Please install Gnu MP, and specify its install prefix in local.sh")
endif()

set(WITH_GMP 1 CACHE INTERNAL "GMP is being used")
set(HAVE_GMP 1 CACHE INTERNAL "GMP is being used")
set(gmp_libname "gmp")

#pointless, since we have workarounds anyway
#message(STATUS "Testing whether gmpxx is available")
#try_compile(HAVE_GMPXX
#    ${PROJECT_BINARY_DIR}/config
#    ${PROJECT_SOURCE_DIR}/config/gmpxx-test.cpp
#    CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${gmp_libname};${gmp_libname}xx")
#
#if (HAVE_GMPXX)
#message(STATUS "Testing whether gmpxx is available -- yes")
#else()
#message(STATUS "Testing whether gmpxx is available -- no")
#endif()
