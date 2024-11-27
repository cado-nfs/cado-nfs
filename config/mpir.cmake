
# You can force a path to mpir.h using the environment variables MPIR, or
# MPIR_INCDIR and MPIR_LIBDIR
string(COMPARE NOTEQUAL "$ENV{MPIR}" "" HAS_MPIR_OVERRIDE)
if (HAS_MPIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPIR} to the search path for MPIR")
    set(MPIR_INCDIR_HINTS "$ENV{MPIR}/include" ${MPIR_INCDIR_HINTS})
    set(MPIR_INCDIR_HINTS "$ENV{MPIR}"         ${MPIR_INCDIR_HINTS})
    set(MPIR_LIBDIR_HINTS "$ENV{MPIR}/lib"     ${MPIR_LIBDIR_HINTS})
    set(MPIR_LIBDIR_HINTS "$ENV{MPIR}/.libs"   ${MPIR_LIBDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPIR_INCDIR}" "" HAS_MPIR_INCDIR_OVERRIDE)
if (HAS_MPIR_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPIR_INCDIR} to the search path for MPIR")
    set(MPIR_INCDIR_HINTS "$ENV{MPIR_INCDIR}" ${MPIR_INCDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPIR_LIBDIR}" "" HAS_MPIR_LIBDIR_OVERRIDE)
if (HAS_MPIR_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPIR_LIBDIR} to the search path for MPIR")
    set(MPIR_LIBDIR_HINTS "$ENV{MPIR_LIBDIR}"     ${MPIR_LIBDIR_HINTS})
endif()

# First try overrides, really. We want cmake to shut up.
if (NOT MPIR_INCDIR)
    find_path   (MPIR_INCDIR mpir.h PATHS ${MPIR_INCDIR_HINTS} DOC "MPIR headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif()
if (NOT MPIR_INCDIR)
    find_path   (MPIR_INCDIR mpir.h HINTS ${MPIR_INCDIR_HINTS} DOC "MPIR headers"
        NO_DEFAULT_PATH
    )
endif()
if (NOT MPIR_INCDIR)
    find_path   (MPIR_INCDIR mpir.h HINTS ${MPIR_INCDIR_HINTS} DOC "MPIR headers")
endif()

find_library(MPIR_LIB    mpir   HINTS ${MPIR_LIBDIR_HINTS} DOC "MPIR library" NO_DEFAULT_PATH)
if(NOT MPIR_LIBDIR)
    find_library(MPIR_LIB    mpir   HINTS ${MPIR_LIBDIR_HINTS} DOC "MPIR library")
endif()

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(MPIR_LIBDIR ${MPIR_LIB} PATH)
message(STATUS "MPIR_INCDIR=${MPIR_INCDIR}")
message(STATUS "MPIR_LIBDIR=${MPIR_LIBDIR}")
if(MPIR_INCDIR)
include_directories(${MPIR_INCDIR})
else()
    message(FATAL_ERROR "mpir.h cannot be found. Please install MPIR, and specify its install prefix in local.sh (optionally, Gnu MP may be used as well)")
endif()
if(MPIR_LIBDIR)
link_directories(${MPIR_LIBDIR})
else()
    message(FATAL_ERROR "mpir.h cannot be found. Please install MPIR, and specify its install prefix in local.sh (optionally, Gnu MP may be used as well)")
endif()


# gmp.h matches too many times in this repository. We can't promise to
# use a placeholder for either mpir.h or gmp.h. On the other hand, if
# mpir was requested, installed in no-compatibility mode, while a
# system-wide gmp exists, we would like mpir to be caught by #include
# "gmp.h". To make sure this happens, we provide a wrapper which
# #includes mpir.h.
# (we could consider doing this only if no gmp.h exists in MPIR_INCDIR)

set(WITH_MPIR 1 CACHE INTERNAL "MPIR is being used")
set(HAVE_MPIR 1 CACHE INTERNAL "MPIR is being used")
file(WRITE ${PROJECT_BINARY_DIR}/gmp.h "#include \"mpir.h\"")
set(gmp_libname "mpir")
