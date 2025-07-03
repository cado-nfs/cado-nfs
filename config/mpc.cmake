# You can force a path to mpc.h using the environment variables MPC, or
# MPC_INCDIR and MPC_LIBDIR
string(COMPARE NOTEQUAL "$ENV{MPC}" "" HAS_MPC_OVERRIDE)
if (HAS_MPC_OVERRIDE)
    message(STATUS "Adding $ENV{MPC} to the search path for Gnu MPC")
    set(MPC_INCDIR_HINTS "$ENV{MPC}/include" ${MPC_INCDIR_HINTS})
    set(MPC_INCDIR_HINTS "$ENV{MPC}"         ${MPC_INCDIR_HINTS})
    set(MPC_LIBDIR_HINTS "$ENV{MPC}/lib"     ${MPC_LIBDIR_HINTS})
    set(MPC_LIBDIR_HINTS "$ENV{MPC}/.libs"   ${MPC_LIBDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPC_INCDIR}" "" HAS_MPC_INCDIR_OVERRIDE)
if (HAS_MPC_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPC_INCDIR} to the search path for Gnu MPC")
    set(MPC_INCDIR_HINTS "$ENV{MPC_INCDIR}" ${MPC_INCDIR_HINTS})
endif()
string(COMPARE NOTEQUAL "$ENV{MPC_LIBDIR}" "" HAS_MPC_LIBDIR_OVERRIDE)
if (HAS_MPC_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{MPC_LIBDIR} to the search path for Gnu MPC")
    set(MPC_LIBDIR_HINTS "$ENV{MPC_LIBDIR}"     ${MPC_LIBDIR_HINTS})
endif()

# First try overrides, really. We want cmake to shut up.
if (NOT MPC_INCDIR)
    find_path   (MPC_INCDIR mpc.h PATHS ${MPC_INCDIR_HINTS} DOC "Gnu MPC headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif()
if (NOT MPC_INCDIR)
    find_path   (MPC_INCDIR mpc.h HINTS ${MPC_INCDIR_HINTS} DOC "Gnu MPC headers"
        NO_DEFAULT_PATH
    )
endif()
if (NOT MPC_INCDIR)
    find_path   (MPC_INCDIR mpc.h HINTS ${MPC_INCDIR_HINTS} DOC "Gnu MPC headers")
endif()

find_library(MPC_LIB    mpc   HINTS ${MPC_LIBDIR_HINTS} DOC "Gnu MPC library" NO_DEFAULT_PATH)
if(NOT MPC_LIBDIR)
    find_library(MPC_LIB    mpc   HINTS ${MPC_LIBDIR_HINTS} DOC "Gnu MPC library")
endif()

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(MPC_LIBDIR ${MPC_LIB} PATH)
message(STATUS "MPC_INCDIR=${MPC_INCDIR}")
message(STATUS "MPC_LIBDIR=${MPC_LIBDIR}")

if(MPC_INCDIR)
    include_directories(${MPC_INCDIR})
endif()
if(MPC_LIBDIR)
    link_directories(${MPC_LIBDIR})
endif()

string(COMPARE NOTEQUAL "${MPC_INCDIR}" MPC_INCDIR-NOTFOUND MPC_INCDIR_OK)
string(COMPARE NOTEQUAL "${MPC_LIB}" MPC_LIB-NOTFOUND MPC_LIBDIR_OK)

get_filename_component(MPC_LIBDIR ${MPC_LIB} PATH)

if(MPC_INCDIR_OK AND MPC_LIBDIR_OK)
    set(WITH_MPC 1 CACHE INTERNAL "MPC is being used")
    set(HAVE_MPC 1 CACHE INTERNAL "MPC is being used")
else()
    message(STATUS "MPC library cannot be found, some code will be disabled")
endif()
