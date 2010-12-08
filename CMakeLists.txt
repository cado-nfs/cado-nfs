
project(CADO_NFS)

# May be overridden by environment.
set(CADO_VERSION "1.0-rc3")

#####################################################################
# cmake version -- set_property requires 2.6. Note that haters of progress
# report things (this include me, casually) have to wait for cmake 2.7, since
# it has only recently become possible to disable them:
# http://public.kitware.com/Bug/view.php?id=8726

# Unfortunately cmake 2.6.0 at least seems to have some bugs with the
# find_program command, which prevent mpi from being properly detected.
# version 2.6.3 is known to work
cmake_minimum_required(VERSION 2.6.3)
if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#####################################################################
# What are we building ?

execute_process(COMMAND scripts/version.sh
        WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE CADO_REV)
message(STATUS "Configuring CADO revision ${CADO_REV}")

#####################################################################
# Get some overrides from the environment.

macro(get_override extname intname)
string(COMPARE NOTEQUAL "$ENV{${extname}}" "" HAS_${intname}_OVERRIDE)
if(HAS_${intname}_OVERRIDE)
    set(${intname} "$ENV{${extname}}")
    message(STATUS "${extname}=${${intname}}")
endif(HAS_${intname}_OVERRIDE)
endmacro(get_override varname)

# Syntax is (user-friendly name) (cmake name)
# Make sure everything is exported from scripts/call_cmake.sh
get_override(PREFIX               CMAKE_INSTALL_PREFIX)
get_override(CFLAGS               CMAKE_C_FLAGS)
get_override(CXXFLAGS             CMAKE_CXX_FLAGS)
get_override(CC                   CMAKE_C_COMPILER)
get_override(CXX                  CMAKE_CXX_COMPILER)
get_override(CADO_VERSION         CADO_VERSION)
get_override(MPI                  MPI)

if(EXISTS ${CADO_NFS_SOURCE_DIR}/files.dist)
    set(EXCLUDE_FROM_RELEASE_BUILD )
else(EXISTS ${CADO_NFS_SOURCE_DIR}/files.dist)
    set(EXCLUDE_FROM_RELEASE_BUILD EXCLUDE_FROM_ALL)
endif(EXISTS ${CADO_NFS_SOURCE_DIR}/files.dist)

#####################################################################
# Add the CFLAGS that are needed in all cases.
# Append them at the *end* of the parameter list so as to allow tweaking
# of e.g. warnings by the user-defined overrides.

if(CMAKE_COMPILER_IS_GNUCC)
    set (CMAKE_C_FLAGS   "-std=c99 -g -W -Wall ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-g -W -Wall ${CMAKE_CXX_FLAGS}")
endif(CMAKE_COMPILER_IS_GNUCC)

#####################################################################
# Set dependencies for external software.

# locates gmp. Obeys the env variables GMP GMP_LIBDIR GMP_INCDIR
include(${CADO_NFS_SOURCE_DIR}/config/gmp.cmake)

#####

set(GF2X_C_FLAGS "${CMAKE_C_FLAGS} -I${GMP_INCDIR}")
set(GF2X_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${GMP_INCDIR}")
set(GF2X_CONFIGURE_FLAGS --disable-shared --disable-dependency-tracking CFLAGS=${GF2X_C_FLAGS} CXXFLAGS=${GF2X_CXX_FLAGS} $ENV{GF2X_CONFIGURE_EXTRA_FLAGS})
get_override(GF2X_CONFIGURE_FLAGS GF2X_CONFIGURE_FLAGS)

#####################################################################
# Obtain information from the host.

include(CheckFunctionExists)

# posix_memalign
CHECK_FUNCTION_EXISTS(posix_memalign      HAVE_POSIX_MEMALIGN)

include(${CADO_NFS_SOURCE_DIR}/config/sse2.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/sse41.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/gcc-assembly.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/curl.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/alloca.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/cabsl.cmake)
include(${CADO_NFS_SOURCE_DIR}/config/log2.cmake)

######################################################################
# Try to detect whether pthreads and/or MPI are available.
# This checks the MPI and PTHREADS environment variables.
# As a result:
# - If pthreads are available:
#   * WITH_PTHREADS is defined to 1 as a cmake variable.
#   * -DWITH_PTHREADS is added for the compilation of every file (because
#     it can't hurt).
#   In fact, in recent versions, pthread support is mandatory. Still,
#   this piece of cmake checks the avaibility of the function
#   pthread_barrier_wait, which is optional.
# - If mpi is available:
#   * WITH_MPI is defined to 1 as a cmake variable.
#   * _NOTHING_ from the build environment is changed by default, because
#     we do not want everything to be compiled with mpicc. Thus, relevant
#     code parts should check for WITH_MPI in the CMakeLists.txt file,
#     and, if relevant, use the MPI_C_COMPILER, MPI_CXX_COMPILER, MPIEXEC
#     variables (which are guaranteed to be defined).
#   * The cado_mpi_config.h file is possibly modified in order to trigger
#     a rebuild of all mpi sources whenever the mpi library gets changed.

include (${CADO_NFS_SOURCE_DIR}/config/pthreads.cmake)
include (${CADO_NFS_SOURCE_DIR}/config/mpi.cmake)

configure_file(cado_mpi_config_h.in cado_mpi_config.h ESCAPE_QUOTES @ONLY)

#####################################################################
# OK, it's freakin ugly. Don't look.

# there would have been ways to do it without a generated helper script,
# I guess -- unfortunately I couldn't find one ;-((
configure_file(config/ccwrap.pl.in ccwrap.pl ESCAPE_QUOTES @ONLY)
set(CMAKE_C_COMPILER "${CADO_NFS_BINARY_DIR}/ccwrap.pl")

configure_file(config/cxxwrap.pl.in cxxwrap.pl ESCAPE_QUOTES @ONLY)
set(CMAKE_CXX_COMPILER "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")

set(CMAKE_C_LINKER_PREFERENCE "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")
set(CMAKE_CXX_LINKER_PREFERENCE "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")

#####################################################################
# now the real stuff.

# First, create the config file gathering everything we've got so far.
configure_file(cado_config_h.in cado_config.h ESCAPE_QUOTES @ONLY)

include_directories (${GMP_INCLUDE_PATH})
include_directories (${CADO_NFS_SOURCE_DIR})
include_directories (${CADO_NFS_SOURCE_DIR}/utils)
include_directories (${CADO_NFS_BINARY_DIR})

# link_directories (${CADO_NFS_BINARY_DIR}/utils)
# link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
# link_directories (${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x)

add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/gf2x
    COMMAND echo Configuring gf2x with options ${GF2X_CONFIGURE_FLAGS}
    COMMAND cd gf2x && ${CADO_NFS_SOURCE_DIR}/gf2x/configure ${GF2X_CONFIGURE_FLAGS}
    COMMENT "Configuring gf2x"
)

add_custom_target(gf2x-build ALL
    DEPENDS ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    WORKING_DIRECTORY gf2x
    COMMAND ${CMAKE_BUILD_TOOL}
    COMMENT "Building gf2x"
)

# set(GF2X_LIB_FILE  ${CADO_NFS_BINARY_DIR}/gf2x/.libs/libgf2x.a)
# set(GF2X-FFT_LIB_FILE ${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs/libgf2x-fft.a)
# set(CANTOR128_LIB_FILE ${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs/libcantor128.a)

# add_subdirectory (gf2x-fft)
add_subdirectory (utils)
add_subdirectory (polyselect)
add_subdirectory (sieve)
add_subdirectory (sqrt)
add_subdirectory (linalg)
add_subdirectory (filter)
add_subdirectory (linalg/bwc)

set(package_basename cado-nfs-${CADO_VERSION})

add_custom_target(dist
    COMMAND scripts/dist.sh ${package_basename}
    WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
)

install(PROGRAMS cadofactor.pl DESTINATION bin)

set(example_subdir share/${package_basename}/examples)

install(FILES params/params.c59 DESTINATION ${example_subdir})
install(FILES params/params.c79 DESTINATION ${example_subdir})

configure_file(
    ${CADO_NFS_SOURCE_DIR}/factor.sh
    ${CADO_NFS_BINARY_DIR}/factor.sh
    ESCAPE_QUOTES
    @ONLY)

configure_file(
    ${CADO_NFS_SOURCE_DIR}/scripts/fixup_params_file.sh.in
    ${CADO_NFS_BINARY_DIR}/fixup_params_file.sh
    ESCAPE_QUOTES
    @ONLY)

install(PROGRAMS ${CADO_NFS_BINARY_DIR}/factor.sh DESTINATION bin)
install(PROGRAMS ${CADO_NFS_BINARY_DIR}/fixup_params_file.sh DESTINATION bin)
