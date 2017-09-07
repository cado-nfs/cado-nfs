
project(CADO_NFS)

set(CADO_VERSION_MAJOR  3)
set(CADO_VERSION_MINOR  0)
set(CADO_VERSION_PATCHLEVEL 0)
set(CADO_VERSION_STRING  "${CADO_VERSION_MAJOR}.${CADO_VERSION_MINOR}.${CADO_VERSION_PATCHLEVEL}")

# CADO_DIST_ARCHIVE_NAME may be overridden by environment.
set(CADO_DIST_ARCHIVE_NAME "cado-nfs-${CADO_VERSION_STRING}")

#####################################################################

# If you change the min version, also change scripts/call_cmake.sh
cmake_minimum_required(VERSION 2.8.11)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
    # VERSION_GREATER_EQUAL is newish
    if (NOT CMAKE_VERSION VERSION_LESS 3.1)
        # does this get in the way with cmake tests ?? We used to have
        # CMP0054 OLD at some point.
        cmake_policy(SET CMP0054 NEW)
    endif()
endif()

if (MINGW)
    # we could do this as well in other cases.
    # http://public.kitware.com/Bug/view.php?id=8726
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set(CMAKE_COLOR_MAKEFILE OFF)
endif()

#####################################################################
# What are we building ?

execute_process(COMMAND scripts/version.sh
        WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE CADO_REV)
message(STATUS "Configuring CADO revision ${CADO_REV}")

#####################################################################
# Allow building unit test programs with "make test"
enable_testing()
# Also allow testing with "make check", but building all required
# dependencies (las for the sieve test, etc.) as well. We also add the
# "verbose" command line option here to get more meaningful output

# Note the $(ARGS) here. It is expanded by make, not by cmake, so that
# the arguments which are present on the "make check" command line are
# obeyed. This has the consequence that in order to pass a $, you need to
# pass $$ ; well, plus quoting. So that gives, e.g.:
#        make check ARGS="-R ^foo\$\$"
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V $(ARGS))

#####################################################################
# Get some overrides from the environment.

macro(get_override extname intname)
string(COMPARE NOTEQUAL "$ENV{${extname}}" "" HAS_${intname}_OVERRIDE)
if(HAS_${intname}_OVERRIDE)
    set(${intname} "$ENV{${extname}}" CACHE STRING "" FORCE)
    message(STATUS "${extname}=${${intname}} (from environment)")
endif()
endmacro()

# Syntax is (user-friendly name) (cmake name)
# Make sure everything is exported from scripts/call_cmake.sh
get_override(PREFIX               CMAKE_INSTALL_PREFIX)
get_override(CFLAGS               CMAKE_C_FLAGS)
get_override(CXXFLAGS             CMAKE_CXX_FLAGS)
get_override(MAKE                 CMAKE_MAKE_PROGRAM)
get_override(CC                   CMAKE_C_COMPILER)
get_override(CXX                  CMAKE_CXX_COMPILER)
get_override(CADO_DIST_ARCHIVE_NAME         CADO_DIST_ARCHIVE_NAME)
get_override(MPI                  MPI)
get_override(LDFLAGS              CMAKE_EXE_LINKER_FLAGS)
get_override(LDFLAGS              CMAKE_SHARED_LINKER_FLAGS)
get_override(FLAGS_SIZE           CMAKE_FLAGS_SIZE)
get_override(BWC_GF2_ARITHMETIC_BACKENDS BWC_GF2_ARITHMETIC_BACKENDS)
get_override(BWC_GF2_MATMUL_BACKENDS    BWC_GF2_MATMUL_BACKENDS)
get_override(BWC_GFP_ARITHMETIC_BACKENDS BWC_GFP_ARITHMETIC_BACKENDS)
get_override(BWC_GFP_MATMUL_BACKENDS    BWC_GFP_MATMUL_BACKENDS)
get_override(BWC_EXTRA_BACKENDS    BWC_EXTRA_BACKENDS)
get_override(CHECKS_EXPENSIVE   CHECKS_EXPENSIVE)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
	# sigh...
	# https://cmake.org/Bug/view.php?id=10342
	message(STATUS "Using ranlib -c for Darwin platforms")
	SET(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
	SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
endif()

if(EXISTS ${CADO_NFS_SOURCE_DIR}/files.dist)
    set(EXCLUDE_FROM_RELEASE_BUILD )
else()
    set(EXCLUDE_FROM_RELEASE_BUILD EXCLUDE_FROM_ALL)
endif()


set(package_basename ${CADO_DIST_ARCHIVE_NAME})
set(BINSUFFIX bin CACHE STRING "subpath below PREFIX for storing the two main scripts")
set(LIBSUFFIX lib/${package_basename} CACHE STRING "subpath below PREFIX for storing libraries and binaries")
set(DATASUFFIX share/${package_basename} CACHE STRING "subpath below PREFIX for storing parameters")

#####################################################################
# Add the CFLAGS that are needed in all cases.
# Append them at the *end* of the parameter list so as to allow tweaking
# of e.g. warnings by the user-defined overrides.
set (CMAKE_C_FLAGS      "${CMAKE_C_FLAGS} ${CMAKE_FLAGS_SIZE}")
set (CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${CMAKE_FLAGS_SIZE}")

# It's redundant with something that CMake prints the first time it
# detects the compiler. Alas, it's within a cached if, so we often don't
# see it.

# Hmmm, how delicious. cmake-2.8.7 does not have
# CMAKE_CXX_COMPILER_VERSION
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# oh my... this is really ugly. I'd love to have more usable control
# structures in cmake...
if(CMAKE_COMPILER_IS_GNUCC)
    # gcc defines __STRICT_ANSI__ with -std=c99, which under MinGW
    # causes many prototypes and macros not to be defined. Using gnu99
    # prevents the definition of __STRICT_ANSI__
    #
    # Ditto for C++, with the c++98/gnu++98 subtlety.
    set (CMAKE_C_FLAGS "-g -W -Wall ${CMAKE_C_FLAGS}")
    if (MINGW)
        set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    else()
        set (CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
    endif()
    # Someday we'll enable c11. This will provide us aligned_alloc.
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "-g -W -Wall ${CMAKE_CXX_FLAGS}")
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
        # only gcc-4.8 is supposed to be c++11-complete. But for what we
        # use, gcc-4.7 will be ok.
        message(STATUS "Enabling C++11 support")
        set(HAVE_CXX11 1)
    endif()
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7 AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
        # see comment in utils/smallset.h
        message(STATUS "Forcing -fabi-version=6 for g++-4.[789]")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=6")
    endif()
    if(HAVE_CXX11)
        if(MINGW)
            set (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
        else()
            set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
        endif()
        if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix")
        endif()
    else()
        # It seems that even g++-4.6.3, as shipped in ubuntu 12.04, is
        # sufficient for the little we use, provided we use the -std=c++0x
        # flag (-std=c++11 replaced it with gcc-4.7.0).
        if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.6.3)
            set (CMAKE_CXX_FLAGS "-std=c++0x -g -W -Wall ${CMAKE_CXX_FLAGS}")
        else()
            if(MINGW)
                set (CMAKE_CXX_FLAGS "-std=gnu++98 ${CMAKE_CXX_FLAGS}")
            else()
                set (CMAKE_CXX_FLAGS "-std=c++98 ${CMAKE_CXX_FLAGS}")
            endif()
        endif()
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    message(STATUS "Applying flags for icc")
    set (CMAKE_C_FLAGS   "-g -W -Wall -wd13200,13203 ${CMAKE_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "-g -W -Wall -wd13200,13203 ${CMAKE_CXX_FLAGS}")
    set (CMAKE_C_FLAGS   "-std=c99 ${CMAKE_C_FLAGS}")

    # only icc>=15 is full c++11 compliant, it seems. But 14 will do.
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
        set (CMAKE_CXX_FLAGS   "-std=c++11 ${CMAKE_CXX_FLAGS}")
        set(HAVE_CXX11 1)
        message(STATUS "Enabling C++11 support")
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
        message(STATUS "Not enabling C++11 support for ICC 14, see bug #21416")
    else()
        # don't set standard flags explicitly. The default seems fine.
        # set (CMAKE_CXX_FLAGS   "-std=gnu++98 -g -W -Wall ${CMAKE_CXX_FLAGS}")
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_HOST_SYSTEM_NAME
        STREQUAL "Darwin" AND CMAKE_C_COMPILER_ID MATCHES "Clang" )
    set (CMAKE_C_FLAGS   "-std=c99 -g -W -Wall ${CMAKE_C_FLAGS}")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR CMAKE_HOST_SYSTEM_NAME
        STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    # Not the slightest idea when apple's clang became c++11-compliant.
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.0.0")
        set (CMAKE_CXX_FLAGS   "-std=c++11 ${CMAKE_CXX_FLAGS}")
        set (HAVE_CXX11 1)
        message(STATUS "Enabling C++11 support")
    endif()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "^Clang")
    # According to http://clang.llvm.org/cxx_status.html
    # Caveat: Apple's clang identifies itself as plain Clang on some
    # occasions... I don't exactly know when and how, but there's some
    # potential for misconfiguration.
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
        set (CMAKE_CXX_FLAGS   "-std=c++11 ${CMAKE_CXX_FLAGS}")
        set (HAVE_CXX11 1)
        message(STATUS "Enabling C++11 support")
        set (CMAKE_CXX_FLAGS "-std=c++11 -g -W -Wall ${CMAKE_CXX_FLAGS}")
    endif()
endif()
message(STATUS "CFLAGS=${CMAKE_C_FLAGS} (synthetized)")
message(STATUS "CXXFLAGS=${CMAKE_CXX_FLAGS} (synthetized)")

if(NOT HAVE_CXX11)
    message(FATAL_ERROR "C++11 is now required to build cado-nfs")
endif()

#####################################################################
# Set dependencies for external software.
include(${CADO_NFS_SOURCE_DIR}/config/utilities.cmake)

set(some_mpir_variable "$ENV{MPIR}$ENV{MPIR_INCDIR}$ENV{MPIR_LIBDIR}")
if(some_mpir_variable)
# locates gmp. Obeys the env variables GMP GMP_LIBDIR GMP_INCDIR
include(${CADO_NFS_SOURCE_DIR}/config/mpir.cmake)
string_join(GF2X_C_FLAGS " " ${CMAKE_C_FLAGS} "-I${MPIR_INCDIR}" "-I${CADO_NFS_SOURCE_DIR}")
string_join(GF2X_CXX_FLAGS " " ${CMAKE_CXX_FLAGS} "-I${MPIR_INCDIR}" "-I${CADO_NFS_SOURCE_DIR}")
string_join(GF2X_LDFLAGS " " ${CMAKE_SHARED_LINKER_FLAGS} "-Wl,-rpath,${MPIR_LIBDIR}" "-L${MPIR_LIBDIR}")
else()
# locates gmp. Obeys the env variables GMP GMP_LIBDIR GMP_INCDIR
include(${CADO_NFS_SOURCE_DIR}/config/gmp.cmake)
string_join(GF2X_C_FLAGS " " ${CMAKE_C_FLAGS} "-I${GMP_INCDIR}" "-I${CADO_NFS_SOURCE_DIR}")
string_join(GF2X_CXX_FLAGS " " ${CMAKE_CXX_FLAGS} "-I${GMP_INCDIR}" "-I${CADO_NFS_SOURCE_DIR}")
string_join(GF2X_LDFLAGS " " ${CMAKE_SHARED_LINKER_FLAGS} "-Wl,-rpath,${GMP_LIBDIR}" "-L${GMP_LIBDIR}")
endif()

#####
include(CheckTypeSize)
CHECK_TYPE_SIZE("unsigned long int" ulong_size)
CHECK_TYPE_SIZE("unsigned long long int" ulonglong_size)
math (EXPR ULONG_BITS "8*${ulong_size}")
math (EXPR ULONGLONG_BITS "8*${ulonglong_size}")
message(STATUS "bits in unsigned long: ${ULONG_BITS}")
message(STATUS "bits in unsigned long long: ${ULONGLONG_BITS}")

include(${CADO_NFS_SOURCE_DIR}/config/gmp_randstate.cmake)

set(GF2X_CONFIGURE_FLAGS --disable-shared --disable-dependency-tracking --enable-fft-interface CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${GF2X_C_FLAGS} CXXFLAGS=${GF2X_CXX_FLAGS} LDFLAGS=${GF2X_LDFLAGS} $ENV{GF2X_CONFIGURE_EXTRA_FLAGS})

#####################################################################
# Obtain information from the host.

if (DEFINED ENV{NO_SSE})
    list(APPEND GF2X_CONFIGURE_FLAGS --disable-sse2)
    message(STATUS "As per $NO_SSE, disabling all simd intrinsics")
elseif (ULONG_BITS EQUAL 32)
    list(APPEND GF2X_CONFIGURE_FLAGS --disable-sse2)
    message(STATUS "32-bit platform, lazily disabling all simd intrinsics")
else()
    include("${CADO_NFS_SOURCE_DIR}/config/sse2.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/sse3.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/ssse3.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/sse41.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/avx.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/avx2.cmake")
    include("${CADO_NFS_SOURCE_DIR}/config/pclmul.cmake")
endif()

get_override(GF2X_CONFIGURE_FLAGS GF2X_CONFIGURE_FLAGS)

if (NOT DEFINED ENV{NO_INLINE_ASSEMBLY})
  include("${CADO_NFS_SOURCE_DIR}/config/gcc-inline-assembly.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/gcc-bugs.cmake")
if (NOT DEFINED ENV{NO_GAS_ASSEMBLY})
  include("${CADO_NFS_SOURCE_DIR}/config/gas-assembly.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/curl.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/alloca.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/math.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/posix.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/rusage_thread.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/sync_fetch.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/asprintf.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/mingw.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/test-rand.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/odd_address_io_int.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/strlcpy.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/linux.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/hwloc.cmake")
IF (NOT DEFINED ENV{NO_PYTHON_CHECK})
  include("${CADO_NFS_SOURCE_DIR}/config/python.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/perl.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/shared_ptr.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/cxxabi.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/regex.cmake")

set(CADO_SHARED)

if("$ENV{ENABLE_SHARED}")
    # cmake scripts here and there use the cmake variable, not the
    # environment variable.
    set(ENABLE_SHARED 1)
    if(MINGW)
        message(FATAL_ERROR "ENABLE_SHARED is not supported with MINGW")
    endif()
    set(CADO_SHARED SHARED)
    message(STATUS "Building with dynamic libraries (ENABLE_SHARED=$ENV{ENABLE_SHARED})")
    string_join(rpath ":"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/utils"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/polyselect"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/linalg/bwc"
    )
    set(CMAKE_INSTALL_RPATH ${rpath})
else()
    message(STATUS "Building with static libraries (ENABLE_SHARED=\"$ENV{ENABLE_SHARED}\")")
endif()




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
include (${CADO_NFS_SOURCE_DIR}/config/openmp.cmake)

message(STATUS "Using libs: ${math_libs} (math) ${pthread_libs} (pthreads)")

configure_file(cado_mpi_config_h.in cado_mpi_config.h ESCAPE_QUOTES @ONLY)

#####################################################################
# OK, it's freakin ugly. Don't look.

# there would have been ways to do it without a generated helper script,
# I guess -- unfortunately I couldn't find one ;-((
configure_file(config/ccwrap.pl.in ccwrap.pl ESCAPE_QUOTES @ONLY)
configure_file(config/cxxwrap.pl.in cxxwrap.pl ESCAPE_QUOTES @ONLY)

set(CMAKE_C_COMPILER "${CADO_NFS_BINARY_DIR}/ccwrap.pl")
set(CMAKE_CXX_COMPILER "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")

# It's weird. This thing does not seem to apply, or at least not always. I've
# cured the problem by tweaking ccwrap, and let it decide whether it's linking
# or not. If yes, then it uses g++ rather than gcc (or mpic++ vs mpicc)
set(CMAKE_C_LINKER_PREFERENCE "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")
set(CMAKE_CXX_LINKER_PREFERENCE "${CADO_NFS_BINARY_DIR}/cxxwrap.pl")

#####################################################################
# now the real stuff.

# First, create the config file gathering everything we've got so far.
configure_file(cado_config_h.in cado_config.h ESCAPE_QUOTES @ONLY)

include_directories (${CADO_NFS_SOURCE_DIR})
include_directories (${CADO_NFS_SOURCE_DIR}/utils)
include_directories (${CADO_NFS_BINARY_DIR})

# link_directories (${CADO_NFS_BINARY_DIR}/utils)
# link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x)


############################################################



add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/gf2x
    COMMAND echo Configuring gf2x with options ${GF2X_CONFIGURE_FLAGS}
    COMMAND cd gf2x && ${CADO_NFS_SOURCE_DIR}/gf2x/configure ${GF2X_CONFIGURE_FLAGS} < /dev/null
    COMMENT "Configuring gf2x"
)

add_custom_target(gf2x-build ALL
    DEPENDS ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    WORKING_DIRECTORY gf2x
    COMMAND "${CMAKE_MAKE_PROGRAM}"
    COMMENT "Building gf2x"
)

# set(GF2X_LIB_FILE  ${CADO_NFS_BINARY_DIR}/gf2x/.libs/libgf2x.a)
# set(GF2X-FFT_LIB_FILE ${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs/libgf2x-fft.a)
# set(CANTOR128_LIB_FILE ${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs/libcantor128.a)


# add_subdirectory (gf2x-fft)
add_subdirectory (utils)
add_subdirectory (polyselect)
add_subdirectory (numbertheory)
add_subdirectory (sieve)
add_subdirectory (sqrt)
add_subdirectory (linalg)
add_subdirectory (filter)
add_subdirectory (linalg/bwc)
add_subdirectory (misc)
add_subdirectory (scripts/cadofactor)
add_subdirectory (tests EXCLUDE_FROM_ALL)
#add_subdirectory (nfs-hd)

add_custom_target(dist
    COMMAND scripts/dist.sh ${package_basename}
    WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
)

# The trailing / is very important here. Without it, we would create a
# parameters/ directory in ${PREFIX}/share/cado-nfs-X.Y/
install(DIRECTORY parameters/ DESTINATION ${DATASUFFIX})

configure_file(cado-nfs.py cado-nfs.py ESCAPE_QUOTES @ONLY)
configure_file(cado-nfs-client.py cado-nfs-client.py ESCAPE_QUOTES @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cado-nfs.py DESTINATION ${BINSUFFIX})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cado-nfs-client.py DESTINATION ${BINSUFFIX})

# Some scripts are handy to have also in the bin directory.
foreach(script split.py descent.py)
    configure_file(
        ${CADO_NFS_SOURCE_DIR}/scripts/${script}
        ${CADO_NFS_BINARY_DIR}/scripts/${script}
        COPYONLY)
endforeach(script)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/scripts/descent.py DESTINATION ${LIBSUFFIX}/scripts)


message(STATUS "Main scripts will go to ${CMAKE_INSTALL_PREFIX}/${BINSUFFIX}")
message(STATUS "Libraries and binaries will go to ${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}")
message(STATUS "Data will go to ${CMAKE_INSTALL_PREFIX}/${DATASUFFIX}")

IF (DEFINED ENV{CMAKE_DUMP_VARIABLES})
  include("${CADO_NFS_SOURCE_DIR}/config/variables.cmake")
ENDIF()

if(MINGW)
    message(STATUS "Please look up dev_docs/howto-MinGW.txt in the git tree for our mingw testing configuration, and frequently encountered errors")
endif()
