#####################################################################
# First get some overrides from the environment. Maybe it's a rather bad
# idea, but cado-nfs likes to pass things to cmake as environment
# variable. We arrange so that the environment variables take precedence
# over the cmake command-line.  Therefore the code below is called
# _before_ project(). This touches the variable cache. Since we FORCE the
# variable into the cache, it's stronger than cmake command-line
# overrides.

macro(get_override extname intname)
if(DEFINED ENV{${extname}})
    set(HAS_${intname}_OVERRIDE 1)
    set(${intname} "$ENV{${extname}}" CACHE STRING "" FORCE)
    message(STATUS "${extname}=${${intname}} (from environment)")
else()
    set(HAS_${intname}_OVERRIDE 0)
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
get_override(TIMEOUT_SCALE   TIMEOUT_SCALE)
get_override(VALGRIND   VALGRIND)

if(VALGRIND)
    set(USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT 1)
endif()



# If you change the min version, also change scripts/call_cmake.sh
# It must be called before the project() command.
cmake_minimum_required(VERSION 3.18.0)

# if we want to use JOB_SERVER_AWARE, we need 3.28.0
# cmake_minimum_required(VERSION 3.28.0)

project(CADO_NFS)

# must be done early, but after project()
include(config/xcode.cmake)

set(CADO_VERSION_MAJOR  3)
set(CADO_VERSION_MINOR  0)
set(CADO_VERSION_PATCHLEVEL 0)
set(CADO_VERSION_STRING  "${CADO_VERSION_MAJOR}.${CADO_VERSION_MINOR}.${CADO_VERSION_PATCHLEVEL}")

# CADO_DIST_ARCHIVE_NAME may be overridden by environment.
set(CADO_DIST_ARCHIVE_NAME "cado-nfs-${CADO_VERSION_STRING}")

#####################################################################

cmake_policy(SET CMP0003 NEW)

# does this get in the way with cmake tests ?? We used to have
# CMP0054 OLD at some point.
cmake_policy(SET CMP0054 NEW)

# IN_LIST
cmake_policy(SET CMP0057 NEW)

# Here's how we set cmake policies that are unknown to old cmake
# versions.
#
# VERSION_GREATER_EQUAL is for cmake 3.7+
# if (NOT CMAKE_VERSION VERSION_LESS 3.8)
#    cmake_policy(SET CMP1234 NEW)
# endif()

if (MINGW)
    # we could do this as well in other cases.
    # http://public.kitware.com/Bug/view.php?id=8726
    set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
    set(CMAKE_COLOR_MAKEFILE OFF)
endif()

set(CMAKE_REQUIRED_QUIET 0)

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

if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
# Note the $(ARGS) here. It is expanded by make, not by cmake, so that
# the arguments which are present on the "make check" command line are
# obeyed. This has the consequence that in order to pass a $, you need to
# pass $$ ; well, plus quoting. So that gives, e.g.:
#        make check ARGS="-R ^foo\$\$"
add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} -V $(ARGS)
    # JOB_SERVER_AWARE true
)
endif()

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

# It's redundant with something that CMake prints the first time it
# detects the compiler. Alas, it's within a cached if, so we often don't
# see it.

set (CMAKE_C_FLAGS      "${CMAKE_C_FLAGS} ${CMAKE_FLAGS_SIZE}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    # gcc defines __STRICT_ANSI__ with -std=c99, which under MinGW
    # causes many prototypes and macros not to be defined. Using gnu99
    # prevents the definition of __STRICT_ANSI__
    #
    # Ditto for C++, with the c++98/gnu++98 subtlety.
    set (CMAKE_C_FLAGS "-Wall -Wextra ${CMAKE_C_FLAGS}")
    if (MINGW)
        set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    else()
        set (CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
    endif()
    # Someday we'll enable c11. This will provide us aligned_alloc.
elseif(CMAKE_C_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Applying flags for icc (icx)")
    # icc defaults to -fp-model=fast, which is slightly getting in the
    # way. We rarely depend on this in cado-nfs, so let's stick to usual
    # behavior.
    set (CMAKE_C_FLAGS   "-fp-model=precise -std=c99 ${CMAKE_C_FLAGS}")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
    message(STATUS "Applying flags for icc")
    set (CMAKE_C_FLAGS   "-W -Wall -wd13200,13203 ${CMAKE_C_FLAGS}")
    set (CMAKE_C_FLAGS   "-std=c99 ${CMAKE_C_FLAGS}")
elseif(CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_HOST_SYSTEM_NAME
        STREQUAL "Darwin" AND CMAKE_C_COMPILER_ID MATCHES "Clang" )
    set (CMAKE_C_FLAGS   "-std=c99 -W -Wall ${CMAKE_C_FLAGS}")
elseif(CMAKE_C_COMPILER_ID MATCHES "^Clang")
    set (CMAKE_C_FLAGS   "-std=c99 -W -Wall ${CMAKE_C_FLAGS}")
else()
    message(WARNING "C compiler not recognized. Setting -std=c99 as a wild guess")
    set (CMAKE_C_FLAGS   "-std=c99 -W -Wall ${CMAKE_C_FLAGS}")
endif()
message(STATUS "CFLAGS=${CMAKE_C_FLAGS} (synthetized)")



set (CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${CMAKE_FLAGS_SIZE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "-Wall -Wextra ${CMAKE_CXX_FLAGS}")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
        # according to https://gcc.gnu.org/projects/cxx-status.html
        # most of the support is in g++-10. Some important things that we
        # probably want to use, such as using enum, are in g++-11.
        message(FATAL_ERROR "cado-nfs requires g++ version 10 or greater")
    endif()
    # I observe
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=109717 with g++13
    # and not g++14 ; it's a false positive with libfmt.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0 AND
            CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 13.0)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
    endif()

    set (CMAKE_CXX_FLAGS "-std=c++20 -Wno-c++11-compat ${CMAKE_CXX_FLAGS}")
    message(STATUS "Enabling C++20 support")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix")
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0 AND
            CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        # https://gcc.gnu.org/ml/gcc/2017-05/msg00073.html
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
    endif()
    set (CMAKE_CXX_FLAGS "-export-dynamic ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_C_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Applying flags for icc (icx)")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "2023")
        message(FATAL_ERROR "cado-nfs requires icx version 2023 or greater")
    endif()
    set (CMAKE_CXX_FLAGS "-fp-model=precise -std=c++20 ${CMAKE_CXX_FLAGS}")
    message(STATUS "Enabling C++20 support")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
    message(FATAL_ERROR "cado-nfs requires icx and does not work with older icc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" OR CMAKE_HOST_SYSTEM_NAME
        STREQUAL "Darwin" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    # 16 is probably a bit excessive. Vanilla clang-14 works, so it might
    # be that apple clang 14 also works. But I have no way to test.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0.0")
        message(FATAL_ERROR "cado-nfs requires Apple Clang version 16 or greater")
    endif()
    set (CMAKE_CXX_FLAGS   "-std=c++20 ${CMAKE_CXX_FLAGS}")
    message(STATUS "Enabling C++20 support")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^Clang")
    set (CMAKE_CXX_FLAGS "-W -Wall ${CMAKE_CXX_FLAGS}")
    # According to http://clang.llvm.org/cxx_status.html
    # deciding when to enable a given feature set is always some hit and
    # miss, because it takes time (like, a multiple version timespan) for
    # support to be considered complete.
    # Caveat: Apple's clang identifies itself as plain Clang on some
    # occasions... I don't exactly know when and how, but there's some
    # potential for misconfiguration.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "12")
        message(FATAL_ERROR "cado-nfs requires Clang version 12 or greater")
    endif()
    set (CMAKE_CXX_FLAGS   "-std=c++20 ${CMAKE_CXX_FLAGS}")
    message(STATUS "Enabling C++20 support")
else()
    message(WARNING "C++ compiler not recognized. Setting -std=c++20 as a wild guess")
    set (CMAKE_CXX_FLAGS   "-std=c++20 ${CMAKE_CXX_FLAGS}")
endif()
message(STATUS "CXXFLAGS=${CMAKE_CXX_FLAGS} (synthetized)")



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

set(OLD_CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET 1)
include(${CADO_NFS_SOURCE_DIR}/config/check_types.cmake)
set(CMAKE_REQUIRED_QUIET "${OLD_CMAKE_REQUIRED_QUIET}")

include(${CADO_NFS_SOURCE_DIR}/config/pagesize.cmake)

include(${CADO_NFS_SOURCE_DIR}/config/gmp_randstate.cmake)

#####################################################################
# Obtain information from the host.

include(config/check_cpu_feature.cmake)

set(GF2X_CONFIGURE_FLAGS)

if (DEFINED ENV{NO_SSE})
    list(APPEND GF2X_CONFIGURE_FLAGS --disable-sse2)
    message(STATUS "As per $NO_SSE, disabling all simd intrinsics")
elseif (ULONG_BITS EQUAL 32)
    list(APPEND GF2X_CONFIGURE_FLAGS --disable-sse2)
    check_cpu_feature(mmx)
    message(STATUS "32-bit platform, lazily disabling all Intel simd intrinsics")
else()
    # note that even 64-bit mean that we must test mmx (think arm64 !)
    check_cpu_feature(mmx)
    check_cpu_feature(sse2)
    check_cpu_feature(sse3 PREREQUISITES HAVE_SSE2)
    check_cpu_feature(ssse3 PREREQUISITES HAVE_SSE3)
    check_cpu_feature(sse4.1 PREREQUISITES HAVE_SSSE3)
    check_cpu_feature(popcnt PREREQUISITES HAVE_SSSE3)
    check_cpu_feature(avx PREREQUISITES HAVE_SSE41)
    check_cpu_feature(pclmul PREREQUISITES HAVE_SSE41)
    set(b USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT)
    check_cpu_feature(avx2 PREREQUISITES HAVE_AVX BLOCKERS ${b})
    check_cpu_feature(avx512f PREREQUISITES HAVE_AVX BLOCKERS ${b})
    include(config/bug30057.cmake)
    check_cpu_feature(avx512vl PREREQUISITES HAVE_AVX512F BLOCKERS ${b})
    include(config/bug30027.cmake)
    check_cpu_feature(avx512dq PREREQUISITES HAVE_AVX512F BLOCKERS ${b})

    if (HAVE_SSE2 AND USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT)
        set(fallback "westmere")
        if (NOT HAVE_AVX OR NOT HAVE_PCLMUL OR NOT HAVE_POPCNT)
            set(fallback "x86-64")
        endif()
        message(STATUS "As per USE_ONLY_ASSEMBLY_INSTRUCTIONS_THAT_VALGRIND_KNOWS_ABOUT, restricting gf2x build to old enough microarchitecture (${fallback})")
        set(GF2X_C_FLAGS "${GF2X_C_FLAGS} -march=${fallback}")
    endif()
endif()

if (DEFINED ENV{NO_NEON})
    message(STATUS "As per $NO_NEON, disabling all ARM NEON intrinsics")
else()
    check_cpu_feature(neon COMPILER_FLAG -mfpu=neon DISPLAYNAME "ARM NEON")
endif()

if (NOT DEFINED ENV{NO_INLINE_ASSEMBLY})
  include("${CADO_NFS_SOURCE_DIR}/config/gcc-inline-assembly.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/gcc-bugs.cmake")
if (NOT DEFINED ENV{NO_GAS_ASSEMBLY})
  include("${CADO_NFS_SOURCE_DIR}/config/gas-assembly.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/glibc.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/musl.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/math.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/math-special.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/posix.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/bsd.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/rusage_thread.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/clock_gettime.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/sync_fetch.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/asprintf.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/mingw.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/test-rand.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/strlcpy.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/linux.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/hwloc.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/fmt.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/gf2x.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/gmpecm.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/mpfr.cmake")
if(HAVE_MPFR)
include("${CADO_NFS_SOURCE_DIR}/config/mpc.cmake")
endif()
IF (NOT DEFINED ENV{NO_PYTHON_CHECK})
  include("${CADO_NFS_SOURCE_DIR}/config/python.cmake")
endif()
include("${CADO_NFS_SOURCE_DIR}/config/perl.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/cxxabi.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/regex.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/jevents.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/aligned_alloc.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/alignas.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/int128.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/m4ri.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/detect-static-analysis.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/detect-coverage.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/execinfo.cmake")
CHECK_INCLUDE_FILES (sys/endian.h HAVE_SYS_ENDIAN_H)
include("${CADO_NFS_SOURCE_DIR}/config/bswap32.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/compression.cmake")
include("${CADO_NFS_SOURCE_DIR}/config/gdb.cmake")

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_C_FLAGS "-Werror=all -Werror=extra ${CMAKE_C_FLAGS}")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "-Werror=extra -Werror=all ${CMAKE_CXX_FLAGS}")
endif()

set(CADO_SHARED)

if("$ENV{ENABLE_SHARED}" OR BUILD_SHARED_LIBS)
    # cmake scripts here and there use the cmake variable, not the
    # environment variable.
    set(ENABLE_SHARED 1)
    set(CADO_NFS_USE_SHARED_LIBS 1)
    if(MINGW)
        message(FATAL_ERROR "ENABLE_SHARED is not supported with MINGW")
    endif()
    set(CADO_SHARED SHARED)
    message(STATUS "Building with dynamic libraries (ENABLE_SHARED=$ENV{ENABLE_SHARED})")

    # https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    string_join(rpath ":"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/utils"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/polyselect"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/linalg/bblas"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/linalg/bwc"
        "${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/gf2x"
    )
    set(CMAKE_INSTALL_RPATH ${rpath})
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    message(STATUS "Building with static libraries (ENABLE_SHARED=\"$ENV{ENABLE_SHARED}\")")
endif()

if(NOT ENABLE_SHARED)
    list(APPEND GF2X_CONFIGURE_FLAGS --disable-shared)
endif()

list(APPEND GF2X_CONFIGURE_FLAGS --disable-dependency-tracking)
list(APPEND GF2X_CONFIGURE_FLAGS --prefix ${CMAKE_INSTALL_PREFIX} --libdir ${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}/gf2x)
list(APPEND GF2X_CONFIGURE_FLAGS CC=${CMAKE_C_COMPILER})
list(APPEND GF2X_CONFIGURE_FLAGS CXX=${CMAKE_CXX_COMPILER})
list(APPEND GF2X_CONFIGURE_FLAGS CFLAGS=${GF2X_C_FLAGS})
list(APPEND GF2X_CONFIGURE_FLAGS CXXFLAGS=${GF2X_CXX_FLAGS})
list(APPEND GF2X_CONFIGURE_FLAGS LDFLAGS=${GF2X_LDFLAGS})
list(APPEND GF2X_CONFIGURE_FLAGS $ENV{GF2X_CONFIGURE_EXTRA_FLAGS})

get_override(GF2X_CONFIGURE_FLAGS GF2X_CONFIGURE_FLAGS)


######################################################################
# Try to detect whether pthreads and/or MPI are available.
# As a result:
# - pthreads are mandatory, but we also check the avaibility of the
#   function pthread_barrier_wait, which is optional.
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

# Included after pthreads.cmake as the tests need pthread to detect newer
# versions of gmp-ecm (the gmpecm.cmake script relies on ${pthread_libs} begin
# set)
include("${CADO_NFS_SOURCE_DIR}/config/gmpecm.cmake")

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

if (NOT "$ENV{IWYU}" MATCHES "^(0|NO|no|OFF|off|)$")
    message(STATUS "Enabling iwyu when compiling c++ code")
    set(findprog_flags)
    if(ENV{IWYU} MATCHES "/")
        # If MPI contains a /, then we assume it should be a path
        list(APPEND findprog_flags HINTS "$ENV{IWYU}" "$ENV{IWYU}/bin")
    endif()
    find_program(iwyu_path NAMES include-what-you-use iwyu ${findprog_flags})
    if(iwyu_path)
        message(STATUS "Enabling iwyu when compiling c++ code -- Success")
        set(CMAKE_C_INCLUDE_WHAT_YOU_USE "${iwyu_path}")
        set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${iwyu_path}")
        if(NOT "$ENV{IWYU_FILTER}" MATCHES "^(0|NO|no|OFF|off|)$")
            configure_file(config/iwyu.pl.in "${CADO_NFS_BINARY_DIR}/iwyu.pl")
            set(CMAKE_C_INCLUDE_WHAT_YOU_USE "${CADO_NFS_BINARY_DIR}/iwyu.pl")
            set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${CADO_NFS_BINARY_DIR}/iwyu.pl")
        endif()
        # oddly enough, the following does not work. (with or without a
        # CMAKE_ prefix)
        # set_property(GLOBAL PROPERTY CMAKE_C_INCLUDE_WHAT_YOU_USE ${iwyu_path})
        # set_property(GLOBAL PROPERTY CMAKE_CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
    else()
        message(WARNING "Enabling iwyu when compiling c++ code -- Failed, no binary found")
    endif()
endif()

message(STATUS "CFLAGS=${CMAKE_C_FLAGS} (synthetized, after feature checks)")
message(STATUS "CXXFLAGS=${CMAKE_CXX_FLAGS} (synthetized, after feature checks)")

#####################################################################
# now the real stuff.

# First, create the config file gathering everything we've got so far.
configure_file(cado_config_h.in cado_config.h ESCAPE_QUOTES @ONLY)

include_directories (${CADO_NFS_SOURCE_DIR})
include_directories (${CADO_NFS_SOURCE_DIR}/utils)
if(NOT HAVE_FMT)
include_directories (${CADO_NFS_SOURCE_DIR}/utils/embedded)
endif()
include_directories (${CADO_NFS_BINARY_DIR})

# link_directories (${CADO_NFS_BINARY_DIR}/utils)
# link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)

if(NOT HAVE_GF2X)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x)

add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/gf2x
    COMMAND echo Configuring gf2x with options ${GF2X_CONFIGURE_FLAGS}
    COMMAND cd gf2x && "${CADO_NFS_SOURCE_DIR}/gf2x/configure" ${GF2X_CONFIGURE_FLAGS} < /dev/null
    COMMENT "Configuring gf2x"
)

add_custom_target(gf2x-build ALL
    DEPENDS ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    WORKING_DIRECTORY gf2x
    COMMAND ${CMAKE_MAKE_PROGRAM}
    COMMENT "Building gf2x"
)

install(CODE "execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} install WORKING_DIRECTORY gf2x)")
endif()


add_subdirectory (utils)
add_subdirectory (polyselect)
add_subdirectory (sieve)
add_subdirectory (sqrt)
add_subdirectory (linalg)
add_subdirectory (filter)
add_subdirectory (linalg/bwc)
add_subdirectory (misc)
add_subdirectory (scripts)
add_subdirectory (tests EXCLUDE_FROM_ALL)

add_custom_target(dist
    COMMAND scripts/dist.sh ${package_basename}
    WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
)

# The trailing / is very important here. Without it, we would create a
# parameters/ directory in ${PREFIX}/share/cado-nfs-X.Y/
install(DIRECTORY parameters/ DESTINATION ${DATASUFFIX})

configure_file(cado-nfs.py cado-nfs.py ESCAPE_QUOTES @ONLY)

# file(CONFIGURE  is only for cmake 3.18 or later :-( ; this should work
# just the same
configure_file(config/source-location.txt.in source-location.txt ESCAPE_QUOTES @ONLY)
# file(CONFIGURE OUTPUT source-location.txt CONTENT "@CADO_NFS_SOURCE_DIR@" @ONLY)
configure_file(cado-nfs-client.py cado-nfs-client.py ESCAPE_QUOTES @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cado-nfs.py DESTINATION ${BINSUFFIX})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cado-nfs-client.py DESTINATION ${BINSUFFIX})


message(STATUS "Main scripts will go to ${CMAKE_INSTALL_PREFIX}/${BINSUFFIX}")
message(STATUS "Libraries and binaries will go to ${CMAKE_INSTALL_PREFIX}/${LIBSUFFIX}")
message(STATUS "Data will go to ${CMAKE_INSTALL_PREFIX}/${DATASUFFIX}")

IF (DEFINED ENV{CMAKE_DUMP_VARIABLES})
  include("${CADO_NFS_SOURCE_DIR}/config/variables.cmake")
ENDIF()

if(MINGW)
    message(STATUS "Please look up dev_docs/howto-MinGW.txt in the git tree for our mingw testing configuration, and frequently encountered errors")
endif()
