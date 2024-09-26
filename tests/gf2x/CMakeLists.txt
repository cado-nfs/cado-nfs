
# note that cmake sees gf2x as an external library, so adding gf2x with
# target_link_libraries cannot have the effect that the target depends on
# building gf2x. That dependency (on gf2x-build, in fact) has to be added
# by hand.

if(DEFINED ENV{CHECKS_EXPENSIVE})
cado_define_test(TEST_NAME test_gf2x
    PROGRAM ${CMAKE_MAKE_PROGRAM} check
    WORKING_DIRECTORY ${CADO_NFS_BINARY_DIR}/gf2x
    TARGET_DEPENDENCIES gf2x-build
    TIMEOUT 90)
endif()

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft/mpfq)

cado_define_test(test_cantor.c
    ARGUMENTS 100000
    LIBRARIES gf2x-fft gf2x tests ${gmp_libname}
    TARGET_DEPENDENCIES gf2x-build)
# This **MUST** match the #define in gf2x/Makefile.am
set_source_files_properties(test_cantor.c COMPILE_FLAGS
        "-DCANTOR_BASE_FIELD_SIZE=128")
