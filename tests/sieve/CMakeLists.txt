# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions

# Build makefb and las before running check
include_directories(${PROJECT_SOURCE_DIR}/sieve)

# Get paths to the makefb and las binaries and store in variables
# makefb_path and las_path, resp.
# for cmake pre-3.0, we did:
# get_property(makefb_path TARGET makefb PROPERTY LOCATION)

# now it seems to be forbidden with cmake 3.0+ (see cmake --help-policy
# CMP0026). The following is probably safe enough.
set(makefb_path         ${CADO_NFS_BINARY_DIR}/sieve/makefb)
set(las_path            ${CADO_NFS_BINARY_DIR}/sieve/las)
set(las_tracek_path     ${CADO_NFS_BINARY_DIR}/sieve/las_tracek)
set(las_descent_path    ${CADO_NFS_BINARY_DIR}/sieve/las_descent)
set(dupsup_path         ${CADO_NFS_BINARY_DIR}/sieve/dupsup)

set(F9_roots "${CMAKE_CURRENT_BINARY_DIR}/F9.roots")

cado_define_test(F9_makefbtest
        PROGRAM ${makefb_path} -poly
        ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/F9.poly -lim 1200000 -maxbits 10 -out ${F9_roots}
        TARGET_DEPENDENCIES makefb)

# Before running F9_sievetest, remove the file with the captured checksums
# We have to make this a "test" so that ctest runs it, even though it does
# not actually test anything
cado_define_test(F9_sievetest_cleanup 
        PROGRAM ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/F9.checksums)

cado_define_test(F9_sievetest
        SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh
        ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/F9.checksums
        TARGET_DEPENDENCIES las
        TEST_DEPENDENCIES F9_sievetest_cleanup F9_makefbtest
        )

cado_define_test(F9_batchtest
        SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/F9.batchtest.sh
        ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/F9.checksums
        TARGET_DEPENDENCIES las
        TEST_DEPENDENCIES F9_sievetest_cleanup F9_makefbtest
        )


# We make these depend on sievetest, for checking the checksum file.
# Since sievetest depends on makefbtest, we don't have to specify it
# again.
cado_divert_test(F9_sievetest unsieve1 APPEND_ARGUMENTS -unsievethresh 1)

cado_define_test(F9_cofactest
    SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/F9.cofactest.sh
        ${F9_sievetest_ARGUMENTS}
        -file-cofact ${CMAKE_CURRENT_SOURCE_DIR}/F9.st
    TARGET_DEPENDENCIES las
    TEST_DEPENDENCIES F9_sievetest)

# The TRACE_K test does only a single special-q, so we don't compare output
# relations or checksums
cado_define_test(F9_tracektest
    SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh
        ${F9_roots} ${las_tracek_path} ${CADO_NFS_SOURCE_DIR} -traceab 8517,584707
    TARGET_DEPENDENCIES las_tracek
    TEST_DEPENDENCIES F9_makefbtest)

# Try sieving I (see bug 18814)
cado_define_test(sievetest_I
    PROGRAM
        ${las_path} -poly ${CMAKE_CURRENT_SOURCE_DIR}/c120.poly -I 12
        -fb1 ${CMAKE_CURRENT_SOURCE_DIR}/c120.roots1.gz -lim0 10000
        -lim1 10000 -lpb0 26 -lpb1 26 -mfb0 28 -mfb1 28 -q0 15350028 -q1
        15350029 -t 2 -powlim0 8000 -powlim1 8000
    TARGET_DEPENDENCIES las)

# Try a too-large lambda value
cado_define_test(F9_lambda0test
    PROGRAM env lambda0=5 ${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    -regex "lambda>.* does not make sense"
    TARGET_DEPENDENCIES las
    TEST_DEPENDENCIES F9_makefbtest)
cado_define_test(F9_lambda1test
    PROGRAM env lambda1=5 ${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    -regex "lambda>.* does not make sense"
    TARGET_DEPENDENCIES las
    TEST_DEPENDENCIES F9_makefbtest)

# Try the todo-list mode of las
cado_define_test(F9_todotest
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.todotest.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    ${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.todo
    TARGET_DEPENDENCIES las
    TEST_DEPENDENCIES F9_makefbtest)

# Try the descent 
cado_define_test(F9_descenttest
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.descenttest.sh
    ${F9_roots} ${las_descent_path} ${CADO_NFS_SOURCE_DIR}
    ${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.hint 
    TARGET_DEPENDENCIES las_descent
    TEST_DEPENDENCIES F9_makefbtest)

# Try duplicate online suppression
cado_define_test(F9_dupsuptest
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.dupsup.sh
    ${F9_roots} ${las_path} ${dupsup_path} ${CADO_NFS_SOURCE_DIR}
    TARGET_DEPENDENCIES dupsup
    TEST_DEPENDENCIES F9_makefbtest)

# Try fake relation generation
cado_define_test(F9_fakereltest
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.fakerels.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    TARGET_DEPENDENCIES fake_rels
    TEST_DEPENDENCIES F9_makefbtest)

cado_define_test(test_memset_things.cpp LIBRARIES las-norms utils)

# Try with default TRIALDIV_MAXLEN, whatever it may be
cado_define_test(test_trialdiv.c ${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c LIBRARIES utils tests trialdiv)

macro (add_trialdiv_test MAXLEN)
    cado_define_test(test_trialdiv_${MAXLEN} SOURCES test_trialdiv.c ${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c LIBRARIES utils tests)
    set_property(TARGET test_trialdiv_${MAXLEN}
                PROPERTY COMPILE_DEFINITIONS TRIALDIV_MAXLEN=${MAXLEN})
endmacro()
# Try specific TRIALDIV_MAXLEN from 1 to 6
add_trialdiv_test(1)
add_trialdiv_test(2)
add_trialdiv_test(3)
add_trialdiv_test(4)
add_trialdiv_test(5)
add_trialdiv_test(6)

cado_define_test(test_init_norms_bucket_region.cpp
    LIBRARIES las-norms utils tests)

add_subdirectory (ecm)
add_subdirectory (strategies)
