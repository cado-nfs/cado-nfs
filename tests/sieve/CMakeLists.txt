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


set(F9_roots ${CMAKE_CURRENT_BINARY_DIR}/F9.roots)
set(F9_poly ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/F9.poly)
set(F9_checksums ${CMAKE_CURRENT_BINARY_DIR}/F9.checksums)
set(F9_st ${CMAKE_CURRENT_SOURCE_DIR}/F9.st)

set(F9_baseargs
        lim0=2300000 lpb0=26 mfb0=52 lambda0=2.1
        lim1=1200000 lpb1=26 mfb1=52 lambda1=2.2
        maxbits=10
        I=12 q0=1200000 q1=1200200
        poly=${F9_poly}
        fb=${F9_roots}
        fbc=.
        CHECKSUM_FILE=${F9_checksums}
        LAS_BINARY=${las_path}
)

# Before running F9_sievetest, remove the file with the captured checksums
# We have to make this a "test" so that ctest runs it, even though it does
# not actually test anything
cado_define_test(F9_sievetest_cleanup
        PROGRAM ${CMAKE_COMMAND} -E remove -f ${F9_checksums})

cado_define_test(F9_makefb
        PROGRAM
            ${makefb_path}
            -poly ${F9_poly} -lim 1200000 -maxbits 10
            -out ${F9_roots}
        TARGET_DEPENDENCIES makefb)

cado_define_test(F9
        ENVIRONMENT ${F9_baseargs}
        SCRIPT
            "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v -v --adjust-strategy 0
        TARGET_DEPENDENCIES las
        TEST_DEPENDENCIES F9_sievetest_cleanup F9_makefb
        NO_DEFAULT_RUN
        )

cado_divert_test(F9 sieve
        APPEND_ENVIRONMENT
            REL_COUNT=431
            REFERENCE_SHA1=e04ae591795ff30af703a2235a4f9fd09d17b380
            REFERENCE_REVISION=ede066464ba58b2d6127fa256cef01706282e0ed
        )

cado_divert_test(F9 sieve_bug_21480a
        APPEND_ENVIRONMENT
            q0=1200109 rho=116642 REL_COUNT=57 CHECKSUM_FILE=
            LAS_BINARY=${las_tracek_path}
        APPEND_ARGUMENTS
            -traceab -1098913,1678890
        TARGET_DEPENDENCIES
            las_tracek
        )

cado_divert_test(F9 sieve_bug_21480b
        APPEND_ENVIRONMENT
            q0=1200007 rho=554209 REL_COUNT=54 CHECKSUM_FILE=
            LAS_BINARY=${las_tracek_path}
        APPEND_ARGUMENTS
            -traceab 190203,834245
        TARGET_DEPENDENCIES
            las_tracek
        )

cado_divert_test(F9 batch
        APPEND_ENVIRONMENT
            lpb1=24 mfb1=72 lambda1=2.9 lambda0=2.2
            batch0=. batch1=.
            REL_COUNT=272
            REFERENCE_SHA1=3aba8ade73e9e47dd76b0f85a2c681408f34a189
            REFERENCE_REVISION=5d9b68511102a957f941b272eb6d753526bcc3b1
        APPEND_ARGUMENTS -batch
        )

if (0)
cado_divert_test(F9 sublat
        APPEND_ENVIRONMENT
            REFERENCE_SHA1=090c7a9172ba1e70a459df8b0d7e024c4289d093
            REFERENCE_REVISION=604a709c968e05ec09423dba1a977e100f038a90
        APPEND_ARGUMENTS -sublat 3
        )
endif()

# We make these depend on sievetest, for checking the checksum file.
# Since sievetest depends on makefb, we don't have to specify it
# again.
cado_divert_test(F9_sieve unsieve1 APPEND_ARGUMENTS -unsievethresh 1)
cado_divert_test(F9_sieve cofactest APPEND_ARGUMENTS -file-cofact ${F9_st})


# Do a diversion for a single special-q. This is taking the same data set
# as above, but with the real las binary.
cado_divert_test(F9 quick
        APPEND_ENVIRONMENT
            q0=1200007 rho=554209      # force (q0, rho) to be used.
            REFERENCE_SHA1=0ac3255855fe39fe729b771147a36dd72fc60d1b
            REFERENCE_REVISION=9253f54766c5438fe3198053f9edd811ee254a17
            CHECKSUM_FILE=${F9_checksums}_quick
        NO_DEFAULT_RUN)

# The TRACE_K test does only a single special-q.
# We redefine the test completely, since the binary changes.
cado_divert_test(F9_quick tracek
        APPEND_ENVIRONMENT
            LAS_BINARY=${las_tracek_path}
        # we need to replace the full argument list, because the binary
        # is not the same.
        APPEND_ARGUMENTS
            -traceab -8517,584707
        TARGET_DEPENDENCIES las_tracek
        TEST_DEPENDENCIES F9_makefb F9_sieve
        )

# Try a too-large lambda value. We could divert on any test we like, but
# it's mildly quicker to do it on the short test.
cado_divert_test(F9_quick lambda0
        APPEND_ENVIRONMENT lambda0=5 "REGEX=lambda>.*does not make sense")
cado_divert_test(F9_quick lambda1
        APPEND_ENVIRONMENT lambda1=5 "REGEX=lambda>.*does not make sense")


set(c120_poly ${CMAKE_CURRENT_SOURCE_DIR}/c120.poly)
set(c120_fb ${CMAKE_CURRENT_SOURCE_DIR}/c120.roots1.gz)
set(c120_baseargs
        LAS_BINARY=${las_path}
        lim0=10000 lpb0=26 mfb0=28 powlim0=8000
        lim1=10000 lpb1=26 mfb1=28 powlim1=8000
        I=12 q0=15350029 rho=14463285
        poly=${c120_poly}
        fb=${c120_fb}
        t=2)

# Try sieving I (see bug 18814)
cado_define_test(sievetest_I
        SCRIPT
            env ${c120_baseargs}
            REFERENCE_SHA1=da39a3ee5e6b4b0d3255bfef95601890afd80709
            REFERENCE_REVISION=88ddb4f779defe692aed14a79ef8c952a4f68efa
            "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
        TARGET_DEPENDENCIES las
        )

# Try the todo-list mode of las
cado_divert_test(F9 todo
        APPEND_ENVIRONMENT
            REFERENCE_SHA1=451406185eca205ee28990c7381f85a1c4942c5f
            REFERENCE_REVISION=88ddb4f779defe692aed14a79ef8c952a4f68efa
            CHECKSUM_FILE=${F9_checksums}_other
            todo=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.todo
        )

# The descent also goes with the typical sievetest.sh
cado_divert_test(F9 descent
        APPEND_ENVIRONMENT
            q0=235086167414177 rho=46506392392611
            REFERENCE_SHA1=4d64b69f37672b2109441c52f0a57f16f08d813d
            REFERENCE_REVISION=5d9b68511102a957f941b272eb6d753526bcc3b1
            descent_hint=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.hint
            ncurves0=10 ncurves1=10
            CHECKSUM_FILE=
            LAS_BINARY=${las_descent_path}
        APPEND_ARGUMENTS
            -allow-largesq
        TARGET_DEPENDENCIES las_descent
        TEST_DEPENDENCIES F9_makefb
        )

# The underlying script changes, here. We have to do something slightly
# different.
cado_define_test(F9_dupsup
        ENVIRONMENT ${F9_baseargs}
            lpb0=28 mfb0=56 lambda0=2.1
            lpb1=28 mfb1=56 lambda1=2.1
            q0=250000000 nq=1
            ncurves0=10 ncurves1=10
            CHECKSUM_FILE=
            DUPSUP_BINARY=${dupsup_path}
            DUPE_COUNT=12
        SCRIPT
            "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
            ${CMAKE_CURRENT_SOURCE_DIR}/F9.dupsup.sh
            -v -v --adjust-strategy 0
        TARGET_DEPENDENCIES las dupsup
        TEST_DEPENDENCIES F9_sievetest_cleanup F9_makefb F9_sieve
        )

# Try fake relation generation
cado_define_test(F9_fakereltest
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.fakerels.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    TARGET_DEPENDENCIES fake_rels freerel
    TEST_DEPENDENCIES F9_makefb)

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

cado_define_test(test_logapprox.cpp 
    ${PROJECT_SOURCE_DIR}/sieve/logapprox.cpp
    LIBRARIES utils tests)

# This one is a different test, meant to check speed and correctness, and
# yet allow a run to simulate real polynomials because well, I can't be
# sure with the current code... */
cado_define_test(test_init_norms.cpp
    ${PROJECT_SOURCE_DIR}/sieve/las-norms.cpp
    ${PROJECT_SOURCE_DIR}/sieve/logapprox.cpp
    ${PROJECT_SOURCE_DIR}/sieve/las-siever-config.cpp
    ${PROJECT_SOURCE_DIR}/sieve/las-todo-entry.cpp
    ${PROJECT_SOURCE_DIR}/sieve/las-qlattice.cpp
    NO_DEFAULT_RUN
    LIBRARIES utils tests
)
set_target_properties(test_init_norms PROPERTIES 
    COMPILE_FLAGS "${OpenMP_C_FLAGS}"
    LINK_FLAGS "${OpenMP_CXX_FLAGS}")


cado_divert_test(test_init_norms c120_basic
    ARGUMENTS
        -poly ${CMAKE_CURRENT_SOURCE_DIR}/c120.poly
        -I 12 -lpb0 26 -lpb1 26
        -mfb0 28 -mfb1 28
        -q0 15350000 -q1 15360000
        -random-sample 10
        -random-seed 0
        -sqside 1
        -t 2
        --hush-max-jitter 1
)

cado_divert_test(test_init_norms c120_torture
    ARGUMENTS
        -poly ${CMAKE_CURRENT_SOURCE_DIR}/c120.poly
        -I 12 -lpb0 26 -lpb1 26
        -mfb0 28 -mfb1 28
        -t 2
        --hush-max-jitter 1
        --abort-on-jitter 1,5
        --nfills-speed-test 0
    NO_DEFAULT_RUN
)

# a couple of hard cases encountered during testing.
cado_divert_test(test_init_norms_c120_torture 1 APPEND_ARGUMENTS
    --norm-sides 0 -sqside 1 -q0 15359089 -rho 4697045 -check-bucket 74)
cado_divert_test(test_init_norms_c120_torture 2 APPEND_ARGUMENTS
    --norm-sides 0 -sqside 1 -q0 15351793 -rho 6002387 -check-bucket 30)
cado_divert_test(test_init_norms_c120_torture 3 APPEND_ARGUMENTS
    --norm-sides 0 -sqside 1 -q0 15356527 -rho 2085523 -check-bucket 87)
cado_divert_test(test_init_norms_c120_torture 4 APPEND_ARGUMENTS
    --norm-sides 0 -sqside 1 -q0 15355589 -rho 4854502 -check-bucket 67)

cado_divert_test(test_init_norms rsa220
    ARGUMENTS
        -poly ${CMAKE_CURRENT_SOURCE_DIR}/rsa220.poly
        -I 16 -lpb0 34 -lpb1 34 -mfb0 68 -mfb1 102
        -q0 4000000000
        -q1 4001000000
        -random-sample 10
        -random-seed 0
        -sqside 1
        -t 2
        --hush-max-jitter 1
        --abort-on-jitter 1,5
)

add_subdirectory (ecm)
add_subdirectory (strategies)
