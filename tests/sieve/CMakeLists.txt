# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions

# Build makefb and las before running check
include_directories(.)
include_directories(${PROJECT_SOURCE_DIR}/utils)
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
        maxbits=10
        poly=${F9_poly}
        fb1=${F9_roots}
        fbc=.
        LAS_BINARY=${las_path}
        lim0=2300000
        lim1=1200000
)

set(c90_21505_poly ${CMAKE_CURRENT_SOURCE_DIR}/c90_21505.poly)
set(c90_21505_roots ${CMAKE_CURRENT_BINARY_DIR}/c90_21505.roots.gz)

cado_define_test(c90_21505_makefb
        PROGRAM
            ${makefb_path}
            -poly ${c90_21505_poly} -lim 811066 -maxbits 11
            -out ${c90_21505_roots}
        TARGET_DEPENDENCIES makefb
        FIXTURES_SETUP c90_21505.fb
        )

cado_define_test(c90_21505_sieve
        ENVIRONMENT
            las_I=11 poly=${c90_21505_poly} 
            lim0=404327 lim1=811066 lpb0=23 lpb1=23 mfb0=46 mfb1=46 ncurves0=9
            ncurves1=10 fb1=${c90_21505_roots} sqside=1
            q0=375923 rho=119201
            CHECKSUM_FILE=
            LAS_BINARY=${las_path}
            "REGEX=^-18320,3:.*74f"
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v -v --adjust-strategy 0
        TARGET_DEPENDENCIES las
        FIXTURES_REQUIRED c90_21505.fb
        )


set(c90_21553_poly ${CMAKE_CURRENT_SOURCE_DIR}/c90_21553.poly)
set(c90_21553_roots ${CMAKE_CURRENT_BINARY_DIR}/c90_21553.roots.gz)

cado_define_test(c90_21553_makefb
        PROGRAM
            ${makefb_path}
            -poly ${c90_21553_poly} -lim 811066 -maxbits 11
            -out ${c90_21553_roots}
        TARGET_DEPENDENCIES makefb
        FIXTURES_SETUP c90_21553.fb
        )

cado_define_test(c90_21553_sieve
        ENVIRONMENT
            las_I=11 poly=${c90_21553_poly} 
            lim0=404327 lim1=811066 lpb0=23 lpb1=23 mfb0=46 mfb1=46 ncurves0=9
            ncurves1=10 fb1=${c90_21553_roots} sqside=1
            q0=327127 rho=22938
            CHECKSUM_FILE=
            LAS_BINARY=${las_path}
            "REGEX=^22938,1:.*,b,.*"
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v -v --adjust-strategy 0
        TARGET_DEPENDENCIES las
        FIXTURES_REQUIRED c90_21553.fb
        )

set(c95_21684_poly ${CMAKE_CURRENT_SOURCE_DIR}/c95_21684.poly)
# the only needed thing is test_init_norms at the end of the file. The
# excerpt commented out below gives material to do a more thorough test.

# set(c95_21684_roots0 ${CMAKE_CURRENT_BINARY_DIR}/c95_21684.roots0.gz)
# set(c95_21684_roots1 ${CMAKE_CURRENT_BINARY_DIR}/c95_21684.roots1.gz)
# 
# cado_define_test(c95_21684_makefb0
#         PROGRAM
#             ${makefb_path}
#             -poly ${c95_21684_poly}
#             -lim 600000
#             -maxbits 13
#             -out ${c95_21684_roots0}
#             -side 0
#         TARGET_DEPENDENCIES makefb
#         )
# cado_define_test(c95_21684_makefb1
#         PROGRAM
#             ${makefb_path}
#             -poly ${c95_21684_poly}
#             -lim 850000
#             -maxbits 13
#             -out ${c95_21684_roots1}
#             -side 1
#         TARGET_DEPENDENCIES makefb
#         )
# 
# cado_define_test(c95_21684_sieve
#         ENVIRONMENT
#             las_I=13 poly=${c95_21684_poly} 
#             lim0=600000 lim1=850000 lpb0=24 lpb1=25 mfb0=48 mfb1=50 ncurves0=11
#             ncurves1=14 
#             fb0=${c95_21684_roots0}
#             fb1=${c95_21684_roots1}
#             sqside=1
#             q0=310049 rho=270571
#             CHECKSUM_FILE=
#             LAS_BINARY=${las_path}
#             "REGEX=reports"
#         PROVIDE_TEMPORARY_WDIR
#         SCRIPT
#             ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
#             -v -v --adjust-strategy 0
#         TARGET_DEPENDENCIES las
#         TEST_DEPENDENCIES c95_21684_makefb0 c95_21684_makefb1
#         )
# 

# Before running F9_sievetest, remove the file with the captured checksums
# We have to make this a "test" so that ctest runs it, even though it does
# not actually test anything
cado_define_test(F9_sievetest_cleanup
        PROGRAM ${CMAKE_COMMAND} -E remove -f
        ${F9_checksums} ${F9_checksums}_quick
        FIXTURES_SETUP F9.sievetest_cleanup
    )

cado_define_test(F9_makefb
        PROGRAM
            ${makefb_path}
            -poly ${F9_poly} -lim 1200000 -maxbits 10
            -out ${F9_roots}
        TARGET_DEPENDENCIES makefb
        FIXTURES_SETUP F9.fb
    )

# Most F9 tests easily take between 5 and 10 seconds on a laptop,
# and this can grow to several tens of seconds on a container that runs
# on a busy machine.
cado_define_test(F9
        ENVIRONMENT ${F9_baseargs}
            las_I=12
            lpb0=26 mfb0=52 lambda0=2.1
            lpb1=26 mfb1=52 lambda1=2.2
            q0=1200000 q1=1200200
            CHECKSUM_FILE=${F9_checksums}
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v --adjust-strategy 0
        TARGET_DEPENDENCIES las
        FIXTURES_REQUIRED F9.fb F9.sievetest_cleanup
        NO_DEFAULT_RUN
        TIMEOUT 90
        )

cado_divert_test(F9 sieve
        APPEND_ENVIRONMENT
            REL_COUNT=432
            REFERENCE_SHA1=ffe3f2036e7b80b2a7a56bc397f665506430f386
            REFERENCE_REVISION=2c5ac2cff
        )

cado_define_test(bogus_fbc_is_consistent
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/sieve/inspect-fbc-file.pl -fbc
    ${CMAKE_CURRENT_SOURCE_DIR}/bogus${OS_PAGE_SIZE}.fbc)

add_custom_target(copy_bogus_fbc ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/bogus${OS_PAGE_SIZE}.fbc
    ${CMAKE_CURRENT_BINARY_DIR}/bogus-copy.fbc)

cado_divert_test(F9_sieve reuse_foreign_fbc
        APPEND_ENVIRONMENT
        fbc=${CMAKE_CURRENT_BINARY_DIR}/bogus-copy.fbc
        TEST_DEPENDENCIES
        bogus_fbc_is_consistent
        TARGET_DEPENDENCIES
        copy_bogus_fbc
)

# override all arguments for this one.
cado_divert_test(F9 sieve_with_hints
        ENVIRONMENT
            ${F9_baseargs}
            todo=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.with_hints.todo
            hint_table=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.with_hints.hints
            REL_COUNT=153
            REFERENCE_SHA1=e80b973209f78c7df4541956e15814e04d699932
            REFERENCE_REVISION=21aa5e0327129475c37784b849cddc303d05f526
        )


cado_divert_test(F9 sieve_recompute_ssdpos
    APPEND_ARGUMENTS -B 12 -nq 1
    APPEND_ENVIRONMENT
        REL_COUNT=54
        CHECKSUM_FILE=
)

# go wild for the under-allocation. If this is valgrind-ok, then
# surely #21513 is fixed. While we're at it, check all kind of possible
# under-allocation recovers (the base test below serves only as a stem).
#
# Note that among the things we want to guard against, we have leaks:
# exception-unsafe code during the exception windup definitely causes
# leaks. So if those tests complain with valgrind, there's probably some
# of those.
#
# Note that bug #21515 had a specific expensive test, which is now covered
# by the _1long test below.
cado_divert_test(F9 sieve_check_overfull
        NO_DEFAULT_RUN
        APPEND_ENVIRONMENT
            q0=1200109 rho=116642
            bkthresh1=100000
            REL_COUNT=57
            REFERENCE_SHA1=2b973ed88e146aee756d9a27e260157179defc4d
            REFERENCE_REVISION=f06a6284d2cdcb669e300cd0dda349b349761578
            CHECKSUM_FILE=
        )
cado_divert_test(F9_sieve_check_overfull 1short
    APPEND_ENVIRONMENT bkmult=1.1,1s:0.5)
cado_divert_test(F9_sieve_check_overfull 1short_twosteps
    APPEND_ARGUMENTS -B 10
    APPEND_ENVIRONMENT bkmult=1.1,1s:0.5)
cado_divert_test(F9_sieve_check_overfull 1long
    APPEND_ENVIRONMENT bkmult=1.1,1l:0.5)
# working with the _2short overruns is a bit subtle. Because if we go
# wild with these, chances are that we will be *VERY* likely to write
# past the blank zone of 1MB, and then right into the red zone.
# Especially because the case below has only one level-2 bucket anyway
# (it is so contrived). So we are quite likely to get a plain, unpolite,
# Segmentation Fault. The check below has been tuned so that the overrun
# mechanism is indeed triggered, and so should run successfully.
cado_divert_test(F9_sieve_check_overfull 2short
    APPEND_ENVIRONMENT bkmult=1.1,2s:0.48)

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
            batchfile0=. batchfile1=.
            REL_COUNT=272
            REFERENCE_SHA1=3aba8ade73e9e47dd76b0f85a2c681408f34a189
            REFERENCE_REVISION=5d9b68511102a957f941b272eb6d753526bcc3b1
        APPEND_ARGUMENTS -batch
        TIMEOUT 180
        )
cado_define_test(bogus_batch_is_consistent
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/sieve/inspect-batch-file.pl -batch
    ${CMAKE_CURRENT_SOURCE_DIR}/bogus.batch0)
cado_divert_test(F9_batch try_to_reuse_foreign
    APPEND_ENVIRONMENT
    batchfile0=${CMAKE_CURRENT_SOURCE_DIR}/bogus.batch0)
set_tests_properties(F9_batch_try_to_reuse_foreign PROPERTIES WILL_FAIL 1)


cado_divert_test(F9 sublat
        APPEND_ENVIRONMENT
            las_I=10
            REL_COUNT=246
            REFERENCE_SHA1=f35a27aa4d96918a087cf8169071985c3992e659
            REFERENCE_REVISION=2c5ac2cff
        APPEND_ARGUMENTS -sublat 3
        TIMEOUT 150
        )

cado_divert_test(F9 skipped
        APPEND_ENVIRONMENT
            las_I=10
            REL_COUNT=60
        APPEND_ARGUMENTS -skipped 100
        )

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
        TEST_DEPENDENCIES F9_sieve
        FIXTURES_REQUIRED F9.fb
        )

# Try a too-large lambda value. We could divert on any test we like, but
# it's mildly quicker to do it on the short test.
cado_divert_test(F9_quick lambda0
        APPEND_ENVIRONMENT lambda0=5 "REGEX=lambda>.*does not make sense"
        TIMEOUT 180
    )
cado_divert_test(F9_quick lambda1
        APPEND_ENVIRONMENT lambda1=5 "REGEX=lambda>.*does not make sense"
        TIMEOUT 180
    )

set(c120_poly ${CMAKE_CURRENT_SOURCE_DIR}/c120.poly)
set(c120_fb ${CMAKE_CURRENT_SOURCE_DIR}/c120.roots1.gz)
set(c120_base_env
        LAS_BINARY=${las_path}
        lim0=10000 lpb0=26 mfb0=28 powlim0=8000
        lim1=10000 lpb1=26 mfb1=28 powlim1=8000
        las_I=12 q0=15350029 rho=14463285
        poly=${c120_poly}
        fb1=${c120_fb}
        t=2)

# Try sieving I (see bug 18814)
cado_define_test(sievetest_I
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            env ${c120_base_env}
            REFERENCE_SHA1=da39a3ee5e6b4b0d3255bfef95601890afd80709
            REFERENCE_REVISION=88ddb4f779defe692aed14a79ef8c952a4f68efa
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
        TARGET_DEPENDENCIES las
        )

# Try bucket-sieving of projective roots
set(c30_30012_base_env
        las_I=9
        lim0=30000
        lim1=30000
        lpb0=17
        lpb1=17
        mfb0=18
        mfb1=18
        poly=${CMAKE_CURRENT_SOURCE_DIR}/c30_30012.poly
        fb1=${CMAKE_CURRENT_SOURCE_DIR}/c30_30012.roots1.gz
        t=auto
        adjust_strategy=2
)

cado_define_test(c30_30012
    ENVIRONMENT
        ${c30_30012_base_env}
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
    NO_DEFAULT_RUN
    )

# This one doesn't test much. It only verifies that a projective root
# gets bucket sieved. If this fails, it probably means that this
# q-lattice ends up behaving differently, and does not have I=256 as per
# adjust-strategy 2.
cado_divert_test(c30_30012 1
    APPEND_ENVIRONMENT
        LAS_BINARY=${las_tracek_path}
        q0=52919 rho=35788 traceij=1,0 "TRACE_REGEX=Pushed.hit.*p=347"
    TARGET_DEPENDENCIES las_tracek
)

# This one is the real business about issue 30012
# Normally, we should sieve the given prime within
# fill_in_buckets_toplevel
cado_divert_test(c30_30012 2a
    APPEND_ENVIRONMENT
        LAS_BINARY=${las_tracek_path}
        q0=52919 rho=35788 traceij=0,347 "TRACE_REGEX=Pushed.hit.*p=347"
    TARGET_DEPENDENCIES las_tracek
)
cado_divert_test(c30_30012 2b
    APPEND_ENVIRONMENT
        LAS_BINARY=${las_tracek_path}
        q0=52919 rho=35788 traceij=1,694 "TRACE_REGEX=Pushed.hit.*p=347"
    TARGET_DEPENDENCIES las_tracek
)

# This is for the 2-level sieving case.
# Normally, we should sieve the given prime within
# fill_in_buckets_lowlevel (this works because we've restricted to very
# small sieving regions with B=10
cado_divert_test(c30_30012 3
    APPEND_ENVIRONMENT
        LAS_BINARY=${las_tracek_path}
        q0=50707 rho=50345
        las_A=23 las_I= las_B=10 bkthresh1=2048     # kill I explicitly
        traceij=1,1153 "TRACE_REGEX=Pushed.hit.*p=1153"
    TARGET_DEPENDENCIES las_tracek
)

# Try the todo-list mode of las
cado_divert_test(F9 todo
        APPEND_ENVIRONMENT
            REFERENCE_SHA1=39f6c2f198f1bbd7a65558c28f8ed696ebe2d1f1
            REFERENCE_REVISION=2c5ac2cff
            REL_COUNT=426
            CHECKSUM_FILE=
            todo=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.todo
        )

# The descent also goes with the typical sievetest.sh
cado_divert_test(F9 descent
        APPEND_ENVIRONMENT
            q0=235086167414177 rho=46506392392611
            REFERENCE_SHA1=4d64b69f37672b2109441c52f0a57f16f08d813d
            REFERENCE_REVISION=5d9b68511102a957f941b272eb6d753526bcc3b1
            descent_hint_table=${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.hint
            ncurves0=10 ncurves1=10
            CHECKSUM_FILE=
            LAS_BINARY=${las_descent_path}
        APPEND_ARGUMENTS
            -allow-largesq
        TARGET_DEPENDENCIES las_descent
        FIXTURES_REQUIRED F9.fb
        )

# The underlying script changes, here. We have to do something slightly
# different.
cado_define_test(F9_dupsup
        ENVIRONMENT ${F9_baseargs}
            las_I=12
            lpb0=28 mfb0=56 lambda0=2.1
            lpb1=28 mfb1=56 lambda1=2.1
            q0=250000000 nq=1
            ncurves0=10 ncurves1=10
            CHECKSUM_FILE=
            DUPSUP_BINARY=${dupsup_path}
            DUPE_COUNT=13
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/F9.dupsup.sh
            -v -v --adjust-strategy 0
        TARGET_DEPENDENCIES las dupsup
        FIXTURES_REQUIRED F9.fb F9.sievetest_cleanup
        TIMEOUT 60
        )

# Try fake relation generation
cado_define_test(F9_fakereltest
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/F9.fakerels.sh
    ${F9_roots} ${las_path} ${CADO_NFS_SOURCE_DIR}
    TARGET_DEPENDENCIES fake_rels freerel las
    FIXTURES_REQUIRED F9.fb
    TIMEOUT 120
)
# This really taking super long with valgrind in particular, so let's
# bump the timeout. 60 seconds is generally ok for normal runs, but 60
# times 8 for valgrind isn't.

# Check a bug in fake_rels (see 9d52c08c)
cado_define_test(c60_fakereltest
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/c60.fakerels.sh
    TARGET_DEPENDENCIES fake_rels freerel
    TIMEOUT 60
)

# Try with default TRIALDIV_MAXLEN, whatever it may be
# tests in debug mode tend to take very long when used with GLIBCXX_DEBUG
cado_define_test(test_trialdiv.cpp
    $<TARGET_OBJECTS:trialdiv>
    LIBRARIES utils tests
    TIMEOUT 30
)

macro (add_trialdiv_test MAXLEN)
    cado_define_test(test_trialdiv_${MAXLEN} SOURCES test_trialdiv.cpp
        ${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.cpp LIBRARIES utils tests)
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

cado_define_test(test_trialdiv_stdin
    SCRIPT
    ${CMAKE_CURRENT_BINARY_DIR}/test_trialdiv_5 -seed 1 -i 1000
    STDIN ${CMAKE_CURRENT_SOURCE_DIR}/test_trialdiv_example_stdin.txt
    EXPECT_SHA1 62ce4b010a4ecdaa4bacf977687748b5f82dc4c4
    TARGET_DEPENDENCIES
    test_trialdiv_5
    )

cado_define_test(test_logapprox.cpp 
    ${PROJECT_SOURCE_DIR}/sieve/logapprox.cpp
    LIBRARIES utils tests)

# This one is a different test, meant to check speed and correctness, and
# yet allow a run to simulate real polynomials because well, I can't be
# sure with the current code... */
cado_define_test(test_init_norms.cpp
    ${PROJECT_SOURCE_DIR}/sieve/las-multiobj-globals.cpp
    $<TARGET_OBJECTS:las_core_b0>
    $<TARGET_OBJECTS:las_core_a>
    $<TARGET_OBJECTS:las_core_0>
    NO_DEFAULT_RUN
    LIBRARIES utils tests
)

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
    AVOID_CONCURRENT 2
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
    AVOID_CONCURRENT 2
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
    AVOID_CONCURRENT 2
)

cado_divert_test(test_init_norms bug21518
    ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/c160_21518.poly
    -I 14 -lpb0 30 -lpb1 30 -mfb0 60 -mfb1 90
    -q0 20552291 -rho 9632932
        -sqside 1
        --hush-max-jitter 1
        --abort-on-jitter 1,5
        -check-bucket 361
)

cado_divert_test(test_init_norms bug21684a
    ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/c95_21684.poly
    -I 13 -lpb0 24 -lpb1 25 -mfb0 48 -mfb1 50
    -q0 310049 -rho 270571
    -sqside 1
    --hush-max-jitter 1
    # 21684 has two algebraic sides, so it's fine if we have some extra
    # jitter.
    --abort-on-jitter 2,5
    -check-bucket 505
)

cado_divert_test(test_init_norms bug21684b
    ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/c95_21684.poly
    -I 13 -lpb0 24 -lpb1 25 -mfb0 48 -mfb1 50
    -q0 310087 -rho 178382
    -sqside 1
    --hush-max-jitter 1
    # 21684 has two algebraic sides, so it's fine if we have some extra
    # jitter.
    --abort-on-jitter 2,5
    -check-bucket 4
)

cado_divert_test(test_init_norms bug30107
    ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/c97_30107.poly
    -A 28
    -sqside 1 -q0 2031811 -rho 1838948
    -lim1 2097152 -lpb1 21 -mfb1 45 -powlim 80
    -lim0 1024 -lpb0 10 -mfb0 10
    -norm-sides 1
    --check-bucket 127
    --abort-on-jitter 0,5 -B 16
)

# we have no pretense to have optimized code for these cases, so don't
# bother.
if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
cado_define_test(test-smallsieve.cpp 
                ${PROJECT_SOURCE_DIR}/sieve/las-where-am-i-prod.cpp
                LIBRARIES utils ${gmplib}
    ARGUMENTS -q -C -A 18 -I 14 -B 14)
cado_divert_test(test-smallsieve line_fragments1 # bucket region = half line
    ARGUMENTS -q -C -A 18 -I 14 -B 13)
cado_divert_test(test-smallsieve line_fragments2 # bucket region = quarter line
    ARGUMENTS -q -C -A 18 -I 15 -B 13)
endif()

if(HAVE_JEVENTS)
  set(jevents_lib "-ljevents")
else()
  set(jevents_lib "")
endif()
cado_define_test(test-sieve2357.cpp
                 ${PROJECT_SOURCE_DIR}/sieve/las-sieve2357.cpp
                 ${PROJECT_SOURCE_DIR}/sieve/las-where-am-i-prod.cpp
                 LIBRARIES utils tests ${jevents_lib})

cado_define_test(torture-redc.cpp LIBRARIES utils)
set_property(TARGET torture-redc
             PROPERTY COMPILE_DEFINITIONS ALIGN_LOOP_32=1)

cado_define_test(torture-redc+asm_cmov SOURCES torture-redc.cpp
    LIBRARIES utils tests
    NO_DEFAULT_RUN)
set_property(TARGET torture-redc+asm_cmov
             PROPERTY COMPILE_DEFINITIONS ALIGN_LOOP_32=1
                      LAS_ARITH_REDC_USE_ASM=1)

cado_define_test(test_fb_root.cpp
    $<TARGET_OBJECTS:las_core_a>
    LIBRARIES utils tests)

cado_define_test(test_batchinv.cpp LIBRARIES utils tests)

cado_define_test(test_invmod_redc_32.cpp LIBRARIES utils tests)

cado_define_test(test-reduce-plattice.cpp LIBRARIES utils ARGUMENTS
    -I 256
    -B 100000
    -ntests 100000
    -T
    )
cado_divert_test(test-reduce-plattice seed APPEND_ARGUMENTS -seed 1)

cado_define_test(test-inline-polynomials-and-expressions
    PROGRAM
            ${makefb_path}
            -poly "inline-poly://skew=1/poly0=x^3+2/poly1=x-2^43/n=2^128+1"
            -lim 200
            -t 2
    AVOID_CONCURRENT 2
    TARGET_DEPENDENCIES
        makefb
)

# test for las with galois
set(galois_poly ${CMAKE_CURRENT_SOURCE_DIR}/galois.poly)
set(galois_roots ${CMAKE_CURRENT_SOURCE_DIR}/galois.roots.gz)

cado_define_test(galois_sieve_abnormal_orbits
        ENVIRONMENT
            las_I=10 poly=${galois_poly}
            lim0=200 lpb0=9 mfb0=10 fb0=${galois_roots}
            q0=101 q1=102 galois=_y
            lim1=200 lpb1=9 mfb1=10
            # Note: the last 3 args (lim1, lpb1 and mfb1) are not needed but
            # sievetest.sh complains if they are not present
            LAS_BINARY=${las_path}
        REL_COUNT=71
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v -v
        TARGET_DEPENDENCIES las
        )

# test for las with only one side
set(d20_poly ${CMAKE_CURRENT_SOURCE_DIR}/d20.poly)
set(d20_roots ${CMAKE_CURRENT_SOURCE_DIR}/d20.roots.gz)

cado_define_test(d20_one_sided_las_sieve
        ENVIRONMENT
            las_I=14 poly=${d20_poly}
            lim0=800 lpb0=10 mfb0=25 ncurves0=2 fb0=${d20_roots}
            q0=900 nq=5 galois=_y
            lim1=800 lpb1=10 mfb1=25
            # Note: the last 3 args (lim1, lpb1 and mfb1) are not needed but
            # sievetest.sh complains if they are not present
            LAS_BINARY=${las_path}
        REL_COUNT=29
        REFERENCE_SHA1=d30c636a8cf76377441e6f1573f053fd7aa900ee
        REFERENCE_REVISION=xxxxxxxxxxxxxxxxxxxxxxxxxx
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/sievetest.sh
            -v -v
        TARGET_DEPENDENCIES las
        )

cado_divert_test(d20_one_sided_las_sieve batch
        APPEND_ARGUMENTS
            -batch
        APPEND_ENVIRONMENT
            REL_COUNT=20
            REFERENCE_SHA1=0124eff1f091882573748c94f505e9f906b8c9f7
            REFERENCE_REVISION=xxxxxxxxxxxxxxxxxxxxxxxxxx
        )

cado_divert_test(d20_one_sided_las_sieve dup
        APPEND_ARGUMENTS
            -dup -dup-qmin 800
        APPEND_ENVIRONMENT
            REL_COUNT=27
            REFERENCE_SHA1=f9d0c3b58bdfcb51d598180cb01895b2c35931f7
            REFERENCE_REVISION=xxxxxxxxxxxxxxxxxxxxxxxxxx
        )

add_subdirectory (ecm)
add_subdirectory (strategies)
if(strategies_variant STREQUAL "strategies_old")
add_subdirectory (${strategies_variant})
endif()
