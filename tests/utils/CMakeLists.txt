cado_define_test(test_strlcat.c LIBRARIES utils)
cado_define_test(test_strlcpy.c LIBRARIES utils)
cado_define_test(test_double_poly.c LIBRARIES utils tests)
cado_define_test(test_gmp_aux.c LIBRARIES utils tests)
cado_define_test(test_gcd.c LIBRARIES utils tests)
cado_define_test(test_getprime.c LIBRARIES utils)
cado_divert_test(test_getprime 1 APPEND_ARGUMENTS -seek 0 -bound 100 -count 25)
cado_divert_test(test_getprime 2 APPEND_ARGUMENTS -seek 2 -bound 100 -count 25)
cado_divert_test(test_getprime 3 APPEND_ARGUMENTS -seek 3 -bound 100 -count 24)
cado_divert_test(test_getprime 4 APPEND_ARGUMENTS -seek 50 -bound 100 -count 10)
cado_divert_test(test_getprime 5 APPEND_ARGUMENTS -bound 10000000 -count 664579)
cado_divert_test(test_getprime 6 APPEND_ARGUMENTS -seek 5000000 -bound 10000000 -count 316066)
cado_define_test(test_memusage.c LIBRARIES utils)
cado_define_test(test_roots_mod.c LIBRARIES utils ARGUMENTS 100 1050)
cado_define_test(test_omega.c LIBRARIES utils)
cado_define_test(test_timing.c LIBRARIES utils)
cado_define_test(test_rootfinder.c LIBRARIES utils tests)
cado_define_test(test_mod_ul.c LIBRARIES utils tests ARGUMENTS 1000)
cado_define_test(test_modredc_ul.c LIBRARIES utils tests ARGUMENTS 1000)
# These two seem to time out easily in coverage mode
cado_define_test(test_modredc_15ul.c LIBRARIES utils tests ARGUMENTS 1000
    TIMEOUT 60
)
cado_define_test(test_modredc_2ul2.c LIBRARIES utils tests ARGUMENTS 1000
    TIMEOUT 120
)
cado_define_test(test_sm_utils.cpp LIBRARIES utils ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/test_sm_utils.data)
cado_define_test(test_mpz_poly.cpp LIBRARIES utils tests AVOID_CONCURRENT 4 TIMEOUT 60)
cado_define_test(test_mpz_poly_parsing.cpp LIBRARIES utils tests)
cado_define_test(test_json_parsing.cpp LIBRARIES utils tests)
cado_define_test(test_cxx_mpz_mat.cpp LIBRARIES utils tests)
cado_define_test(test_cxx_mpz_poly.cpp LIBRARIES utils tests)
cado_define_test(test_cado_poly.cpp LIBRARIES utils tests)
cado_divert_test(test_cado_poly mnfs5 APPEND_ARGUMENTS
    ${CADO_NFS_SOURCE_DIR}/tests/misc/test_renumber.data/mnfs5.poly)
cado_divert_test(test_cado_poly c105_30003 APPEND_ARGUMENTS
    ${CADO_NFS_SOURCE_DIR}/tests/misc/test_renumber.data/c105_30003.poly)
cado_divert_test(test_cado_poly buggy EXPECT_FAIL APPEND_ARGUMENTS
    ${CADO_NFS_SOURCE_DIR}/tests/utils/buggy.poly)
cado_define_test(test_relation.cpp LIBRARIES utils tests)
cado_define_test(test_cache.c LIBRARIES utils)
cado_define_test(test_usp.c LIBRARIES utils tests)
cado_define_test(test_bit_vector.c LIBRARIES utils tests)
cado_define_test(test_intinv.c LIBRARIES utils tests)
cado_define_test(test_modul_poly.c LIBRARIES utils tests)
cado_define_test(test_modredcul_batchinv.c LIBRARIES utils tests)
cado_define_test(test_modredcul_batchinv_ul.c LIBRARIES utils tests)
cado_define_test(test_dllist.c LIBRARIES tests)
cado_define_test(test_threadpool.cpp LIBRARIES utils tests)
# note: these tests used to be disabled...
cado_define_test(test_verbose.c LIBRARIES utils tests) # Should be run under Helgrind or other thread checkers
cado_define_test(test_modredc2ul2_Q_to_Fp.c LIBRARIES utils tests)
cado_define_test(test_mpz_vector.c LIBRARIES utils)
cado_define_test(test_mpz_poly_bivariate.cpp LIBRARIES utils tests)
cado_define_test(test_smallset.cpp LIBRARIES)
cado_define_test(test_mpz_mat.c LIBRARIES utils)
cado_define_test(test_gpf.c LIBRARIES utils)
cado_define_test(test-sort.cpp ARGUMENTS --log-bigsize 14 --max-span 128 LIBRARIES ${gmp_libname})
cado_define_test(test_u64arith.cpp LIBRARIES utils tests)
cado_define_test(test_mod.cpp LIBRARIES utils tests)
cado_define_test(test_modint.cpp LIBRARIES utils tests)
cado_define_test(test_istream_matcher.cpp LIBRARIES utils tests)
cado_define_test(test_mmappable_vector.cpp LIBRARIES utils
    PRECOMMAND "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --arg --tmpdir
    TIMEOUT 60
    )
cado_define_test(test_sha1.c LIBRARIES utils ARGUMENTS --test)
cado_define_test(test_fstream_maybe_compressed.cpp LIBRARIES utils
    NO_DEFAULT_RUN
    PRECOMMAND "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --arg --wdir)
cado_divert_test(test_fstream_maybe_compressed gz ARGUMENTS test.gz)
cado_divert_test(test_fstream_maybe_compressed bz2 ARGUMENTS test.bz2)
cado_divert_test(test_fstream_maybe_compressed none ARGUMENTS test)
cado_define_test(test_min_max_heap.cpp LIBRARIES tests)
cado_define_test(test_random_distributions.cpp LIBRARIES utils tests)

cado_define_test(
    SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/test-badideals.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.small-example
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.rsa220
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.p1024
    PRECOMMAND
        ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env wdir --
    TARGET_DEPENDENCIES
        numbertheory_tool
    )

cado_define_test(test_subdivide_primes_interval.cpp LIBRARIES utils)
cado_define_test(test_galois_action.cpp LIBRARIES utils tests)

add_subdirectory(numbertheory)

if(SAGE)
    cado_define_test(
        SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/test_explain_relations.sh
        PRECOMMAND "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env wdir --
        NO_DEFAULT_RUN
        TARGET_DEPENDENCIES
        makefb
        freerel
        las
        dup1
        dup2
        explain_indexed_relation
    )

    # note that the sage preparser takes an absurd amount of time here,
    # so it's important to split into very small data sets.
    cado_divert_test(test_explain_relations
        c30
        -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_explain_relations.data/c30.poly
        -I 9
        -q0 56230 -q1 56240
        -lim0 30000 -lpb0 17 -mfb0 18
        -lim1 30000 -lpb1 17 -mfb1 18
    )
    cado_divert_test(test_explain_relations_c30 dl APPEND_ARGUMENTS -dl)
endif()
