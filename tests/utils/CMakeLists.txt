cado_define_test(test_strlcat.c LIBRARIES utils)
cado_define_test(test_strlcpy.c LIBRARIES utils)
cado_define_test(test_double_poly.c LIBRARIES utils tests)
cado_define_test(test_polynomial.cpp LIBRARIES utils tests)
cado_define_test(test_number_literal.cpp LIBRARIES utils)
if(HAVE_MPFR)
    target_link_libraries(test_polynomial mpfr)
endif()
cado_define_test(test_cado_expression_parser.cpp LIBRARIES utils tests)
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
cado_define_test(test_trial_division.cpp LIBRARIES utils tests ARGUMENTS -q)
cado_define_test(test_rootfinder.c LIBRARIES utils tests)

add_subdirectory(arith)
add_subdirectory(arithxx)

cado_define_test(test_sm_utils.cpp LIBRARIES utils ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/test_sm_utils.data)
cado_divert_test(test_sm_utils failure EXPECT_FAIL ARGUMENTS ${CMAKE_CURRENT_SOURCE_DIR}/test_sm_utils_failure.data)
cado_define_test(test_mpz_poly.cpp
    LIBRARIES utils tests
    AVOID_CONCURRENT 4
    ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
    TIMEOUT 60)
cado_define_test(test_mpz_poly_parsing.cpp LIBRARIES utils tests)
cado_define_test(test_json_parsing.cpp LIBRARIES utils tests)
if(HAVE_MPFR)
cado_define_test(test_cxx_mpfr.cpp LIBRARIES utils tests)
endif()
if(HAVE_MPC)
cado_define_test(test_cxx_mpc.cpp LIBRARIES utils tests)
endif()
cado_define_test(test_cado_math_aux.cpp LIBRARIES utils tests)
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
cado_define_test(test_dllist.c LIBRARIES tests)
cado_define_test(test_threadpool.cpp LIBRARIES utils tests)
# note: these tests used to be disabled...
cado_define_test(test_verbose.c LIBRARIES utils tests) # Should be run under Helgrind or other thread checkers
cado_define_test(test_mpz_vector.c LIBRARIES utils)
cado_define_test(test_mpz_poly_bivariate.cpp ARGUMENTS -iter 100 LIBRARIES utils tests TIMEOUT 60)
cado_define_test(test_smallset.cpp LIBRARIES)
cado_define_test(test_mpz_mat.c LIBRARIES utils)
cado_define_test(test_gpf.c LIBRARIES utils)
cado_define_test(test-sort.cpp ARGUMENTS --log-bigsize 14 --max-span 128
    LIBRARIES utils)
cado_define_test(test_u64arith.cpp LIBRARIES utils tests)
cado_define_test(test_multityped_array.cpp LIBRARIES utils)
cado_define_test(test_mmappable_vector.cpp LIBRARIES utils
    PROVIDE_TEMPORARY_WDIR
    TIMEOUT 60
    )
cado_define_test(test_sha1.c LIBRARIES utils ARGUMENTS --test)
cado_define_test(test_fstream_maybe_compressed.cpp LIBRARIES utils
    PROVIDE_TEMPORARY_WDIR
    NO_DEFAULT_RUN)
cado_define_test(test_runtime_numeric_cast.cpp LIBRARIES utils)
macro(add_test_fstream_maybe_compressed var ext)
if(${var})
cado_divert_test(test_fstream_maybe_compressed ${ext} ARGUMENTS test.${ext})
endif()
endmacro()
add_test_fstream_maybe_compressed(HAVE_GZIP gz)
add_test_fstream_maybe_compressed(HAVE_BZIP2 bz2)
add_test_fstream_maybe_compressed(HAVE_ZSTD zstd)
add_test_fstream_maybe_compressed(HAVE_XZ xz)
add_test_fstream_maybe_compressed(HAVE_LZMA lzma)
cado_divert_test(test_fstream_maybe_compressed none ARGUMENTS test)
cado_define_test(test_min_max_heap.cpp LIBRARIES tests)
cado_define_test(test_random_distributions.cpp LIBRARIES utils tests)
cado_define_test(test_decomposed_path.cpp LIBRARIES utils tests)
cado_define_test(test_prime_count.cpp LIBRARIES utils tests)
cado_define_test(test_integer_partitions.cpp LIBRARIES utils)
cado_define_test(test_mpz_conversions.cpp LIBRARIES utils tests)

cado_define_test(
    SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/test-badideals.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.small-example
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.rsa220
        ${CMAKE_CURRENT_SOURCE_DIR}/numbertheory/badideals.p1024
    PROVIDE_TEMPORARY_WDIR
    TARGET_DEPENDENCIES
        numbertheory_tool
    )

cado_define_test(test_subdivide_primes_interval.cpp LIBRARIES utils)
cado_define_test(test_galois_action.cpp LIBRARIES utils tests)
cado_define_test(test_imaginary_quadratic_class_groups.cpp LIBRARIES utils tests)

cado_define_test(test_cachesize
    PROGRAM
    ${CADO_NFS_BINARY_DIR}/utils/cachesize
    TARGET_DEPENDENCIES cachesize
    # only make sure that the code reaches the end. This code is
    # certainly not robust enough to have a requirement on its output.I
    # wish the situation were a bit better.
    SHA1_ON_REGEXP_LINES "^#define L1_CACHE_SIZE [0-9]*$"
    EXPECT_SHA1 835f61ab6ba177ef4091c3231b0f30eaaad021b4
    SHA1_ON_SED_OUTPUT "s/^\\(#define L1_CACHE_SIZE\\).*/\\1/"
)

add_subdirectory(numbertheory)

if(SAGE)
    cado_define_test(
        SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/test_explain_relations_with_sagemath.sh
        PROVIDE_TEMPORARY_WDIR
        NO_DEFAULT_RUN
        TARGET_DEPENDENCIES
            makefb
            freerel
            las
            dup1
            dup2
            explain_indexed_relation
        TIMEOUT 60
    )
    add_dependencies(all_sagemath_test_dependencies
        test_explain_relations_with_sagemath_dependencies)

    # note that the sage preparser takes an absurd amount of time here,
    # so it's important to split into very small data sets.
    cado_divert_test(test_explain_relations_with_sagemath
        c30
        -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_explain_relations_with_sagemath.data/c30.poly
        -I 9
        -q0 56230 -q1 56240
        -lim0 30000 -lpb0 17 -mfb0 18
        -lim1 30000 -lpb1 17 -mfb1 18
    )
    cado_divert_test(test_explain_relations_with_sagemath_c30 dl APPEND_ARGUMENTS -dl)
endif()

