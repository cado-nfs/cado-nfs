cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testcheck_rels.sh
    ${CADO_NFS_BINARY_DIR}/misc/check_rels
    TARGET_DEPENDENCIES check_rels)

if (HAVE_GMPECM)
cado_define_test(descent_init_Fp_jl
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/misc/descent_init_Fp -jl
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p30.poly
    -mt 2 -minB1 200 -mineff 1000 -maxeff 100000 -side 1 -lpb 22
    -seed 42 701173953068971112417987441927 128476523876523762325
    TARGET_DEPENDENCIES descent_init_Fp)
endif()

# test renumber
cado_define_test(test_renumber_on_the_fly_mnfsdl
    PROGRAM ${PROJECT_BINARY_DIR}/misc/debug_renumber
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.data/mnfs5.poly
    -lpbs 11,10,10,10,10
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber)

cado_define_test(test_renumber_on_the_fly_snfs
    PROGRAM ${PROJECT_BINARY_DIR}/misc/debug_renumber
    -poly ${PROJECT_SOURCE_DIR}/parameters/polynomials/F9.poly
    -lpbs 10,10
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber)

cado_define_test(test_renumber_on_the_fly_gnfs
    PROGRAM ${PROJECT_BINARY_DIR}/misc/debug_renumber
    -poly ${PROJECT_SOURCE_DIR}/parameters/polynomials/rsa768.poly
    -lpbs 10,10
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber)

# This is one of the failure conditions that we see in bug 30003
# times out easily in coverage mode
cado_define_test(test_renumber_on_the_fly_gnfs2_30003
    PROGRAM ${PROJECT_BINARY_DIR}/misc/debug_renumber
    -poly ${PROJECT_SOURCE_DIR}/parameters/polynomials/rsa768.poly
    -lpbs 20,21
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber
    TIMEOUT 60
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
)

# keep a file that is saved by freerel and read by debug_renumber.
cado_define_test(test_renumber_file
    NO_DEFAULT_RUN
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.sh
    -b ${CADO_NFS_BINARY_DIR} 
    PROVIDE_TEMPORARY_WDIR
    TARGET_DEPENDENCIES freerel debug_renumber)

# Test MNFS with 5 polys (interesting primes are 2 and 3)
cado_divert_test(test_renumber_file mnfs5
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.data/mnfs5.poly
    -lpbs 11,10,10,10,10
    )

if(DEFINED ENV{CHECKS_EXPENSIVE})
cado_define_test(test_renumber_on_the_fly_c105_30003
    PROGRAM ${PROJECT_BINARY_DIR}/misc/debug_renumber
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.data/c105_30003.poly
    -lpbs 22,23
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber
    TIMEOUT 60
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
)

cado_divert_test(test_renumber_file c105_30003
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.data/c105_30003.poly
    -lpbs 22,23
    TIMEOUT 60
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
    )
endif()

if (HAVE_GMPECM)
    cado_define_test(misc_factor
        SCRIPT
        ${CADO_NFS_BINARY_DIR}/misc/factor -v -v -N
        SHA1_ON_REGEXP_LINES "^factorization: "
        TARGET_DEPENDENCIES factor
        NO_DEFAULT_RUN
    )

    # p is easily found, and the remaining cofactor happens to be a
    # square
    cado_divert_test(misc_factor pq2r2 APPEND_ARGUMENTS
        9448184470610043110564621588359899555130833989152313
        EXPECT_SHA1 8c21546928dc23c5eb0b51a2f176499cd7ab334c
    )

    cado_divert_test(misc_factor complete APPEND_ARGUMENTS
        66045451192627703404231139322032983287173736630736
        EXPECT_SHA1 f03dea0ac1a2de340cafc423885ddcbd94bc55df
    )

    # C77 contains two large factors, by using -effort 0 we ensure that they
    # will not be found and the binary will return 1
    cado_divert_test(misc_factor incomplete APPEND_ARGUMENTS
        84884147409828091725676728670213067387206838101828807864190286991865870575397
        -effort 0
        EXPECT_FAIL
    )

    cado_divert_test(misc_factor error1 EXPECT_FAIL APPEND_ARGUMENTS 0)
    cado_divert_test(misc_factor error2 EXPECT_FAIL APPEND_ARGUMENTS
        65537 trialdiv_bound=1
    )
    cado_divert_test(misc_factor error3 EXPECT_FAIL APPEND_ARGUMENTS
        65537 --isprime-niter 0
    )
    cado_divert_test(misc_factor error4 EXPECT_FAIL ARGUMENTS)
    cado_divert_test(misc_factor error5 EXPECT_FAIL ARGUMENTS squabble=42)

    set(C214 1491853389933411857714186950405505504385922663380730145523965911837516584504855532767352946927077785080156547980196625376697559081935835336314328338395762630753432863463017078387556027886346533171763392847320002862)
    # 65537 is a bogus hint, but the code should happily discard it
    set(HINTS_C214
        2103814835105734276966733904285860990868274236590576967285421378423293661987882101115537118852349489076143360401417774069682013332768445918577301453043,121,544478094571057,65537)
    cado_divert_test(misc_factor complete_with_hints
        APPEND_ARGUMENTS
        ${C214} -hints ${HINTS_C214}
        EXPECT_SHA1 4fce16df8636756a244aa4928c25d01878fff6e2
    )
endif()
