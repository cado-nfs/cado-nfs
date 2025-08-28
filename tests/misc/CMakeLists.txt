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

cado_define_test(misc_cl_structure
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/misc/cl_structure_from_group_order -v -v
    SHA1_ON_REGEXP_LINES "^structure: "
    TARGET_DEPENDENCIES cl_structure_from_group_order
    NO_DEFAULT_RUN
)

set(D20_POLY_1 "inline-poly://n=-13784472188539234591/skew=3712744562.0/poly0=13784472188539234591,0,1")
cado_divert_test(misc_cl_structure d20_no_pSylow
    APPEND_ARGUMENTS
    -poly "${D20_POLY_1}" -order "2*617*2341309"
    EXPECT_SHA1 ad426daf55ed3ffdf96058602abb25c95c364476
)

set(D20_POLY_2 "inline-poly://n=-31415926535897932387/skew=5604991216.0/poly0=31415926535897932387,0,1")
cado_divert_test(misc_cl_structure d20_2and3Sylow
    APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order "2^2*3^4*7*13*53*541"
    EXPECT_SHA1 dd24813b875f503cf03f783f6ed8ac9d5aa6d902
)

# Test with an unnecessary factor in the group order
# The hash should be the same as above
cado_divert_test(misc_cl_structure d20_unnecessary_factor
    APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order "2^2,3^4,7,13,53,101^2,541"
    EXPECT_SHA1 dd24813b875f503cf03f783f6ed8ac9d5aa6d902
    EXPECT_FAIL
)

# Use -pSylow-bound 10 to disable the computation of the 3-Sylow (3^4 > 10)
cado_divert_test(misc_cl_structure d20_should_fail
    APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order "2,2,3,3,3,3,7,13,53,541" -pSylow-bound 10
    EXPECT_FAIL
)

# Test with a prime p < log^2(|D|) dividing the conductor
set(D20_POLY_3 "inline-poly://n=-24341965603669336362/skew=4933757756.0/poly0=24341965603669336362,0,1")
cado_divert_test(misc_cl_structure p_div_conductor
    APPEND_ARGUMENTS
    -poly "${D20_POLY_3}" -order "2^10,3^4,113,599"
    EXPECT_SHA1 b4f3f5f52b5d382ec2a1799d07ff4cb8c2528276
)

# Use only prime ideals below 100 (-B 100) to shorten the tests (the results is
# still correct as it is enough to generate the complete class group)
cado_divert_test(misc_cl_structure d130_2Sylow
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/d130.poly -order 2^9,3,5,26984470254026794913,37322833896097955713,5947362367424588311471 -B 100
    EXPECT_SHA1 bf6bbfdda8e54627fe491d7f7c418d856dcdb19c
)

cado_divert_test(misc_cl_structure error1 EXPECT_FAIL ARGUMENTS)
cado_divert_test(misc_cl_structure error2 EXPECT_FAIL ARGUMENTS squabble=42)
# fails because discriminant is positive
cado_divert_test(misc_cl_structure error3 EXPECT_FAIL APPEND_ARGUMENTS
    -poly "inline-poly://n=3/poly0=-3,0,1" -order 2
)
# fails due to missing args
cado_divert_test(misc_cl_structure error4 EXPECT_FAIL APPEND_ARGUMENTS
    -poly "${D20_POLY_2}"
)
# fails because values in -exps should be positive
cado_divert_test(misc_cl_structure error5 EXPECT_FAIL APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order 2,3,3^-1
)
# fails because values in -exps should be positive
cado_divert_test(misc_cl_structure error6 EXPECT_FAIL APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order 2^2,3^0
)
# fails because factors are missing from the given group order
cado_divert_test(misc_cl_structure error7 EXPECT_FAIL APPEND_ARGUMENTS
    -poly "${D20_POLY_2}" -order 2^2,3^4
)
