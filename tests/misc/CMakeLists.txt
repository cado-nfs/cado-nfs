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
)

# keep a file that is saved by freerel and read by debug_renumber.
cado_define_test(test_renumber_file
    NO_DEFAULT_RUN
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.sh
    -b ${CADO_NFS_BINARY_DIR} 
    PRECOMMAND
    "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
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
    -lpbs 25,26
    -check -build -quiet
    TARGET_DEPENDENCIES debug_renumber
    TIMEOUT 240
)

cado_divert_test(test_renumber_file c105_30003
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/test_renumber.data/c105_30003.poly
    -lpbs 25,26
    TIMEOUT 240
    )
endif()
