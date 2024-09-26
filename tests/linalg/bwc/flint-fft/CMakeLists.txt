
# The flint-fft library is only used if GF(p) lowlevel arithmetic is
# enabled, and we never worked on that in the 32-bit context. So we need
# to disable the flint-fft-only checks if BWC_GFP_ARITHMETIC_BACKENDS is
# empty.

# note that these tests are randomized, and takes variable time.
if(BWC_GFP_ARITHMETIC_BACKENDS)
    include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/flint-fft)
    include_directories(${PROJECT_SOURCE_DIR}/utils)
    cado_define_test(test_transform_length.cpp LIBRARIES flint-fft utils)
    set_target_properties(test_transform_length PROPERTIES
        COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
        LINK_FLAGS "${OpenMP_CXX_FLAGS}")

    cado_define_test(test-flint.c
        NO_DEFAULT_RUN
        LIBRARIES flint-fft utils ${gmp_libname})
    set_target_properties(test-flint PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS}"
        LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")
    # The mul0 test should activate the MFA code. Presently it's buggy with
    # openmp
    cado_divert_test(test-flint mul0 mul0 -s 20000)
    cado_divert_test(test-flint mul mul)
    # cado_divert_test(test-flint mulmod mulmod)
    cado_divert_test(test-flint mul_fppol mul_fppol)
    cado_divert_test(test-flint mp_fppol mp_fppol)

    if(DEFINED ENV{CHECKS_EXPENSIVE} AND MAGMA)
        add_executable(test-flint-debug EXCLUDE_FROM_ALL test-flint.c)
        target_link_libraries(test-flint-debug flint-fft-debug utils ${gmp_libname})
        set_target_properties(test-flint-debug PROPERTIES
            COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS} -DDEBUG_FFT"
            LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")
        cado_define_test(test-flint-debug-magma
            SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env TMPDIR
                ${CMAKE_CURRENT_SOURCE_DIR}/test-flint-debug-magma.sh
                ${MAGMA}
                ${CMAKE_CURRENT_SOURCE_DIR}/test-flint.m
                ${CMAKE_CURRENT_BINARY_DIR}/test-flint-debug
            TARGET_DEPENDENCIES test-flint-debug
            NO_DEFAULT_RUN
            # let's hope that even with the test randomness, all
            # executions will remain comfortably below this timeout.
            TIMEOUT 60
        )
    cado_divert_test(test-flint-debug-magma mul0 APPEND_ARGUMENTS mul0 -s 20000)
    cado_divert_test(test-flint-debug-magma mul APPEND_ARGUMENTS mul)
    # cado_divert_test(test-flint-debug-magma mulmod APPEND_ARGUMENTS mulmod)
    cado_divert_test(test-flint-debug-magma mul_fppol APPEND_ARGUMENTS mul_fppol)
    cado_divert_test(test-flint-debug-magma mp_fppol APPEND_ARGUMENTS mp_fppol)
    endif()

endif()
