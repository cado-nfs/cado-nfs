
# The flint-fft library is only used if GF(p) lowlevel arithmetic is
# enabled, and we never worked on that in the 32-bit context. So we need
# to disable the flint-fft-only checks if BWC_GFP_ARITHMETIC_BACKENDS is
# empty.

# note that these tests are randomized, and takes variable time.
if(BWC_GFP_ARITHMETIC_BACKENDS)
    include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/flint-fft)
    include_directories(${PROJECT_SOURCE_DIR}/utils)
    cado_define_test(test_transform_length.cpp LIBRARIES flint-fft utils)
    mark_targets_as_openmp(test_transform_length)

    cado_define_test(test-flint.c
        NO_DEFAULT_RUN
        LIBRARIES flint-fft utils ${gmp_libname})
    mark_targets_as_openmp(test-flint)
    mark_targets_as_mpi(test-flint)
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
        mark_targets_as_mpi(test-flint-debug)
        mark_targets_as_openmp(test-flint-debug)
        set_property(TARGET test-flint-debug
            APPEND PROPERTY COMPILE_DEFINITIONS DEBUG_FFT)
        cado_define_test(test-flint-debug-magma
            PROVIDE_TEMPORARY_WDIR
            SCRIPT
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
    # randomness really has an impact on the length of these tests, so
    # let's us seeds.
    cado_divert_test(test-flint-debug-magma mul_fppol1 APPEND_ARGUMENTS mul_fppol -seed 19365)
    cado_divert_test(test-flint-debug-magma mul_fppol2 APPEND_ARGUMENTS mul_fppol -seed 17994)
    cado_divert_test(test-flint-debug-magma mp_fppol1 APPEND_ARGUMENTS mp_fppol -seed 9254)
    cado_divert_test(test-flint-debug-magma mp_fppol2 APPEND_ARGUMENTS mp_fppol -seed 3939)
    endif()

endif()
