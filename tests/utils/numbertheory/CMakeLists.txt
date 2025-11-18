
include_directories(${CADO_NFS_SOURCE_DIR}/numbertheory)

add_executable(test-ant test-ant.cpp)
target_link_libraries(test-ant utils)

if (SAGE)
    macro(numbertheory_sagemath_test name)
        cado_divert_test(${name} with_sagemath
            PROVIDE_TEMPORARY_WDIR
            PREPEND_PRECOMMAND
            ${CADO_NFS_SOURCE_DIR}/tests/pipe_to_sage.sh
        )

        add_dependencies(all_sagemath_test_dependencies
            ${name}_with_sagemath_dependencies)
    endmacro()
endif()

cado_define_test(
    numbertheory
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    TARGET_DEPENDENCIES
    test-ant
    NO_DEFAULT_RUN
    )

cado_divert_test(numbertheory p_maximal_order_batch
    APPEND_ARGUMENTS
        -test p-maximal-order-batch
        -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    )

cado_divert_test(numbertheory factorization_of_prime_batch
    APPEND_ARGUMENTS
        -test factorization-of-prime-batch
        -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    )

cado_divert_test(numbertheory valuations_of_ideal_batch
    APPEND_ARGUMENTS
        -test valuations-of-ideal-batch
        -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    )

cado_divert_test(numbertheory maximal_order_single
    APPEND_ARGUMENTS
        -test maximal-order
        -poly "10*x^3-x-2" -bound 20
    )

cado_divert_test(numbertheory p_maximal_order_single
    APPEND_ARGUMENTS
        -test p-maximal-order
        -poly "10*x^3-x-2" -prime 5
    )

cado_divert_test(numbertheory factorization_of_prime_single
    APPEND_ARGUMENTS
        -test factorization-of-prime
        -poly "10*x^3-x-2" -prime 5
    )

cado_divert_test(numbertheory valuations_of_ideal_single
    APPEND_ARGUMENTS
        -test valuations-of-ideal
        -poly "10*x^3-x-2" -prime 5 -elements "5,320*x^2-210*x+3"
    )

if (SAGE)
    numbertheory_sagemath_test(numbertheory_maximal_order_single)
    numbertheory_sagemath_test(numbertheory_p_maximal_order_single)
    numbertheory_sagemath_test(numbertheory_factorization_of_prime_single)
    numbertheory_sagemath_test(numbertheory_valuations_of_ideal_single)
endif()

cado_divert_test(numbertheory object_interface
    APPEND_ARGUMENTS
    -test nt-object-interface
    NO_DEFAULT_RUN
    )

cado_divert_test(numbertheory_object_interface test1 APPEND_ARGUMENTS -poly "10*x^3-x-2")
cado_divert_test(numbertheory_object_interface test2 APPEND_ARGUMENTS -poly "8*x^3 + 5*x + 7")

if (SAGE)
    numbertheory_sagemath_test(numbertheory_object_interface_test1)
    numbertheory_sagemath_test(numbertheory_object_interface_test2)
endif()

find_program(MAGMA magma HINTS ENV MAGMA)
if(MAGMA AND DEFINED ENV{CHECKS_EXPENSIVE})
cado_define_test(
    numbertheory_p_maximal_order_expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test p-maximal-order-batch -batch <(${MAGMA} -b -s /dev/null -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
cado_define_test(
    numbertheory_factorization_of_prime_expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test factorization-of-prime-batch -batch <(${MAGMA} -b -s /dev/null -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
cado_define_test(
    numbertheory_valuations_of_ideal_expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test valuations-of-ideal-batch -batch <(${MAGMA} -b -s /dev/null -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
endif()
