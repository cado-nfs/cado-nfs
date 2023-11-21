
include_directories(${CADO_NFS_SOURCE_DIR}/numbertheory)

add_executable(test-ant test-ant.cpp)
target_link_libraries(test-ant utils)

cado_define_test(
    numbertheory-p-maximal-order
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test p-maximal-order-batch
    -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    TARGET_DEPENDENCIES
    test-ant
    )
cado_define_test(
    numbertheory-factorization-of-prime
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test factorization-of-prime-batch
    -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    TARGET_DEPENDENCIES
    test-ant
    )

cado_define_test(
    numbertheory-valuations-of-ideal
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test valuations-of-ideal-batch
    -batch ${CMAKE_CURRENT_SOURCE_DIR}/test.data
    TARGET_DEPENDENCIES
    test-ant
    )

find_program(MAGMA magma HINTS ENV MAGMA)
if(MAGMA AND DEFINED ENV{CHECKS_EXPENSIVE})
cado_define_test(
    numbertheory-p-maximal-order-expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test p-maximal-order-batch -batch <(${MAGMA} -b -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
cado_define_test(
    numbertheory-factorization-of-prime-expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test factorization-of-prime-batch -batch <(${MAGMA} -b -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
cado_define_test(
    numbertheory-valuations-of-ideal-expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test valuations-of-ideal-batch -batch <(${MAGMA} -b -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    TIMEOUT 60
    )
endif()
