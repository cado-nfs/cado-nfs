
include_directories(${CADO_NFS_SOURCE_DIR}/linalg/bblas)

add_executable(test_bblas EXCLUDE_FROM_ALL
        test_bblas.cpp
        test_bblas_base.cpp
        test_bblas_level2.cpp
        test_bblas_level3.cpp
        test_bblas_level4_pluq.cpp
        test_bblas_level4_gauss.cpp
        test_bblas_level5.cpp
    )
target_link_libraries(test_bblas bblas utils)

cado_define_test(test_bblas_level2 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l2
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_bblas_level3 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l3
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_bblas_level4 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l4
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_bblas_level5 PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test_bblas --tests l5
    --fast n=1000
    TARGET_DEPENDENCIES test_bblas)

cado_define_test(test_gauss.c
    LIBRARIES bblas ${gmp_libname})
