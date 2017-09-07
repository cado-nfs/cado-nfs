
include_directories(${CADO_NFS_SOURCE_DIR}/linalg)

cado_define_test(test_bitlinalg_matops.cpp 
        ARGUMENTS --test-fast 1000 --n 1000
        LIBRARIES bitlinalg utils)
cado_define_test(test_bitlinalg_gauss.c
    LIBRARIES bitlinalg ${gmp_libname})


add_subdirectory (bwc)
