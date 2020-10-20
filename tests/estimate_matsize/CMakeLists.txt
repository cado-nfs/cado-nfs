
cado_define_test(estimate_matsize
    SCRIPT
    ${CADO_NFS_SOURCE_DIR}/scripts/estimate_matsize.sh
    -params ${CMAKE_CURRENT_SOURCE_DIR}/params.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/p30.poly
    ENVIRONMENT
    CADO_BUILD=${PROJECT_BINARY_DIR}
    TARGET_DEPENDENCIES makefb las fake_rels purge merge-dl
)
