cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testcheck_rels.sh
    ${CADO_NFS_BINARY_DIR}/misc/check_rels
    TARGET_DEPENDENCIES check_rels)

if (HAVE_GMPECM)
cado_define_test(descent_init_Fp_jl
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/misc/descent_init_Fp -jl
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p30.poly
    -mt 2 -minB1 200 -mineff 1000 -maxeff 100000 -side 1 -target 22
    -seed 42 701173953068971112417987441927 128476523876523762325
    TARGET_DEPENDENCIES descent_init_Fp)
endif()
