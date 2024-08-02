cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testsm.sh
    ${CADO_NFS_BINARY_DIR}/filter/sm
    PRECOMMAND
    "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
    TARGET_DEPENDENCIES sm)

set(sage_opt)
if(SAGE)
    set(sage_opt _with_sagemath)
endif()

cado_define_test(test_dup2${sage_opt}
   NO_DEFAULT_RUN
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_dup2.sh
   -b ${CADO_NFS_BINARY_DIR} 
   PRECOMMAND
   "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh"
   --env WORKDIR
   TARGET_DEPENDENCIES freerel debug_renumber dup2)

if(DEFINED ENV{CHECKS_EXPENSIVE})
cado_divert_test(test_dup2${sage_opt} c105_30003
    APPEND_ARGUMENTS
    -poly ${CADO_NFS_SOURCE_DIR}/tests/misc/test_renumber.data/c105_30003.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/c105_30003.rels.gz
    -lpbs 25,26
    TIMEOUT 120
    )
endif()

cado_divert_test(test_dup2${sage_opt} c105_30004
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.rels.gz
    -lpbs 14,14
    -dl
    -S daa017e17f50607820e7896dcc3a05370364b3b4
    -R 2959973f1
    )
