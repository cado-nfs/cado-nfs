cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testsm.sh
    ${CADO_NFS_BINARY_DIR}/filter/sm
    PRECOMMAND
    "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh" --env WORKDIR
    TARGET_DEPENDENCIES sm)

cado_define_test(test_dup2
   NO_DEFAULT_RUN
   SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_dup2.sh
   -b ${CADO_NFS_BINARY_DIR} 
   PRECOMMAND
   "${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh"
   --env WORKDIR
   TARGET_DEPENDENCIES freerel debug_renumber dup2)

if(DEFINED ENV{CHECKS_EXPENSIVE})
cado_divert_test(test_dup2 c105_30003
    APPEND_ARGUMENTS
    -poly ${CADO_NFS_SOURCE_DIR}/tests/misc/test_renumber.data/c105_30003.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/c105_30003.rels.gz
    -lpbs 25,26
    TIMEOUT 120
    )
endif()

cado_divert_test(test_dup2 c105_30004
    APPEND_ARGUMENTS
    -poly ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.poly
    -rels ${CMAKE_CURRENT_SOURCE_DIR}/p40_30004.rels.gz
    -lpbs 14,14
    -dl
    -S c8c5a0aa4750790e8e5e150bf13de7b1e77ade01
    -R b40784b0e
    )
# b40784b0e is old working commit 122e3868, plus this test
