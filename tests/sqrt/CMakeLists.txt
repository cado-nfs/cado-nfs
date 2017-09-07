if (HAVE_CABSL AND HAVE_CLOG)
cado_define_test(test_crtalgsqrt
  SCRIPT ${CADO_NFS_SOURCE_DIR}/tests/sqrt/test_crtalgsqrt.sh
    ${CADO_NFS_BINARY_DIR}/sqrt/crtalgsqrt
    ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.dep.000
    ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.ratdep.000
    ${CADO_NFS_SOURCE_DIR}/tests/sqrt/c30.poly
  TARGET_DEPENDENCIES crtalgsqrt)
endif()
