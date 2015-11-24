
# doctests in cadocommand don't seem to work.
set(PYTHON_DOCTEST_SOURCES
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadoparams.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadoprograms.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadotask.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/workunit.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/wudb.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/wuserver.py
    ${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/toplevel.py
)

foreach (doctest_source ${PYTHON_DOCTEST_SOURCES})
  # Get filename without extension, e.g., cadoparams
  get_filename_component (basename ${doctest_source} NAME_WE)
  # Add doctest of Python sources, e.g., doctest cadoparams.py under test
  # name cadoparams
  cado_define_test(TEST_NAME test_python_${basename}
      PROGRAM python3 -m doctest ${doctest_source}
      TEST_DEPENDENCIES)
  set_tests_properties(test_python_${basename} PROPERTIES RUN_SERIAL 1)
endforeach (doctest_source)
