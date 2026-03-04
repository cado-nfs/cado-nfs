
# This is git grep -l '>>>' scripts/cadofactor/
set(PYTHON_DOCTEST_SOURCES
    scripts/cadofactor/cadoparams.py
    scripts/cadofactor/cadoprograms.py
    scripts/cadofactor/cadotask.py
    scripts/cadofactor/cadoutils.py
    scripts/cadofactor/workunit.py
    scripts/cadofactor/wudb.py
    scripts/cadofactor/toplevel.py
)

foreach (doctest_source ${PYTHON_DOCTEST_SOURCES})
  # Get filename without extension, e.g., cadoparams
  get_filename_component (basename ${doctest_source} NAME_WE)
  # Add doctest of Python sources, e.g., doctest cadoparams.py under test
  # name cadoparams
  cado_define_test(TEST_NAME test_python_${basename}
      ENVIRONMENT PYTHONPATH=${CADO_NFS_SOURCE_DIR}/scripts
      PRECOMMAND ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env TMPDIR
      PROGRAM
      python3 -m doctest ${CADO_NFS_SOURCE_DIR}/${doctest_source}
      )
  set_tests_properties(test_python_${basename} PROPERTIES RUN_SERIAL 1)
endforeach (doctest_source)

cado_define_test(test_python_wudb_cmdline
    PROVIDE_TEMPORARY_WDIR
    ENVIRONMENT PYTHONPATH=${CADO_NFS_SOURCE_DIR}/scripts
    SCRIPT
    "${CMAKE_CURRENT_SOURCE_DIR}/test_python_wudb_cmdline.sh" "${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/"
    TIMEOUT 60
)

if(DEFINED ENV{CHECKS_EXPENSIVE})
    cado_define_test(test_python_stress_test
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
        "${CMAKE_CURRENT_SOURCE_DIR}/test_python_stress_test.sh"
        AVOID_CONCURRENT
        TIMEOUT 240
        TARGET_DEPENDENCIES full_F7_dependencies
        NO_DEFAULT_RUN
    )
    cado_divert_test(test_python_stress_test sqlite3
        APPEND_ARGUMENTS --db-backend sqlite3
    )
    if(HAVE_PYTHON_MYSQL_CONNECTOR)
        cado_divert_test(test_python_stress_test mysql
            APPEND_ARGUMENTS --db-backend mysql
        )
    endif()

    # this one isn't really expensive, but it does pull
    # full_F7_dependencies, which belong to an expensive test
    cado_define_test(test_python_resubmit
        PROVIDE_TEMPORARY_WDIR
        SCRIPT
        "${CMAKE_CURRENT_SOURCE_DIR}/test_python_resubmit.sh"
        AVOID_CONCURRENT
        TARGET_DEPENDENCIES full_F7_dependencies
        TIMEOUT 90
    )
endif()
