function(cado_divert_test TEST_BASE DIVERSION)
    # This function is appropriate when a test defined by
    # cado_define_test (possibly with NO_DEFAULT_RUN) needs to be
    # called several times with various parameter sets.
    # The different ways to call this function are:
    # 
    # optional arguments:
    #           TEST_DEPENDENCIES
    #           TARGET_DEPENDENCIES
    #           ARGUMENTS
    #           APPEND_ARGUMENTS
    #           NO_DEFAULT_RUN (useful if we want to divert further !)
    #
    # example use:
    #
    #     cado_define_test(TEST_NAME foo .....)
    #
    #     cado_divert_test(foo 1 --fast)
    #     cado_divert_test(foo 2 --thorough)
    #     cado_divert_test(foo extra --fallback TEST_DEPENDENCIES test3)
    #     cado_divert_test(foo extra TEST_DEPENDENCIES test3
    #                                ARGUMENTS --fallback)
    #
    # this defines tests foo_1, foo_2, and foo_extra, with the specific
    # arguments given. The two last syntaxes are equivalent.
    #
    # NOTE that if the base test name does indeed correspond to a test
    # (that is, it was not defined with NO_DEFAULT_RUN), then the
    # diverted test is made to depend on this base test.
    set(APPEND_ARGUMENTS)
    set(APPEND_ENVIRONMENT)
    set(PREPEND_PRECOMMAND)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    foreach(x
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
        )
        set(${x})
    endforeach()
    set(current ARGUMENTS)
    foreach(x ${ARGN})
        if (x STREQUAL "APPEND_ARGUMENTS")
            SET(current APPEND_ARGUMENTS)
        elseif (x STREQUAL "APPEND_ENVIRONMENT")
            SET(current APPEND_ENVIRONMENT)
        elseif (x STREQUAL "PREPEND_PRECOMMAND")
            SET(current PREPEND_PRECOMMAND)
        elseif (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        elseif (x IN_LIST CADO_NFS_TEST_KEYWORDS_USUAL)
            set(current ${x})
        elseif (x IN_LIST CADO_NFS_TEST_KEYWORDS_MULTI)
            set(current ${x})
        elseif (x IN_LIST CADO_NFS_TEST_KEYWORDS_FLAGS)
            set(current ${x})
            set(${x} 1)
        elseif (current IN_LIST CADO_NFS_TEST_KEYWORDS_LAST_WINS)
            set(${current} ${x})
        else()
            list(APPEND ${current} ${x})
        endif()
    endforeach(x)

    set(TEST_NAME "${TEST_BASE}_${DIVERSION}")

    cado_nfs_test_compute_lengths(APPEND_ARGUMENTS APPEND_ENVIRONMENT PREPEND_PRECOMMAND)

    cado_nfs_test_check_argument_rules()

    cado_nfs_test_propagate_all()
    cado_nfs_test_print_debug()
    cado_epilogue_create_test()
endfunction()
