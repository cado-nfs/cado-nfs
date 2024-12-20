macro(cado_nfs_test_get_TEST_NAME)
    # direct the implicit parameters somewhere
    if(nIMPLICIT GREATER 0) 
        if (((nSOURCES GREATER 0) OR (nSCRIPT GREATER 0)) AND (nTEST_NAME GREATER 0))
            message(FATAL_ERROR "bad syntax with implicit parameter list and both (SOURCES or SCRIPT) and TEST_NAME defined")
        elseif (((nSOURCES GREATER 0) OR (nSCRIPT GREATER 0)) AND (nTEST_NAME EQUAL 0))
            set(TEST_NAME ${IMPLICIT})
            list(LENGTH TEST_NAME nTEST_NAME)
        # from then on we know that SOURCES and SCRIPT are empty
        else()
            set(SOURCES ${IMPLICIT})
            list(LENGTH SOURCES nSOURCES)
        endif()
    endif()
    if(nSCRIPT GREATER 1)
        if (nARGUMENTS EQUAL 0)
            # Then use the tail as an argument list
            LIST(GET SCRIPT 0 x)
            LIST(REMOVE_AT SCRIPT 0)
            SET(ARGUMENTS ${SCRIPT})
            SET(SCRIPT ${x})
            list(LENGTH SCRIPT nSCRIPT)
        else()
            message(WARNING "Cannot have both SCRIPT and ARGUMENTS")
            LIST(GET SCRIPT 0 x)
            SET(SCRIPT ${x})
            list(LENGTH SCRIPT nSCRIPT)
        endif()
    endif()
    if(nSCRIPT GREATER 0)
        if (nSOURCES GREATER 0)
            message(FATAL_ERROR "SCRIPT and SOURCES incompatible")
        endif()
        if (nLIBRARIES GREATER 0)
            message(FATAL_ERROR "SCRIPT and LIBRARIES incompatible")
        endif()
    endif()

    if(nTEST_NAME EQUAL 0)
        # # message(STATUS "trying to find a test name")
        # # message(STATUS "SOURCES = ${SOURCES}")
        # # message(STATUS "SCRIPT = ${SCRIPT}")
        if(nSOURCES GREATER 0)
            # define the test name as the base name without extension of the
            # first source file.
            LIST(GET SOURCES 0 x)
            get_filename_component (TEST_NAME ${x} NAME_WE)
        elseif(nSCRIPT GREATER 0)
            # define the test name as the base name without extension of the
            # first source file.
            LIST(GET SCRIPT 0 x)
            get_filename_component (TEST_NAME ${x} NAME_WE)
        else()
            message(FATAL_ERROR "cannot find a name for test")
        endif()
    endif()

    if(NOT TEST_NAME)
        message(FATAL_ERROR "could not find a test name")
    endif()
endmacro()

macro(cado_nfs_test_handle_script_or_binary)
    # If we have a binary to build, then the meaning of
    # TARGET_DEPENDENCIES is most probably that we want the binary target
    # itself to depend on them. But bear in mind that the main use for
    # TARGET_DEPENDENCIES is probably for scripts anyway.
    if (nSCRIPT EQUAL 0)
        add_executable(${TEST_NAME} ${SOURCES})
        target_link_libraries(${TEST_NAME} ${LIBRARIES})
        set(EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})
        list(APPEND TARGET_DEPENDENCIES ${TEST_NAME})
        list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    else()
        if (${SCRIPT} MATCHES "\\.sh$")
            set(EXECUTABLE env bash ${SCRIPT})
        else()
            set(EXECUTABLE ${SCRIPT})
        endif()
    endif()
endmacro()

function(cado_define_test)
    # This function is appropriate for a test which consists in running a
    # binary. The different ways to call this function are:
    # 
    # cado_define_test(SOURCES foo.c bar.c
    #                        [TEST_NAME x]
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # cado_define_test(foo.c bar.c
    #                        [TEST_NAME x]
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # cado_define_test(foo
    #                        SOURCES foo.c bar.c
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # 
    # cado_define_test(foo
    #                        SCRIPT test.sh
    #                        [WORKING_DIRECTORY x]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # 
    # we define the "implicit" parameters as those occuring before any
    # tag of the form [SOURCES|LIBRARIES|...|...] (e.g foo.c bar.c in the
    # second example above).
    #
    # The test is named in either of the following ways:
    #   - after the TEST_NAME parameter if one is given,
    #   - as the implicit parameter if there is an explicit SOURCES or
    #     SCRIPT parameter and implicit parameters (there must not be
    #     more than one, then).
    #   - after the basename (with extension stripped) of the first source
    #     file or the SCRIPT parameter otherwise.
    #
    # The source files to be compiled are those defined by the SOURCES
    # parameter, or the implicit parameters otherwise.
    #
    # If the test consists in running a script (or any arbitrary external
    # command, really), then an explicit SCRIPT parameter must be
    # present, and explicit or implicit SOURCES parameters are forbidden.
    # If the SCRIPT parameter consists of more than one parameter and no
    # explicit ARGUMENTS parameter list is given, then the tail of the
    # SCRIPT parameter list is taken as an argument list.
    #
    # the PROGRAM tag has the exact same meaning as SCRIPT (handling is
    # identical)
    #
    # LIBRARIES, TARGET_DEPENDENCIES, TEST_DEPENDENCIES should be self
    # explanatory. LIBRARIES is incompatible with SCRIPT. Note that
    # TARGET_DEPENDENCIES is mostly useful for scripts (it ensures that
    # the targets have been built before running the script), but if
    # instead source files are given, the built executable is made to
    # depend on TARGET_DEPENDENCIES.
    #
    # ARGUMENTS defines which arguments are to be passed to the binary.
    # Optional.
    #
    # NO_DEFAULT_RUN indicates that no test is really defined, only a
    # baseline for further calls to cado_divert_test. (note that despite
    # the fact that no test is created in cmake's mind, we do need a test
    # name for further reference in cado_divert_test).
    #
    # Giving ARGUMENTS together with NO_DEFAULT_RUN just specifies a
    # common set of arguments to be passed to all diversions, provided
    # these are defined with APPEND_ARGUMENTS

    # Some additional features have been added after this function was
    # documented. We have ENVIRONMENT, PRECOMMAND, as well as their
    # APPEND_ friends.

    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(IMPLICIT)
    set(LIBRARIES)
    set(SOURCES)
    set(SCRIPT)
    set(TEST_NAME)
    foreach(x
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
        )
        set(${x})
    endforeach()
    set(current IMPLICIT)
    set(TIMEOUT 20)
    foreach(x ${ARGN})
        if (x STREQUAL "TEST_NAME")
            SET(current TEST_NAME)
        elseif (x STREQUAL "LIBRARIES")
            SET(current LIBRARIES)
        elseif (x STREQUAL "SOURCES")
            SET(current SOURCES)
        elseif ((x STREQUAL "SCRIPT") OR (x STREQUAL "PROGRAM"))
            SET(current SCRIPT)
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

    cado_nfs_test_compute_lengths(TEST_NAME SOURCES SCRIPT LIBRARIES IMPLICIT)

    cado_nfs_test_get_TEST_NAME()

    cado_nfs_test_check_argument_rules()

    cado_nfs_test_handle_script_or_binary()

    cado_nfs_test_print_debug()
    cado_epilogue_create_test()
endfunction()
