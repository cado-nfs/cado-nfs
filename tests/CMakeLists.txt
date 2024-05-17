add_custom_target(all_test_dependencies)

# cmake understands tests and targets as being two different things. It
# is perfectly common to have both a target and a test having the same
# name X. Yet testing X will not trigger building X, and that is
# admittedly unfortunate.

# INTERNALS:
#
# we define several macros for dealing with this in a better way.
#
#  For each test X, we define another test builddep_X.
#
#  The test builddep_X which only takes the necessary steps for building
#  the dependencies of the test X. Those dependencies form a target
#  named X_dependencies.
#
#  The target X_dependencies may perhaps include a target named X.
#
#  The test X is made to depend on builddep_X
#
#  The target all_test_dependencies is made to depend on X_dependencies
#
#  Tests may have other dependency relationships, which should be set
#  in a manner which is compatible with our dependency X <- builddep_X.
#  This is taken care of by update_test_build_depends

# WHICH CMAKE COMMANDS SHOULD BE USED ?
#
# only the following:
#
#       cado_define_test()
#       cado_divert_test()
#
# those are documented below within their source definitions.
#
# Note that for test X, there are some use cases where it is nice to know
# about X_dependencies (for adding extra dependencies late), or the
# targets X (for setting compile time definitions).

# set(CADO_NFS_LOCK_TOOL)
# find_program(HAVE_LOCKF_TOOL NAMES flock DOC "locking tool lockf")
# if (HAVE_LOCKF_TOOL)
#     set(CADO_NFS_LOCK_TOOL ${HAVE_LOCKF_TOOL})
# else()
#     find_program(HAVE_LCKDO_TOOL NAMES flock DOC "locking tool lockf")
#     if (HAVE_LCKDO_TOOL)
#         set(CADO_NFS_LOCK_TOOL ${HAVE_LCKDO_TOOL} -w)
#     endif()
# endif()
# 
math(EXPR test_index 0)

# openmp jobs really don't like competing with eachother, so we must be
# cautious. Setting OMP_DYNAMIC=true seems to alleviate this problem
# significantly
math(EXPR test_concurrency_tolerance 2)

if(TIMEOUT_SCALE)
    message(STATUS "TIMEOUT_SCALE is set to ${TIMEOUT_SCALE}")
endif()

set(CADO_NFS_TEST_KEYWORDS_LAST_WINS TIMEOUT AVOID_CONCURRENT)
set(CADO_NFS_TEST_KEYWORDS_FIRST_WINS_NOISY WORKING_DIRECTORY TEST_NAME)

macro(cado_set_test_timeout TEST_NAME TIMEOUT)
    if(TIMEOUT_SCALE)
        MATH(EXPR scaled_timeout "${TIMEOUT}*${TIMEOUT_SCALE}")
        set_property(TEST ${TEST_NAME} PROPERTY TIMEOUT
            ${scaled_timeout})
    else()
        set_property(TEST ${TEST_NAME} PROPERTY TIMEOUT ${TIMEOUT})
    endif()
endmacro()

macro(cado_define_test)
    # This macro is appropriate for a test which consists in running a
    # binary. The different ways to call this macro are:
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

    # Some additional features have been added after this macro was
    # documented. We have ENVIRONMENT, PRECOMMAND, as well as their
    # APPEND_ friends.

    set(IMPLICIT)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(LIBRARIES)
    set(SOURCES)
    set(SCRIPT)
    set(ARGUMENTS)
    set(ENVIRONMENT)
    set(PRECOMMAND)
    set(TEST_NAME)
    set(WORKING_DIRECTORY)
    set(NO_DEFAULT_RUN)
    set(AVOID_CONCURRENT)
    set(EXPECT_FAIL)
    set(current IMPLICIT)
    set(TIMEOUT 20)
    foreach(x ${ARGN})
        if (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        elseif (x STREQUAL "LIBRARIES")
            SET(current LIBRARIES)
        elseif (x STREQUAL "SOURCES")
            SET(current SOURCES)
        elseif ((x STREQUAL "SCRIPT") OR (x STREQUAL "PROGRAM"))
            SET(current SCRIPT)
        elseif (x STREQUAL "ARGUMENTS")
            SET(current ARGUMENTS)
        elseif (x STREQUAL "ENVIRONMENT")
            SET(current ENVIRONMENT)
        elseif (x STREQUAL "PRECOMMAND")
            SET(current PRECOMMAND)
        elseif (x STREQUAL "TEST_NAME")
            SET(current TEST_NAME)
        elseif (x STREQUAL "WORKING_DIRECTORY")
            SET(current WORKING_DIRECTORY)
        elseif (x STREQUAL "NO_DEFAULT_RUN")
            SET(current NO_DEFAULT_RUN)
            # special, because it's a flag, really
            set(NO_DEFAULT_RUN 1)
        elseif (x STREQUAL "AVOID_CONCURRENT")
            SET(current AVOID_CONCURRENT)
            set(AVOID_CONCURRENT 1)
        elseif (x STREQUAL "EXPECT_FAIL")
            SET(current EXPECT_FAIL)
            set(EXPECT_FAIL 1)
        elseif (x STREQUAL "TIMEOUT")
            SET(current TIMEOUT)
        elseif (current IN_LIST CADO_NFS_TEST_KEYWORDS_LAST_WINS)
            set(${current} ${x})
        else()
            list(APPEND ${current} ${x})
        endif()
    endforeach(x)
    # do some checking.
    #message(STATUS "IMPLICIT is ${IMPLICIT}")
    #message(STATUS "TEST_DEPENDENCIES is ${TEST_DEPENDENCIES}")
    #message(STATUS "TARGET_DEPENDENCIES is ${TARGET_DEPENDENCIES}")
    #message(STATUS "LIBRARIES is ${LIBRARIES}")
    #message(STATUS "SOURCES is ${SOURCES}")
    #message(STATUS "SCRIPT is ${SCRIPT}")
    #message(STATUS "ARGUMENTS is ${ARGUMENTS}")
    #message(STATUS "ENVIRONMENT is ${ENVIRONMENT}")
    #message(STATUS "PRECOMMAND is ${PRECOMMAND}")
    #message(STATUS "TEST_NAME is ${TEST_NAME}")
    #message(STATUS "WORKING_DIRECTORY is ${WORKING_DIRECTORY}")
    #message(STATUS "NO_DEFAULT_RUN is ${NO_DEFAULT_RUN}")
    list(LENGTH TEST_NAME nTEST_NAME)
    list(LENGTH WORKING_DIRECTORY nWORKING_DIRECTORY)
    list(LENGTH TIMEOUT nTIMEOUT)
    list(LENGTH SOURCES nSOURCES)
    list(LENGTH SCRIPT nSCRIPT)
    list(LENGTH LIBRARIES nLIBRARIES)
    list(LENGTH IMPLICIT nIMPLICIT)
    list(LENGTH NO_DEFAULT_RUN nNO_DEFAULT_RUN)
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
    list(LENGTH ENVIRONMENT nENVIRONMENT)
    list(LENGTH PRECOMMAND nPRECOMMAND)
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
    foreach(x ${CADO_NFS_TEST_KEYWORDS_FIRST_WINS_NOISY})
        if(n${x} GREATER 1)
            message(WARNING "too many ${x} directives for ${TEST_NAME} ${${x}}, retaining only first")
            LIST(GET ${x} 0 y)
            SET(${x} ${y})
        endif()
    endforeach()

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
    if(nNO_DEFAULT_RUN GREATER 1)
        message(FATAL_ERROR "discarded arguments after NO_DEFAULT_RUN for test ${TEST_NAME}")
        set(NO_DEFAULT_RUN 1)
    endif()

    # If we have a binary to build, then the meaning of
    # TARGET_DEPENDENCIES is most probably that we want the binary target
    # itself to depend on them. But bear in mind that the main use for
    # TARGET_DEPENDENCIES is probably for scripts anyway.
    if (nSCRIPT EQUAL 0)
        add_executable(${TEST_NAME} ${SOURCES})
        target_link_libraries(${TEST_NAME} ${LIBRARIES})
        if(nTARGET_DEPENDENCIES GREATER 0)
            add_dependencies(${TEST_NAME} ${TARGET_DEPENDENCIES})
        endif()
        set(TARGET_DEPENDENCIES ${TEST_NAME})
        list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    endif()

    if(NOT TEST_NAME)
        message(FATAL_ERROR "could not find a test name")
        endif()
    # even script tests may define target dependencies, even though it's
    # probably of little use.
    add_custom_target(${TEST_NAME}_dependencies)
    add_dependencies(all_test_dependencies ${TEST_NAME}_dependencies)

    if(nTARGET_DEPENDENCIES GREATER 0)
        add_dependencies(${TEST_NAME}_dependencies ${TARGET_DEPENDENCIES})
        # if(CADO_NFS_LOCK_TOOL)
        # add_test(builddep_${TEST_NAME} ${CADO_NFS_LOCK_TOOL} ${PROJECT_BINARY_DIR}/build-lock ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)
        # else()
        add_test(builddep_${TEST_NAME} ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)
        # endif()
    endif()

    if (nSCRIPT EQUAL 0)
        set(${TEST_NAME}_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})
    else()
        set(binary ${SCRIPT})
        if (${SCRIPT} MATCHES "\\.sh$")
            set(binary env bash ${SCRIPT})
        endif()
        set(${TEST_NAME}_EXECUTABLE ${binary})
    endif()

    # save for use by diverted tests
    set(${TEST_NAME}_ARGUMENTS ${ARGUMENTS})
    set(${TEST_NAME}_ENVIRONMENT ${ENVIRONMENT})
    set(${TEST_NAME}_PRECOMMAND ${PRECOMMAND})
    set(${TEST_NAME}_AVOID_CONCURRENT ${AVOID_CONCURRENT})
    set(${TEST_NAME}_EXPECT_FAIL ${EXPECT_FAIL})
    set(${TEST_NAME}_TIMEOUT ${TIMEOUT})

    if(NOT NO_DEFAULT_RUN)
        set(test_env env)
        foreach(v CADO_NFS_SOURCE_DIR CADO_NFS_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR CMAKE_CURRENT_BINARY_DIR PROJECT_BINARY_DIR)
            set(test_env ${test_env} ${v}=${${v}})
        endforeach()
        if (nENVIRONMENT GREATER 0)
            set(test_env ${test_env} ${${TEST_NAME}_ENVIRONMENT})
        endif()
        MATH(EXPR test_index "${test_index}+1")
        set(A2F)
        if(EXPECT_FAIL)
            list(APPEND A2F ${CADO_NFS_SOURCE_DIR}/tests/abort-to-fail.sh)
        endif()
        # remove the pre-call to ${${TEST_NAME}_LOCKF} (right before
        # A2F), because it doesn't play nicely with TIMEOUT. It's likely
        # that using RESOURCE_LOCK below is a better idea.
        set(add_test_args NAME ${TEST_NAME} COMMAND ${test_env} ${A2F} ${${TEST_NAME}_PRECOMMAND} ${CADO_NFS_SOURCE_DIR}/tests/test.sh ${${TEST_NAME}_EXECUTABLE} ${ARGUMENTS})
        if(nWORKING_DIRECTORY GREATER 0)
            list(APPEND add_test_args WORKING_DIRECTORY ${WORKING_DIRECTORY})
        endif()
        add_test(${add_test_args})
        if(AVOID_CONCURRENT)
            set(${TEST_NAME}_LOCKED_RESOURCES)
            set(n ${AVOID_CONCURRENT})
            if(n GREATER ${test_concurrency_tolerance})
                set(n ${test_concurrency_tolerance})
            endif()
            foreach(j RANGE 1 ${n})
                MATH(EXPR testmod "(${test_index} + ${j}) % ${test_concurrency_tolerance}")
                LIST(APPEND ${TEST_NAME}_LOCKED_RESOURCES lock.${testmod})
            endforeach()
            if(n GREATER_EQUAL 2)
                message(STATUS "Test ${TEST_NAME} locks resources ${${TEST_NAME}_LOCKED_RESOURCES}")
            endif()
            set_property(TEST ${TEST_NAME} PROPERTY RESOURCE_LOCK
                ${${TEST_NAME}_LOCKED_RESOURCES})
            # set_property(TEST ${TEST_NAME} PROPERTY RESOURCE_LOCK lock.${testmod})
        endif()
        # if(CADO_NFS_LOCK_TOOL AND AVOID_CONCURRENT)
        #     MATH(EXPR testmod "${test_index} % ${test_concurrency_tolerance}")
        #     set(${TEST_NAME}_LOCKF ${CADO_NFS_LOCK_TOOL} ${PROJECT_BINARY_DIR}/lock.${testmod})
        # endif()
        if(nTIMEOUT GREATER 0)
            cado_set_test_timeout(${TEST_NAME} ${TIMEOUT})
        endif()
        set(${TEST_NAME}_TEST_EXISTS 1)
        if(nTARGET_DEPENDENCIES GREATER 0)
            # seems that we need set_tests_properties to take only a
            # ;-separated list here.
            set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
            string_join(all_deps ";" ${deps})
            set_tests_properties (${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
        endif()
        set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 125)
        if (EXPECT_FAIL)
            set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL 1)
        endif()
    else()
        # There is a special case. Even a NO_DEFAULT_RUN test may have
        # test dependencies, because that's a handy way to have all
        # inherited tests depend on these.
        set(${TEST_NAME}_TEST_DEPENDENCIES ${TEST_DEPENDENCIES})
    endif()
        if(nTIMEOUT GREATER 0)
            # message(STATUS "Timeout for ${TEST_NAME} is ${TIMEOUT}")
        endif()
endmacro()

macro(cado_divert_test TEST_BASE DIVERSION)
    # This macro is appropriate when a test defined by
    # cado_define_test (possibly with NO_DEFAULT_RUN) needs to be
    # called several times with various parameter sets.
    # The different ways to call this macro are:
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
    set(ARGUMENTS)
    set(ENVIRONMENT)
    set(PRECOMMAND)
    set(APPEND_ARGUMENTS)
    set(APPEND_ENVIRONMENT)
    set(PREPEND_PRECOMMAND)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(NO_DEFAULT_RUN)
    set(AVOID_CONCURRENT)
    set(EXPECT_FAIL)
    set(current ARGUMENTS)
    set(TIMEOUT)
    foreach(x ${ARGN})
        if (x STREQUAL "ARGUMENTS")
            SET(current ARGUMENTS)
        elseif (x STREQUAL "ENVIRONMENT")
            SET(current ENVIRONMENT)
        elseif (x STREQUAL "PRECOMMAND")
            SET(current PRECOMMAND)
        elseif (x STREQUAL "APPEND_ARGUMENTS")
            SET(current APPEND_ARGUMENTS)
        elseif (x STREQUAL "APPEND_ENVIRONMENT")
            SET(current APPEND_ENVIRONMENT)
        elseif (x STREQUAL "PREPEND_PRECOMMAND")
            SET(current PREPEND_PRECOMMAND)
        elseif (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        elseif (x STREQUAL "NO_DEFAULT_RUN")
            SET(current NO_DEFAULT_RUN)
            # special, because it's a flag, really
            set(NO_DEFAULT_RUN 1)
        elseif (x STREQUAL "AVOID_CONCURRENT")
            SET(current AVOID_CONCURRENT)
            set(AVOID_CONCURRENT 1)
        elseif (x STREQUAL "TIMEOUT")
            SET(current TIMEOUT)
        elseif (x STREQUAL "EXPECT_FAIL")
            SET(current EXPECT_FAIL)
            set(EXPECT_FAIL 1)
        elseif (current IN_LIST CADO_NFS_TEST_KEYWORDS_LAST_WINS)
            set(${current} ${x})
        else()
            list(APPEND ${current} ${x})
        endif()
    endforeach(x)
    # do some checking.
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
    list(LENGTH ENVIRONMENT nENVIRONMENT)
    list(LENGTH PRECOMMAND nPRECOMMAND)
    list(LENGTH TIMEOUT nTIMEOUT)
    list(LENGTH APPEND_ARGUMENTS nAPPEND_ARGUMENTS)
    list(LENGTH APPEND_ENVIRONMENT nAPPEND_ENVIRONMENT)
    list(LENGTH PREPEND_PRECOMMAND nPREPEND_PRECOMMAND)
    list(LENGTH NO_DEFAULT_RUN nNO_DEFAULT_RUN)
    if ((nAPPEND_ARGUMENTS GREATER 0) AND (nARGUMENTS GREATER 0))
        message(FATAL_ERROR "ARGUMENTS and APPEND_ARGUMENTS are incompatible")
    endif()
    if ((nAPPEND_ENVIRONMENT GREATER 0) AND (nENVIRONMENT GREATER 0))
        message(FATAL_ERROR "ENVIRONMENT and APPEND_ENVIRONMENT are incompatible")
    endif()
    if ((nPREPEND_PRECOMMAND GREATER 0) AND (nPRECOMMAND GREATER 0))
        message(FATAL_ERROR "PRECOMMAND and PREPEND_PRECOMMAND are incompatible")
    endif()

    set(TEST_NAME "${TEST_BASE}_${DIVERSION}")

    if (NOT nARGUMENTS GREATER 0)
        set(ARGUMENTS ${${TEST_BASE}_ARGUMENTS} ${APPEND_ARGUMENTS})
        list(LENGTH ARGUMENTS nARGUMENTS)
    endif()

    if (NOT nENVIRONMENT GREATER 0)
        set(ENVIRONMENT ${${TEST_BASE}_ENVIRONMENT} ${APPEND_ENVIRONMENT})
        list(LENGTH ENVIRONMENT nENVIRONMENT)
    endif()

    if (NOT nPRECOMMAND GREATER 0)
        set(PRECOMMAND ${PREPEND_PRECOMMAND} ${${TEST_BASE}_PRECOMMAND})
        list(LENGTH PRECOMMAND nPRECOMMAND)
    endif()

    if (NOT AVOID_CONCURRENT GREATER 0)
        set(AVOID_CONCURRENT ${${TEST_BASE}_AVOID_CONCURRENT})
    endif()

    if (NOT EXPECT_FAIL GREATER 0)
        set(EXPECT_FAIL ${${TEST_BASE}_EXPECT_FAIL})
    endif()

    if (NOT nTIMEOUT GREATER 0)
        set(TIMEOUT ${${TEST_BASE}_TIMEOUT})
        # message(STATUS "Timeout for ${TEST_NAME} is inherited from ${TEST_BASE}: ${TIMEOUT}")
        list(LENGTH TIMEOUT nTIMEOUT)
    endif()

    # we want to define builddep_X unconditionally, depending at least on
    # the dependencies of the base test.
    add_custom_target(${TEST_NAME}_dependencies)
    add_dependencies(all_test_dependencies ${TEST_NAME}_dependencies)
    add_dependencies(${TEST_NAME}_dependencies ${TEST_BASE}_dependencies ${TARGET_DEPENDENCIES})
    add_test(builddep_${TEST_NAME} ${CMAKE_COMMAND} --build
        ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)


    set(${TEST_NAME}_EXECUTABLE ${${TEST_BASE}_EXECUTABLE})
    set(${TEST_NAME}_ARGUMENTS ${ARGUMENTS})
    set(${TEST_NAME}_ENVIRONMENT ${ENVIRONMENT})
    set(${TEST_NAME}_PRECOMMAND ${PRECOMMAND})
    set(${TEST_NAME}_AVOID_CONCURRENT ${AVOID_CONCURRENT})
    set(${TEST_NAME}_EXPECT_FAIL ${EXPECT_FAIL})
    set(${TEST_NAME}_TIMEOUT ${TIMEOUT})

    if(NOT NO_DEFAULT_RUN)
        set(test_env env)
        foreach(v CADO_NFS_SOURCE_DIR CADO_NFS_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR CMAKE_CURRENT_BINARY_DIR PROJECT_BINARY_DIR)
            set(test_env ${test_env} ${v}=${${v}})
        endforeach()
        if (nENVIRONMENT GREATER 0)
            set(test_env ${test_env} ${${TEST_NAME}_ENVIRONMENT})
        endif()
        MATH(EXPR test_index "${test_index}+1")
        # if(CADO_NFS_LOCK_TOOL AND AVOID_CONCURRENT)
        #     MATH(EXPR testmod "${test_index} % ${test_concurrency_tolerance}")
        #     set(${TEST_NAME}_LOCKF ${CADO_NFS_LOCK_TOOL} ${PROJECT_BINARY_DIR}/lock.${testmod})
        # endif()
        set(A2F)
        if(EXPECT_FAIL)
            list(APPEND A2F ${CADO_NFS_SOURCE_DIR}/tests/abort-to-fail.sh)
        endif()
        # dropping ${${TEST_NAME}_LOCKF} (see in cado_define_test)
        add_test(NAME ${TEST_NAME} COMMAND ${test_env} ${A2F} ${${TEST_NAME}_PRECOMMAND} ${CADO_NFS_SOURCE_DIR}/tests/test.sh ${${TEST_BASE}_EXECUTABLE} ${ARGUMENTS})
        if(AVOID_CONCURRENT)
            set(${TEST_NAME}_LOCKED_RESOURCES)
            set(n ${AVOID_CONCURRENT})
            if(n GREATER ${test_concurrency_tolerance})
                set(n ${test_concurrency_tolerance})
            endif()
            foreach(j RANGE 1 ${n})
                MATH(EXPR testmod "(${test_index} + ${j}) % ${test_concurrency_tolerance}")
                LIST(APPEND ${TEST_NAME}_LOCKED_RESOURCES lock.${testmod})
            endforeach()
            if(n GREATER_EQUAL 2)
                message(STATUS "Test ${TEST_NAME} locks resources ${${TEST_NAME}_LOCKED_RESOURCES}")
            endif()
            set_property(TEST ${TEST_NAME} PROPERTY RESOURCE_LOCK
                ${${TEST_NAME}_LOCKED_RESOURCES})
            # MATH(EXPR testmod "${test_index} % ${test_concurrency_tolerance}")
            # set_property(TEST ${TEST_NAME} PROPERTY RESOURCE_LOCK lock.${testmod})
        endif()
        if(nTIMEOUT GREATER 0)
            cado_set_test_timeout(${TEST_NAME} ${TIMEOUT})
        endif()
        set(${TEST_NAME}_TEST_EXISTS 1)
        # prepare a ;-separated list for set_tests_properties
        # note that we do not need to make the test depend on
        # builddep_${TEST_BASE}, since we have made ${TEST_NAME}_dependencies
        # depend on ${TEST_BASE}_dependencies. It's in fact relieving, given
        # that builddep_${TEST_BASE} does not necessarily exist.
        set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
        if(${TEST_BASE}_TEST_EXISTS)
            list(APPEND deps ${TEST_BASE})
        else()
            # ok, TEST_BASE does not exist, but we've still built a list
            # of its test dependencies that we want to obey.
            list(APPEND deps ${${TEST_BASE}_TEST_DEPENDENCIES})
        endif()
        string_join(all_deps ";" ${deps})
        set_tests_properties(${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
        set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 125)
        if (EXPECT_FAIL)
            set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL 1)
        endif()
    else()
        set(deps ${TEST_DEPENDENCIES})
        if(${TEST_BASE}_TEST_EXISTS)
            list(APPEND deps ${TEST_BASE})
        else()
            list(APPEND deps ${${TEST_BASE}_TEST_DEPENDENCIES})
        endif()
        set(${TEST_NAME}_TEST_DEPENDENCIES ${deps})
    endif()
endmacro()

add_custom_target(full_dependencies_base)
add_dependencies(full_dependencies_base
        polyselect polyselect_ropt polyselect3
        freerel makefb
        las
        dup1 dup2 purge
        )

# MPI binaries are not properly launched by cado-nfs.py (see #21819)
if(NOT MINGW AND NOT HAVE_MPI)

cado_define_test(full_c30
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/cado-nfs.py 999073468111577057576445816581
        --server-threads 2
    TARGET_DEPENDENCIES
        full_dependencies_base
        merge replay
        characters
        bwc_full_gf2
        sqrt
    AVOID_CONCURRENT
    TIMEOUT 360
        )
    if(pz IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
        # no need to have both the JL and non-JL in normal tests, I
        # think, given that there's always at least one machine that does
        # the expensive tests.
        cado_define_test(full_p30
            SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            -dlp -ell 445257213376519703258944813
            -t 2 -dlp-no-keep
            target=92800609832959449330691138186
            3037544709654617415632521514287
            TARGET_DEPENDENCIES full_dependencies_base
            bwc_full_gfp
            lingen_pz sm sm_simple las_descent reconstructlog-dl
            merge-dl replay-dl numbertheory_tool
            AVOID_CONCURRENT
            TIMEOUT 360
        )
        # The following case was hitting the bug #21707
        set(full_p30_JL_params
            fastSM=true
            jlpoly=true tasks.polyselect.bound=3 tasks.polyselect.modm=5
            tasks.polyselect.degree=3
            tasks.reconstructlog.checkdlp=false
            tasks.lim0=4000 tasks.lim1=600
            tasks.lpb0=12
            tasks.lpb1=11
            tasks.descent.lpb1=11
            tasks.descent.init_lpb=20
            tasks.sieve.mfb0=12
            tasks.sieve.mfb1=22
            tasks.sieve.lambda0=1.1
            tasks.sieve.lambda1=2.1
            tasks.qmin=700
            tasks.sieve.qrange=100
            tasks.sieve.sqside=0
            tasks.sieve.rels_wanted=1500
        )
        set(full_p30_JL_deps
            full_dependencies_base
            bwc_full_gfp
            lingen_pz sm sm_simple las_descent reconstructlog-dl
            dlpolyselect merge-dl replay-dl numbertheory_tool
            )
        if(HAVE_GMPECM)
            set(full_p30_JL_deps ${full_p30_JL_deps} descent_init_Fp)
            set(full_p30_JL_params ${full_p30_JL_params} target=2,3)
        else()
            set(full_p30_JL_params ${full_p30_JL_params} tasks.reconstructlog.checkdlp=false)
        endif()
        cado_define_test(full_p30_JL
            SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            -dlp -ell 350586976534485556208993720963
            -t 2 -dlp-no-keep
            701173953068971112417987441927
            ${full_p30_JL_params}
            TARGET_DEPENDENCIES
            ${full_p30_JL_deps}
            AVOID_CONCURRENT
            TIMEOUT 360
            )
    endif()

    # Test importing an SNFS polynomial with no skew: line
    cado_define_test(full_F7
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/cado-nfs.py
        --server-threads 2
        -s 1
        340282366920938463463374607431768211457
        ${CADO_NFS_SOURCE_DIR}/parameters/factor/parameters.F7
        tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/F7.poly
        slaves.hostnames=localhost
    TARGET_DEPENDENCIES
        full_dependencies_base
        skewness
        merge replay
        characters
        bwc_full_gf2
        sqrt
    AVOID_CONCURRENT
    TIMEOUT 300
        )
endif()
# MPI binaries are not properly launched by cado-nfs.py (see #21819)
if(NOT HAVE_MPI AND DEFINED ENV{CHECKS_EXPENSIVE})
    # we used to test separately with and without a user-supplied
    # workdir. We no longer do that, given that we rely first and
    # foremost on the shell wrapper to be a better safety net and catch
    # SIGHUP yet still clean up after itself. And we want to do that in
    # all cases, and that means that we're passing a user-supplied
    # workdir.
    cado_define_test(full_c59
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg --workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            90377629292003121684002147101760858109247336549001090677693
            --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT
        TIMEOUT 360
            )
    cado_define_test(full_c65_128
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg --workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            22381151342911794563007770625342242843893411711318155026160427871
            tasks.linalg.m=128 tasks.linalg.n=128 --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT
        TIMEOUT 360
            )
    cado_define_test(full_c59_nonlinear
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            --server-threads 2
            71641520761751435455133616475667090434063332228247871795429
            tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/c59_nonlinear.poly
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT
        TIMEOUT 360
        )
    if (BWC_GFP_ARITHMETIC_BACKENDS AND pz IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
#        this test is probably not very useful, given that the other
#        (quicker) tests cover everything covered by this test.
#        cado_define_test(full_p59
#            SCRIPT
#               ${CADO_NFS_BINARY_DIR}/cado-nfs.py
#               43341748620473677010074177283795146221310971425909898235183
#               -dlp -t 2 -s auto -dlpnokeep
#                ${CADO_NFS_SOURCE_DIR}
#            TARGET_DEPENDENCIES full_p30_dependencies
#        )
        cado_define_test(full_gfp3
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --env t
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp3.sh
                ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
            TARGET_DEPENDENCIES full_p30_dependencies
            lingen_pz sm reconstructlog-dl merge-dl replay-dl numbertheory_tool
        AVOID_CONCURRENT
        TIMEOUT 360
            )
        cado_define_test(full_gfp2
            SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg --workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py 100000000000000000039 -dlp -gfpext 2 -t 2 -ell 164354743277891 -dlp-no-keep
            TARGET_DEPENDENCIES full_p30_dependencies
            filter_galois polyselect_gfpn
        AVOID_CONCURRENT
        TIMEOUT 360
            )
# see bug https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/21767
#         cado_define_test(full_gfp2_2_2g
#             SCRIPT
#                 ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp2_2_2g.sh
#                 ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
#             TARGET_DEPENDENCIES full_c30_dependencies
#             plingen_pz sm reconstructlog-dl merge-dl replay-dl numbertheory_tool
# 	    filter_galois)
    endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library (tests STATIC test_iter.c tests_common.c)
target_link_libraries (tests ${gmp_libname})

cado_define_test(test_tests_common.c NO_DEFAULT_RUN LIBRARIES tests)
cado_divert_test(test_tests_common 1)
cado_divert_test(test_tests_common 2 -seed)
cado_divert_test(test_tests_common 3 -seed 1)
set_tests_properties(test_tests_common_3 PROPERTIES PASS_REGULAR_EXPRESSION "Using random seed=1")
cado_divert_test(test_tests_common 4 -seed a)
cado_divert_test(test_tests_common 5 -seed "")
set_tests_properties(test_tests_common_2 test_tests_common_4 test_tests_common_5 PROPERTIES WILL_FAIL 1)
cado_divert_test(test_tests_common 6 -iter 10)
set_tests_properties(test_tests_common_6 PROPERTIES PASS_REGULAR_EXPRESSION "Using 10 iterations")

cado_define_test(test_iceildiv.c)

# used by some tests.
find_program(MAGMA magma HINTS ENV MAGMA)
if(MAGMA)
    message(STATUS "magma program found, some bwc checks will run Magma code")
endif()
find_program(SAGE sage HINTS ENV SAGE)
if(SAGE)
    message(STATUS "sage program found, some bwc checks will run SageMath code")
endif()

add_subdirectory (scripts)
add_subdirectory (utils)
add_subdirectory (sieve)
add_subdirectory (polyselect)
add_subdirectory (filter)
add_subdirectory (estimate_matsize)
add_subdirectory (linalg)
if(NOT HAVE_GF2X)
add_subdirectory (gf2x)
endif()
add_subdirectory (misc)
add_subdirectory (sqrt)
