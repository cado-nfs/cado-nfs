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

if(TIMEOUT_SCALE)
    message(STATUS "TIMEOUT_SCALE is set to ${TIMEOUT_SCALE}")
endif()

set(CADO_NFS_TEST_KEYWORDS_LAST_WINS TIMEOUT AVOID_CONCURRENT)
set(CADO_NFS_TEST_KEYWORDS_FIRST_WINS_NOISY WORKING_DIRECTORY TEST_NAME)

function(cado_set_test_timeout TEST_NAME TIMEOUT)
    if(TIMEOUT_SCALE)
        MATH(EXPR scaled_timeout "${TIMEOUT}*${TIMEOUT_SCALE}")
        set_property(TEST ${TEST_NAME} PROPERTY TIMEOUT
            ${scaled_timeout})
    else()
        set_property(TEST ${TEST_NAME} PROPERTY TIMEOUT ${TIMEOUT})
    endif()
endfunction()


set(CADO_NFS_TEST_KEYWORDS_USUAL
    # these args all thare the same kind of treatment in cado-nfs tests:
    # only one argument are accepted, and they get propagated to diverted
    # tests.
    EXPECT_SHA1
    SHA1_ON_REGEXP_LINES
    STDIN
    TIMEOUT
    WORKING_DIRECTORY
)

set(CADO_NFS_TEST_KEYWORDS_MULTI
    # pretty much the same, but multi args are ok (and expected)
    ENVIRONMENT
    PRECOMMAND
    ARGUMENTS
    FIXTURES_SETUP
    FIXTURES_CLEANUP
    FIXTURES_REQUIRED
)

set(CADO_NFS_TEST_KEYWORDS_FLAGS
    NO_DEFAULT_RUN
    AVOID_CONCURRENT
    EXPECT_FAIL
    PROVIDE_TEMPORARY_WDIR
)

set(CADO_NFS_TEST_KEYWORDS_DO_NOT_PROPAGATE
    NO_DEFAULT_RUN
    FIXTURES_SETUP
    FIXTURES_CLEANUP
)

set(CADO_NFS_TEST_KEYWORDS_SPECIAL_PROPAGATE_RULES
    ENVIRONMENT
    PRECOMMAND
    ARGUMENTS
    EXECUTABLE
    TARGET_DEPENDENCIES
)

list(APPEND CADO_NFS_TEST_KEYWORDS_FIRST_WINS_NOISY
    ${CADO_NFS_TEST_KEYWORDS_USUAL})


macro(cado_nfs_test_propagate_ARGUMENTS)
    # message(STATUS "propagating ARGUMENTS for ${TEST_NAME} ; ${TEST_BASE}_EXECUTABLE=${${TEST_BASE}_EXECUTABLE}; ${TEST_BASE}_ARGUMENTS=${${TEST_BASE}_ARGUMENTS}; APPEND_ARGUMENTS=${APPEND_ARGUMENTS}")
    if (NOT nARGUMENTS GREATER 0)
        set(ARGUMENTS ${${TEST_BASE}_ARGUMENTS} ${APPEND_ARGUMENTS})
        list(LENGTH ARGUMENTS nARGUMENTS)
    elseif (nAPPEND_ARGUMENTS)
        message(FATAL_ERROR "ARGUMENTS and APPEND_ARGUMENTS are incompatible")
    endif()
endmacro()

macro(cado_nfs_test_propagate_ENVIRONMENT)
    if (NOT nENVIRONMENT GREATER 0)
        set(ENVIRONMENT ${${TEST_BASE}_ENVIRONMENT} ${APPEND_ENVIRONMENT})
        list(LENGTH ENVIRONMENT nENVIRONMENT)
    elseif (nAPPEND_ENVIRONMENT)
        message(FATAL_ERROR "ENVIRONMENT and APPEND_ENVIRONMENT are incompatible")
    endif()
endmacro()

macro(cado_nfs_test_propagate_PRECOMMAND)
    if (NOT nPRECOMMAND GREATER 0)
        set(PRECOMMAND ${PREPEND_PRECOMMAND} ${${TEST_BASE}_PRECOMMAND})
        list(LENGTH PRECOMMAND nPRECOMMAND)
    elseif (nPREPEND_PRECOMMAND)
        message(FATAL_ERROR "PRECOMMAND and PREPEND_PRECOMMAND are incompatible")
    endif()
endmacro()

macro(cado_nfs_test_propagate_EXECUTABLE)
    set(EXECUTABLE ${${TEST_BASE}_EXECUTABLE})
endmacro()

macro(cado_nfs_test_propagate_TARGET_DEPENDENCIES)
    # we want to define builddep_X unconditionally, depending at least on
    # the dependencies of the base test.
    # message(STATUS "-------------- propagating TARGET_DEPENDENCIES from ${TEST_BASE} to ${TEST_NAME}")
    list(APPEND TARGET_DEPENDENCIES ${TEST_BASE}_dependencies)
endmacro()


macro(cado_nfs_test_process_EXPECT_SHA1)
    list(APPEND wrapper_args --expect-sha1 ${EXPECT_SHA1})
endmacro()

macro(cado_nfs_test_process_SHA1_ON_REGEXP_LINES)
    # if(NOT EXPECT_SHA1)
    #     message(FATAL_ERROR "SHA1_ON_REGEXP_LINES ignored if EXPECT_SHA1 is not set")
    # endif()
    list(APPEND wrapper_args --filter-output "${SHA1_ON_REGEXP_LINES}")
endmacro()

macro(cado_nfs_test_process_STDIN)
    list(APPEND wrapper_args --stdin ${STDIN})
endmacro()

function(cado_nfs_test_process_AVOID_CONCURRENT)
		if(AVOID_CONCURRENT EQUAL 1)
			message(STATUS "Test ${TEST_NAME} runs serially")
			set_property(TEST ${TEST_NAME} PROPERTY RUN_SERIAL)
		else()
			set_property(TEST ${TEST_NAME} PROPERTY PROCESSORS ${AVOID_CONCURRENT})
		endif()
endfunction()

function(cado_nfs_test_process_TIMEOUT)
    cado_set_test_timeout(${TEST_NAME} ${TIMEOUT})
endfunction()

macro(cado_nfs_test_process_EXPECT_FAIL)
    # This gets called once before add_test, and once after
    if(${TEST_NAME}_TEST_EXISTS)
        set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL 1)
    else()
        LIST(APPEND wrapper_args --abort-to-fail)
    endif()
endmacro()

macro(cado_nfs_test_process_PROVIDE_TEMPORARY_WDIR)
    set(provide_wdir ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env wdir)
endmacro()

macro(cado_nfs_test_process_keyword ARG)
    if(n${ARG} GREATER 0)
        cmake_language(CALL cado_nfs_test_process_${ARG} ${ARGN})
    endif()
endmacro()

macro(cado_nfs_test_compute_lengths)
    foreach(x
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
            TARGET_DEPENDENCIES
            ${ARGN}
        )
        list(LENGTH ${x} n${x})
    endforeach()
endmacro()

macro(cado_nfs_test_propagate_all)
    foreach(x
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
            EXECUTABLE
            TARGET_DEPENDENCIES
        )
        if (${x} IN_LIST CADO_NFS_TEST_KEYWORDS_DO_NOT_PROPAGATE)
        elseif (${x} IN_LIST
                CADO_NFS_TEST_KEYWORDS_SPECIAL_PROPAGATE_RULES)
            # message(STATUS "### ${TEST_BASE} --> ${TEST_NAME} --> cado_nfs_test_propagate_${x}")
            cmake_language(CALL cado_nfs_test_propagate_${x})
            list(LENGTH ${x} n${x})
        else()
            if (NOT n${x} GREATER 0 AND ${TEST_BASE}_${x})
                # message(STATUS "propagating ${x} from ${TEST_BASE} to ${TEST_NAME}")
                set(${x} ${${TEST_BASE}_${x}})
                list(LENGTH ${x} n${x})
            endif()
        endif()
    endforeach()
endmacro()

macro(cado_nfs_test_process_environment)
    set(test_env)
    foreach(v
            CADO_NFS_SOURCE_DIR
            CADO_NFS_BINARY_DIR
            CMAKE_CURRENT_SOURCE_DIR
            CMAKE_CURRENT_BINARY_DIR
            PROJECT_BINARY_DIR
            MAGMA
            SAGE
        )
        if(${v})
            list(APPEND test_env ${v}=${${v}})
        endif()
    endforeach()
    list(APPEND test_env ${ENVIRONMENT})
endmacro()

function(cado_nfs_test_print_debug)
    foreach(x
            EXECUTABLE
            IMPLICIT
            TEST_DEPENDENCIES
            TARGET_DEPENDENCIES
            LIBRARIES
            SOURCES
            SCRIPT
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
        )
        LIST(LENGTH ${x} n)
        # if (${n})
        #     message(STATUS "${TEST_NAME}: ${x}=${${x}}")
        # endif()
    endforeach()
endfunction()

macro(cado_nfs_test_forward_to_test_property)
    foreach(p FIXTURES_CLEANUP FIXTURES_SETUP FIXTURES_REQUIRED)
        if(n${p})
            set_property(TEST ${TEST_NAME} PROPERTY ${p} ${${p}})
        endif()
    endforeach()
endmacro()

macro(cado_nfs_test_check_argument_rules)
    foreach(x ${CADO_NFS_TEST_KEYWORDS_FIRST_WINS_NOISY})
        if(n${x} GREATER 1)
            message(WARNING "too many ${x} directives for ${TEST_NAME} ${${x}}, retaining only first")
            LIST(GET ${x} 0 y)
            SET(${x} ${y})
        endif()
    endforeach()

    foreach(x ${CADO_NFS_TEST_KEYWORDS_FLAGS})
        if(n${x} GREATER 1)
            message(FATAL_ERROR "discarded arguments after ${x} for test ${TEST_NAME}")
            set(${x} 1)
        endif()
    endforeach()
endmacro()


macro(cado_epilogue_create_test)
    add_custom_target(${TEST_NAME}_dependencies)
    add_dependencies(all_test_dependencies ${TEST_NAME}_dependencies)

    if(nTARGET_DEPENDENCIES GREATER 0)
        add_dependencies(${TEST_NAME}_dependencies ${TARGET_DEPENDENCIES})
        add_test(builddep_${TEST_NAME}
            ${CMAKE_COMMAND}
            --build ${CMAKE_BINARY_DIR}
            --target ${TEST_NAME}_dependencies)
        list(APPEND TEST_DEPENDENCIES builddep_${TEST_NAME})
    endif()

    # save for use by diverted tests
    foreach(x
            ${CADO_NFS_TEST_KEYWORDS_USUAL}
            ${CADO_NFS_TEST_KEYWORDS_MULTI}
            ${CADO_NFS_TEST_KEYWORDS_FLAGS}
            EXECUTABLE
        )
        set(${TEST_NAME}_${x} ${${x}} PARENT_SCOPE)
    endforeach()

    set(${TEST_NAME}_TEST_DEPENDENCIES ${TEST_DEPENDENCIES})
    set(${TEST_NAME}_TEST_DEPENDENCIES ${TEST_DEPENDENCIES} PARENT_SCOPE)

    if(NOT NO_DEFAULT_RUN)
        cado_nfs_test_process_environment()
        MATH(EXPR test_index "${test_index}+1")

        set(wrapper_args)
        set(provide_wdir)

        cado_nfs_test_process_keyword(EXPECT_SHA1)
        cado_nfs_test_process_keyword(SHA1_ON_REGEXP_LINES)
        cado_nfs_test_process_keyword(STDIN)
        cado_nfs_test_process_keyword(EXPECT_FAIL)
        cado_nfs_test_process_keyword(PROVIDE_TEMPORARY_WDIR)

        set(add_test_args
            NAME ${TEST_NAME}
            COMMAND
                # ${CMAKE_CTEST_COMMAND}
                # --build-and-test    "${CMAKE_SOURCE_DIR}"
                #                     "${CMAKE_BINARY_DIR}"
                # --build-generator   "${CMAKE_GENERATOR}"
                # --build-makeprogram "${CMAKE_MAKE_PROGRAM}"
                # --build-target ${TEST_NAME}_dependencies
                # --test-command
                env ${test_env}
                ${provide_wdir}
                ${PRECOMMAND}
                ${CADO_NFS_SOURCE_DIR}/tests/test.sh
                ${wrapper_args}
                --
                ${EXECUTABLE}
                ${ARGUMENTS}
            )
        if(nWORKING_DIRECTORY GREATER 0)
            list(APPEND add_test_args WORKING_DIRECTORY ${WORKING_DIRECTORY})
        endif()
        add_test(${add_test_args})

        cado_nfs_test_process_keyword(AVOID_CONCURRENT)
        cado_nfs_test_process_keyword(TIMEOUT)

        set(${TEST_NAME}_TEST_EXISTS 1)
        set(${TEST_NAME}_TEST_EXISTS 1 PARENT_SCOPE)

        set_property(TEST ${TEST_NAME}
            PROPERTY DEPENDS ${${TEST_NAME}_TEST_DEPENDENCIES})

        cado_nfs_test_forward_to_test_property()

        set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 125)
        cado_nfs_test_process_keyword(EXPECT_FAIL)
    endif()
endmacro()
