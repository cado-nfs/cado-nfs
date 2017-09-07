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
    set(IMPLICIT)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(LIBRARIES)
    set(SOURCES)
    set(SCRIPT)
    set(ARGUMENTS)
    set(TEST_NAME)
    set(WORKING_DIRECTORY)
    set(NO_DEFAULT_RUN)
    set(current IMPLICIT)
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
        elseif (x STREQUAL "TEST_NAME")
            SET(current TEST_NAME)
        elseif (x STREQUAL "WORKING_DIRECTORY")
            SET(current WORKING_DIRECTORY)
        elseif (x STREQUAL "NO_DEFAULT_RUN")
            SET(current NO_DEFAULT_RUN)
            # special, because it's a flag, really
            set(NO_DEFAULT_RUN 1)
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
    #message(STATUS "TEST_NAME is ${TEST_NAME}")
    #message(STATUS "WORKING_DIRECTORY is ${WORKING_DIRECTORY}")
    #message(STATUS "NO_DEFAULT_RUN is ${NO_DEFAULT_RUN}")
    list(LENGTH TEST_NAME nTEST_NAME)
    list(LENGTH WORKING_DIRECTORY nWORKING_DIRECTORY)
    list(LENGTH SOURCES nSOURCES)
    list(LENGTH SCRIPT nSCRIPT)
    list(LENGTH LIBRARIES nLIBRARIES)
    list(LENGTH IMPLICIT nIMPLICIT)
    list(LENGTH NO_DEFAULT_RUN nNO_DEFAULT_RUN)
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
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
            message(WARNING "too many script names, retaining only first")
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
    if(nWORKING_DIRECTORY GREATER 1)
        message(WARNING "too many test names ${WORKING_DIRECTORY}, retaining only first")
        LIST(GET WORKING_DIRECTORY 0 x)
        SET(WORKING_DIRECTORY ${x})
    endif()
    if(nTEST_NAME GREATER 1)
        message(WARNING "too many test names ${TEST_NAME}, retaining only first")
        LIST(GET TEST_NAME 0 x)
        SET(TEST_NAME ${x})
    elseif(nTEST_NAME EQUAL 0)
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
        add_test(builddep_${TEST_NAME} ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)
    endif()

    if (nSCRIPT EQUAL 0)
        set(${TEST_NAME}_EXECUTABLE ${TEST_NAME})
    else()
        set(binary ${SCRIPT})
        if (${SCRIPT} MATCHES "\\.sh$")
            set(binary env bash ${SCRIPT})
        endif()
        set(${TEST_NAME}_EXECUTABLE ${binary})
    endif()

    set(${TEST_NAME}_ARGUMENTS ${ARGUMENTS})

    if(NOT NO_DEFAULT_RUN)
        set(add_test_args NAME ${TEST_NAME} COMMAND ${${TEST_NAME}_EXECUTABLE} ${ARGUMENTS})
        if(nWORKING_DIRECTORY GREATER 0)
            list(APPEND add_test_args WORKING_DIRECTORY ${WORKING_DIRECTORY})
        endif()
        add_test(${add_test_args})
        set(${TEST_NAME}_TEST_EXISTS 1)
        if(nTARGET_DEPENDENCIES GREATER 0)
            # seems that we need set_tests_properties to take only a
            # ;-separated list here.
            set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
            string_join(all_deps ";" ${deps})
            set_tests_properties (${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
        endif()
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
    set(APPEND_ARGUMENTS)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(NO_DEFAULT_RUN)
    set(current ARGUMENTS)
    foreach(x ${ARGN})
        if (x STREQUAL "ARGUMENTS")
            SET(current ARGUMENTS)
        elseif (x STREQUAL "APPEND_ARGUMENTS")
            SET(current APPEND_ARGUMENTS)
        elseif (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        elseif (x STREQUAL "NO_DEFAULT_RUN")
            SET(current NO_DEFAULT_RUN)
            # special, because it's a flag, really
            set(NO_DEFAULT_RUN 1)
        else()
            list(APPEND ${current} ${x})
        endif()
    endforeach(x)
    # do some checking.
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
    list(LENGTH APPEND_ARGUMENTS nAPPEND_ARGUMENTS)
    list(LENGTH NO_DEFAULT_RUN nNO_DEFAULT_RUN)
    if ((nAPPEND_ARGUMENTS GREATER 0) AND (nARGUMENTS GREATER 0))
        message(FATAL_ERROR "ARGUMENTS and APPEND_ARGUMENTS are incompatible")
    endif()

    set(TEST_NAME "${TEST_BASE}_${DIVERSION}")

    if (nAPPEND_ARGUMENTS GREATER 0)
        set(ARGUMENTS ${${TEST_BASE}_ARGUMENTS} ${APPEND_ARGUMENTS})
        list(LENGTH ARGUMENTS nARGUMENTS)
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

    if(NOT NO_DEFAULT_RUN)
        add_test(${TEST_NAME} ${${TEST_BASE}_EXECUTABLE} ${ARGUMENTS})
        # prepare a ;-separated list for set_tests_properties
        # note that we do not need to make the test depend on
        # builddep_${TEST_BASE}, since we have made ${TEST_NAME}_dependencies
        # depend on ${TEST_BASE}_dependencies. It's in fact relieving, given
        # that builddep_${TEST_BASE} does not necessarily exist.
        set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
        if(${TEST_BASE}_TEST_EXISTS)
            list(APPEND deps ${TEST_BASE})
        endif()
        string_join(all_deps ";" ${deps})
        set_tests_properties (${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
    endif()
endmacro()


if(NOT MINGW)
cado_define_test(full_c30
    SCRIPT
        ${CADO_NFS_SOURCE_DIR}/cado-nfs.py 999073468111577057576445816581
        --server-threads 2
    TARGET_DEPENDENCIES
        polyselect polyselect_ropt
        freerel makefb
        las
        dup1 dup2 purge merge replay
        characters
        bwc_full
        sqrt)
    if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        # see linalg/bwc/CMakeLists.txt for an explanation
        cado_define_test(full_p30
            SCRIPT ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
            -dlp -ell 445257213376519703258944813
            -t 2 -dlp-no-keep
            target=92800609832959449330691138186
            3037544709654617415632521514287
            TARGET_DEPENDENCIES full_c30_dependencies debug_renumber
            plingen_pz sm sm_simple las_descent reconstructlog-dl
            merge-dl replay-dl badideals)
    endif()
endif()
if(DEFINED ENV{CHECKS_EXPENSIVE})
    cado_define_test(full_c59
        SCRIPT
            ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
            90377629292003121684002147101760858109247336549001090677693
            --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
            )
    cado_define_test(full_c59workdir
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg workdir
            ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
	    90377629292003121684002147101760858109247336549001090677693
        TARGET_DEPENDENCIES
	    full_c30_dependencies
            )
    cado_define_test(full_c65_128
        SCRIPT
            ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
            22381151342911794563007770625342242843893411711318155026160427871
            tasks.linalg.m=128 tasks.linalg.n=128 --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
            )
    cado_define_test(nonlinear_c59
        SCRIPT
            ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
            --server-threads 2
            71641520761751435455133616475667090434063332228247871795429
            tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/c59_nonlinear.poly
        TARGET_DEPENDENCIES
            full_c30_dependencies
        )
        if (BWC_GFP_ARITHMETIC_BACKENDS)
#        this test is probably not very useful, given that the other
#        (quicker) tests cover everything covered by this test.
#        cado_define_test(full_p59
#            SCRIPT
#               ${CADO_NFS_SOURCE_DIR}/cado-nfs.py
#               43341748620473677010074177283795146221310971425909898235183
#               -dlp -t 2 -s auto -dlpnokeep
#                ${CADO_NFS_SOURCE_DIR}
#            TARGET_DEPENDENCIES full_c30_dependencies debug_renumber
#            plingen_pz sm reconstructlog-dl merge-dl replay-dl badideals
#        )
        cado_define_test(test_gfp3
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp3.sh
                ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
            TARGET_DEPENDENCIES full_c30_dependencies debug_renumber
            plingen_pz sm reconstructlog-dl merge-dl replay-dl badideals)
        cado_define_test(test_gfp2
            SCRIPT ${CADO_NFS_SOURCE_DIR}/cado-nfs.py 100000000000000000039 -dlp -gfpext 2 -t 2 -ell 164354743277891 -dlp-no-keep
            TARGET_DEPENDENCIES full_c30_dependencies debug_renumber
            plingen_pz sm reconstructlog-dl merge-dl replay-dl badideals
            filter_galois polyselect_gfpn)
        cado_define_test(test_gfp2_2_2g
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp2_2_2g.sh
                ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
            TARGET_DEPENDENCIES full_c30_dependencies debug_renumber
            plingen_pz sm reconstructlog-dl merge-dl replay-dl badideals
	    filter_galois)
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

add_subdirectory (scripts)
add_subdirectory (utils)
add_subdirectory (sieve)
add_subdirectory (polyselect)
add_subdirectory (numbertheory)
add_subdirectory (filter)
add_subdirectory (linalg)
add_subdirectory (gf2x)
add_subdirectory (misc)
add_subdirectory (sqrt)
#add_subdirectory (nfs-hd)
