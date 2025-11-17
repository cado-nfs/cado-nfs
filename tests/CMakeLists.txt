add_custom_target(all_test_dependencies)

# used by some tests.
find_program(MAGMA magma HINTS ENV MAGMA)
if(MAGMA)
    message(STATUS "magma program found, some bwc checks will run Magma code")
endif()
find_program(SAGE sage HINTS ENV SAGE)
if(SAGE)
    # ugly workaround for $30126
    if (EXISTS /proc/cpuinfo AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        file(STRINGS /proc/cpuinfo has_adx_flags REGEX "^flags.* adx ")
        if(NOT has_adx_flags)
            message(STATUS "sage program found, but this linux machine does not have adx isntructions, so we can't use sagemath. See #30126")
            set(SAGE NOTFOUND)
        endif()
    endif()
endif()
if(SAGE)
    message(STATUS "sage program found, some bwc checks will run SageMath code")
endif()
add_custom_target(all_magma_test_dependencies)
add_custom_target(all_sagemath_test_dependencies)


include(cado_test.cmake)
include(cado_define_test.cmake)
include(cado_divert_test.cmake)
include(predefined_primes.cmake)




add_custom_target(full_dependencies_base)
add_dependencies(full_dependencies_base
        polyselect polyselect_ropt polyselect3
        freerel makefb
        las
        dup1 dup2 purge
        )

set(AVOID_CONCURRENT_UNDER_COVERAGE)
if(COVERAGE_TEST)
    set(AVOID_CONCURRENT_UNDER_COVERAGE AVOID_CONCURRENT 16)
endif()

# MPI binaries are not properly launched by cado-nfs.py (see #21819)
if(NOT MINGW AND NOT HAVE_MPI)

cado_define_test(full_c30
    SCRIPT
    ${CADO_NFS_BINARY_DIR}/cado-nfs.py 999073468111577057576445816581
        slaves.downloadretry=1
        --server-threads 2
    TARGET_DEPENDENCIES
        full_dependencies_base
        merge replay
        characters
        bwc_full_gf2
        sqrt
    AVOID_CONCURRENT 16
    ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
    TIMEOUT 360
        )
    cado_divert_test(full_c30 pq2
        ARGUMENTS
        1000000008500000023430000020691
        slaves.downloadretry=1
        --server-threads 2
        SHA1_ON_REGEXP_LINES "^[0-9 ]*$"
        EXPECT_SHA1 a56b28a490e2a8f240f6102c695456e59d119a9e
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
            slaves.downloadretry=1
            target=92800609832959449330691138186
            3037544709654617415632521514287
            TARGET_DEPENDENCIES full_dependencies_base
            bwc_full_gfp
            lingen_pz sm sm_simple las_descent reconstructlog-dl
            merge-dl replay-dl numbertheory_tool
            AVOID_CONCURRENT 16
            ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
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
            tasks.sieve.qrange=400
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
            slaves.downloadretry=1
            701173953068971112417987441927
            ${full_p30_JL_params}
            TARGET_DEPENDENCIES
            ${full_p30_JL_deps}
            AVOID_CONCURRENT 16
            ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
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
        slaves.downloadretry=1
    TARGET_DEPENDENCIES
        full_dependencies_base
        skewness
        merge replay
        characters
        bwc_full_gf2
        sqrt
        score
    AVOID_CONCURRENT 16
    ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
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
            slaves.downloadretry=1
            --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT 16
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
        TIMEOUT 360
            )
    cado_define_test(full_c65_128
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg --workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            22381151342911794563007770625342242843893411711318155026160427871
            slaves.downloadretry=1
            tasks.linalg.m=128 tasks.linalg.n=128 --server-threads 2
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT 16
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
        TIMEOUT 360
            )
    cado_define_test(full_c59_nonlinear
        SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py
            --server-threads 2
            71641520761751435455133616475667090434063332228247871795429
            slaves.downloadretry=1
            tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/c59_nonlinear.poly
        TARGET_DEPENDENCIES
            full_c30_dependencies
        AVOID_CONCURRENT 16
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
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
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp3.sh
            PROVIDE_TEMPORARY_WDIR
            TARGET_DEPENDENCIES full_p30_dependencies
            lingen_pz sm reconstructlog-dl merge-dl replay-dl numbertheory_tool
        AVOID_CONCURRENT 16
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
        TIMEOUT 360
            )
        cado_define_test(full_gfp2
            SCRIPT
            ${CMAKE_CURRENT_SOURCE_DIR}/provide-wdir.sh  --arg --workdir
            ${CADO_NFS_BINARY_DIR}/cado-nfs.py 100000000000000000039
            slaves.downloadretry=1
            -dlp -gfpext 2 -t 2 -ell 164354743277891 -dlp-no-keep
            TARGET_DEPENDENCIES full_p30_dependencies
            filter_galois polyselect_gfpn
        AVOID_CONCURRENT 16
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
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

# let's make this one an exception. As regards tests, I think it's okay
# to have it static
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

if (HAVE_OPENMP)
    add_executable(omp_get_max_threads omp_get_max_threads.c)
    mark_targets_as_openmp(omp_get_max_threads)
    # virtually all targets in cado-nfs depend on utils, and we want to
    # be sure that in all situations, the check environment will have
    # this binary available after the build step. Having utils depend on
    # it sounds like a semi-reasonable thing to do.
    add_dependencies(utils omp_get_max_threads)
endif()

if (HAVE_HWLOC)
    add_executable(hwloc_cado_helper hwloc_cado_helper.c)
    target_link_libraries(hwloc_cado_helper hwloc)
    mark_targets_as_openmp(hwloc_cado_helper)
    # same trick
    add_dependencies(utils hwloc_cado_helper)
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
