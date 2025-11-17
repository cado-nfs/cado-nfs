# NOTE: most GF(p) checks here depend on lingen_pz ; this is artifical,
# and caused by the fact that bwc.pl is for the moment unable to call
# anything else than lingen_pz

# some tests below use the bwc include files.
include_directories(${PROJECT_SOURCE_DIR}/linalg)
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc)
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/linalg/bblas)

# add_bwc_test is used both for mod2 and modp tests, and therefore we
# depend on bwc_full, not the finer grained bwc_full_gf2 and bwc_full_gfp
macro(add_bwc_test nm script)
    set(mpiconfs ${${nm}_mpi_tests})
    if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
        set(mpiconfs)
    endif()
    foreach(mpiconf nompi ${mpiconfs})
        set(suffix _${mpiconf})
        if(suffix STREQUAL _nompi)
            set(suffix)
            set(default_timeout 50)
        else()
            # Yes it is quite weird, but when we run the test suite with
            # mpi, we are often overcommitting the resource that we're
            # using, and it ends up degrading performance quite a bit.
            set(default_timeout 180)
        endif()
        # place the TIMEOUT stanza early on, so that ARGN can
        # override it.
        cado_define_test(test_bwc_${nm}${suffix}
            AVOID_CONCURRENT 16
            ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16
            TARGET_DEPENDENCIES bwc_full random_matrix
            # We have cases where _many_ bwc tests end up taking more
            # than the default 20 seconds. Maybe it's better to add
            # the cutoff.
            TIMEOUT ${default_timeout}
            SCRIPT
            ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg mats --other --arg wdir --
            ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
            ${CADO_NFS_SOURCE_DIR}/tests/linalg/bwc/${script}
            bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc pre_wipe=1 ${ARGN}
        )
        # these verbose_flags don't play with the idea of adding a
        # divert_bwc_test macro.
        # -- verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch,^bwc-loading-mksol-files,^bwc-cache-major-info,^bwc-cpubinding,^perl-checks,^perl-sections
    endforeach()
endmacro()

macro(divert_bwc_test nm variant)
    set(mpiconfs ${${nm}_mpi_tests})
    if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
        set(mpiconfs)
    endif()
    foreach(mpiconf nompi ${mpiconfs})
        set(suffix _${mpiconf})
        if(suffix STREQUAL _nompi)
            set(suffix)
        endif()
        cado_divert_test(test_bwc_${nm}${suffix} ${variant} ${ARGN})
    endforeach()
endmacro()

if(HAVE_PERL_DIGEST_MD5)
    macro(add_lingen_binary_test nm)
        set(mpiconfs ${lingen_binary_${nm}_mpi_tests})
        if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
            set(mpiconfs)
        endif()
        foreach(mpiconf nompi ${mpiconfs})
            set(suffix _${mpiconf})
            if(suffix STREQUAL _nompi)
                set(suffix)
            endif()
            # place the TIMEOUT stanza early on, so that ARGN can
            # override it.
            #
            # AVOID_CONCURRENT 4 is very extreme. As it turns out, this
            # test is very sensitive to the load average of the host (of
            # the container). The other options are the (even more
            # extreme) RUN_SERIAL, or possibly refuse to run the test if
            # the load average is too large (but compared to what, it's
            # not completely clear).
            cado_define_test(test_bwc_lingen_binary_${nm}${suffix}
                TARGET_DEPENDENCIES lingen_b64
                AVOID_CONCURRENT 4
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
                TIMEOUT 60
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-lingen.sh
                ${ARGN}
                TARGET_DEPENDENCIES lingen_b64
                )
        endforeach()
    endmacro()

    # the sha1sums here are from the time were we were transposing the
    # output. Now we no longer do that.
    # add_lingen_binary_test(test0 1dea0d8879d125aa280e26b8573f932f6bccea75,d0e758e619c8e81b3e28bbd4fe69c32826e9a4bd 64 64 160 1 lingen_mpi_threshold=10)

    # todo presumably, on 32 bit, we get something different for the sha1sum.

    # the sha1sum 8938... is for the old binary. The newer one gets d0e7...
    # and it seems to be a question of varying normalization.

    # the block-based version gets 5917486028cd96a5d0c55c7caef89cc48b70d1b1

    # square2 is too large for this tiny test
    set(lingen_binary_test0_mpi_tests mpi_square1)
    add_lingen_binary_test(test0
        expect_sha1_F=d0e758e619c8e81b3e28bbd4fe69c32826e9a4bd,893801e0727d840353a35b2623430fb758e4d379,5917486028cd96a5d0c55c7caef89cc48b70d1b1
        m=64 n=64 sequence_length=480 seed=1
        tuning_thresholds=recursive:64,notiming:0,collective:128
        TIMEOUT 300
        )

    # We have significant difficulties with tuning for the binary case.
    # This test takes absurdly long.
    # 
    # 3fd4e3439ad7f76145329bb631c78ccfbf939a90 is the new hash obtained
    # by the block-based lingen basecase
    if(DEFINED ENV{CHECKS_EXPENSIVE})
        # square2 is too large for this tiny test
        set(lingen_binary_test2_mpi_tests mpi_square1)
        add_lingen_binary_test(test2
            expect_sha1_F=fefcefb595751b7120d593d6f10f0a5bb5466dfb,3fd4e3439ad7f76145329bb631c78ccfbf939a90
            m=128 n=128 sequence_length=480 seed=1
            tuning_thresholds=recursive:128,notiming:0,collective:128
            )
    endif()

    # The test below fails on 32bit. This should be investigated.
    # add_lingen_binary_test(test1
    # expect_sha1_F=30b90451fb828f08caa160f96c585c4b1c0951a8 m=64 n=64 sequence_length=200 seed=1  -- lingen_mpi_threshold=10)
endif()

# Now we only ever test p1, p3, and pz.
set(BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS)
foreach(layer ${BWC_GFP_ARITHMETIC_BACKENDS})
    parse_arith_layer(${layer} dummy)
    set(layer ${canonical_arith_layer})

    string(REGEX MATCH "^(p1|p3|pz)$" t "${layer}")
    if(t)
        list(APPEND BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS ${t})
    else()
        message(STATUS "Note: tests for ${layer} are skipped, because we restrict the routine tests")
    endif()
endforeach()
        
macro(create_magma_sage_diversions t)
    # the magma code doesn't like 32-bits
    if(MAGMA AND ULONG_BITS EQUAL 64)
        divert_bwc_test(${t} with_magma APPEND_ARGUMENTS magma=${MAGMA}
            TIMEOUT 180
        )
        add_dependencies(all_magma_test_dependencies test_bwc_${t}_with_magma_dependencies)
    endif()
    if(SAGE)
        divert_bwc_test(${t} with_sagemath APPEND_ARGUMENTS wordsize=${ULONG_BITS} sage=${SAGE}
            TIMEOUT 180
        )
        add_dependencies(all_sagemath_test_dependencies test_bwc_${t}_with_sagemath_dependencies)
    endif()
endmacro()

foreach(layer ${BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS})
    set(extra)
    string(REGEX MATCH "p([0-9]+)" t "${layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
        set(p "${prime_${width}}")
    elseif(layer STREQUAL "pz")
        set(t TRUE)
        set(width "z")
        set(extra "lingen_program=lingen_pz")
        set(p "${prime_${width}}")
    else()
        message(FATAL_ERROR "bad arithmetic backend ${layer}")
    endif()
    if(layer STREQUAL "pz")
    elseif("pz" IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
        set(t modp_inhomogeneous_${layer})
        if(layer STREQUAL "p1" OR layer STREQUAL "p13")
            set(${t}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
        endif()

        add_bwc_test(${t}
            bwc-ptrace.sh prime=${p} random_matrix_size=500 nrhs=2 m=6 n=4 seed=1 ${extra}
            TARGET_DEPENDENCIES
            lingen_pz lingen_${layer})
        cado_set_test_timeout(test_bwc_${t} 180)

        if(layer STREQUAL "p1")
            if(DEFINED ENV{CHECKS_EXPENSIVE})
                # it's quick, but let's avoid redundant tests
                set(t modp_inhomogeneous_minimal_3rhs_mn6)
                add_bwc_test(${t} bwc-ptrace.sh prime=1009
                    random_matrix_size=50 nrhs=3 m=6 n=6 seed=1
                    TARGET_DEPENDENCIES
                    lingen_pz lingen_${layer}
                )
                create_magma_sage_diversions(${t})
            endif()

            set(t modp_inhomogeneous_minimal_2rhs_mn4)
            set(${t}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
            add_bwc_test(${t} bwc-ptrace.sh
                prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1
                TARGET_DEPENDENCIES
                lingen_pz lingen_${layer})

            set(t modp_inhomogeneous_left_minimal_2rhs_mn4)
            set(${t}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
            add_bwc_test(${t} bwc-ptrace.sh
                nullspace=LEFT
                prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1
                TARGET_DEPENDENCIES
                lingen_pz lingen_${layer})
                create_magma_sage_diversions(${t})

            if(DEFINED ENV{CHECKS_EXPENSIVE})
                set(t modp_inhomogeneous_minimal_2rhs_mn4_expensive)
                add_bwc_test(${t} bwc-ptrace.sh
                    prime=1237940039285380274899124357
                    random_matrix_size=200 nrhs=2 m=4 n=4 seed=1
                    TARGET_DEPENDENCIES
                    lingen_pz lingen_${layer})
                create_magma_sage_diversions(${t})
            endif()

            if(DEFINED ENV{CHECKS_EXPENSIVE})
                set(t modp_inhomogeneous_minimal_1rhs_mn1)
                add_bwc_test(${t} bwc-ptrace.sh
                    prime=1009 random_matrix_size=50 nrhs=1 m=1 n=1 seed=1
                TARGET_DEPENDENCIES
                    lingen_pz lingen_${layer})
                create_magma_sage_diversions(${t})
            endif()

            # This one is a bit of a lie, to be honest. Depending on the
            # seed, it sometimes crashes (it is wordsize dependent).
            # seed=8 crashes on 64-bit. seed=1 crashes on 32-bit. Anyway,
            # it's fixable. The only thing is that we reach a situation
            # where transpose(x), transpose(x)*M, and the next iterates,
            # do not constrain the resulting vector to a dimension 0
            # subspace. So what we expect to be zero is actually of
            # dimension maybe 1 in that case. This is detected in
            # gather.cpp, and yields an error message.
            set(t modp_homogeneous_minimal_mn4)
            add_bwc_test(${t} bwc-ptrace.sh
                prime=1009 random_matrix_size=50 m=4 n=4 seed=2
                random_matrix_minkernel=1
                TARGET_DEPENDENCIES
                lingen_pz lingen_${layer})
            create_magma_sage_diversions(${t})

            ## this one isn't really useful
            # set(t modp_inhomogeneous_minimal_2rhs_mn2)
            # add_bwc_test(${t} bwc-ptrace.sh
            #   prime=1009 random_matrix_size=50 nrhs=2 m=2 n=2
            #   seed=1
            #   TARGET_DEPENDENCIES
            #   lingen_pz lingen_${layer})
            # create_magma_sage_diversions(${t})
        endif()

        set(t modp_secure_continuation_${layer})
        if(layer STREQUAL "p1" OR layer STREQUAL "p13")
            set(${t}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
        endif()
        add_bwc_test(${t} bwc-test-secure-continuation.sh prime=${p} random_matrix_size=500 nrhs=4 m=8 n=8 seed=1 interval=32 ${extra})
    endif()

    # Also add some lingen tests if it so occurs that the current
    # backend can process them
    # Each test has the following arguments
    # lingen test name
    # reference sha1 sum for generator
    # m
    # n
    # length of sequence to compute
    # p
    # seed for random generation (entries are (2^k*(seed/1000) mod q) + (3^k*(seed%1000) mod q) mod p, with q=2^20-3).
    # arguments to be passed to lingen (including mpi arguments)

    # Note that the mpi and lingen_mpi_threshold arguments here are used
    # only in MPI context, and ignored otherwise (that is, the test runs
    # without mpi). This is not the same for lingen_mpi_threshold, which
    # in reality is a bit of a misnomer. This flag is one of the means by
    # which the automatic tuning of lingen can be forced to use a given
    # threshold, and that also applies to the non-mpi case.  In MPI
    # context, this triggers the mpi portions of the code (actually the
    # test gets run twice, once with the threads, once without. Both have
    # to succeed, of course).

    # Note that we're only enabling the p1 tests because the checksums
    # below were obtained for p=1009 specifically.
    if(layer STREQUAL "p1")
        # The --ascii here is annoying. It gets in the way of our
        # automatic guess of the input sequence length.
        #
        # We should be able to drop it. However we have an
        # annoying situation, since the sha1sum is computed on the
        # generator as it is written on disk, of course, and that follows
        # --ascii too.
        set(mpiconfs mpi_square1 mpi_square2)
        if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
            set(mpiconfs)
        endif()
        foreach(mpiconf nompi ${mpiconfs})
            set(suffix _${mpiconf})
            if(suffix STREQUAL _nompi)
                set(suffix)
            endif()
            cado_define_test(test_bwc_lingen_${layer}_test0${suffix}
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                2b6428178054fc0d3931dcdc9257ccf3bdea62f9 4 2 36 1009 17042
                mpi_magic=${mpiconf} thr=2x2 --ascii
                tuning_thresholds=collective:64,recursive:32,notiming:0
                ${extra}
                TARGET_DEPENDENCIES lingen_${layer}
                AVOID_CONCURRENT 4
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
                TIMEOUT 60
                )
            cado_divert_test(test_bwc_lingen_${layer}_test0${suffix} frag
                APPEND_ARGUMENTS io_block_size=1024)
            cado_define_test(test_bwc_lingen_${layer}_test1${suffix}
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
                thr=2x2 --ascii
                tuning_thresholds=collective:64,recursive:32,notiming:0
                ${extra}
                TARGET_DEPENDENCIES lingen_${layer}
                AVOID_CONCURRENT 4
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
                TIMEOUT 60
                )
            cado_define_test(test_bwc_lingen_${layer}_test2${suffix}
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                b448493e04299d8d54be6cc1273c43293520ac9e 3 3 1000 1009 1111
                thr=3x3 --ascii
                tuning_thresholds=collective:64,recursive:32,notiming:0
                ${extra}
                TARGET_DEPENDENCIES lingen_${layer}
                AVOID_CONCURRENT 9
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=9 OMP_THREAD_LIMIT=9
                TIMEOUT 60
                )
            # We really _can't_ have ascii with the checkpoint test. ascii
            # modifies the anticipation of the schedule, which in turns means
            # that the list of cached schedules is modified in the course of
            # the execution, and that eventually goes to the checkpoint file
            # as well, wreaking obvious havoc.
            cado_define_test(test_bwc_lingen_${layer}_test2_checkpoints${suffix}
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-lingen-checkpoints.sh
                c82b888f2ba844f6db1e30e091444d1023052f87 3 3 1000 1009 1111
                thr=3x3
                checkpoint_threshold=40
                tuning_thresholds=notiming:0,recursive:100,collective:400
                ${extra}
                TARGET_DEPENDENCIES
                lingen_${layer}
                lingen_verify_checkpoints_${layer}
                AVOID_CONCURRENT 9
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=9 OMP_THREAD_LIMIT=9
                TIMEOUT 180
                )
        endforeach()
    endif()

    if(layer STREQUAL "pz" OR (ULONG_BITS EQUAL 64 AND layer STREQUAL "p3") OR (ULONG_BITS EQUAL 32 AND layer STREQUAL "p6"))
        set(mpiconfs mpi_square1 mpi_square2)
        if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
            set(mpiconfs)
        endif()
        foreach(mpiconf nompi ${mpiconfs})
            set(suffix _${mpiconf})
            if(suffix STREQUAL _nompi)
                set(suffix)
            endif()
            cado_define_test(test_bwc_lingen_${layer}_test4${suffix}
                PROVIDE_TEMPORARY_WDIR
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                # For some reason, the 32-bit test returns an equivalent,
                # but slightly different matrix. We should really
                # normalize the output, but presently we don't. Let's
                # take that digest as valid.
                00a481d0e0e2ee7403fa234b7f8b818dc2ec8db4,2ac9961e5a13652247660678b5fabd2782ab740a 3 3 1000
                1532495540865888858358347027150309183618739122183602191 1111
                thr=3x3
                tuning_thresholds=notiming:0,recursive:160,collective:320
                lingen_program=lingen_${layer}
                tuning_thresholds=collective:64,recursive:32,notiming:0
                ${extra}
                TARGET_DEPENDENCIES lingen_${layer}
                AVOID_CONCURRENT 9
                ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=9 OMP_THREAD_LIMIT=9
                TIMEOUT 90
            )
        endforeach()
    endif()

    if(layer STREQUAL "p1")
        cado_define_test(test_bwc_lingen_21744
            AVOID_CONCURRENT 4
            ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
            PROGRAM ${CADO_NFS_BINARY_DIR}/linalg/bwc/lingen_p1
            ARGUMENTS
            --afile ${CMAKE_CURRENT_SOURCE_DIR}/A.21744 ffile=/dev/null m=4 n=2
            prime=485263
                tuning_thresholds=collective:64,recursive:32,notiming:0
            TARGET_DEPENDENCIES lingen_p1
            TIMEOUT 60)

        cado_define_test(test_bwc_lingen_21744_bis
            AVOID_CONCURRENT
            PROGRAM ${CADO_NFS_BINARY_DIR}/linalg/bwc/lingen_p1
            ARGUMENTS
            --random-input-with-length 20   m=4 n=2  prime=1009 --seed 1377
                tuning_thresholds=collective:64,recursive:32,notiming:0
            TARGET_DEPENDENCIES lingen_p1
            TIMEOUT 60)
    endif()

    cado_define_test(TEST_NAME test_matpoly_${layer}
        SOURCES test_matpoly.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_qcode_prime.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_expected_pi_length.cpp
        ARGUMENTS --prime ${p}
        LIBRARIES lingen_${layer}_support
        AVOID_CONCURRENT 4
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
        )

    set(our_compile_flags)
    parse_arith_layer(${layer} our_compile_flags)
    set(layer ${canonical_arith_layer})

    add_executable(time_matpoly_ft_parallel_${layer} EXCLUDE_FROM_ALL
        time_matpoly_ft_parallel.cpp)
    target_link_libraries(time_matpoly_ft_parallel_${layer}
        lingen_${layer}_support)

    add_executable(test_lingen_mini_${layer} EXCLUDE_FROM_ALL
        test_lingen_mini.cpp
    )
    target_link_libraries(test_lingen_mini_${layer}
        lingen_${layer}_core)

    mark_targets_as_mpi(
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer}
        test_lingen_mini_${layer})

    mark_targets_as_openmp(
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer}
        test_lingen_mini_${layer})

    set_property(
        TARGET
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer}
        test_lingen_mini_${layer}
        APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})

endforeach()

foreach(layer b64)
    # The binary matpoly tests are a nightmare because of openmp. Setting
    # OMP_WAIT_POLICY=passive GOMP_SPINCOUNT=0 (as we do in ci.sh) is a
    # good way to alleviate the slowdown. OMP_PROC_BIND is another, but
    # we can't afford it since it impacts other tests negatively.
    #
    # So even though on a normal machine these tests take roughly one
    # second, the only option that we have is to increase the timeout
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
    include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
    include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
    cado_define_test(TEST_NAME test_matpoly_${layer}
        SOURCES test_matpoly.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_qcode_binary.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_expected_pi_length.cpp
        LIBRARIES lingen_${layer}_support
        AVOID_CONCURRENT 4
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
        TIMEOUT 150
    )

    set(our_compile_flags)
    parse_arith_layer(${layer} our_compile_flags)
    list(APPEND our_compile_flags -DLINGEN_BINARY)

    # (unused for p=2 anyway)
    # ARGUMENTS --prime 2

    add_executable(time_matpoly_ft_parallel_${layer} EXCLUDE_FROM_ALL
        time_matpoly_ft_parallel.cpp)
    target_link_libraries(time_matpoly_ft_parallel_${layer}
        lingen_${layer}_support)

    mark_targets_as_mpi(
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer})

    mark_targets_as_openmp(
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer})

    set_property(
        TARGET
        test_matpoly_${layer}
        time_matpoly_ft_parallel_${layer}
        APPEND PROPERTY COMPILE_OPTIONS ${our_compile_flags})


    cado_define_test(TEST_NAME test_matpoly_basecase_${layer}
        PROGRAM
        ${CMAKE_CURRENT_BINARY_DIR}/test_matpoly_${layer}
        ARGUMENTS --test-basecase m=128 n=64 len1=256
        TARGET_DEPENDENCIES test_matpoly_${layer}
        AVOID_CONCURRENT 4
        ENVIRONMENT OMP_DYNAMIC=true OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
        TIMEOUT 180
        )
endforeach()

add_executable(check_mpi_problems EXCLUDE_FROM_ALL
        ${CADO_NFS_SOURCE_DIR}/linalg/bwc/check_mpi_problems.c)

mark_targets_as_mpi(check_mpi_problems)
set_property(TARGET check_mpi_problems
    APPEND PROPERTY COMPILE_DEFINITIONS WANT_MAIN)


set(mod2_plain_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
add_bwc_test(mod2_plain bwc-ptrace.sh prime=2 nullspace=left
    random_matrix_size=800 m=64 n=64 seed=1 thr=2x2
    TARGET_DEPENDENCIES mf_scan2
    AVOID_CONCURRENT 4
    )
create_magma_sage_diversions(mod2_plain)
if(SAGE)
    cado_set_test_timeout(test_bwc_mod2_plain_with_sagemath 480)
endif()

add_bwc_test(mod2_padded bwc-ptrace.sh prime=2 nullspace=left
    random_matrix_size=700,600 m=64 n=64 seed=1 thr=2x2
    TARGET_DEPENDENCIES mf_scan2
    AVOID_CONCURRENT 4
    )

add_bwc_test(mod2_secure_continuation bwc-test-secure-continuation.sh
    prime=2 nullspace=left random_matrix_size=800 m=64 n=64 seed=1 thr=2x2 interval=64
    TARGET_DEPENDENCIES mf_scan2
    TIMEOUT 60
    AVOID_CONCURRENT 4
)

add_bwc_test(mod2_secure_consistent-defaults
    bwc-test-secure-consistent-defaults.sh
    prime=2 nullspace=left random_matrix_size=800 m=64 n=64 seed=1 thr=2x2 interval=64
    TARGET_DEPENDENCIES mf_scan2
    TIMEOUT 60
    AVOID_CONCURRENT 4
)

# Note: because of stop_at_step here, we're not verifying what we're
# doing, which is a bit stupid.
add_bwc_test(mod2_rolling bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=2000
    interval=10 keep_rolling_checkpoints=2 checkpoint_precious=40 m=64 n=64
    seed=1 thr=2x2 stop_at_step=lingen
    AVOID_CONCURRENT 4
    TARGET_DEPENDENCIES mf_scan2)

add_bwc_test(mod2_8k bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=8000 m=128 n=128 seed=1 thr=2
    TARGET_DEPENDENCIES mf_scan2
    AVOID_CONCURRENT 4
    TIMEOUT 180
)

if (COMPILE_MPFQ_BINARY_FIELD_m128)
    add_bwc_test(mod2_8k_simd128 bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=8000 m=128 n=128 seed=1 thr=2 simd=128
    TARGET_DEPENDENCIES mf_scan2
    AVOID_CONCURRENT 4
    TIMEOUT 90
)
endif()

set(mod2_oddthr mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
add_bwc_test(mod2_oddthr bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=800 m=128 n=128 seed=1 thr=3
    AVOID_CONCURRENT 4
    TARGET_DEPENDENCIES mf_scan2)

add_bwc_test(mod2_interleaving bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=800 m=128 n=128
    interleaving=1 seed=1 thr=2
    AVOID_CONCURRENT 4
    TARGET_DEPENDENCIES mf_scan2)
 
add_bwc_test(mod2_m_neq_n bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=400 m=128 n=64
    seed=1 thr=5
    AVOID_CONCURRENT 6
    TARGET_DEPENDENCIES mf_scan2)

if(HAVE_HWLOC)
    cado_define_test(test_cpubinding.cpp
        ARGUMENTS
        --cpubinding ${CMAKE_CURRENT_SOURCE_DIR}/cpubind_config_file_for_tests.conf
        LIBRARIES utils bwc_base)
endif()

#########################################################################

add_executable(test-parallelizing_info test-parallelizing_info.cpp)
target_link_libraries(test-parallelizing_info utils bwc_base)
mark_targets_as_mpi(test-parallelizing_info)
set(mpiconfs mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
    set(mpiconfs)
endif()
foreach(mpiconf nompi ${mpiconfs})
    set(suffix _${mpiconf})
    if(suffix STREQUAL _nompi)
        set(suffix)
    endif()
    cado_define_test(test-parallelizing_info${suffix} PROGRAM
        ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
        ${PROJECT_BINARY_DIR}/linalg/bwc/bwc.pl
        :mpirun
        -- ${CMAKE_CURRENT_BINARY_DIR}/test-parallelizing_info -v
        thr=2x2
        AVOID_CONCURRENT 4
        TARGET_DEPENDENCIES test-parallelizing_info
        )
endforeach()

#########################################################################

set(mpiconfs mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
    set(mpiconfs)
endif()
foreach(mpiconf nompi ${mpiconfs})
    set(suffix _${mpiconf})
    if(suffix STREQUAL _nompi)
        set(suffix)
    endif()
    cado_define_test(test_bwc_dispatch_matmul${suffix} SCRIPT
        ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
        ${CMAKE_CURRENT_SOURCE_DIR}/bwc-dispatch-matmul.sh
        density=2
        thr=2x2
        AVOID_CONCURRENT 4
        TARGET_DEPENDENCIES bwc_matmul_gf2 dispatch mf_scan2 mf_bal
        random_matrix
        )
endforeach()

#########################################################################

if (WITH_MPI)
add_definitions(-DWITH_MPI)
endif()

#########################################################################

add_executable(spmv_test EXCLUDE_FROM_ALL spmv_test.cpp)
mark_targets_as_mpi(spmv_test)
target_link_libraries(spmv_test bwc_base)

add_executable(short_matmul short_matmul.c)
target_link_libraries(short_matmul bwc_base)

set(test_spmv_binary_b64_alltests YES)
set(test_spmv_binary_b128_alltests YES)
set(test_spmv_binary_b256_alltests YES)
set(test_spmv_prime_p3_alltests YES)

macro(create_test_spmv test_spmv_base)
    if(HAVE_PERL_DIGEST_MD5)
        if (test_spmv_base MATCHES "^prime")
            set(deps bwc_matmul_gf2)
        else()
            set(deps bwc_matmul_gfp)
        endif()
        set(mpiconfs)
        if(${test_spmv_${test_spmv_base}_alltests})
            set(mpiconfs mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
        endif()
        if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
            set(mpiconfs)
        endif()
        foreach(mpiconf nompi ${mpiconfs})
            set(suffix _${mpiconf})
            if(suffix STREQUAL _nompi)
                set(suffix)
            endif()
            cado_define_test(
                TEST_NAME test_spmv_${test_spmv_base}${suffix}
                TIMEOUT 40
                NO_DEFAULT_RUN
                TARGET_DEPENDENCIES
                random_matrix spmv_test short_matmul mf_scan2 mf_bal bwc_matmul_gf2
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-spmv.sh
                --matrix-size 1000 --density 10
                ${ARGN}
                )


            cado_divert_test(test_spmv_${test_spmv_base}${suffix} 11 APPEND_ARGUMENTS thr=1x1)
            if(${test_spmv_${test_spmv_base}_alltests})
                if(DEFINED ENV{CHECKS_EXPENSIVE})
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix}
                        12
                        APPEND_ARGUMENTS thr=1x2
                        AVOID_CONCURRENT 2
                        )
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix}
                        22
                        APPEND_ARGUMENTS thr=2x2
                        AVOID_CONCURRENT 4
                        )
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix}
                        13
                        APPEND_ARGUMENTS thr=1x3
                        AVOID_CONCURRENT 3
                        )
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix}
                        32
                        APPEND_ARGUMENTS thr=3x2
                        AVOID_CONCURRENT 6
                        )
                endif()
                cado_divert_test(test_spmv_${test_spmv_base}${suffix}
                    23
                    APPEND_ARGUMENTS thr=2x3
                    AVOID_CONCURRENT 6
                    )
            endif()
        endforeach()
    else()
        # Then this macro will expand to nothing, that's life.
    endif()
endmacro()

#########################################################################

set(mpiconfs mpi_rect1 mpi_rect2)
if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
    set(mpiconfs)
endif()
foreach(mpiconf nompi ${mpiconfs})
    set(suffix _${mpiconf})
    if(suffix STREQUAL _nompi)
        set(suffix)
    endif()
    cado_define_test(test_random_matrix${suffix}
        SCRIPT
        ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
        ${PROJECT_BINARY_DIR}/linalg/bwc/bwc.pl :mpirun
        --
        ${CADO_NFS_BINARY_DIR}/linalg/bwc/krylov --random_matrix 20000,density=100 m=64 n=64 ys=0..64 prime=2 interval=100 end=200
        TARGET_DEPENDENCIES bwc_full
        TIMEOUT 60
        )
    cado_divert_test(test_random_matrix${suffix} thr22 APPEND_ARGUMENTS thr=2x2
        TIMEOUT 180)
endforeach()


foreach(arith_layer ${BWC_GF2_ARITHMETIC_BACKENDS})
    parse_arith_layer(${arith_layer} dummy)
    set(arith_layer ${canonical_arith_layer})

    string_join(mm_backends " " ${BWC_GF2_MATMUL_BACKENDS})

    create_test_spmv(binary-${arith_layer}
        --backends "${mm_backends}"
        --arith-layer ${arith_layer}
        TIMEOUT 40)
endforeach()

foreach(arith_layer ${BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS})
    set(extra)
    string(REGEX MATCH "p([0-9]+)" t "${arith_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    endif()
    string_join(mm_backends " " ${BWC_GFP_MATMUL_BACKENDS})
    create_test_spmv(prime-${arith_layer}
        --backends "${mm_backends}"
        --arith-layer ${arith_layer}
        --prime "${prime_${width}}")
endforeach()

cado_define_test(TEST_NAME spmv_bucket_30075
    PROVIDE_TEMPORARY_WDIR
    PROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/bug30075.sh
    TARGET_DEPENDENCIES
    spmv_test bwc_matmul_gf2 mf_bal mf_scan2
    AVOID_CONCURRENT 16
    TIMEOUT 40
)


#########################################################################

cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-matcache.sh
    --matrix-size 1000 --density 10
    TARGET_DEPENDENCIES
    random_matrix build_matcache
    bench_matcache
    bwc_matmul_gf2
)

# This test is not very important. And honestly. since we don't use its
# output, it way very well be total garbage. The whole point is that
# bench_polmatmul contains important test code that should perhaps
# persist, and therefore we want this code to at least keep basic
# functionality.
cado_define_test(PROGRAM
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bench_polmatmul
    5000 4 4
    TARGET_DEPENDENCIES bench_polmatmul TIMEOUT 40)

#########################################################################

cado_define_test(TEST_NAME bwc_staged_krylov
    PROVIDE_TEMPORARY_WDIR
    PROGRAM
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/krylov thr=2x2 nullspace=left
    interval=200  mn=64 prime=2  ys=0..64 start=0 end=200
    skip_online_checks=1 rebuild_cache=1 seed=1
    sequential_cache_build=1 random_matrix=10000,density=150,seed=1
    no_save_cache=1 verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch
    TARGET_DEPENDENCIES bwc_full_gf2
    TIMEOUT 180
)

# Tests which use bwc_base will need bwc_config.h
include_directories(${PROJECT_BINARY_DIR}/linalg/bwc)

#########################################################################

add_executable(test_hello test_hello.cpp)
target_link_libraries(test_hello bwc_base)
mark_targets_as_mpi(test_hello)
set(mpiconfs mpi_rect1 mpi_square1)
if(DEFINED ENV{CHECKS_EXPENSIVE})
    set(mpiconfs ${mpiconfs} mpi_rect2 mpi_square2)
endif()
if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
    set(mpiconfs)
endif()
foreach(mpiconf nompi ${mpiconfs})
    set(suffix _${mpiconf})
    if(suffix STREQUAL _nompi)
        set(suffix)
    endif()
    cado_define_test(TEST_NAME test_hello${suffix}
        TARGET_DEPENDENCIES test_hello
        PROGRAM
        ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
        ${PROJECT_BINARY_DIR}/linalg/bwc/bwc.pl
        :mpirun
        -- ${CMAKE_CURRENT_BINARY_DIR}/test_hello -v
        )
    if(DEFINED ENV{CHECKS_EXPENSIVE})
        cado_divert_test(test_hello${suffix} thr12 APPEND_ARGUMENTS thr=1x2)
        cado_divert_test(test_hello${suffix} thr22 APPEND_ARGUMENTS thr=2x2)
    endif()
    cado_divert_test(test_hello${suffix} thr23 APPEND_ARGUMENTS thr=2x3)
endforeach()

cado_define_test(test_bl_extraction.cpp
    ${PROJECT_SOURCE_DIR}/linalg/bwc/blocklanczos_extraction.cpp
    LIBRARIES bwc_base bitlinalg)
mark_targets_as_mpi(test_bl_extraction)

#########################################################################

cado_define_test(TEST_NAME test-blocklanczos
    PROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/test-blocklanczos.sh
        --density 10
    NO_DEFAULT_RUN
    AVOID_CONCURRENT 16
    TARGET_DEPENDENCIES mf_scan2 mf_bal random_matrix blocklanczos
    bwc_full_gf2 short_matmul)
cado_divert_test(test-blocklanczos right
    APPEND_ARGUMENTS --matrix-size "1000,2000" --nullspace RIGHT)
cado_divert_test(test-blocklanczos left
    APPEND_ARGUMENTS --matrix-size "2000,1000" --nullspace left)
cado_divert_test(test-blocklanczos left-thr21
    APPEND_ARGUMENTS --matrix-size "2000,1000" --nullspace left thr=2x1)
cado_divert_test(test-blocklanczos left-thr12
    APPEND_ARGUMENTS --matrix-size "90,69" --nullspace left thr=1x2)
cado_divert_test(test-blocklanczos left-thr23
    APPEND_ARGUMENTS --matrix-size "2500,1000" --nullspace left thr=2x3)

cado_define_test(TEST_NAME test-blocklanczos-staged
    AVOID_CONCURRENT 16
    PROVIDE_TEMPORARY_WDIR
    PROGRAM
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bwc.pl :mpirun --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/blocklanczos
    random_matrix=1200,1000,80,seed=1 m=64 n=64 ys=0..64
    no_save_cache=1 interval=10 nullspace=left
    TARGET_DEPENDENCIES test-blocklanczos_dependencies)

#########################################################################
### if(HAVE_PERL_DIGEST_MD5)
### # This is testing the old binary.
### cado_define_test(TEST_NAME test-lingen-crc
###     AVOID_CONCURRENT
###     PROVIDE_TEMPORARY_WDIR
###     PROGRAM
###     ${CADO_NFS_SOURCE_DIR}/tests/linalg/bwc/test-lingen-crc.sh
###     bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc
###     seed=0 m=64 n=64 random_stem=30000 sequence_length=200
###     # the two crcs below are produced by marginally different versions
###     # (76fbd81d8 and the one after 38c969cfe), but both seem to factor
###     # numbers correctly, so I'm not too worried.
###     expect_crc_pi=c29af350,af9c33a0 expect_crc_F=04f4e167,3aa047a4
###     TARGET_DEPENDENCIES lingen)
### endif()
### 
#########################################################################

cado_define_test(test-arith-modp.cpp LIBRARIES tests utils)

cado_define_test(test_matrix_file.cpp
    ${CADO_NFS_SOURCE_DIR}/linalg/bwc/matrix_file.cpp
    $<TARGET_OBJECTS:bwc_parallelization>
    NO_DEFAULT_RUN
    LIBRARIES utils)

mark_targets_as_mpi(test_matrix_file)
mark_targets_as_openmp(test_matrix_file)

cado_define_test(test_matrix_file_script
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
${CMAKE_CURRENT_SOURCE_DIR}/test_matrix_file.sh
TARGET_DEPENDENCIES test_matrix_file
NO_DEFAULT_RUN)

cado_divert_test(test_matrix_file_script sample1 APPEND_ARGUMENTS sample1)
cado_divert_test(test_matrix_file_script sample2 APPEND_ARGUMENTS sample2)


cado_define_test(test_linalg_example
    SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/test_lingen_example.sh
    PROVIDE_TEMPORARY_WDIR
    NO_DEFAULT_RUN
)

if("p1" IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
cado_divert_test(test_linalg_example p1_30105
    APPEND_ARGUMENTS
    layer=p1
    tuning_thresholds=recursive:10,flint:10,notiming:0
    ${CMAKE_CURRENT_SOURCE_DIR}/testcase_30105.E
    TARGET_DEPENDENCIES
    test_lingen_mini_p1
    lingen_verify_checkpoints_p1
)
endif()

cado_define_test(test_bwc_filenames.cpp
    ${PROJECT_SOURCE_DIR}/linalg/bwc/bwc_filenames.cpp
    LIBRARIES utils)

add_subdirectory (flint-fft)

#########################################################################

cado_define_test(linalg_determinant
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bwc.pl
        :determinant thr=2 nullspace=right
        matrix=${CMAKE_CURRENT_SOURCE_DIR}/d20.sparse.bin
    SHA1_ON_REGEXP_LINES "^determinant modulo p: "
    TARGET_DEPENDENCIES bwc_full_gfp
    NO_DEFAULT_RUN
)

cado_divert_test(linalg_determinant p4294967291
    APPEND_ARGUMENTS prime=4294967291
    EXPECT_SHA1 68a613f825efc6d86039ed8fbc8d221ce4bab632
    NO_DEFAULT_RUN
)

cado_divert_test(linalg_determinant_p4294967291 n1 APPEND_ARGUMENTS n=1 m=2)
cado_divert_test(linalg_determinant_p4294967291 n3 APPEND_ARGUMENTS n=3 m=3)
cado_divert_test(linalg_determinant_p4294967291 n5 APPEND_ARGUMENTS n=5 m=10)

cado_divert_test(linalg_determinant p2341309
    APPEND_ARGUMENTS prime=2341309
    EXPECT_FAIL
    EXPECT_SHA1 8db6dd9ba108474d478ddd27daf698847f12cd6f
    NO_DEFAULT_RUN
)

cado_divert_test(linalg_determinant_p2341309 n1 APPEND_ARGUMENTS n=1 m=2)
cado_divert_test(linalg_determinant_p2341309 n3 APPEND_ARGUMENTS n=3 m=6)
cado_divert_test(linalg_determinant_p2341309 n4 APPEND_ARGUMENTS n=4 m=4)

cado_divert_test(linalg_determinant plarge
    APPEND_ARGUMENTS prime=340282365653287863235145205935065006031 n=2 m=4
    EXPECT_SHA1 fdab6722a96703c0463e5e1a407e60988d66a039
)

cado_define_test(linalg_determinant_nonsquare
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bwc.pl
        :determinant n=1 m=2 prime=101 thr=2 nullspace=right random_matrix=50,51
    EXPECT_FAIL
    TARGET_DEPENDENCIES bwc_full_gfp
)

cado_define_test(linalg_determinant_p2
    PROVIDE_TEMPORARY_WDIR
    SCRIPT
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bwc.pl
        :determinant n=64 m=128 prime=2 thr=2 random_matrix=100
    EXPECT_FAIL
    TARGET_DEPENDENCIES bwc_full_gfp
)
