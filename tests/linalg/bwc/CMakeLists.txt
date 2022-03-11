# NOTE: most GF(p) checks here depend on lingen_pz ; this is artifical,
# and caused by the fact that bwc.pl is for the moment unable to call
# anything else than lingen_pz

# some tests below use the bwc include files.
include_directories(${PROJECT_SOURCE_DIR}/linalg)
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc)
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq)
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
        endif()
        cado_define_test(test_bwc_${nm}${suffix} SCRIPT
            ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg mats --other --arg wdir --
            ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
            ${CADO_NFS_SOURCE_DIR}/tests/linalg/bwc/${script}
            bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc pre_wipe=1 ${ARGN} -- verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch,^bwc-loading-mksol-files,^bwc-cache-major-info,^bwc-cpubinding,^perl-checks,^perl-sections
            AVOID_CONCURRENT
            TARGET_DEPENDENCIES bwc_full)
    endforeach()
endmacro()


# All this is 64-bit only !

set(prime_1 281474976710677)
set(prime_2 5192296858534827628530496329220121)
set(prime_3 95780971304118053647396689196894323976171195136475563)
set(prime_4 1766847064778384329583297500742918515827483896875618958121606201292619891)
set(prime_5 32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358053)
set(prime_6 601226901190101306339707032778070279008174732520529886901066488712245510429339761526706943586500787976175353983)
set(prime_7 11090678776483259438313656736572334813745748301503266300681918322458485231222502492159897624416558312389564843845614287315896632389)
set(prime_8 204586912993508866875824356051724947013540127877691549342705710506008362275292159680204380770369009821930417757972504438076078534117837065833032974919)
set(prime_9 3773962424821541352241554580988268890916921220416440428376206300245624162392148852086126725177658767541468375030763844899770584629924792632561434251432696043649395327187)
set(prime_10 909932137091507349600288597675321757101414873334402438858626164916349345767541752041377306747887895885200013437238135716251354881646154492844839284172487068480415175470480234519706599)
set(prime_11 13753158638690205675902021824448304531266335651977023929726467446560173439346691276391932466456818176470980294225128109021325275194117892575170183673185378686943300130203536714365357041913407373)
set(prime_12 4692851634592474092451166992324491990345415997425726384327421710769161180436738237983736300974130547508887839887328010265354129712831899610758743135488356905756675150570418835123110229290807866956335696997123808253)
set(prime_13 13388377234648697797204242513192222495613229645153417754162592753075059328303771597899086676966128161635082877836594479456533503960668428615439212786511751016880713962371421583171013773630512114726152981840500396109472830638409781678764122812210613)
set(prime_14 116794188855301759465775830943452977444472544566396791289770302684456866852120153216374564917242878827570194742827987846047833882294285805095893735571230206011354707561346289783366234291185550813375522907780819050448884878716880881681486735868973620798833004312957)
set(prime_15 358752939420752931158375023855513400064244447993094387657025186121127116148346374410346539594710018174137751708688854198856062241744412981520402102281094231132870913674050365678978481894462185188236318219972031917211474544119229541236960242093722745601167293308448934364898080327)
set(prime_z 69617318994479297159441705409245167921344429126717528237597542082203295398081625160307507496908132931192662194421301381083506846944815643283884602656894137393981852330936660004926669193667)

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
            cado_define_test(test_bwc_lingen_binary_${nm}${suffix} SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-lingen.sh
                ${ARGN}
                TARGET_DEPENDENCIES lingen_u64k1
                AVOID_CONCURRENT
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
        tuning_thresholds=recursive:16,notiming:16,collective:16
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
            tuning_thresholds=recursive:128,notiming:16,collective:128
            )
    endif()

    # The test below fails on 32bit. This should be investigated.
    # add_lingen_binary_test(test1
    # expect_sha1_F=30b90451fb828f08caa160f96c585c4b1c0951a8 m=64 n=64 sequence_length=200 seed=1  -- lingen_mpi_threshold=10)
endif()

# Now we only ever test p_1, p_3, and pz.
set(BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS)
foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})
    string(REGEX MATCH "^(p_1|p_3|pz)$" t "${gfp_layer}")
    if(t)
        list(APPEND BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS ${t})
    else()
        message(STATUS "Note: tests for ${gfp_layer} are skipped, because we restrict the routine tests")
    endif()
endforeach()
        

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS})
    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
        set(p "${prime_${width}}")
    elseif(gfp_layer STREQUAL "pz")
        set(t TRUE)
        set(width "z")
        set(extra "lingen_program=lingen_pz")
        set(p "${prime_${width}}")
    else()
        message(FATAL_ERROR "bad arithmetic backend ${gfp_layer}")
    endif()
    if(gfp_layer STREQUAL "pz")
    elseif("pz" IN_LIST BWC_GFP_ARITHMETIC_BACKENDS)
        if(gfp_layer STREQUAL "p_1" OR gfp_layer STREQUAL "p_13")
            set(modp_inhomogeneous_${gfp_layer}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
        endif()
        add_bwc_test(modp_inhomogeneous_${gfp_layer} bwc-ptrace.sh prime=${p} random_matrix_size=500 nrhs=2 m=6 n=4 seed=1 ${extra})
        add_dependencies(test_bwc_modp_inhomogeneous_${gfp_layer}_dependencies
            lingen_pz lingen_${gfp_layer})

        if(gfp_layer STREQUAL "p_1")
            set(mextra)
            # Only try to enable magma if CHECKS_EXPENSIVE is set.
            if(DEFINED ENV{CHECKS_EXPENSIVE} AND MAGMA)
                set(mextra magma=${MAGMA})
            endif()
            if(DEFINED ENV{CHECKS_EXPENSIVE})
                # it's quick, but let's avoid redundant tests
                add_bwc_test(modp_inhomogeneous_minimal_3rhs_mn6 bwc-ptrace.sh prime=1009
                    random_matrix_size=50 nrhs=3 m=6 n=6 seed=1 ${mextra})
                add_dependencies(test_bwc_modp_inhomogeneous_minimal_3rhs_mn6_dependencies
                    lingen_pz lingen_${gfp_layer})
            endif()

            set(modp_inhomogeneous_minimal_2rhs_mn4_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
            add_bwc_test(modp_inhomogeneous_minimal_2rhs_mn4 bwc-ptrace.sh
                prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1 ${mextra})
            add_dependencies(test_bwc_modp_inhomogeneous_minimal_2rhs_mn4_dependencies
                lingen_pz lingen_${gfp_layer})

            if(DEFINED ENV{CHECKS_EXPENSIVE})
                add_bwc_test(modp_inhomogeneous_minimal_2rhs_mn4_expensive bwc-ptrace.sh
                    prime=1237940039285380274899124357 random_matrix_size=200 nrhs=2 m=4 n=4 seed=1 ${mextra})
                add_dependencies(test_bwc_modp_inhomogeneous_minimal_2rhs_mn4_expensive_dependencies
                    lingen_pz lingen_${gfp_layer})
            endif()

            if(DEFINED ENV{CHECKS_EXPENSIVE})
                add_bwc_test(modp_inhomogeneous_minimal_1rhs_mn1 bwc-ptrace.sh
                    prime=1009 random_matrix_size=50 nrhs=1 m=1 n=1 seed=1 ${mextra})
                add_dependencies(test_bwc_modp_inhomogeneous_minimal_1rhs_mn1_dependencies
                    lingen_pz lingen_${gfp_layer})
            endif()

            add_bwc_test(modp_homogeneous_minimal_mn4 bwc-ptrace.sh
                prime=1009 random_matrix_size=50 m=4 n=4 seed=1
                random_matrix_minkernel=1 ${mextra})
            add_dependencies(test_bwc_modp_homogeneous_minimal_mn4_dependencies
                lingen_pz lingen_${gfp_layer})
            #        add_bwc_test(modp_inhomogeneous_minimal_2rhs_mn2 bwc-ptrace.sh
            #            prime=1009 random_matrix_size=50 nrhs=2 m=2 n=2 seed=1 ${mextra})
            #        add_dependencies(modp-inhomogeneous-minimal-2rhs-mn2_dependencies
            #        lingen_pz lingen_${gfp_layer})
        endif()
        if(gfp_layer STREQUAL "p_1" OR gfp_layer STREQUAL "p_13")
            set(modp_secure_continuation_${gfp_layer}_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
        endif()
        add_bwc_test(modp_secure_continuation_${gfp_layer} bwc-test-secure-continuation.sh prime=${p} random_matrix_size=500 nrhs=4 m=8 n=8 seed=1 interval=32 ${extra})
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

    # Note that we're only enabling the p_1 tests because the checksums
    # below were obtained for p=1009 specifically.
    if(gfp_layer STREQUAL "p_1")
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
            cado_define_test(test_bwc_lingen_${gfp_layer}_test0${suffix}
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env WDIR
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                2b6428178054fc0d3931dcdc9257ccf3bdea62f9 4 2 36 1009 17042
                mpi_magic=${mpiconf} thr=2x2 --ascii
                ${extra}
                TARGET_DEPENDENCIES lingen_${gfp_layer}
                AVOID_CONCURRENT
                )
            cado_divert_test(test_bwc_lingen_${gfp_layer}_test0${suffix} frag
                APPEND_ARGUMENTS io_block_size=1024)
            cado_define_test(test_bwc_lingen_${gfp_layer}_test1${suffix}
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env WDIR
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
                thr=2x2 --ascii
                ${extra}
                TARGET_DEPENDENCIES lingen_${gfp_layer}
                AVOID_CONCURRENT
                )
            cado_define_test(test_bwc_lingen_${gfp_layer}_test2${suffix}
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env WDIR
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                b448493e04299d8d54be6cc1273c43293520ac9e 3 3 1000 1009 1111
                thr=3x3 --ascii
                ${extra}
                TARGET_DEPENDENCIES lingen_${gfp_layer}
                AVOID_CONCURRENT
                )
            # We really _can't_ have ascii with the checkpoint test. ascii
            # modifies the anticipation of the schedule, which in turns means
            # that the list of cached schedules is modified in the course of
            # the execution, and that eventually goes to the checkpoint file
            # as well, wreaking obvious havoc.
            cado_define_test(test_bwc_lingen_${gfp_layer}_test2_checkpoints${suffix}
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env WDIR
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-lingen-checkpoints.sh
                c82b888f2ba844f6db1e30e091444d1023052f87 3 3 1000 1009 1111
                thr=3x3
                checkpoint_threshold=40
                ${extra}
                TARGET_DEPENDENCIES
                lingen_${gfp_layer}
                lingen_verify_checkpoints_${gfp_layer}
                AVOID_CONCURRENT
                )
        endforeach()
    endif()

    if(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_3")
        set(mpiconfs mpi_square1 mpi_square2)
        if (NOT HAVE_MPI OR CMAKE_VERSION VERSION_LESS 3.9)
            set(mpiconfs)
        endif()
        foreach(mpiconf nompi ${mpiconfs})
            set(suffix _${mpiconf})
            if(suffix STREQUAL _nompi)
                set(suffix)
            endif()
            cado_define_test(test_bwc_lingen_${gfp_layer}_test4${suffix}
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --env WDIR
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
                00a481d0e0e2ee7403fa234b7f8b818dc2ec8db4 3 3 1000
                1532495540865888858358347027150309183618739122183602191 1111
                thr=3x3
                lingen_program=lingen_${gfp_layer}
                ${extra}
                TARGET_DEPENDENCIES lingen_${gfp_layer}
                AVOID_CONCURRENT
            )
        endforeach()
    endif()

    if(gfp_layer STREQUAL "p_1")
        cado_define_test(test_bwc_lingen_21744
            AVOID_CONCURRENT
            PROGRAM ${CADO_NFS_BINARY_DIR}/linalg/bwc/lingen_p_1
            ARGUMENTS
            --afile ${CMAKE_CURRENT_SOURCE_DIR}/A.21744 ffile=/dev/null m=4 n=2
            prime=485263
            TARGET_DEPENDENCIES lingen_p_1)

        cado_define_test(test_bwc_lingen_21744_bis
            AVOID_CONCURRENT
            PROGRAM ${CADO_NFS_BINARY_DIR}/linalg/bwc/lingen_p_1
            ARGUMENTS
            --random-input-with-length 20   m=4 n=2  prime=1009 --seed 1377
            TARGET_DEPENDENCIES lingen_p_1)
    endif()

    cado_define_test(TEST_NAME test_matpoly_${gfp_layer}
        SOURCES test_matpoly.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_qcode_prime.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_expected_pi_length.cpp
        ARGUMENTS --prime ${p}
        LIBRARIES lingen_${gfp_layer}_support)

    set_target_properties(test_matpoly_${gfp_layer} PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${gfp_layer}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")

    add_executable(time_matpoly_ft_parallel_${gfp_layer} EXCLUDE_FROM_ALL
        time_matpoly_ft_parallel.cpp)
    target_link_libraries(time_matpoly_ft_parallel_${gfp_layer}
        lingen_${gfp_layer}_support)

    set_target_properties(time_matpoly_ft_parallel_${gfp_layer} PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${gfp_layer}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")
endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

foreach(gfp_layer u64k1)
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/fft/.libs)
    link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
    include_directories(${CADO_NFS_BINARY_DIR}/gf2x/fft)
    include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
    cado_define_test(TEST_NAME test_matpoly_${gfp_layer}
        SOURCES test_matpoly.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_qcode_binary.cpp
        ${PROJECT_SOURCE_DIR}/linalg/bwc/lingen_expected_pi_length.cpp
        LIBRARIES lingen_${gfp_layer}_support)
    # (unused for p=2 anyway)
    # ARGUMENTS --prime 2
    set_target_properties(test_matpoly_${gfp_layer} PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${gfp_layer}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")
    add_executable(time_matpoly_ft_parallel_${gfp_layer} EXCLUDE_FROM_ALL
        time_matpoly_ft_parallel.cpp)
    target_link_libraries(time_matpoly_ft_parallel_${gfp_layer}
        lingen_${gfp_layer}_support)
    set_target_properties(time_matpoly_ft_parallel_${gfp_layer} PROPERTIES
        COMPILE_FLAGS "--mpi ${OpenMP_CXX_FLAGS} -DSELECT_MPFQ_LAYER_${gfp_layer}"
        LINK_FLAGS "--mpi ${OpenMP_CXX_FLAGS}")

    cado_define_test(TEST_NAME test_matpoly_basecase_${gfp_layer}
        PROGRAM
        ${CMAKE_CURRENT_BINARY_DIR}/test_matpoly_${gfp_layer}
        ARGUMENTS --test-basecase m=128 n=64 len1=256
        TARGET_DEPENDENCIES test_matpoly_${gfp_layer}
        )
endforeach()

add_executable(check_mpi_problems EXCLUDE_FROM_ALL
        ${CADO_NFS_SOURCE_DIR}/linalg/bwc/check_mpi_problems.c)
set_target_properties(check_mpi_problems PROPERTIES
    COMPILE_FLAGS "--mpi ${OpenMP_C_FLAGS} -DWANT_MAIN"
    LINK_FLAGS "--mpi ${OpenMP_C_FLAGS}")


set(mod2_plain_mpi_tests mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
add_bwc_test(mod2_plain bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=800 m=64 n=64 seed=1 thr=2x2)
add_dependencies(test_bwc_mod2_plain_dependencies mf_scan)
add_bwc_test(mod2_padded bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=700,600 m=64 n=64 seed=1 thr=2x2)
add_dependencies(test_bwc_mod2_padded_dependencies mf_scan)

add_bwc_test(mod2_secure_continuation bwc-test-secure-continuation.sh
    prime=2 nullspace=left random_matrix_size=800 m=64 n=64 seed=1 thr=2x2 interval=64)
add_dependencies(test_bwc_mod2_secure_continuation_dependencies mf_scan)
add_bwc_test(mod2_secure_consistent-defaults
    bwc-test-secure-consistent-defaults.sh
    prime=2 nullspace=left random_matrix_size=800 m=64 n=64 seed=1 thr=2x2 interval=64)
add_dependencies(test_bwc_mod2_secure_continuation_dependencies mf_scan)

# Note: because of stop_at_step here, we're not verifying what we're
# doing, which is a bit stupid.
add_bwc_test(mod2_rolling bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=2000
    interval=10 keep_rolling_checkpoints=2 checkpoint_precious=40 m=64 n=64
    seed=1 thr=2x2 stop_at_step=lingen)
add_dependencies(test_bwc_mod2_rolling_dependencies mf_scan)

add_bwc_test(mod2_8k bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=8000 m=128 n=128 seed=1 thr=2)
add_dependencies(test_bwc_mod2_8k_dependencies mf_scan)

if (COMPILE_MPFQ_BINARY_FIELD_m128)
add_bwc_test(mod2_8k_simd128 bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=8000 m=128 n=128 seed=1 thr=2 simd=128)
add_dependencies(test_bwc_mod2_8k_simd128_dependencies mf_scan)
endif()

set(mod2_oddthr mpi_rect1 mpi_rect2 mpi_square1 mpi_square2)
add_bwc_test(mod2_oddthr bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=800 m=128 n=128 seed=1 thr=3)
add_dependencies(test_bwc_mod2_oddthr_dependencies mf_scan)

add_bwc_test(mod2_interleaving bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=800 m=128 n=128
    interleaving=1 seed=1 thr=2)
add_dependencies(test_bwc_mod2_oddthr_dependencies mf_scan)
 
add_bwc_test(mod2_m_neq_n bwc-ptrace.sh prime=2 nullspace=left random_matrix_size=400 m=128 n=64
    seed=1 thr=5)
add_dependencies(test_bwc_mod2_m_neq_n_dependencies mf_scan)

add_subdirectory(mpfq)


if(HAVE_CXX11 AND HAVE_HWLOC)
    cado_define_test(test_cpubinding.c
        ARGUMENTS
        --cpubinding ${CMAKE_CURRENT_SOURCE_DIR}/cpubind_config_file_for_tests.conf
        LIBRARIES utils bwc_base)
endif()

#########################################################################

add_executable(test-parallelizing_info test-parallelizing_info.cpp)
target_link_libraries(test-parallelizing_info utils bwc_base)
set_target_properties(test-parallelizing_info PROPERTIES
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")
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
        TARGET_DEPENDENCIES bwc_matmul_gf2 dispatch mf_scan mf_bal
        random_matrix
        )
endforeach()

#########################################################################

if (WITH_MPI)
add_definitions(-DWITH_MPI)
endif()

#########################################################################

add_executable(spmv_test EXCLUDE_FROM_ALL spmv_test.c)
set_target_properties(spmv_test PROPERTIES
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
target_link_libraries(spmv_test bwc_base)

add_executable(short_matmul short_matmul.c)
target_link_libraries(short_matmul bwc_base)

set(test_spmv_binary_u64k1_alltests YES)
set(test_spmv_binary_u64k4_alltests YES)
set(test_spmv_binary_m128_alltests YES)
set(test_spmv_prime_p_3_alltests YES)

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
                SCRIPT
                ${CADO_NFS_SOURCE_DIR}/tests/do_with_mpi.sh "${mpiconf}"
                ${CMAKE_CURRENT_SOURCE_DIR}/test-spmv.sh
                --matrix-size 1000 --density 10
                ${ARGN}
                NO_DEFAULT_RUN
                TARGET_DEPENDENCIES
                random_matrix spmv_test short_matmul mf_scan mf_bal bwc_matmul_gf2)


            cado_divert_test(test_spmv_${test_spmv_base}${suffix} 11 APPEND_ARGUMENTS thr=1x1)
            if(${test_spmv_${test_spmv_base}_alltests})
                if(DEFINED ENV{CHECKS_EXPENSIVE})
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix} 12 APPEND_ARGUMENTS thr=1x2)
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix} 22 APPEND_ARGUMENTS thr=2x2)
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix} 13 APPEND_ARGUMENTS thr=1x3)
                    cado_divert_test(test_spmv_${test_spmv_base}${suffix} 32 APPEND_ARGUMENTS thr=3x2)
                endif()
                cado_divert_test(test_spmv_${test_spmv_base}${suffix} 23 APPEND_ARGUMENTS thr=2x3)
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
        )
    cado_divert_test(test_random_matrix${suffix} thr22 APPEND_ARGUMENTS thr=2x2)
endforeach()


foreach(arith_layer ${BWC_GF2_ARITHMETIC_BACKENDS})
    if (arith_layer STREQUAL "m128" AND NOT HAVE_SSE41)
    else()
        string_join(mm_backends " " ${BWC_GF2_MATMUL_BACKENDS})

        create_test_spmv(binary-${arith_layer}
            --backends "${mm_backends}"
            --arith-layer ${arith_layer})
    endif()
endforeach()

foreach(arith_layer ${BWC_GFP_ARITHMETIC_BACKENDS_FOR_TESTS})
    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${arith_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    endif()
    string_join(mm_backends " " ${BWC_GFP_MATMUL_BACKENDS})
    create_test_spmv(prime-${arith_layer}
        --backends "${mm_backends}"
        --arith-layer ${arith_layer}
        --prime "${prime_${width}}")
endforeach()

#########################################################################

cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-matcache.sh
    --matrix-size 1000 --density 10
    TARGET_DEPENDENCIES
    random_matrix build_matcache bench_matcache)

#########################################################################

cado_define_test(TEST_NAME bwc_staged_krylov PROGRAM
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/krylov thr=2x2 nullspace=left
    interval=200  mn=64 prime=2  ys=0..64 start=0 end=200
    skip_online_checks=1 rebuild_cache=1 seed=1
    sequential_cache_build=1 random_matrix=10000,density=150,seed=1
    no_save_cache=1 verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch
    TARGET_DEPENDENCIES bwc_full_gf2)

# Tests which use bwc_base will need bwc_config.h
include_directories(${PROJECT_BINARY_DIR}/linalg/bwc)

#########################################################################

add_executable(test_hello test_hello.c)
target_link_libraries(test_hello bwc_base)
set_target_properties(test_hello PROPERTIES
    LINKER_LANGUAGE CXX
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")
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

cado_define_test(test_bl_extraction.cpp LIBRARIES bwc_base bitlinalg)
set_target_properties(test_bl_extraction PROPERTIES
                LINKER_LANGUAGE CXX
                COMPILE_FLAGS "--mpi"
                LINK_FLAGS "--mpi")

#########################################################################

cado_define_test(TEST_NAME test-blocklanczos
    PROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/test-blocklanczos.sh
        --density 10
        --bindir ${CADO_NFS_BINARY_DIR}/linalg/bwc
    NO_DEFAULT_RUN
    TARGET_DEPENDENCIES mf_scan mf_bal random_matrix blocklanczos
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
    AVOID_CONCURRENT
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
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
###     PROGRAM
###     ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
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

add_subdirectory (flint-fft)
