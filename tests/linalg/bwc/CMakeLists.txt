# NOTE: most GF(p) checks here depend on plingen_pz ; this is artifical,
# and caused by the fact that bwc.pl is for the moment unable to call
# anything else than plingen_pz

macro(add_bwc_test nm script)
    cado_define_test(test-bwc-${nm} SCRIPT
        ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg mats --arg wdir --
        ${CADO_NFS_SOURCE_DIR}/linalg/bwc/${script}
        bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc pre_wipe=1 ${ARGN} -- verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch,^bwc-loading-mksol-files,^bwc-cache-major-info,^bwc-cpubinding,^perl-checks,^perl-sections
        TARGET_DEPENDENCIES bwc_full)
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
macro(add_lingen_binary_test testname)
    cado_define_test(test-bwc-lingen-binary-${testname} SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/test-lingen.sh
        ${CADO_NFS_BINARY_DIR} ${ARGN}
        TARGET_DEPENDENCIES lingen)
endmacro()

# the sha1sums here are from the time were we were transposing the
# output. Now we no longer do that.
# add_lingen_binary_test(test0 1dea0d8879d125aa280e26b8573f932f6bccea75,d0e758e619c8e81b3e28bbd4fe69c32826e9a4bd 64 64 160 1 lingen_threshold=10)

# todo presumably, on 32 bit, we get something different for the sha1sum.
add_lingen_binary_test(test0 893801e0727d840353a35b2623430fb758e4d379 64 64 160 1 lingen_threshold=10)

# The test below fails on 32bit. This should be investigated.
# add_lingen_binary_test(test1 30b90451fb828f08caa160f96c585c4b1c0951a8 64 64 200 1 lingen_threshold=10)
endif()

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    else()
        if(gfp_layer STREQUAL "pz")
            set(t TRUE)
            set(width "z")
            set(extra "plingen_program=plingen_pz")
        endif()
    endif()
    set(p "${prime_${width}}")
    if(gfp_layer STREQUAL "pz")
    else()
        add_bwc_test(modp-inhomogeneous-${gfp_layer} bwc-ptrace.sh prime=${p} random_matrix_size=500 nrhs=4 m=8 n=8 seed=1 ${extra})
        add_dependencies(test-bwc-modp-inhomogeneous-${gfp_layer}_dependencies
            plingen_pz plingen_${gfp_layer})
        if(gfp_layer STREQUAL "p_1")
            set(mextra)
            # Only try to enable magma if CHECKS_EXPENSIVE is set.
            if(DEFINED ENV{CHECKS_EXPENSIVE})
                find_program(MAGMA magma HINTS ENV MAGMA)
                if (MAGMA)
                    set(mextra magma=${MAGMA})
                endif()
            endif()
            add_bwc_test(modp-inhomogeneous-minimal-3rhs-mn6 bwc-ptrace.sh prime=1009
                random_matrix_size=50 nrhs=3 m=6 n=6 seed=1 ${mextra})
            add_dependencies(test-bwc-modp-inhomogeneous-minimal-3rhs-mn6_dependencies
                plingen_pz plingen_${gfp_layer})

            add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn4 bwc-ptrace.sh
                prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1 ${mextra})
            add_dependencies(test-bwc-modp-inhomogeneous-minimal-2rhs-mn4_dependencies
                plingen_pz plingen_${gfp_layer})

            if (HAVE_MPI)
                add_bwc_test(mpi1-modp-inhomogeneous-minimal-2rhs-mn4 bwc-ptrace.sh
                    prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1
                    mpi=2x3 lingen_mpi=2x2 ${mextra})
                add_dependencies(test-bwc-mpi1-modp-inhomogeneous-minimal-2rhs-mn4_dependencies
                    plingen_pz plingen_${gfp_layer})
                add_bwc_test(mpi2-modp-inhomogeneous-minimal-2rhs-mn4 bwc-ptrace.sh
                    prime=1009 random_matrix_size=50 nrhs=2 m=4 n=4 seed=1
                    mpi=3x2 lingen_mpi=2x2 ${mextra})
                add_dependencies(test-bwc-mpi2-modp-inhomogeneous-minimal-2rhs-mn4_dependencies
                    plingen_pz plingen_${gfp_layer})
            endif()

            if(DEFINED ENV{CHECKS_EXPENSIVE})
                add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn4-expensive bwc-ptrace.sh
                    prime=1237940039285380274899124357 random_matrix_size=200 nrhs=2 m=4 n=4 seed=1 ${mextra})
                add_dependencies(test-bwc-modp-inhomogeneous-minimal-2rhs-mn4-expensive_dependencies
                    plingen_pz plingen_${gfp_layer})
            endif()

            add_bwc_test(modp-inhomogeneous-minimal-1rhs-mn1 bwc-ptrace.sh
                prime=1009 random_matrix_size=50 nrhs=1 m=1 n=1 seed=1 ${mextra})
            add_dependencies(test-bwc-modp-inhomogeneous-minimal-1rhs-mn1_dependencies
                plingen_pz plingen_${gfp_layer})
            #        add_bwc_test(modp-homogeneous-minimal-mn4 bwc-ptrace.sh
            #            prime=1009 random_matrix_size=50 m=4 n=4 seed=1 ${mextra})
            #        add_dependencies(modp-homogeneous-minimal-mn4_dependencies
            #        plingen_pz plingen_${gfp_layer})
            #        add_bwc_test(modp-inhomogeneous-minimal-2rhs-mn2 bwc-ptrace.sh
            #            prime=1009 random_matrix_size=50 nrhs=2 m=2 n=2 seed=1 ${mextra})
            #        add_dependencies(modp-inhomogeneous-minimal-2rhs-mn2_dependencies
            #        plingen_pz plingen_${gfp_layer})
        endif()
    endif()

    # Also add some plingen tests if it so occurs that the current
    # backend can process them
    # Each test has the following arguments
    # plingen test name
    # reference sha1 sum for generator
    # m
    # n
    # length of sequence to compute
    # p
    # seed for random generation (entries are (2^k*(seed/1000) mod q) + (3^k*(seed%1000) mod q) mod p, with q=2^20-3).
    # arguments to be passed to plingen (including mpi arguments)

    # we create a phony first test, so as to limit the risk of having
    # several tests trigger the same build dependency on
    # plingen_${gfp_layer} at the same time. Of course, the bare test
    # here does not have the require arguments, so it will fail (we avoid
    # NO_DEFAULT_RUN in this case, precisely because we want a phony
    # test).
    cado_define_test(test-bwc-plingen-${gfp_layer} SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/test-plingen.sh
        TARGET_DEPENDENCIES plingen_${gfp_layer})
    set_tests_properties(test-bwc-plingen-${gfp_layer} PROPERTIES WILL_FAIL 1)

    # Note that the thr argument here is used only in MPI
    # context, and ignored otherwise. In MPI context, this triggers the
    # mpi portions of the code (actually the test gets run twice, once
    # with the threads, once without. Both have to succeed, of course).
    if(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_1")
        cado_divert_test(test-bwc-plingen-${gfp_layer} test0
            ${CADO_NFS_BINARY_DIR}
            2b6428178054fc0d3931dcdc9257ccf3bdea62f9 4 2 36 1009 17042
            lingen_threshold=10 lingen_mpi_threshold=10 thr=2x2 --ascii
            ${extra})
        cado_divert_test(test-bwc-plingen-${gfp_layer} test1
            ${CADO_NFS_BINARY_DIR}
            f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
            lingen_threshold=10 lingen_mpi_threshold=10 thr=2x2 --ascii
            ${extra})
        cado_divert_test(test-bwc-plingen-${gfp_layer} test2
            ${CADO_NFS_BINARY_DIR}
            f9bc59f0caaca0fd523ef6a020a6f49fa97f38e9 4 2 200 1009 88888
            lingen_threshold=10 lingen_mpi_threshold=40 thr=2x2 --ascii
            ${extra})
        cado_divert_test(test-bwc-plingen-${gfp_layer} test3
            ${CADO_NFS_BINARY_DIR}
            b448493e04299d8d54be6cc1273c43293520ac9e 3 3 1000 1009 1111
            lingen_threshold=10 lingen_mpi_threshold=40 thr=3x3 --ascii
            ${extra})
    endif()
    if(gfp_layer STREQUAL "pz" OR gfp_layer STREQUAL "p_3")
        cado_divert_test(test-bwc-plingen-${gfp_layer} test4
            ${CADO_NFS_BINARY_DIR}
            00a481d0e0e2ee7403fa234b7f8b818dc2ec8db4 3 3 1000
            1532495540865888858358347027150309183618739122183602191 1111
            lingen_threshold=10 lingen_mpi_threshold=40 thr=3x3
            ${extra})
    endif()
endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})


add_bwc_test(mod2-plain bwc-trace.sh random_matrix_size=800 m=64 n=64 seed=1 thr=2x2)
add_dependencies(test-bwc-mod2-plain_dependencies mf_scan)
add_bwc_test(mod2-padded bwc-trace.sh random_matrix_size=700,600 m=64 n=64 seed=1 thr=2x2)
add_dependencies(test-bwc-mod2-padded_dependencies mf_scan)

# Note: because of stop_at_step here, we're not verifying what we're
# doing, which is a bit stupid.
add_bwc_test(mod2-rolling bwc-trace.sh random_matrix_size=2000
    interval=10 keep_rolling_checkpoints=2 checkpoint_precious=40 m=64 n=64
    seed=1 thr=2x2 stop_at_step=lingen)
add_dependencies(test-bwc-mod2-rolling_dependencies mf_scan)

add_bwc_test(mod2-8k bwc-trace.sh random_matrix_size=8000 m=128 n=128 seed=1 thr=2)
add_dependencies(test-bwc-mod2-8k_dependencies mf_scan)

add_bwc_test(mod2-oddthr bwc-trace.sh random_matrix_size=800 m=128 n=128 seed=1 thr=3)
add_dependencies(test-bwc-mod2-oddthr_dependencies mf_scan)

# temporarily disabling this test, until I get to re-enabling the
# feature.
# add_bwc_test(mod2-interleaving bwc-trace.sh random_matrix_size=800 m=128 n=128
#     interleaving=1 seed=1 thr=2)
# add_dependencies(test-bwc-mod2-oddthr_dependencies mf_scan)
 
add_bwc_test(mod2-m-neq-n bwc-trace.sh random_matrix_size=400 m=128 n=64
    seed=1 thr=5)
add_dependencies(test-bwc-mod2-m-neq-n_dependencies mf_scan)
if (HAVE_MPI)
add_bwc_test(mpi1-mod2-m-neq-n bwc-trace.sh random_matrix_size=400 m=128 n=64
    seed=1 thr=5 mpi=2x3 lingen_mpi=2x2)
add_dependencies(test-bwc-mpi1-mod2-m-neq-n_dependencies mf_scan)
add_bwc_test(mpi2-mod2-m-neq-n bwc-trace.sh random_matrix_size=400 m=128 n=64
    seed=1 thr=5 mpi=3x2 lingen_mpi=2x2)
add_dependencies(test-bwc-mpi2-mod2-m-neq-n_dependencies mf_scan)
endif()

add_subdirectory(mpfq)

# some tests below use the bwc include files.
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc)
include_directories(${PROJECT_SOURCE_DIR}/utils)

if(HAVE_CXX11 AND HAVE_HWLOC)
    cado_define_test(test_cpubinding.c
        ARGUMENTS
        --cpubinding ${CMAKE_CURRENT_SOURCE_DIR}/cpubind_config_file_for_tests.conf
        LIBRARIES utils bwc_base)
endif()

macro(test_dispatch_matmul testname)
    cado_define_test(${testname} SCRIPT
        ${CMAKE_CURRENT_SOURCE_DIR}/bwc-dispatch-matmul.sh
        bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc ${ARGN}
        TARGET_DEPENDENCIES bwc_matmul dispatch mf_scan mf_bal
        random_matrix)
endmacro()

# mingw has antique bash (3.1) with different regexp syntax, and having
# the script below work with both isn't a very funny perspective.
if (NOT MINGW)
test_dispatch_matmul(dispatch-matmul-65536-thr nrows=65536 density=2 thr=2x2)
if (HAVE_MPI)
test_dispatch_matmul(dispatch-matmul-65536-mpi nrows=65536 density=2 mpi=2x2)
endif()
endif()

if (WITH_MPI)
add_definitions(-DWITH_MPI)
endif()

cado_define_test(test_hello.c NO_DEFAULT_RUN LIBRARIES bwc_base)
set_target_properties(test_hello PROPERTIES
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")

add_executable(spmv_test spmv_test.c)
set_target_properties(spmv_test PROPERTIES
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
target_link_libraries(spmv_test bwc_base)

add_executable(short_matmul short_matmul.c)
target_link_libraries(short_matmul bwc_base)

set(test-spmv-binary-u64k1_alltests YES)
set(test-spmv-binary-u64k4_alltests YES)
set(test-spmv-prime-p_3_alltests YES)

macro(create_test_spmv test_spmv_base)
if(HAVE_PERL_DIGEST_MD5)
cado_define_test(
    TEST_NAME test-spmv-${test_spmv_base}
    SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-spmv.sh
    ARGUMENTS
    --matrix-size 1000 --density 10 --bindir ${CADO_NFS_BINARY_DIR}
    ${ARGN}
    NO_DEFAULT_RUN
    TARGET_DEPENDENCIES
    random_matrix spmv_test short_matmul mf_scan mf_bal bwc_matmul)


cado_divert_test(test-spmv-${test_spmv_base} 11 APPEND_ARGUMENTS thr=1x1)
if(${test-spmv-${test_spmv_base}_alltests})
cado_divert_test(test-spmv-${test_spmv_base} 12 APPEND_ARGUMENTS thr=1x2)
cado_divert_test(test-spmv-${test_spmv_base} 22 APPEND_ARGUMENTS thr=2x2)
cado_divert_test(test-spmv-${test_spmv_base} 13 APPEND_ARGUMENTS thr=1x3)
cado_divert_test(test-spmv-${test_spmv_base} 32 APPEND_ARGUMENTS thr=3x2)
cado_divert_test(test-spmv-${test_spmv_base} 23 APPEND_ARGUMENTS thr=2x3)
if(HAVE_MPI)
cado_divert_test(test-spmv-${test_spmv_base} 2122 APPEND_ARGUMENTS mpi=2x1 thr=2x2)
cado_divert_test(test-spmv-${test_spmv_base} 1322 APPEND_ARGUMENTS mpi=1x3 thr=2x2)
cado_divert_test(test-spmv-${test_spmv_base} 3222 APPEND_ARGUMENTS mpi=3x2 thr=2x2)
endif()
endif()
else()
    # Then this macro will expand to nothing, that's life.
endif()
endmacro()

cado_define_test(test_random_matrix
    PROGRAM ${CADO_NFS_BINARY_DIR}/linalg/bwc/krylov --random_matrix 20000,density=100 m=64 n=64 ys=0..64 prime=2 interval=100 end=200
    TARGET_DEPENDENCIES bwc_full
    )
cado_divert_test(test_random_matrix thr22 APPEND_ARGUMENTS thr=2x2)


foreach(arith_layer ${BWC_GF2_ARITHMETIC_BACKENDS})
    string_join(mm_backends " " ${BWC_GF2_MATMUL_BACKENDS})

    create_test_spmv(binary-${arith_layer}
        --backends "${mm_backends}"
        --arith-layer ${arith_layer})
endforeach()

foreach(arith_layer ${BWC_GFP_ARITHMETIC_BACKENDS})
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

cado_define_test(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-matcache.sh
    --matrix-size 1000 --density 10 --bindir
    ${CADO_NFS_BINARY_DIR}/linalg/bwc
    TARGET_DEPENDENCIES
    random_matrix build_matcache bench_matcache)

cado_define_test(TEST_NAME bwc_staged_krylov PROGRAM
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/krylov thr=2x2 nullspace=left
    interval=200  mn=64 prime=2  ys=0..64 start=0 end=200
    skip_online_checks=1 rebuild_cache=1 seed=1
    sequential_cache_build=1 random_matrix=10000,density=150,seed=1
    no_save_cache=1 verbose_flags=^all-bwc-sub-timings,^all-bwc-dispatch
    TARGET_DEPENDENCIES bwc_full)

# Tests which use bwc_base will need nwc_config.h
include_directories(${PROJECT_BINARY_DIR}/linalg/bwc)

cado_divert_test(test_hello 11 ARGUMENTS -v thr=1x1)
cado_divert_test(test_hello 12 ARGUMENTS -v thr=1x2)
cado_divert_test(test_hello 23 ARGUMENTS -v thr=2x3)

if (HAVE_MPI)
    # yes, the mpi= and thr= arguments are given twice here. This is
    # intentional, because our only logic for discovering the scheduling
    # environment is within bwc.pl, so we have to use it.
    cado_define_test(TEST_NAME test_hello_mpi PROGRAM
        ${PROJECT_BINARY_DIR}/linalg/bwc/bwc.pl
        NO_DEFAULT_RUN
        TARGET_DEPENDENCIES test_hello
        TEST_DEPENDENCIES test_hello)
    macro(test_hello_mpi code)
        cado_divert_test(test_hello_mpi ${code} ARGUMENTS
            :mpirun ${ARGN}
            -- ${CMAKE_CURRENT_BINARY_DIR}/test_hello -v ${ARGN}
            TEST_DEPENDENCIES test_hello)
    endmacro()
    test_hello_mpi(2122 mpi=2x1 thr=2x2)
    test_hello_mpi(1322 mpi=1x3 thr=2x2)
    test_hello_mpi(3222 mpi=3x2 thr=2x2)

cado_define_test(test_random_matrix_mpi
    PROGRAM
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/bwc.pl krylov --random_matrix 20000,density=100 m=64 n=64 ys=0..64 prime=2 interval=100 end=200
    TARGET_DEPENDENCIES bwc_full
    TEST_DEPENDENCIES test_random_matrix
    )
    cado_divert_test(test_random_matrix_mpi 22 APPEND_ARGUMENTS mpi=2x2)
    cado_divert_test(test_random_matrix_mpi 23 APPEND_ARGUMENTS mpi=2x3)
endif()

cado_define_test(test_bl_extraction.c LIBRARIES bwc_base bitlinalg)
set_target_properties(test_bl_extraction PROPERTIES
                LINKER_LANGUAGE CXX
                COMPILE_FLAGS "--mpi"
                LINK_FLAGS "--mpi")

cado_define_test(TEST_NAME test-blocklanczos
    PROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/test-blocklanczos.sh
        --density 10
        --bindir ${CADO_NFS_BINARY_DIR}/linalg/bwc
    NO_DEFAULT_RUN
    TARGET_DEPENDENCIES mf_scan mf_bal random_matrix blocklanczos
    bwc_full short_matmul)
cado_divert_test(test-blocklanczos right
    APPEND_ARGUMENTS --matrix-size "1000,2000" --nullspace RIGHT)
cado_divert_test(test-blocklanczos left
    APPEND_ARGUMENTS --matrix-size "2000,1000" --nullspace left)
cado_divert_test(test-blocklanczos left-thr21
    APPEND_ARGUMENTS --matrix-size "2000,1000" --nullspace left thr=2x1)
cado_divert_test(test-blocklanczos left-thr12
    APPEND_ARGUMENTS --matrix-size "90,69" --nullspace left thr=1x2)
cado_divert_test(test-blocklanczos left-thr23
    APPEND_ARGUMENTS --matrix-size "1500,1000" --nullspace left thr=2x3)
if (HAVE_MPI)
    cado_divert_test(test-blocklanczos left-mpi1-thr23
        APPEND_ARGUMENTS --matrix-size "2500,1500" --nullspace left mpi=2x3 thr=2x3)
    cado_divert_test(test-blocklanczos left-mpi2-thr23
        APPEND_ARGUMENTS --matrix-size "2500,1500" --nullspace left mpi=3x2 thr=2x3)
endif()

cado_define_test(TEST_NAME test-blocklanczos-staged
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_BINARY_DIR}/linalg/bwc/blocklanczos
    random_matrix=1200,1000,80,seed=1 m=64 n=64 ys=0..64
    no_save_cache=1 interval=10 nullspace=left
    TARGET_DEPENDENCIES test-blocklanczos_dependencies)

if(HAVE_PERL_DIGEST_MD5)
cado_define_test(TEST_NAME test-lingen-crc
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/tests/provide-wdir.sh --arg wdir --
    ${CADO_NFS_SOURCE_DIR}/tests/linalg/bwc/test-lingen-crc.sh
    bindir=${CADO_NFS_BINARY_DIR}/linalg/bwc
    seed=0 m=64 n=64 random_stem=30000 sequence_length=200
    expect_crc_pi=c29af350 expect_crc_F=04f4e167
    TARGET_DEPENDENCIES lingen)
endif()

cado_define_test(test-arith-modp.cpp LIBRARIES tests utils)
