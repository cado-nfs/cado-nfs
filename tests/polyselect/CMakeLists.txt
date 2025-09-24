# Build polyselect before running check

include_directories(${PROJECT_SOURCE_DIR}/polyselect)

macro(polyselect_test n)
    cado_define_test(test_polsel_${n}
        PROGRAM ${CADO_NFS_BINARY_DIR}/polyselect/polyselect ${ARGN}
        TARGET_DEPENDENCIES polyselect)
endmacro()
set(example_N 90377629292003121684002147101760858109247336549001090677693)
polyselect_test(2 -N ${example_N} degree=2 P=1000 admax=100  nq=8)
polyselect_test(3 -N ${example_N} degree=3 P=1000 admax=500  nq=10)
polyselect_test(4 -N ${example_N} degree=4 P=1000 admax=1000 nq=10)
polyselect_test(5 -N ${example_N} degree=5 P=1000 admax=1000 nq=10)
polyselect_test(6 -N ${example_N} degree=6 P=1000 admax=1000 nq=10)
polyselect_test(7 -N ${example_N} degree=7 P=1000 admax=1000 nq=10 TIMEOUT 60)

if(DEFINED ENV{CHECKS_EXPENSIVE})
# in the following example, there are 6 roots for q=241, but only 0 or 2 roots
# for the other small primes
set(example_N 412023436986659543855531365332575948179811699844327982845455626433876445565248426198098870423161841879261420247188869492560931776375033421130982397485150944909106910269861031862704114880866970564902903653658867433731720813104105190864254793282601391257624033946373269391)
polyselect_test(8 -N ${example_N} degree=6 P=1000000 admin=3720
    admax=3721 nq=7776 TIMEOUT 120
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
)
endif()

## we disable the following test since twocubics is not in "make dist"
## test also the twocubics binary
# cado_define_test(twocubics PROGRAM
#    ${CADO_NFS_BINARY_DIR}/polyselect/twocubics -n 999073468111577057576445816581 -P 500 -incr 30 -admax 1000
#    TARGET_DEPENDENCIES twocubics)

cado_define_test(test_ropt_auxiliary.cpp LIBRARIES utils tests polyselect_common)
cado_define_test(test_rho.c LIBRARIES polyselect_common tests)
# This one is more expensive than a few others, as it takes 5 seconds on
# a laptop, and can easily take tens of seconds on a container that runs
# on a busy machine.
cado_define_test(test_polyselect_shash.cpp
    ${PROJECT_SOURCE_DIR}/polyselect/polyselect_shash.c
    ${PROJECT_SOURCE_DIR}/polyselect/polyselect_poly_header.c
    ${PROJECT_SOURCE_DIR}/polyselect/polyselect_match.c
    LIBRARIES polyselect_common tests
    TIMEOUT 90
)

cado_define_test(test_ropt_adv
    PROGRAM ${CADO_NFS_BINARY_DIR}/polyselect/polyselect_ropt --adv
    TARGET_DEPENDENCIES polyselect_ropt
    NO_DEFAULT_RUN
    TIMEOUT 60
    EXPECT_SHA1 42
    SHA1_ON_REGEXP_LINES "^[^#]"
    )

# set(ropt_sha1_c80 68c14f886966071c2b79983f94f34cfd218d123b)
# set(ropt_sha1_c75 8d08dd3fd67e358e291c0aeb7cc5935682b91758)
# # 38deff1c0c2d303cc87960e4f98d0fbb3a58ef24 is on freebsd. See #30098
# set(ropt_sha1_c70 3164916d4b4bc5512e9ae283d6bb604fe6ad3930,38deff1c0c2d303cc87960e4f98d0fbb3a58ef24)
# set(ropt_sha1_c65 5c7adf825238c1a75fcb33567b2460e2e3b11486)
# set(ropt_sha1_c60 d945f8b4b1f51b07d5efdbcc8f32d74d185f7f76)

# Higher precision skewness computation has changed data a bit.
set(ropt_sha1_c80 abd946d0837e008ca1775f3dc9bebe02b6d9087e)
set(ropt_sha1_c75 0580622f2b8e6ca129ef63d6f5aacb045586ab47)
set(ropt_sha1_c70 c59e7eeb346c88b1374df40a89f8a10e54a7a32a,d335f74d729e599c30d73b81de0a738d89101656)
set(ropt_sha1_c65 ad7ba67781efa987b32e2d9fb891e64dbdeb1a5e)
set(ropt_sha1_c60 847a5d86728e33f191ac7e58fa84fc278ef570d6)


foreach(example c60 c65 c70 c75 c80)
    cado_divert_test(test_ropt_adv ${example}
        APPEND_ARGUMENTS
        -f ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/${example}.poly
        EXPECT_SHA1 ${ropt_sha1_${example}}
    )
endforeach()

cado_divert_test(test_ropt_adv c80msieve
    APPEND_ARGUMENTS
    -fm ${CMAKE_CURRENT_SOURCE_DIR}/c80msieve.sample
    -d 4
    -n 12022203778835263712651529862414424087010841660557141141978630366614660861225497
    -v -v -amin -16 -amax 16 -e 2 2 1 -bmax 10 -cmax 10000
    # EXPECT_SHA1 cf2da7482ab71b79b298204a66fcccdb4593c633
    EXPECT_SHA1 a22b68bd09a10122ea57d996650598e24f5a2501
)

# This should take only a few seconds. Under certain circumstances, it
# seems to take a lot more...
cado_divert_test(test_ropt_adv deg6
    APPEND_ARGUMENTS
    -f ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/rsa768.poly
    -v -v -d 6   -amin -4 -amax 4 -e 1 1 -bmax 6 -cmax 12
    TIMEOUT 240
    # EXPECT_SHA1 d451c895d0fbc5cf1a58253235d05151f15c159b
    EXPECT_SHA1 8f8e34ff85a9369b52bdc961514bfa9da52a6209
    ${AVOID_CONCURRENT_UNDER_COVERAGE}
    )

if(DEFINED ENV{CHECKS_EXPENSIVE})
    # We need this to cover the ropt_quadratic code.
    # Honestly, it's way too long, it takes two minutes or
    # so on my laptop. Furthermore, it looks very odd, since in this
    # example I have the impression that the output of the "optimization"
    # is actually not as good as the input... (and it's also quite
    # bizarre that on the rsa768 polynomial, the output is exactly the
    # same as test_ropt_adv_deg6)
    #
    # note that most of the runtime-affecting parameters in here are
    # actually hardcoded, e.g. in polyselect/ropt_param.c ; there's very
    # little that can be done to change the running time, short of
    # changing the size of the number to be factored.
    cado_divert_test(test_ropt_adv deg5large
        APPEND_ARGUMENTS
        -f ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/rsa155.poly
        -v -v
        -amin -16 -amax 16 -bmax 10 -cmax 10000
        -e 3 5 2 1
        # EXPECT_SHA1 e592a8f38c898163fe0e5995580df6dc2694fcd3
        EXPECT_SHA1 59f634b37e478cf4e57780728c746ef6333ac483
        # this one is for rsa768.poly (degree 6)
        # EXPECT_SHA1 d451c895d0fbc5cf1a58253235d05151f15c159b

        TIMEOUT 360
    )

    # takes less than a minute on my laptop. The ropteffort parameter
    # isn't much tunable, unfortunately. I'm not sure that we really need
    # this addition in order to cover new pieces of code.
    cado_define_test(test_ropt_list_rsa250
        PROGRAM
        ${CADO_NFS_BINARY_DIR}/polyselect/polyselect_ropt
        -inputpolys ${CADO_NFS_SOURCE_DIR}/tests/polyselect/ropt.rsa250.txt
        -Bf 2.100e+09 -Bg 1.800e+09 -area 2.400e+19
        -ropteffort 0.05
        # EXPECT_SHA1 fa313790a88f19182f4c4f7d1dc7739d05adda34
        EXPECT_SHA1 8264a5eb5f2788d68ab87bcdde84c32e306edb44
        SHA1_ON_REGEXP_LINES "^[^#]"
        TARGET_DEPENDENCIES polyselect_ropt
        TIMEOUT 360
        )
endif()

cado_define_test(
    gfpkdlpolyselect_test.c
    table_t_Py_f_deg4_type0_h1_t-200--200.c
    ${PROJECT_SOURCE_DIR}/polyselect/gfpkdlpolyselect.c
    LIBRARIES polyselect_common utils
    NO_DEFAULT_RUN
    )
    
cado_divert_test(gfpkdlpolyselect_test 1 -p 10000000019 -k 2 -label test01)
cado_divert_test(gfpkdlpolyselect_test 2 -p 10000000019 -k 2 -label test02)
cado_divert_test(gfpkdlpolyselect_test 3 -p 10000000019 -k 2 -label test03)
foreach(p 10000000019 10000000033 10000000061 10000000069 10000000097)
    cado_divert_test(gfpkdlpolyselect_test ${p} APPEND_ARGUMENTS -k 2 -p ${p} -label test${p})
endforeach()

cado_define_test(test_polyselect_gfpn
        PROGRAM
        ${CADO_NFS_BINARY_DIR}/polyselect/polyselect_gfpn
        ARGUMENTS
        -n 2 -p
        TARGET_DEPENDENCIES polyselect_gfpn
        NO_DEFAULT_RUN
        )
cado_divert_test(test_polyselect_gfpn 0 APPEND_ARGUMENTS 750955589968013246993551)
cado_divert_test(test_polyselect_gfpn 1 APPEND_ARGUMENTS 325296040168596930710621)
cado_divert_test(test_polyselect_gfpn 2 APPEND_ARGUMENTS 584064835518808197191269)
cado_divert_test(test_polyselect_gfpn 3 APPEND_ARGUMENTS 45086221481123679612257)

include(CMakeLists-nodist.txt OPTIONAL)

cado_define_test(test_sopt PROGRAM
    ${CADO_NFS_BINARY_DIR}/polyselect/sopt
    -sopteffort 2
    NO_DEFAULT_RUN
    SHA1_ON_REGEXP_LINES "^[^#]"
    TARGET_DEPENDENCIES sopt
)

cado_divert_test(test_sopt 1
    APPEND_ARGUMENTS
    -inputpolys ${CMAKE_CURRENT_SOURCE_DIR}/sopt.input1.txt
    # EXPECT_SHA1 465209edaa2c8bbad4e7254a7d6fcf2f121a71b8
    # the sha1 above is when skewness is limited to 10 bits of precision.
    # There is little point in placing this restriction now, so we get
    # this instead.
    EXPECT_SHA1 f30ec6c6f2d6d2c90a6e8e96d5b822ae5821fab5
)
cado_divert_test(test_sopt 2
    APPEND_ARGUMENTS
    -inputpolys ${CMAKE_CURRENT_SOURCE_DIR}/sopt.input2.txt
    # EXPECT_SHA1 a3aadd110944689565adb70ae7456b69d8d97dbc
    # same as above
    EXPECT_SHA1 83a512c121fec35711d03a91901e04f1aeee21b6
)


cado_define_test(polyselect_add_score
    PROGRAM
        ${CADO_NFS_BINARY_DIR}/polyselect/score
        -Bf 1e+04 -Bg 2e+05 -area 5e+10
        ${CADO_NFS_SOURCE_DIR}/tests/misc/c60.poly
    EXPECT_SHA1 27954264d0ab7e6bf8bbc07a8376ff9c63649ba8
    TARGET_DEPENDENCIES score)

cado_define_test(polyselect_add_skewness
    PROGRAM
        ${CADO_NFS_BINARY_DIR}/polyselect/skewness
        ${CADO_NFS_SOURCE_DIR}/tests/misc/c60.poly
    EXPECT_SHA1 70639bbafd1b6d4ac3bbcb4a47f621fa446c267f
    TARGET_DEPENDENCIES skewness)
