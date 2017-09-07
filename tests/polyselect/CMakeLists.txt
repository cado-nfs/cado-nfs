# Build polyselect before running check

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
polyselect_test(7 -N ${example_N} degree=7 P=1000 admax=1000 nq=10)

if(DEFINED ENV{CHECKS_EXPENSIVE})
# in the following example, there are 6 roots for q=241, but only 0 or 2 roots
# for the other small primes
set(example_N 412023436986659543855531365332575948179811699844327982845455626433876445565248426198098870423161841879261420247188869492560931776375033421130982397485150944909106910269861031862704114880866970564902903653658867433731720813104105190864254793282601391257624033946373269391)
polyselect_test(8 -N ${example_N} degree=6 P=1000000 admin=3720 admax=3721 nq=7776)
endif()

## we disable the following test since twocubics is not in "make dist"
## test also the twocubics binary
# cado_define_test(twocubics PROGRAM
#    ${CADO_NFS_BINARY_DIR}/polyselect/twocubics -n 999073468111577057576445816581 -P 500 -incr 30 -admax 1000
#    TARGET_DEPENDENCIES twocubics)

cado_define_test(test_auxiliary.c LIBRARIES utils tests polyselect_common)
cado_define_test(test_rho.c LIBRARIES polyselect_common tests)
