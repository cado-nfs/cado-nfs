
cado_define_test(test_mod.cpp LIBRARIES utils tests)
cado_define_test(test_lucas.cpp LIBRARIES utils tests)
cado_define_test(test_modint.cpp LIBRARIES utils tests)

# Timings are only meaningful when run by hand, see the source; as a
# test, this checks that the old arith layer and arithxx agree.
cado_define_test(bench_arith_vs_arithxx
    SOURCES bench_arith_vs_arithxx.cpp
        bench_arith_vs_arithxx_old_ul.cpp
        bench_arith_vs_arithxx_old_15ul.cpp
        bench_arith_vs_arithxx_old_2ul2.cpp
        bench_arith_vs_arithxx_old_mpz.cpp
    ARGUMENTS -quick
    LIBRARIES utils)
