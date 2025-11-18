#include "cado.h" // IWYU pragma: keep

#include <gmp.h>

#include "bblas_gauss.h"
#include "bblas_mat64.hpp"
#include "bblas_level3a.hpp"      // for mat64_fill_random
#include "bblas_level4.hpp"
#include "bblas_perm_matrix.hpp"
#include "gmp_aux.h"              // for memfill_random
#include "test_bblas_base.hpp"    // for test_bblas_base, test_bblas_base::t...
#include "test_bblas_level4.hpp"
#include "time_bblas_common.hpp"

static int gauss_MN_C(unsigned int bM, unsigned int bN, gmp_randstate_t rstate)
{
    constexpr const unsigned int B = mat64::width;
    unsigned int const M = B * bM;
    unsigned int const N = B * bN;
    mat64 * mm = mat64::alloc(bM * bN);
    memfill_random(mm, bM * bN * sizeof(mat64), rstate);
    int const r = kernel((mp_limb_t*)mm, nullptr, M, N, N/ULONG_BITS, M/ULONG_BITS);
    mat64::free(mm, bM * bN);
    return r;
}


test_bblas_base::tags_t test_bblas_level4::gauss_tags { "gauss", "l4" };
void test_bblas_level4::gauss() {
    mat64 m;
    mat64 e;
    mat64 mm;
    mat64 l, u, p;
    mat64_fill_random(m, rstate);
    mat64_fill_random(e, rstate);
    mat64_fill_random(mm, rstate);
    mat64 m4[4];
    mat64 u4[4];
    memfill_random(m4, 4 * sizeof(mat64), rstate);
    // printf("-- for reference: best matrix mult, 64x64 --\n");
    // TIME1(2, mul_6464_6464, mm, e, m);
    // TIME1(2, mul_N64_T6464, mm, e, m, 64);
    // TIME1(2, gauss_6464_C, mm, e, m);
    // TIME1(2, gauss_6464_imm, mm, e, m);
    // TIME1(2, PLUQ64_inner, nullptr, l, u, m, 0);
    int phi[128];
    {
        perm_matrix p, q;
        mat64 m[4], l[4], u[4];
        memfill_random(m, 4 * sizeof(mat64), rstate);
        perm_matrix_init(p, 128);
        perm_matrix_init(q, 128);
        TIME1(2, PLUQ128, p, l, u, q, m);
        perm_matrix_clear(p);
        perm_matrix_clear(q);
    }
    int n=2;
    {
        cxx_gmp_randstate rstate2;
        TIME1N(2, memfill_random, m4, n*sizeof(mat64), rstate2);
    }
    TIME1N_SPINS(, 2, PLUQ64_n, phi, l, u4, m4, 64*n);
    {
        cxx_gmp_randstate rstate2 = rstate;
        TIME1N_SPINS(memfill_random(m4, n*sizeof(mat64), rstate2), 2, PLUQ64_n, phi, l, u4, m4, 64*n);
    }
    TIME1(2, LUP64_imm, l, u, p, m);
    TIME1(2, full_echelon_6464_imm, mm, e, m);
    TIME1(2, gauss_128128_C, m4);
    TIME1(2, gauss_MN_C, 10, 2, rstate);
    TIME1(2, gauss_MN_C, 2, 10, rstate);
    TIME1(2, gauss_MN_C, 100, 2, rstate);
    TIME1(2, gauss_MN_C, 2, 100, rstate);
    TIME1(2, gauss_MN_C, 1000, 2, rstate);
    TIME1(2, gauss_MN_C, 2, 1000, rstate);
    TIME1(2, gauss_MN_C, 32, 32, rstate);
#ifdef  HAVE_M4RI
    m4ri_plu_tests(64);
    m4ri_plu_tests(128);
    m4ri_plu_tests(256);
    m4ri_plu_tests(512);
    m4ri_plu_tests(1024);
#endif
}

