#include "cado.h" // IWYU pragma: keep
#include <stdint.h>                   // for uint64_t
#ifndef LINGEN_BINARY
#include <cmath>
#include <gmp.h>                      // for mpz_get_d, mpz_fdiv_q_ui, mpz_s...
#include <stddef.h>                   // for size_t
#include "cxx_mpz.hpp"
#endif
#include "macros.h"                   // for ASSERT_ALWAYS
#include "lingen_average_matsize.hpp"
#include "lingen_matpoly_select.hpp"

/*{{{ avg_matsize */
template<bool over_gf2>
double avg_matsize(matpoly::arith_hard *, unsigned int m, unsigned int n, int ascii);

template<>
double avg_matsize<true>(matpoly::arith_hard *, unsigned int m, unsigned int n, int ascii)
{
    ASSERT_ALWAYS(!ascii);
    ASSERT_ALWAYS((m*n) % 64 == 0);
    unsigned int nwords = m * n / 64;
    return nwords * sizeof(uint64_t);
}

#ifndef LINGEN_BINARY
template<>
double avg_matsize<false>(matpoly::arith_hard * ab, unsigned int m, unsigned int n, int ascii)
{
    if (!ascii) {
        /* Easy case first. If we have binary input, then we know a priori
         * that the input data must have size a multiple of the element size.
         */
        return ab->vec_elt_stride(m * n);
    }

    /* Ascii is more complicated. We're necessarily fragile here.
     * However, assuming that each coefficient comes with only one space,
     * and each matrix with an extra space (this is how the GPU program
     * prints data -- not that this ends up having a considerable impact
     * anyway...), we can guess the number of bytes per matrix. */

    /* Formula for the average number of digits of an integer mod p,
     * written in base b:
     *
     * (k-((b^k-1)/(b-1)-1)/p)  with b = Ceiling(Log(p)/Log(b)).
     */
    double avg;
    cxx_mpz a;
    double pd = mpz_get_d(ab->characteristic());
    unsigned long k = ceil(log(pd)/log(10));
    unsigned long b = 10;
    mpz_ui_pow_ui(a, b, k);
    mpz_sub_ui(a, a, 1);
    mpz_fdiv_q_ui(a, a, b-1);
    avg = k - mpz_get_d(a) / pd;
    // printf("Expect roughly %.2f decimal digits for integers mod p.\n", avg);
    double matsize = (avg + 1) * m * n + 1;
    // printf("Expect roughly %.2f bytes for each sequence matrix.\n", matsize);
    return matsize;
}
#endif

double average_matsize(matpoly::arith_hard * ab, unsigned int m, unsigned int n, int ascii)
{
    return avg_matsize<matpoly::over_gf2>(ab, m, n, ascii);
}
/*}}}*/

