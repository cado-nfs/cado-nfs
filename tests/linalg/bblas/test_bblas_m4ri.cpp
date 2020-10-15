#include "cado.h"
#include "bblas.hpp"
#include <cinttypes>
#include "test_bblas_level4.hpp"
#include "time_bblas_common.hpp"

#ifdef  HAVE_M4RI
/* m4ri types are defined in m4ri/misc.h
 * rci_t: rows and column indices
 * wi_t: "indices for array of words that make up a row"
 * word: typical data structure to represent packed bits.
 *
 * rci_t and wi_t are currently int, and will probably remain
 * undistinguishable anyway.
 *
 * word is currently uint64_t, and we are really assuming that this stays
 * so
 */
static_assert(std::is_same<word,uint64_t>::value, "m4ri word MUST be uint64_t");

void my_mzd_randomize(mzd_t * A, gmp_randstate_t rstate)/*{{{*/
{
    for(rci_t i = 0 ; i < A->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) A->rows[i];
        memfill_random(ptr, A->width * sizeof(uint64_t), rstate);
    }
}/*}}}*/

void test_bblas_base::mzd_set_mem(mzd_t * M, const uint64_t * s, unsigned int n)/*{{{*/
{
    ASSERT_ALWAYS(M->width == 1);
    ASSERT_ALWAYS(M->nrows == (rci_t) n);
    for(rci_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        /* back in the days (before version 20111203 at least), m4ri had
         * limbs in the wrong order */
        ptr[0] = s[i];
    }
}/*}}}*/

void test_bblas_base::mzd_set_memT(mzd_t * M, const uint64_t * s, unsigned int n)/*{{{*/
{
    mzd_set_mem(M, s, n);
    mzd_transpose(M, M);
}/*}}}*/

void test_bblas_base::mzd_check_mem(mzd_t * M, uint64_t * s, unsigned int n)/*{{{*/
{
    ASSERT_ALWAYS(M->width == 1);
    ASSERT_ALWAYS(M->nrows >= (rci_t) n);
    for(unsigned int i = 0 ; i < n ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        /* back in the days (before version 20111203 at least), m4ri had
         * limbs in the wrong order */
        if (ptr[0] != s[i]) {
            fprintf(stderr, "Rows %d differ: %016" PRIx64 " != %016" PRIx64 "\n",
                    (int) i, ptr[0], s[i]);
            abort();
        }
    }
}/*}}}*/

#if 0
void mzd_check_memT(mzd_t * M, uint64_t * s, unsigned int n)/*{{{*/
{
    mzd_t * tmp = mzd_transpose(NULL, M);
    mzd_check_mem(M, s, n);
    mzd_free(tmp);
}/*}}}*/
#endif

#endif  /* HAVE_M4RI */





/*************/


#ifdef  HAVE_M4RI
static inline void mzd_mypluq(mzd_t * LU, mzd_t * M, mzp_t * P, mzp_t * Q, int c)
{
    mzd_copy(LU, M);
    mzd_pluq(LU,P,Q,c);
}
static inline void mzd_myechelonize_m4ri(mzd_t * E, mzd_t * M, int full, int k)
{
    mzd_copy(E, M);
    mzd_echelonize_m4ri(E,full,k);
}
static inline void mzd_myechelonize_pluq(mzd_t * E, mzd_t * M, int full)
{
    mzd_copy(E, M);
    mzd_echelonize_pluq(E,full);
}
#endif


#ifdef  HAVE_M4RI
void test_bblas_level4::m4ri_plu_tests(int n __attribute__((unused)))
{
    mzd_t * M;
    mzd_t * LU;
    mzp_t * P, *Q;
    M = mzd_init(n, n);
#if 0
    mzd_set_mem(M, m, n);
    uint64_t * m = new uint64_t[n*n/64];
    memfill_random(m, (64) * sizeof(uint64_t), rstate);
    delete[] m;
#else
    my_mzd_randomize(M, rstate);
#endif
    LU = mzd_init(n, n);
    P = mzp_init(n);
    Q = mzp_init(n);
    TIME1N(2, mzd_mypluq, LU, M, P, Q, 0);
    TIME1N(2, mzd_myechelonize_m4ri, LU, M, 0, 0);
    TIME1N(2, mzd_myechelonize_pluq, LU, M, 0);
    mzd_free(M);
    mzd_free(LU);
    mzp_free(P);
    mzp_free(Q);
}
#endif



