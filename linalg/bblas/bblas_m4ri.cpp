#include "cado.h"

#ifdef  HAVE_M4RI
void my_mzd_randomize(mzd_t * A)/*{{{*/
{
    for(size_t i = 0 ; i < A->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) A->rows[i];
        for(size_t j = 0 ; j < A->width ; j++) {
            ptr[j] = rand64();
        }
    }
}/*}}}*/

void mzd_set_mem(mzd_t * M, const uint64_t * s, unsigned int n)/*{{{*/
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        ptr[0] = sometimes_bitrev(s[i]);
    }
}/*}}}*/

void mzd_set_memT(mzd_t * M, const uint64_t * s, unsigned int n)/*{{{*/
{
    mzd_set_mem(M, s, n);
    mzd_transpose(M, M);
}/*}}}*/

void mzd_check_mem(mzd_t * M, uint64_t * s, unsigned int n)/*{{{*/
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        if (ptr[0] != sometimes_bitrev(s[i])) {
            fprintf(stderr, "Rows %zu differ: %016" PRIx64 " != %016" PRIx64 "\n",
                    i, sometimes_bitrev(ptr[0]), s[i]);
            abort();
        }
    }
}/*}}}*/

void mzd_check_memT(mzd_t * M, uint64_t * s, unsigned int n)/*{{{*/
{
    mzd_t * tmp = mzd_transpose(NULL, M);
    mzd_check_mem(M, s, n);
    mzd_free(tmp);
}/*}}}*/
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


