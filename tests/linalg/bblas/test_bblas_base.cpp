#include "cado.h" // IWYU pragma: keep
#include "bblas_level3a.hpp"
#include "test_bblas_base.hpp"
#include "memory.h"

int test_bblas_base::test_accel;

test_bblas_base::test_bblas_base(unsigned int nmax) : nmax(nmax) {/*{{{*/
    gmp_randinit_default(rstate);
    xr = (uint64_t *) malloc_aligned(nmax * sizeof(uint64_t), 64);
    r =  (uint64_t *) malloc_aligned(nmax * sizeof(uint64_t), 64);
    a =  (uint64_t *) malloc_aligned(nmax * sizeof(uint64_t), 64);
    b =  (uint64_t *) malloc_aligned(nmax * sizeof(uint64_t), 64);
#ifdef  HAVE_M4RI
    R = mzd_init(nmax, 64);
    A = mzd_init(nmax, 64);
    R64 = mzd_init(64, 64);
    A64 = mzd_init(64, 64);
    W = mzd_init(64, 64);
    WT = mzd_init(64, 64);
#endif  /* HAVE_M4RI */

    /* a, w, wt are assumed constant */
    memfill_random(a, (nmax) * sizeof(uint64_t), rstate);
    mat64_fill_random(w, rstate);
    mat64_transpose(wt, w);

#ifdef  HAVE_M4RI
    mzd_set_mem(A, a, nmax);
    if (nmax >= 64) {
        mzd_set_mem(A64, a, 64);
    }
    mzd_set_mem(W, w.data(), 64);
    mzd_set_mem(WT, wt.data(), 64);
#endif  /* HAVE_M4RI */
}/*}}}*/

test_bblas_base::~test_bblas_base() {/*{{{*/
    free_aligned(xr);
    free_aligned(r);
    free_aligned(a);
    free_aligned(b);
#ifdef  HAVE_M4RI
    mzd_free(R);
    mzd_free(A);
    mzd_free(R64);
    mzd_free(A64);
    mzd_free(W);
    mzd_free(WT);
#endif  /* HAVE_M4RI */
    gmp_randclear(rstate);
}/*}}}*/
bool test_bblas_base::matches(std::string const & s, tags_t const & T, bool & match)/*{{{*/
{
    if (s == "all") return match = true;
    for(auto const & t : T) {
        if (s == t)
            return match = true;
    }
    return false;
}/*}}}*/
