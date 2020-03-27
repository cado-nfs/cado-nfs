#include "cado.h"
#include "test_bblas_base.hpp"

int test_bblas_base::test_accel;

test_bblas_base::test_bblas_base(unsigned int nmax) : nmax(nmax) {/*{{{*/
    gmp_randinit_default(rstate);
    xr = new uint64_t[nmax];
    r =  new uint64_t[nmax];
    a =  new uint64_t[nmax];
    b =  new uint64_t[nmax];
    w =  new uint64_t[64];
    wt = new uint64_t[64];
#ifdef  HAVE_M4RI
    R = mzd_init(nmax, 64);
    A = mzd_init(nmax, 64);
    W = mzd_init(64, 64);
    WT = mzd_init(64, 64);
#endif  /* HAVE_M4RI */

    /* a, w, wt are assumed constant */
    memfill_random(a, (nmax) * sizeof(uint64_t), rstate);
    memfill_random(w, (64) * sizeof(uint64_t), rstate);
    mat64_transpose(wt, w);

#ifdef  HAVE_M4RI
    mzd_set_mem(A, a, nmax);
    mzd_set_mem(W, w, 64);
    mzd_set_mem(WT, wt, 64);
#endif  /* HAVE_M4RI */
}/*}}}*/

test_bblas_base::~test_bblas_base() {/*{{{*/
    delete[] xr;
    delete[] r;
    delete[] a;
    delete[] b;
    delete[] w;
    delete[] wt;
#ifdef  HAVE_M4RI
    mzd_free(D->R);
    mzd_free(D->A);
    mzd_free(D->W);
    mzd_free(D->WT);
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
