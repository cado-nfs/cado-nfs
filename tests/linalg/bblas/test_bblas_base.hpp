#ifndef TEST_BBLAS_BASE_HPP_
#define TEST_BBLAS_BASE_HPP_

#include "bblas.hpp"
#include "test_bblas_base.hpp"
#include <gmp.h>
#include <cstdint>
#include <vector>
#include <string>
#include "gmp_aux.h"

static inline uint64_t uint64_random(gmp_randstate_t rstate)
{
    uint64_t a;
    memfill_random(&a, sizeof(uint64_t), rstate);
    return a;
}

struct test_bblas_base {
    gmp_randstate_t rstate;

    unsigned int nmax;

    uint64_t * xr;
    uint64_t * r;
    uint64_t * a;
    uint64_t * b;
    uint64_t * w;
    uint64_t * wt;

#ifdef  HAVE_M4RI
    mzd_t *R;
    mzd_t *A;
    mzd_t *W;
    mzd_t *WT;
#endif  /* HAVE_M4RI */

    typedef std::vector<std::string> tags_t;

    test_bblas_base(unsigned int nmax) : nmax(nmax) {/*{{{*/
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
    ~test_bblas_base() {/*{{{*/
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

    void set_seed(int seed) {
        gmp_randseed_ui(rstate, seed);
    }

    bool matches(std::string const & s, tags_t const & T, bool & match)/*{{{*/
    {
        if (s == "all") return match = true;
        for(auto const & t : T) {
            if (s == t)
                return match = true;
        }
        return false;
    }/*}}}*/
};

#endif	/* TEST_BBLAS_BASE_HPP_ */
