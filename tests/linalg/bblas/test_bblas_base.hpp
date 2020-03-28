#ifndef TEST_BBLAS_BASE_HPP_
#define TEST_BBLAS_BASE_HPP_

#include "bblas.hpp"
#include "test_bblas_base.hpp"
#include <gmp.h>
#include <cstdint>
#include <vector>
#include <string>
#include <set>
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
    mat64 w;
    mat64 wt;

#ifdef  HAVE_M4RI
    mzd_t *R;
    mzd_t *A;
    mzd_t *W;
    mzd_t *WT;
#endif  /* HAVE_M4RI */

    static int test_accel;

    typedef std::vector<std::string> tags_t;

    test_bblas_base(unsigned int nmax);
    ~test_bblas_base();

    void set_seed(int seed) { gmp_randseed_ui(rstate, seed); }

    bool matches(std::string const & s, tags_t const & T, bool & match);
};

#endif	/* TEST_BBLAS_BASE_HPP_ */
