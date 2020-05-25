#ifndef TEST_BBLAS_BASE_HPP_
#define TEST_BBLAS_BASE_HPP_

#include "bblas_mat64.hpp"
#include <gmp.h>
#include <cstdint>
#include <vector>
#include <string>
#include "gmp_aux.h"
#ifdef HAVE_M4RI
/* To test against m4ri routines, include a checkout of
 * https://bitbucket.org/malb/m4ri.git under linalg/m4ri, and run
 * "autoreconf -i" there ; cado-nfs cmake logic then detects it and
 * enables the corresponding code here (and in a few other places in
 * test_bblas). */
#include "m4ri/m4ri.h"
#endif

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
    mzd_t *R64;
    mzd_t *A64;
    mzd_t *W;
    mzd_t *WT;

    /* These functions are defined in test_bblas_m4ri.cpp */
    static void mzd_set_mem(mzd_t * M, const uint64_t * s, unsigned int n);
    static void mzd_set_memT(mzd_t * M, const uint64_t * s, unsigned int n);
    static void mzd_check_mem(mzd_t * M, uint64_t * s, unsigned int n);
#endif  /* HAVE_M4RI */

    static int test_accel;

    typedef std::vector<std::string> tags_t;

    test_bblas_base(unsigned int nmax);
    ~test_bblas_base();

    void set_seed(int seed) { gmp_randseed_ui(rstate, seed); }

    bool matches(std::string const & s, tags_t const & T, bool & match);
};

#endif	/* TEST_BBLAS_BASE_HPP_ */
