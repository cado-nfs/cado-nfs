#ifndef CADO_LINALG_BWC_BLOCKLANCZOS_HPP
#define CADO_LINALG_BWC_BLOCKLANCZOS_HPP

#include <memory>

#include <gmp.h>                 // for gmp_randclear, gmp_randinit_default

#include "gmp_aux.h"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"
#include "bblas_mat64.hpp"
#include "bit_vector.h"

struct blstate {
    /* We'll need several intermediary n*n matrices. These will be
     * allocated everywhere (we set flags to 0, so as to avoid having
     * shared vectors) */
    mat64 L[3];

    std::unique_ptr<arith_generic> A;
    matmul_top_data mmt;

    std::unique_ptr<arith_cross_generic> AxA;

    mmt_vec y, my;

    mmt_vec V[3];
    bit_vector D[3];

    cxx_gmp_randstate rstate;

    /* Here are the semantics of the data fields above.
     *
     * For iteration n, we let n0 = n%3, n1 = (n-1)%3, n2 = (n-2)%3.
     *
     * V[n0] is the V vector which is the input to iteration n. It does
     *       not consist of independent  vectors.
     * D[n0] is the extracted set of columns, computed from V[n0] (and
     *       also from D[n1]. Together, V[n0] and D[n0] allow to compute
     *       the basis of the n-th sub vector space W, although no
     *       explicit data is reserved to its storage.
     * L[n0] is computed from V[n0] and D[n0], and is some sort of local
     *       inverse (iirc, it's W_i^{inv} in most accounts).
     *
     * of course *[n1] and *[n2] are the same for the previous steps.
     *
     * In order to compute the vector for step n + 1, the data to be used
     * is V[*], L[*], and D[n0, n1]
     */
    blstate(parallelizing_info_ptr pi, cxx_param_list & pl);
    blstate(blstate const&) = delete;
    blstate& operator=(blstate const&) = delete;
    ~blstate();
    void set_start();
    void load(unsigned int iter);
    void save(unsigned int iter);
    void save_result(unsigned int iter);
    int operator()(parallelizing_info_ptr pi);
};

#endif	/* LINALG_BWC_BLOCKLANCZOS_HPP_ */
