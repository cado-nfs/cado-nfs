#ifndef PLINGEN_HPP_
#define PLINGEN_HPP_

/* This contains some definitions for lingen mod p.
 *
 * This code is to be understood as scotch tape around an old and quirky
 * implementation.
 */

#include <cstdlib>
#include <vector>

#include "tree_stats.hpp"
#ifndef SELECT_MPFQ_LAYER_u64k1
#include "mpfq_layer.h"
#else
#include "mpfq_fake.hpp"
#endif
#include "select_mpi.h"
#include "lingen_hints.hpp"

/* TODO: Rename ! */
struct bw_dimensions {
    unsigned int m,n,nrhs;
    abfield ab;
};

struct bmstatus {/*{{{*/
    bw_dimensions d;
    unsigned int t;
    std::vector<int> lucky;

    double t_basecase;
    double t_mp;
    double t_mul;
    double t_cp_io;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm com[3]; /* [0]: MPI_COMM_WORLD, reordered.
                        [1]: row-wise
                        [2]: column-wise */

    lingen_hints hints;

    tree_stats stats;

    int depth() const { return stats.non_transition_depth(); }

    bmstatus(unsigned int m, unsigned int n)/*{{{*/
    {
        memset(&d, 0, sizeof(bw_dimensions));
        d.m = m;
        d.n = n;
        lucky.assign(m+n, 0);
    }/*}}}*/
    lingen_call_companion & companion(int depth, size_t L);
    bool recurse(int depth, size_t L) {/*{{{*/
        return companion(depth, L).recurse;
    }/*}}}*/
    bool recurse(size_t L) {/*{{{*/
        return companion(depth(), L).recurse;
    }/*}}}*/
};/*}}}*/

#endif	/* PLINGEN_HPP_ */
