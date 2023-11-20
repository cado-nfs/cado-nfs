#ifndef LINGEN_BMSTATUS_HPP_
#define LINGEN_BMSTATUS_HPP_
// IWYU pragma: no_include <algorithm>
#include <cstring>
#include <vector>
#include <tuple>
#include "lingen_bw_dimensions.hpp"
#include "lingen_hints.hpp"
#include "lingen_call_companion.hpp"
#include "tree_stats.hpp"
#include "select_mpi.h"               // for MPI_Comm
#include "cxx_mpz.hpp"

struct bmstatus {
    bw_dimensions d;
    unsigned int t = 0;
    std::vector<int> lucky;
    std::vector<unsigned int> delta;
    int done = 0;

    double t_basecase = 0;
    double t_mp = 0;
    double t_mul = 0;
    double t_cp_io = 0;

    // unsigned int lingen_threshold;
    // unsigned int lingen_mpi_threshold;
    
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm com[3]; /* [0]: MPI_COMM_WORLD, reordered.
                        [1]: row-wise
                        [2]: column-wise */

    lingen_hints hints;

    tree_stats stats;

    int depth() const { return stats.non_transition_depth(); }

    bmstatus(unsigned int m, unsigned int n, cxx_mpz const & p)
        : d(m, n, p)
    { /*{{{*/
        lucky.assign(m+n, 0);
        delta.assign(d.m + d.n, 0);
    }/*}}}*/
    bmstatus(unsigned int m, unsigned int n, mpz_srcptr p)
        : bmstatus(m, n, cxx_mpz(p)) {}

    void set_t0(unsigned int t0) {
        t = t0;
        delta.assign(d.m + d.n, t);
    }
    /* Attention: reloading a checkpoint invalidates this reference !! */
    lingen_call_companion & companion(int depth, size_t L);
    bool recurse(int depth, size_t L) {/*{{{*/
        return companion(depth, L).recurse();
    }/*}}}*/
    void display_deltas() const;
    bool recurse(size_t L) {/*{{{*/
        return companion(depth(), L).recurse();
    }/*}}}*/
    std::tuple<unsigned int, unsigned int> get_minmax_delta_on_solutions() const;
    unsigned int get_max_delta_on_solutions() const;
};

#endif	/* LINGEN_BMSTATUS_HPP_ */
