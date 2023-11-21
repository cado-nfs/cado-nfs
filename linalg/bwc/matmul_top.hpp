#ifndef MATMUL_TOP_HPP_
#define MATMUL_TOP_HPP_

#include <stddef.h>
#include <stdint.h>

#include <string>
#include <array>
#include <gmp.h>                 // for gmp_randstate_t
struct timing_data;

#include "select_mpi.h"
#include "parallelizing_info.hpp"
#include "matmul.hpp"
#include "params.h"
#include "balancing.hpp"

struct matmul_top_matrix {
    // global stuff.
    //
    // n[0] is the number of rows, includding padding.
    // n[1] is the number of columns, includding padding.
    //
    // n0[] is without padding.
    //
    // Note that a communicator in direction 0 handles in total a dataset
    // of size n[1] (the number of items in a matrix row is the number of
    // columns.
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.

    // this really ends up within the mm field. It's not a complete file
    // name though. We lack the implementation extension, the possible
    // transposition tag, as well as the .bin extension.
    std::string locfile;

    /* These two are global to all threads and jobs (well, each thread
     * has its own pointer though, it's not a shared_malloc. It could be,
     * but it isn't).
     *
     * For random matrices, both strings below are empty.
     */
    std::string mname;
    std::string bname;

    balancing bal;

    /* These are local excerpts of the balancing permutation: arrays of
     * pairs (row index in the (sub)-matrix) ==> (row index in the
     * original matrix), but only when both coordinates of this pair are
     * in the current row and column range. This can be viewed as the set
     * of non-zero positions in the permutation matrix if it were split
     * just like the current matrix is. */
    std::vector<std::array<unsigned int, 2>> perm[2]; // row,col

    /* XXX XXX XXX dangerous: has_perm_local is only a local thing!
     * Some threads/jobs may have has_perm_local(d)==false while others don't!
     */
    inline bool has_perm_local(int d) const { return !perm[d].empty(); }

    /* This is the global version, as prepared by
     * matmul_top_init_prepare_local_permutations
     */
    int has_perm_map[2];
    inline bool has_perm(int d) const { return has_perm_map[d]; }


    matmul_ptr mm;
};

struct matmul_top_data {
    arith_generic * abase;
    parallelizing_info_ptr pi;
    pi_datatype_ptr pitype;
    /* These n[] and n0[] correspond to the dimensions of the product
     *
     * n[0] is matrices[0]->n[0]
     * n[1] is matrices[nmatrices-1]->n[1]
     * n0[0] is matrices[0]->n0[0]
     * n0[1] is matrices[nmatrices-1]->n0[1]
     */
    unsigned int n[2];
    unsigned int n0[2]; // n0: unpadded.

    /* The matrix we are dealing with is
     * matrices[0] * matrices[1] * ... * matrices[nmatrices-1]
     */
    std::vector<matmul_top_matrix> matrices;
    matmul_top_data(
        arith_generic * abase,
        parallelizing_info_ptr pi,
        param_list_ptr pl,
        int optimized_direction);
    ~matmul_top_data();
};

struct mmt_vec;


extern void matmul_top_decl_usage(param_list_ptr pl);
extern void matmul_top_lookup_parameters(param_list_ptr pl);
extern void matmul_top_report(matmul_top_data & mmt, double scale, int full);
extern unsigned int matmul_top_rank_upper_bound(matmul_top_data & mmt);
#if 0
extern void matmul_top_fill_random_source(matmul_top_data & mmt, int d);
#endif
extern void mmt_vec_truncate(matmul_top_data & mmt, mmt_vec & v);
extern void mmt_vec_truncate_above_index(matmul_top_data & mmt, mmt_vec & v, unsigned int idx);
extern void mmt_vec_truncate_below_index(matmul_top_data & mmt, mmt_vec & v, unsigned int idx);
extern void matmul_top_mul_cpu(matmul_top_data & mmt, int midx, int d, mmt_vec & w, mmt_vec const & v);
extern void matmul_top_comm_bench(matmul_top_data & mmt, int d);
extern void matmul_top_mul_comm(mmt_vec & w, mmt_vec & v);

/* v is both input and output. w is temporary */
extern void matmul_top_mul(matmul_top_data & mmt, mmt_vec *w, struct timing_data * tt);

extern void mmt_vec_twist(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_vec_untwist(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_vec_apply_T(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_vec_unapply_T(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_vec_apply_S(matmul_top_data & mmt, int midx, mmt_vec & y);
extern void mmt_vec_unapply_S(matmul_top_data & mmt, int midx, mmt_vec & y);
extern void mmt_vec_apply_P(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_vec_unapply_P(matmul_top_data & mmt, mmt_vec & y);
extern void mmt_apply_identity(mmt_vec & w, mmt_vec const & v);
extern void indices_twist(matmul_top_data & mmt, uint32_t * xs, unsigned int n, int d);

#endif	/* MATMUL_TOP_HPP_ */
