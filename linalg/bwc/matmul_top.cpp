#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>               // for access, unlink, ssize_t, R_OK
#include <sys/stat.h>             // for stat, mkdir
#include <pthread.h>              // for pthread_mutex_lock, pthread_mutex_u...
#include <gmp.h>
#include "async.hpp"              // for timing_next_timer, timing_data (ptr...
#include "balancing_workhorse.hpp"
#include "intersections.h"
#include "bwc_config.h" // IWYU pragma: keep
#include "macros.h"     // ASSERT_ALWAYS // IWYU pragma: keep
#include "matmul.hpp"
#include "matmul_top.hpp"
#include "mf_bal.hpp"
#include "misc.h"
#include "params.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "random_matrix.hpp"
#include "raw_matrix_u32.h"       // for matrix_u32
#include "select_mpi.h"
#include "timing.h"     // wct_seconds
#include "verbose.h"    // CADO_VERBOSE_PRINT_BWC_CACHE_BUILD
#include "matmul_top_comm.hpp"

///////////////////////////////////////////////////////////////////
/* Start with stuff that does not depend on abase at all -- this
 * provides a half-baked interface */

void matmul_top_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "matrix",
            "the matrix file (binary)");
    param_list_decl_usage(pl, "balancing",
            "the matrix balancing file, as computed by mf_bal");
    param_list_decl_usage(pl, "random_matrix",
            "characteristics of a random matrix to be used for staged runs.");
    param_list_decl_usage(pl, "static_random_matrix",
            "(unset or set to something arbitrary): indicate that the matrix is fake, and that there is no need to bother with the generation of vectors");

    param_list_decl_usage(pl, "rebuild_cache",
            "force rebuilding matrix caches");
    param_list_decl_usage(pl, "export_cachelist",
            "print the per-node needed cache files and exit");
    param_list_decl_usage(pl, "save_submatrices",
            "after dispatching, save a copy of the local uncompressed matrix before creating the cache file");
    param_list_decl_usage(pl, "sequential_cache_build",
            "build the cache files sequentially on each node");
    param_list_decl_usage(pl, "sequential_cache_read",
            "read the cache files sequentially on each node");
    param_list_decl_usage(pl, "balancing_options",
            "options to pass to the balancing subprogram (see mf_bal_adjust_from_option_string)");
    param_list_decl_usage(pl, "multi_matrix",
            "whether to chain several matrices (experimental)");
    balancing_decl_usage(pl);
    matmul_decl_usage(pl);
}

void matmul_top_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "matrix");
    param_list_lookup_string(pl, "balancing");
    param_list_lookup_string(pl, "random_matrix");
    param_list_lookup_string(pl, "static_random_matrix");
    param_list_lookup_string(pl, "rebuild_cache");
    param_list_lookup_string(pl, "export_cachelist");
    param_list_lookup_string(pl, "save_submatrices");
    param_list_lookup_string(pl, "sequential_cache_build");
    param_list_lookup_string(pl, "sequential_cache_read");
    param_list_lookup_string(pl, "balancing_options");
    param_list_lookup_string(pl, "multi_matrix");
    balancing_lookup_parameters(pl);
    matmul_lookup_parameters(pl);
}

//////////////////////////////////////////////////////////////////////////

void matmul_top_mul(matmul_top_data & mmt, mmt_vec * v, struct timing_data * tt)/*{{{*/
{
    /* Do all matrices in turn.
     *
     * We represent M as * M0*M1*..*M_{n-1}.
     *
     * The input vector is v[0], and the result is put in v[0] again.
     *
     * The direction in which we apply the product is given by v[0].d.
     * For v[0].d == 0, we do v[0]*M. For v[0].d==1, we do M*v[0]
     *
     * We use temporaries as follows.
     * For v[0].d == 0:
     *  v[1] <- v[0] * M0 ; v[1].d == 1
     *  v[2] <- v[1] * M1 ; v[2].d == 0
     *  if n is odd:
     *  v[n] <- v[n-1] * M_{n-1} ; v[n].d == 1
     *  if n is even:
     *  v[0] <- v[n-1] * M_{n-1} ; v[0].d == 0
     *
     * For v[0].d == 1:
     *  v[1] <- M0 * v[0] ; v[1].d == 0
     *  v[2] <- M1 * v[1] ; v[2].d == 1
     *  if n is odd:
     *  v[n] <- M_{n-1} * v[n-1] ; v[n].d == 0
     *  if n is even:
     *  v[0] <- M_{n-1} * v[n-1] ; v[0].d == 1
     *
     * This has the consequence that v must hold exactly n+(n&1) vectors.
     *
     * Appropriate communication operations are run after each step.
     *
     * If tt is not NULL, it should be a timing_data structure holding
     * exactly 4*n timers (or only 4, conceivably). Timers are switched
     * exactly that many times.
     *
     * If the mmt.pi->interleaving setting is on, we interleave
     * computations and communications. We do 2n flips. Communications
     * are forbidden both before and after calls to this function (in the
     * directly adjacent code fragments before the closest flip() call,
     * that is).
     */

    int d = v[0].d;
    int nmats_odd = mmt.matrices.size() & 1;
    int midx = (d ? (mmt.matrices.size() - 1) : 0);
    for(size_t l = 0 ; l < mmt.matrices.size() ; l++) {
        mmt_vec const & src = v[l];
        bool last = l == (mmt.matrices.size() - 1);
        size_t lnext = last && !nmats_odd ? 0 : (l+1);
        mmt_vec & dst = v[lnext];

        ASSERT_ALWAYS(src.consistency == 2);
        matmul_top_mul_cpu(mmt, midx, d, dst, src);
        ASSERT_ALWAYS(dst.consistency == 0);

        timing_next_timer(tt);
        /* now measuring jitter */
        pi_interleaving_flip(mmt.pi);
        serialize(mmt.pi->m);
        timing_next_timer(tt);

        /* Now we can resume MPI communications. */
        if (last && nmats_odd) {
            ASSERT_ALWAYS(lnext == mmt.matrices.size());
            matmul_top_mul_comm(v[0], dst);
        } else {
            mmt_vec_allreduce(dst);
        }
        timing_next_timer(tt);
        /* now measuring jitter */
        pi_interleaving_flip(mmt.pi);
        serialize(mmt.pi->m);

        timing_next_timer(tt);
        midx += d ? -1 : 1;
    }
    ASSERT_ALWAYS(v[0].consistency == 2);
}
/*}}}*/

#if 0/*{{{*/
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source_generic(matmul_top_data & mmt, size_t stride, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt.wr[d]->v;

    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do mmt_vec_broadcast
    // afterwards...
    if ((v.flags & THREAD_SHARED_VECTOR) == 0 || mmt.pi->wr[d]->trank == 0)
        mpfq_generic_random(stride, v.v, mmt.wr[d]->i1 - mmt.wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    mmt_vec_broadcast(mmt, v, d);
}
#endif/*}}}*/

/**********************************************************************/
/* Utility stuff for applying standard permutations to vectors. */

/* We have four permutations defined in the context of bwc.
 *
 * Sr -- row permutation
 * Sc -- column permutation
 * P  -- shuffled product permutation
 * T  -- decorrelating permutation
 *
 * As a convention, in the code as well as in the documentation below, we
 * freely handle the following equivalent representations of a
 * permutation:
 *  - a function f,
 *  - a list of images [f(i), * 0<=i<n], 
 *  - a list of pairs [<i,f(i)>, * 0<=i<n], 
 *  - a matrix F of size n*n, with only entries at row i and column f(i) being 1
 *
 * We let Mpad be the "padded" matrix, where #rows and #columns are both
 * multiples of nh*nv
 *
 * T is a constant permutation on the columns of the matrix M. It is
 * mostly working as a preconditioner, meant to eliminate some nasty
 * correlation effects between row and column weights. In every respect,
 * the matrix we work with consistently is the matrix Mpad*T. T is computed
 * by balancing_pre_shuffle and balancing_pre_unshuffle.
 *
 * row i of T has its non-zero coefficient in column
 * balancing_pre_shuffle(i)
 *
 * Sr and Sc are defined in the balancing file (as list of images). They
 * are such that the matrix Mtwisted = Sr*Mpad*T*Sc^-1 is well balanced
 * across the different jobs and threads. For square matrices, depending
 * on how the balancing was computed, one may be implicitly defined from
 * the other (but not equal, for the reason below).
 *
 * P is defined only for square matrices. It reflects the fact that
 * although the matrix is "geographically" stored as Mtwisted =
 * Sr*Mpad*T*Sc^-1 in the jobs and threads, the matmul_top_mul code
 * multiplies by a matrix which is:
 *      - for matrix times vector: Sc*Mpad*T*Sc^-1
 *              (that is, v<-v*Transpose(Sc*Mpad*T*Sc^-1) )
 *      - for vector times matrix: P*Sc*Mpad*T*Sc^-1*P^-1
 *
 * Implicit definition of Sr for square matrices
 * =============================================
 *
 * Sc, read as the colperm[] array in the balancing file, is such that column
 * i of the twisted matrix is in fact column Sc(i) in the original matrix,
 * which means that we build Mtwisted as [...]*M*Sc^-1
 *
 * for square matrices, with FLAG_REPLICATE, the Sr permutation is chosen
 * implicitly.
 *
 * here's how we compute the (forward) row permutation in the code (here, xr
 * is the colperm array).
 *              ix = (i * nv + j) * elem;
 *              iy = (j * nh + i) * elem;
 *              for(unsigned int k = 0 ; k < elem ; k++) {
 *                  m->fw_rowperm[xr[iy+k]] = ix+k;
 *
 * We denote by Sr the permutation, implicitly computed here, such that row i
 * of the twisted matrix is row Sr(i) in the original matrix -- so that
 * Mtwisted is in fact Sr*M*Sc^-1. Here, fw_rowperm is the inverse: row i in
 * the original matrix goes to row fw_rowperm[i] in the twisted matrix. So
 * that fw_rowperm == Sr^-1.
 *
 * Our code does fw_rowperm(Sc(P(x)))=x, for P the permutation which sends
 * sub-block nv*i+j to sub-block nh*j+i. Writing as operations on row vectors
 * (as magma does), this gives:
 *      P * Sc * (Sr^-1) = id
 *      Sr = P * Sc
 *
 * So in this case the twisted matrix is in fact Sr*Mpad*Sc^-1, and Sr = P*Sc;
 *
 * Implicit definition of Sc for square matrices
 * =============================================
 *
 * There is no code doing this at the moment, but one could imagine
 * defining this as well. In such a situation, we write down what Sc
 * would be.
 *
 * Matrix is geographically stored as Sr*Mpad*T*Sc^-1 in the jobs and
 * threads, and the matmul_top_mul code multiplies by a matrix which is:
 *      - for matrix times vector: P^-1*Sr*Mpad*T*Sc^-1
 *      - for vector times matrix: Sr*Mpad*T*Sc^-1*P^-1
 * Therefore we want Sc^-1*P^-1 == Sr^-1, whence Sc == P^-1 * Sr
 *
 * Action of P with matmul_top_mul_comm
 * ====================================
 *
 * After v*Mtwisted, reduce_sameside of the result (in direction 1) followed
 * by broadcast (in direction 0) produces a twisted vector. Since v*Mtwisted
 * is in direction 1, then blocks of the vector are to be read in column-major
 * order. If we reduce then broadcast, then sub-block nh*j+i will go to
 * sub-block nv*i+j. This means that v*Mtwisted will be transformed into
 * v*Mtwisted*P^-1
 *
 * After v*Transpose(Mtwisted), reduce_sameside on the result which is in
 * direction 0 transforms it to v*transpose(Mtwisted)*P (here P is
 * transpose of P^-1: we transpose indices blocks in the other
 * direction), which is thus v*transpose(P^-1*Mtwisted).
 *
 * Conclusion: with Mtwisted = Sr*Mpad*Sc^-1, and Sr = P*Sc, when we do
 * matmul_top_mul_cpu followed by matmul_top_mul_comm, we multiply by:
 *      - for matrix times vector: Sc*Mpad*Sc^-1
 *              (that is, v<-v*Transpose(Sc*Mpad*T*Sc^-1) )
 *      - for vector times matrix: P*Sc*Mpad*Sc^-1*P^-1
 *
 * Applying P to vectors
 * =====================
 *
 * We can emulate the multiplications by P and P^-1 with appropriate
 * combinations of apply_identity and matmul_top_mul_comm. So we give a few
 * details.
 *
 * Let v be a distributed vector. With mmt_apply_identity, the output is
 * distributed in the other direction, but not consistent. If we do
 * mmt_vec_reduce_sameside (or mmt_vec_allreduce, which does more), the
 * resulting vector (looking at locally-owned pieces only) is the same as
 * v, in the other direction.
 *
 * And then, matmul_top_mul_comm on this resulting vector does the P
 * action as above. Therefore:
 *
 * For v.d == 0, applying in sequence the functions mmt_apply_identity
 * matmul_top_mul_comm, then we get v*P^-1
 *
 * For v.d == 1, the same sequence produces v*Transpose(P^-1)==v*P
 *
 * Doing the converse is feasible. For v.d==0, if we do instead
 * matmul_top_mul_comm then mmt_apply_identity, we get v*P
 *
 * For v.d==1, if we do matmul_top_mul_comm then mmt_apply_identity, we
 * get v*P^-1
 *
 */

void mmt_apply_identity(mmt_vec & w, mmt_vec const & v)
{
    /* input: fully consistent */
    /* output: inconsistent ! 
     * Need mmt_vec_allreduce or mmt_vec_reduce_sameside, or
     * matmul_top_mul_comm, depending on what we want to do. */
    ASSERT_ALWAYS(v.consistency == 2);
    ASSERT_ALWAYS(w.abase == v.abase);
    ASSERT_ALWAYS(v.d != w.d);
    ASSERT_ALWAYS(v.n == w.n);

    arith_generic * A = v.abase;

    serialize_threads(w.pi->m);
    mmt_full_vec_set_zero(w);
    serialize_threads(w.pi->m);

    unsigned int v_off, w_off;
    unsigned int how_many = intersect_two_intervals(&v_off, &w_off,
            v.i0, v.i1, w.i0, w.i1);

    A->vec_set(A->vec_subvec(w.v, w_off), A->vec_subvec(v.v, v_off), how_many);
    w.consistency = 1;
}

void mmt_vec_apply_or_unapply_P_inner(matmul_top_data & mmt, mmt_vec & y, int apply)
{
    ASSERT_ALWAYS(y.consistency == 2);
    mmt_vec yt(mmt, y.abase, y.pitype, !y.d, 0, y.n);
    if ((apply ^ y.d) == 0) {
        // y.d == 0: get v*P^-1
        // y.d == 1: get v*P
        mmt_apply_identity(yt, y);
        matmul_top_mul_comm(y, yt);
    } else {
        // y.d == 0: get v*P
        // y.d == 1: get v*P^-1
        mmt_vec_downgrade_consistency(y);
        matmul_top_mul_comm(yt, y);
        mmt_apply_identity(y, yt);
        mmt_vec_allreduce(y);
    }
    ASSERT_ALWAYS(y.consistency == 2);
}

void mmt_vec_unapply_P(matmul_top_data & mmt, mmt_vec & y)
{
    mmt_vec_apply_or_unapply_P_inner(mmt, y, 0);
}

void mmt_vec_apply_P(matmul_top_data & mmt, mmt_vec & y)
{
    mmt_vec_apply_or_unapply_P_inner(mmt, y, 1);
}

/* apply == 1 for apply, apply == 0 for unapply *
 * apply == 1 d == 0 Sr defined:   v <- v * Sr
 * apply == 1 d == 0 Sr implicit:  v <- v * Sc
 * apply == 1 d == 1 Sc defined:   v <- v * Sc
 * apply == 1 d == 1 Sc implicit:  v <- v * Sr
 * apply == 0 d == 0 Sr defined:   v <- v * Sr^-1
 * apply == 0 d == 0 Sr implicit:  v <- v * Sc^-1
 * apply == 0 d == 1 Sc defined:   v <- v * Sc^-1
 * apply == 0 d == 1 Sc implicit:  v <- v * Sr^-1
 *
 * See that when, say, Sr is implicitly defined (to P*Sc), this function
 * only applies Sc, not P !
 */
void mmt_vec_apply_or_unapply_S_inner(matmul_top_data & mmt, int midx, mmt_vec & y, int apply)
{
    ASSERT_ALWAYS(y.consistency == 2);
    /* input: fully consistent */
    /* output: fully consistent */
    int d = y.d;
    arith_generic * A = y.abase;

    serialize_threads(y.pi->m);

    /* We'll have two vectors of size n[d], one named y in direction d,
     * and one named yt in direction !d.
     * 
     *
     * In the permutation mmt.perm[d], the pairs (i,j) are such that,
     * given two vectors of size n[d], one named y in direction d,
     * and one named yt in direction !d, we have:
     *  i in [y.i0..y.i1[
     *  j in [yt->i0..yt->i1[
     */
    matmul_top_matrix const & Mloc = mmt.matrices[midx];

    /* For square matrices, we'll use the other permutation transparently
     * with this piece of code. Note though that when we do so, applying
     * the permutation actually goes in the opposite direction. */
    int xd = d;
    if (!Mloc.has_perm(xd) && (Mloc.bal.flags & FLAG_REPLICATE)) {
        ASSERT_ALWAYS(Mloc.n[0] == Mloc.n[1]);
        xd = !d;
    }

    /* it could well be that we have nothing to do */
    if (!Mloc.has_perm(xd)) return;

    auto const & s(Mloc.perm[xd]);

    if ((apply^d^xd) == 0) {
        /*
         * apply == 0 d == 0 Sr defined:   v <- v * Sr^-1
         * apply == 0 d == 1 Sc defined:   v <- v * Sc^-1
         * apply == 1 d == 0 Sr implicit:  v <- v * Sc
         * apply == 1 d == 1 Sc implicit:  v <- v * Sr
         */
        mmt_vec yt(mmt, A, y.pitype, !d, 0, y.n);

        mmt_apply_identity(yt, y);
        mmt_vec_allreduce(yt);
        mmt_full_vec_set_zero(y);
        serialize_threads(y.pi->m);
        for(auto const & uv : s) {
            if (uv[d^xd] < y.i0 || uv[d^xd] >= y.i1)
                continue;
            if (uv[d^xd^1] < yt.i0 || uv[d^xd^1] >= yt.i1)
                continue;
            A->set( A->vec_item( y.v, uv[d^xd]  - y.i0),
                    A->vec_item(yt.v, uv[d^xd^1] - yt.i0));
        }
        y.consistency = 1;
        serialize_threads(y.pi->m);
        mmt_vec_allreduce(y);
    } else {
        /*
         * apply == 1 d == 0 Sr defined:   v <- v * Sr
         * apply == 1 d == 1 Sc defined:   v <- v * Sc
         * apply == 0 d == 0 Sr implicit:  v <- v * Sc^-1
         * apply == 0 d == 1 Sc implicit:  v <- v * Sr^-1
         */
        mmt_vec yt(mmt, A, y.pitype, !d, 0, y.n);
        for(auto const & uv : s) {
            if (uv[d^xd] < y.i0 || uv[d^xd] >= y.i1)
                continue;
            if (uv[d^xd^1] < yt.i0 || uv[d^xd^1] >= yt.i1)
                continue;
            A->set( A->vec_item(yt.v, uv[d^xd^1] - yt.i0),
                    A->vec_item( y.v, uv[d^xd] -  y.i0));
        }
        yt.consistency = 1;
        mmt_vec_allreduce(yt);
        mmt_apply_identity(y, yt);
        mmt_vec_allreduce(y);
    }
    serialize_threads(y.pi->m);
    ASSERT_ALWAYS(y.consistency == 2);
}

/* multiply v by Sr^-1 if v.d == 0, by Sc^-1 if v.d == 1 */
void mmt_vec_unapply_S(matmul_top_data & mmt, int midx, mmt_vec & y)
{
    matmul_top_matrix & Mloc = mmt.matrices[midx];
    mmt_vec_apply_or_unapply_S_inner(mmt, midx, y, 0);
    if ((Mloc.bal.flags & FLAG_REPLICATE) && !Mloc.has_perm(y.d)) {
        if (y.d == 0) {
            /* implicit Sr^-1 is Sc^-1*P^-1 */
            mmt_vec_unapply_P(mmt, y);
        } else {
            /* implicit Sc^-1 is Sr^-1*P */
            mmt_vec_apply_P(mmt, y);
        }
    }
}

/* multiply v by Sr if v.d == 0, by Sc if v.d == 1 */
void mmt_vec_apply_S(matmul_top_data & mmt, int midx, mmt_vec & y)
{
    matmul_top_matrix & Mloc = mmt.matrices[midx];
    if ((Mloc.bal.flags & FLAG_REPLICATE) && !Mloc.has_perm(y.d)) {
        if (y.d == 0) {
            /* implicit Sr is P * Sc */
            mmt_vec_apply_P(mmt, y);
        } else {
            /* implicit Sc is P^-1 * Sr */
            mmt_vec_unapply_P(mmt, y);
        }
    }
    mmt_vec_apply_or_unapply_S_inner(mmt, midx, y, 1);
}

/* for square matrix products, the inner loops do
 *   v <- v * (Sr=P*Sc)*Mpad*Sc^-1*P^-1              (for vector times matrix)
 *   v <- v * Transpose(P^-1*(Sr=P*Sc)*Mpad*Sc^-1)   (for matrix times vector)
 * while for non-square (no FLAG_REPLICATE), it's simply by
 *      Sr*Mpad*Sc^-1
 *
 * therefore the rules for twisting and untwisting are as follows.
 *
 * twisting v.d == 0
 *      we assume we want here to change v to become a good *input* for a
 *      vector times matrix operation. Therefore we do:
 *              v <- v * Sr^-1
 * twisting v.d == 1
 *      v <- v * Sc^-1
 *
 */

void mmt_vec_twist(matmul_top_data & mmt, mmt_vec & y)
{
    mmt_vec_unapply_S(mmt, y.d == 0 ? 0 : (mmt.matrices.size()-1), y);
}

void mmt_vec_untwist(matmul_top_data & mmt, mmt_vec & y)
{
    mmt_vec_apply_S(mmt, y.d == 0 ? 0 : (mmt.matrices.size()-1), y);
}

/* {{{ mmt_vec_{un,}appy_T -- this applies the fixed column
 * permutation which we use unconditionally in bwc to avoid correlation
 * of row and column weights.
 */
    // pshuf indicates two integers a,b such that the COLUMN i of the input
    // matrix is in fact mapped to column a*i+b mod n in the matrix we work
    // with. pshuf_inv indicates the inverse permutation. a and b do
void mmt_vec_apply_or_unapply_T_inner(matmul_top_data & mmt, mmt_vec & y, int apply)
{
    /* apply: coefficient i of the vector goes to coefficient
     * balancing_pre_shuffle[i]
     *
     * For row vectors this means: apply == (v <- v * T)
     * For column vectors this means: apply == (v <- T^-1 * v)
     *
     */
    if (y.d == 0) return;
    matmul_top_matrix & Mloc = mmt.matrices[mmt.matrices.size() - 1];
    ASSERT_ALWAYS(y.consistency == 2);
    serialize_threads(y.pi->m);
    mmt_vec yt(mmt, y.abase, y.pitype, !y.d, 0, y.n);
    for(unsigned int i = y.i0 ; i < y.i1 ; i++) {
        unsigned int j;
        if (apply) {
            j = balancing_pre_shuffle(Mloc.bal, i);
        } else {
            j = balancing_pre_unshuffle(Mloc.bal, i);
        }
        if (j >= yt.i0 && j < yt.i1) {
            y.abase->set(
                    y.abase->vec_item(yt.v, j - yt.i0),
                    y.abase->vec_item(y.v, i - y.i0));
        }
    }
    yt.consistency = 0;
    mmt_vec_allreduce(yt);
    mmt_apply_identity(y, yt);
    mmt_vec_allreduce(y);
}
void mmt_vec_unapply_T(matmul_top_data & mmt, mmt_vec & v)
{
    mmt_vec_apply_or_unapply_T_inner(mmt, v, 0);
}

void mmt_vec_apply_T(matmul_top_data & mmt, mmt_vec & v)
{
    mmt_vec_apply_or_unapply_T_inner(mmt, v, 1);
}
/* }}} */

/* {{{ Application of permutations to indices */

/* Given indices which relate to a vector in direction d, modify them in
 * place so that they correspond to coefficients of the matching twisted vector.
 *
 * So we want: v_i == twist(v)_{twist(i)}
 *
 * for d == 0, twist(v) = v * Sr^-1, so that index i in v goes to
 * position S(i) in in the twisted vector. For implicit Sr, we have to
 * take into account the fact that Sr=P*Sc, hence S(i) is in fact
 * Sc[P[i]].
 *
 * The algorithm is simple: only one thread knows about each single
 * (i,Sr(i)) pair. So this one changes the relevant xs[] values, leaving
 * the other to zero. And then we do a global allreduce() on the xs[]
 * values.
 *
 * (here, "relevant" may mean something which includes the P
 * permutation).
 */

/* returns the two intervals such that for all pairs (i,j) in
 * mmt.perm->x, we have ii[0] <= i < ii[1], and jj[0] <= j < jj[1]
 */
static void get_local_permutations_ranges(matmul_top_data & mmt, int d, unsigned int ii[2], unsigned int jj[2])
{
    int pos[2];

    for(int dir = 0 ; dir < 2 ; dir++)  {
        pi_comm_ptr piwr = mmt.pi->wr[dir];
        pos[dir] = piwr->jrank * piwr->ncores + piwr->trank;
    }

    size_t e = mmt.n[d] / mmt.pi->m->totalsize;
    ii[0] = e *  pos[!d]    * mmt.pi->wr[d]->totalsize;
    ii[1] = e * (pos[!d]+1) * mmt.pi->wr[d]->totalsize;
    jj[0] = e *  pos[d]     * mmt.pi->wr[!d]->totalsize;
    jj[1] = e * (pos[d]+1)  * mmt.pi->wr[!d]->totalsize;
}

void indices_twist(matmul_top_data & mmt, uint32_t * xs, unsigned int n, int d)
{
    int midx = d == 0 ? 0 : (mmt.matrices.size() - 1);
    matmul_top_matrix & Mloc = mmt.matrices[midx];
    /* d == 1: twist(v) = v*Sc^-1
     * d == 0: twist(v) = v*Sr^-1
     */
    unsigned int ii[2], jj[2];

    if (Mloc.has_perm(d)) {
        /* explicit S */
        /* coordinate S[i] in the original vector becomes i in the
         * twisted vector.
         */
        get_local_permutations_ranges(mmt, d, ii, jj);
        unsigned int * r = (unsigned int *) malloc((jj[1] - jj[0]) * sizeof(unsigned int));
        memset(r, 0, (jj[1] - jj[0]) * sizeof(unsigned int));
        for(auto const & uv : Mloc.perm[d]) {
            ASSERT_ALWAYS(uv[0] >= ii[0]);
            ASSERT_ALWAYS(uv[0] <  ii[1]);
            ASSERT_ALWAYS(uv[1] >= jj[0]);
            ASSERT_ALWAYS(uv[1] <  jj[1]);
            // r[uv[0] - ii[0]] = uv[1];
            r[uv[1] - jj[0]] = uv[0];
        }
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = xs[k];
            if (j >= jj[0] && j < jj[1])
                xs[k] = r[j - jj[0]];
            else
                xs[k] = 0;
        }
        free(r);
        pi_allreduce(NULL, xs, n * sizeof(uint32_t), BWC_PI_BYTE, BWC_PI_BXOR, mmt.pi->m);
    } else if (Mloc.has_perm(!d) && (Mloc.bal.flags & FLAG_REPLICATE)) {
        ASSERT_ALWAYS(Mloc.n[0] == Mloc.n[1]);
        /* implicit S -- first we get the bits about the S in the other
         * direction, because the pieces we have are for the other
         * ranges, which is a bit disturbing...
         */
        get_local_permutations_ranges(mmt, !d, ii, jj);
        unsigned int * r = (unsigned int *) malloc((jj[1] - jj[0]) * sizeof(unsigned int));
        memset(r, 0, (jj[1] - jj[0]) * sizeof(unsigned int));
        for(auto const & uv : Mloc.perm[!d]) {
            ASSERT_ALWAYS(uv[0] >= ii[0]);
            ASSERT_ALWAYS(uv[0] <  ii[1]);
            ASSERT_ALWAYS(uv[1] >= jj[0]);
            ASSERT_ALWAYS(uv[1] <  jj[1]);
            r[uv[1] - jj[0]] = uv[0];
        }
        /* nh and nv are the same for all submatrices, really */
        unsigned int nn[2] = { Mloc.bal.nh, Mloc.bal.nv };
        /* for d == 0, we have implicit Sr = P * Sc.
         * for d == 1, we have implicit Sc = P^-1 * Sc
         */
        size_t z = Mloc.bal.trows / (nn[0]*nn[1]);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = xs[k];
            /* for d == 0, index i goes to P^-1[Sc^-1[i]] */
            if (j >= jj[0] && j < jj[1]) {
                unsigned int i = r[j - jj[0]];
                unsigned int qz = i / z;
                unsigned int rz = i % z;
                /* P    sends sub-block nv*i+j to sub-block nh*j+i */
                /* P^-1 sends sub-block nh*j+i to sub-block nv*i+j */
                unsigned int qh = qz / nn[d];
                unsigned int rh = qz % nn[d];
                ASSERT_ALWAYS(qz == qh * nn[d] + rh);
                ASSERT_ALWAYS(i == (qh * nn[d] + rh)*z + rz);
                xs[k] = (rh * nn[!d] + qh)*z + rz;
            } else {
                xs[k] = 0;
            }
        }
        free(r);
        pi_allreduce(NULL, xs, n * sizeof(uint32_t), BWC_PI_BYTE, BWC_PI_BXOR, mmt.pi->m);
    }
}
/* }}} */

/**********************************************************************/
#if 0
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source(matmul_top_data & mmt, int d)
{
    matmul_top_fill_random_source_generic(mmt, mmt.vr->stride, NULL, d);
}
#endif

/* For d == 0: do w = v * M
 * For d == 1: do w = M * v
 *
 * We do not necessarily have v.d == d, although this will admittedly be
 * the case most often:
 * - in block Wiedemann, we have a vector split in some direction (say
 *   d==1 when we wanna solve Mw=0), we compute w=Mv, and then there's
 *   the matmul_top_mul_comm step which moves stuff to v again.
 * - in block Lanczos, say we start from v.d == 1 again. We do w=M*v,
 *   so that w.d==0. But then we want to compute M^T * w, which is w^T *
 *   M. So again, w.d == 0 is appropriate with the second product being
 *   in the direction foo*M.
 * - the only case where this does not necessarily happen so is when we
 *   have several matrices.
 */
void matmul_top_mul_cpu(matmul_top_data & mmt, int midx, int d, mmt_vec & w, mmt_vec const & v)
{
    matmul_top_matrix & Mloc = mmt.matrices[midx];
    ASSERT_ALWAYS(v.consistency == 2);
    ASSERT_ALWAYS(w.abase == v.abase);
    unsigned int di_in  = v.i1 - v.i0;
    unsigned int di_out = w.i1 - w.i0;
    ASSERT_ALWAYS(Mloc.mm->dim[!d] == di_out);
    ASSERT_ALWAYS(Mloc.mm->dim[d] == di_in);

    ASSERT_ALWAYS(w.siblings); /* w must not be shared */

    pi_log_op(mmt.pi->m, "[%s:%d] enter matmul_mul", __func__, __LINE__);

    /* Note that matmul_init copies the calling abase argument to the
     * lower-level mm structure. It can quite probably be qualified as a
     * flaw.
     */
    matmul_mul(Mloc.mm, w.v, v.v, d);
    w.consistency = 0;
}

/* This takes partial results in w, and puts the
 * collected and re-broadcasted results in the areas mmt.wd[d]->v
 *
 * Note that for the shuffled product, this is not equivalent to a trivial
 * operation.
 */
void matmul_top_mul_comm(mmt_vec & v, mmt_vec & w)
{
    /* this takes inconsistent input.
     * XXX if we have fully consistent input, then a reduce() is much
     * undesired !
     */
    ASSERT_ALWAYS(w.consistency != 2);
    pi_log_op(v.pi->m, "[%s:%d] enter mmt_vec_reduce", __func__, __LINE__);
    mmt_vec_reduce(v, w);
    ASSERT_ALWAYS(v.consistency == 1);
    pi_log_op(v.pi->m, "[%s:%d] enter mmt_vec_broadcast", __func__, __LINE__);
    mmt_vec_broadcast(v);
    ASSERT_ALWAYS(v.consistency == 2);

    /* If we have shared input data for the column threads, then we'd
     * better make sure it has arrived completely, because while all
     * threads will need the data, only one is actually importing it.
     */
    if (!v.siblings) {
        pi_log_op(v.pi->wr[v.d], "[%s:%d] serialize threads", __func__, __LINE__);
        serialize_threads(v.pi->wr[v.d]);
    }
}


/* _above and _below functions here do not use mmt, but we activate this
 * same interface nevertheless, for consistency with mmt_vec_truncate */
void mmt_vec_truncate_above_index(matmul_top_data & mmt MAYBE_UNUSED, mmt_vec & v, unsigned int idx)
{
    if (idx <= v.i0) idx = v.i0;
    if (v.i0 <= idx && idx < v.i1) {
        if (v.siblings) {
            v.abase->vec_set_zero(
                    v.abase->vec_subvec(v.v, idx - v.i0),
                    v.i1 - idx);
        } else {
            serialize_threads(v.pi->wr[v.d]);
            if (v.pi->wr[v.d]->trank == 0)
                v.abase->vec_set_zero(
                        v.abase->vec_subvec(v.v, idx - v.i0),
                        v.i1 - idx);
            serialize_threads(v.pi->wr[v.d]);
        }
    }
}

void mmt_vec_truncate_below_index(matmul_top_data & mmt MAYBE_UNUSED, mmt_vec & v, unsigned int idx)
{
    if (idx >= v.i1) idx = v.i1;
    if (v.i0 <= idx && idx < v.i1) {
        if (v.siblings) {
            v.abase->vec_set_zero(v.v, idx - v.i0);
        } else {
            serialize_threads(v.pi->wr[v.d]);
            if (v.pi->wr[v.d]->trank == 0)
                v.abase->vec_set_zero(v.v, idx - v.i0);
            serialize_threads(v.pi->wr[v.d]);
        }
    }
}

void mmt_vec_truncate(matmul_top_data & mmt, mmt_vec & v)
{
    mmt_vec_truncate_above_index(mmt, v, mmt.n0[v.d]);
}

/**********************************************************************/
static void matmul_top_read_submatrix(matmul_top_data & mmt, int midx, param_list_ptr pl, int optimized_direction);

/* see matmul_top2.cpp */
extern char * matrix_list_get_item(param_list_ptr pl, const char * key, int midx);

static std::string matrix_get_derived_cache_subdir(std::string const & matrixname, parallelizing_info_ptr pi)
{
    /* input is NULL in the case of random matrices */
    if (matrixname.empty()) return std::string();
    unsigned int nh = pi->wr[1]->totalsize;
    unsigned int nv = pi->wr[0]->totalsize;
    char * copy = strdup(matrixname.c_str());
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int rc = asprintf(&t, "%s.%ux%u", t, nh, nv);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    std::string s(t);
    free(t);
    return s;
}

static void matrix_create_derived_cache_subdir(std::string const & matrixname, parallelizing_info_ptr pi)
{
    std::string d = matrix_get_derived_cache_subdir(matrixname, pi);
    struct stat sbuf[1];
    int rc = stat(d.c_str(), sbuf);
    if (rc < 0 && errno == ENOENT) {
        rc = mkdir(d.c_str(), 0777);
        if (rc < 0 && errno != EEXIST) {
            fprintf(stderr, "mkdir(%s): %s\n", d.c_str(), strerror(errno));
        }
    }
}

/* return an allocated string with the name of a balancing file for this
 * matrix and this mpi/thr split.
 */
static std::string matrix_get_derived_balancing_filename(std::string const & matrixname, parallelizing_info_ptr pi)
{
    /* input is NULL in the case of random matrices */
    if (matrixname.empty()) return std::string();
    std::string dn = matrix_get_derived_cache_subdir(matrixname, pi);
    char * t;
    int rc = asprintf(&t, "%s/%s.bin", dn.c_str(), dn.c_str());
    ASSERT_ALWAYS(rc >=0);
    std::string s(t);
    free(t);
    return s;
}

static std::string matrix_get_derived_cache_filename_stem(std::string const & matrixname, parallelizing_info_ptr pi, uint32_t checksum)
{
    /* input is NULL in the case of random matrices */
    if (matrixname.empty()) return std::string();
    char * copy = strdup(matrixname.c_str());
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int pos[2];
    for(int d = 0 ; d < 2 ; d++)  {
        pi_comm_ptr wr = pi->wr[d];
        pos[d] = wr->jrank * wr->ncores + wr->trank;
    }
    std::string dn = matrix_get_derived_cache_subdir(matrixname, pi);
    int rc = asprintf(&t, "%s/%s.%08" PRIx32 ".h%d.v%d",
            dn.c_str(), dn.c_str(), checksum, pos[1], pos[0]);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    std::string s(t);
    free(t);
    return s;
}

static std::string matrix_get_derived_submatrix_filename(std::string const & matrixname, parallelizing_info_ptr pi)
{
    /* input is empty in the case of random matrices */
    if (matrixname.empty()) return std::string();
    char * copy = strdup(matrixname.c_str());
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int pos[2];
    for(int d = 0 ; d < 2 ; d++)  {
        pi_comm_ptr wr = pi->wr[d];
        pos[d] = wr->jrank * wr->ncores + wr->trank;
    }
    std::string dn = matrix_get_derived_cache_subdir(matrixname, pi);
    int rc = asprintf(&t, "%s/%s.h%d.v%d.bin",
            dn.c_str(), dn.c_str(), pos[1], pos[0]);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    std::string s(t);
    free(t);
    return s;
}

static void matmul_top_init_fill_balancing_header(matmul_top_data & mmt, int i, matmul_top_matrix & Mloc, param_list_ptr pl)
{
    parallelizing_info_ptr pi = mmt.pi;

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        if (Mloc.mname.empty()) {
            random_matrix_fill_fake_balancing_header(Mloc.bal, pi, param_list_lookup_string(pl, "random_matrix"));
        } else {
            if (access(Mloc.bname.c_str(), R_OK) != 0) {
                if (errno == ENOENT) {
                    printf("Creating balancing file %s\n", Mloc.bname.c_str());
                    struct mf_bal_args mba = {
                        .rwfile = NULL,
                        .cwfile = NULL,
                        .mfile = Mloc.mname.empty() ? NULL : Mloc.mname.c_str(),
                        .bfile = Mloc.bname.empty() ? NULL : Mloc.bname.c_str(),
                        .quiet = 0,
                        .nh = (int) pi->wr[1]->totalsize,
                        .nv = (int) pi->wr[0]->totalsize,
                        .withcoeffs = !mmt.abase->is_characteristic_two(),
                        .rectangular = 0,
                        .skip_decorrelating_permutation = 0,
                        .do_perm = { mf_bal_args::MF_BAL_PERM_AUTO, mf_bal_args::MF_BAL_PERM_AUTO },
                    };
                    mf_bal_adjust_from_option_string(&mba, param_list_lookup_string(pl, "balancing_options"));
                    /* withcoeffs being a switch for param_list, it is
                     * clobbered by the configure_switch mechanism */
                    mba.withcoeffs = !mmt.abase->is_characteristic_two();
                    matrix_create_derived_cache_subdir(Mloc.mname, mmt.pi);

                    mf_bal(&mba);
                } else {
                    fprintf(stderr, "Cannot access balancing file %s: %s\n", Mloc.bname.c_str(), strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }
            balancing_read_header(Mloc.bal, Mloc.bname);
        }
    }
    pi_bcast(&Mloc.bal, sizeof(balancing), BWC_PI_BYTE, 0, 0, mmt.pi->m);

    /* check that balancing dimensions are compatible with our run */
    int ok = 1;
    ok = ok && mmt.pi->wr[0]->totalsize == Mloc.bal.nv;
    ok = ok && mmt.pi->wr[1]->totalsize == Mloc.bal.nh;
    if (ok) return;

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        fprintf(stderr, "Matrix %d, %s: balancing file %s"
                " has dimensions %ux%u,"
                " this conflicts with the current run,"
                " which expected_last_iteration dimensions (%ux%u)x(%ux%u).\n",
                i, Mloc.mname.c_str(), Mloc.bname.c_str(),
                Mloc.bal.nh, Mloc.bal.nv,
                mmt.pi->wr[1]->njobs,
                mmt.pi->wr[1]->ncores,
                mmt.pi->wr[0]->njobs,
                mmt.pi->wr[0]->ncores);
    }
    serialize(mmt.pi->m);
    exit(1);
}


static void matmul_top_init_prepare_local_permutations(matmul_top_data & mmt, int i)
{
    matmul_top_matrix & Mloc = mmt.matrices[i];
    /* Here, we get a copy of the rowperm and colperm.
     *
     * For each (job,thread), two pairs of intervals are defined.
     *
     * for row indices: [i0[0], i1[0][ = [mrow->i0, mrow->i1[ and [Xi0[0], Xi1[0][
     * for col indices: [i0[1], i1[1][ = [mcol->i0, mcol->i1[ and [Xi0[1], Xi1[1][
     *
     * The parts we get are:
     *      [i, rowperm[i]] for    i in [mrow->i0, mrow->i1[
     *                         and rowperm[i] in [X0, X1[
     *      [j, colperm[j]] for    j in [mcol->i0, mcol->i1[
     *                         and colperm[i] in [Y0, Y1[
     * where X0 is what would be mcol->i0 if the matrix as many columns as
     * rows, and ditto for X1,Y0,Y1.
     */

    unsigned int rowperm_items=0;
    unsigned int colperm_items=0;

    /* Define a complete structure for the balancing which is shared
     * among threads, so that we'll be able to access it from all threads
     * simultaneously. We will put things in bal_tmp.rowperm and
     * bal_tmp.colperm, but beyond that, the header part will be wrong
     * at non-root nodes.
     */
    balancing * pbal_tmp = (balancing *) shared_malloc_set_zero(mmt.pi->m, sizeof(balancing));
    balancing & bal_tmp = * pbal_tmp;

    if (mmt.pi->m->jrank == 0 && mmt.pi->m->trank == 0) {
        if (!Mloc.bname.empty())
            balancing_read(bal_tmp, Mloc.bname);
        /* It's fine if we have nothing. This just means that we'll have
         * no balancing to deal with (this occurs only for matrices
         * which are generated at random on the fly). */
        rowperm_items = bal_tmp.rowperm != NULL ? bal_tmp.trows : 0;
        colperm_items = bal_tmp.colperm != NULL ? bal_tmp.tcols : 0;
    }
    pi_bcast(&rowperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt.pi->m);
    pi_bcast(&colperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt.pi->m);

    if (mmt.pi->m->trank == 0) {
        if (rowperm_items) {
            ASSERT_ALWAYS(rowperm_items == Mloc.bal.trows);
            if (mmt.pi->m->jrank != 0)
                bal_tmp.rowperm = (uint32_t *) malloc(Mloc.bal.trows * sizeof(uint32_t));
            MPI_Bcast(bal_tmp.rowperm, Mloc.bal.trows * sizeof(uint32_t), MPI_BYTE, 0, mmt.pi->m->pals);
        }
        if (colperm_items) {
            ASSERT_ALWAYS(colperm_items == Mloc.bal.tcols);
            if (mmt.pi->m->jrank != 0)
                bal_tmp.colperm = (uint32_t *) malloc(Mloc.bal.tcols * sizeof(uint32_t));
            MPI_Bcast(bal_tmp.colperm, Mloc.bal.tcols * sizeof(uint32_t), MPI_BYTE, 0, mmt.pi->m->pals);
        }
    }
    serialize_threads(mmt.pi->m);      /* important ! */

    uint32_t * balperm[2] = { bal_tmp.rowperm, bal_tmp.colperm };
    for(int d = 0 ; d < 2 ; d++)  {
        unsigned int ii[2];
        unsigned int jj[2];
        get_local_permutations_ranges(mmt, d, ii, jj);

        if (!balperm[d]) continue;

        static pthread_mutex_t pp = PTHREAD_MUTEX_INITIALIZER;

        pthread_mutex_lock(&pp);
        Mloc.perm[d].reserve(ii[1] - ii[0]);
        pthread_mutex_unlock(&pp);

        /* now create the really local permutation */
        for(unsigned int i = ii[0] ; i < ii[1] ; i++) {
            unsigned int j = balperm[d][i];
            if (j >= jj[0] && j < jj[1])
                Mloc.perm[d].push_back({i, j});
        }
#if 0
        const char * text[2] = { "left", "right", };
        printf("[%s] J%uT%u does %zu/%u permutation pairs for %s vectors\n",
                mmt.pi->nodenumber_s,
                mmt.pi->m->jrank, mmt.pi->m->trank,
                Mloc.perm[d]->n, d ? Mloc.bal.tcols : Mloc.bal.trows,
                text[d]);
#endif
    }

    serialize_threads(mmt.pi->m);      /* important ! */

    if (mmt.pi->m->trank == 0) {
        if (bal_tmp.colperm) free(bal_tmp.colperm);
        if (bal_tmp.rowperm) free(bal_tmp.rowperm);
    }

    shared_free(mmt.pi->m, pbal_tmp);

    /* Do this so that has_perm makes sense globally */
    Mloc.has_perm_map[0] = Mloc.has_perm_local(0);
    Mloc.has_perm_map[1] = Mloc.has_perm_local(1);
    pi_allreduce(NULL, Mloc.has_perm_map, 2, BWC_PI_INT, BWC_PI_MAX, mmt.pi->m);
    if (mmt.pi->m->trank == 0) {
        MPI_Allreduce(MPI_IN_PLACE, Mloc.has_perm_map, 2, MPI_INT, MPI_MAX, mmt.pi->m->pals);
    }
    pi_bcast(Mloc.has_perm_map, 2, BWC_PI_INT, 0, 0, mmt.pi->m);
}

matmul_top_data::matmul_top_data(
        arith_generic * abase,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        param_list_ptr pl,
        int optimized_direction)
    : abase(abase)
    , pi(pi)
    , pitype(pi_alloc_arith_datatype(pi, abase))
{
    matmul_top_data & mmt = *this;

    int nbals = param_list_get_list_count(pl, "balancing");
    int multimat = 0;
    int nmatrices = param_list_lookup_string(pl, "matrix") != NULL;
    param_list_parse_int(pl, "multi_matrix", &multimat);
    if (multimat)
        nmatrices = param_list_get_list_count(pl, "matrix");
    const char * random_description = param_list_lookup_string(pl, "random_matrix");
    const char * static_random_matrix = param_list_lookup_string(pl, "static_random_matrix");


    if (random_description || static_random_matrix) {
        if (nbals || nmatrices) {
            fprintf(stderr, "random_matrix is incompatible with balancing= and matrix=\n");
            exit(EXIT_FAILURE);
        }
    } else if (nbals && !nmatrices) {
        fprintf(stderr, "missing parameter matrix=\n");
        exit(EXIT_FAILURE);
    } else if (!nbals && nmatrices) {
        /* nbals == 0 is a hint towards taking the default balancing file
         * names, that's it */
    } else if (nbals != nmatrices) {
        fprintf(stderr, "balancing= and matrix= have inconsistent number of items\n");
        exit(EXIT_FAILURE);
    }

    if (random_description)
        nmatrices = 1;
    if (static_random_matrix)
        nmatrices = 1;

    serialize_threads(mmt.pi->m);

    /* The initialization goes through several passes */
    for(int i = 0 ; i < nmatrices ; i++) {
        matmul_top_matrix Mloc;
        if (multimat) {
            Mloc.mname = matrix_list_get_item(pl, "matrix", i);
            Mloc.bname = matrix_list_get_item(pl, "balancing", i);
        } else {
            const char * t;
            t = param_list_lookup_string(pl, "matrix");
            Mloc.mname = t ? t : "";
            t = param_list_lookup_string(pl, "balancing");
            Mloc.bname = t ? t : "";
        }
        if (static_random_matrix) {
            ASSERT_ALWAYS(i == 0);
            Mloc.mname = strdup(static_random_matrix);
        }
        if (Mloc.bname.empty()) {
            /* returns NULL is mname is NULL */
            Mloc.bname = matrix_get_derived_balancing_filename(Mloc.mname, mmt.pi);
        }
        /* At this point mname and bname are either NULL or freshly
         * allocated */
        ASSERT_ALWAYS((!Mloc.bname.empty()) == !random_description);

        matmul_top_init_fill_balancing_header(mmt, i, Mloc, pl);

        Mloc.n[0] = Mloc.bal.trows;
        Mloc.n[1] = Mloc.bal.tcols;
        Mloc.n0[0] = Mloc.bal.nrows;
        Mloc.n0[1] = Mloc.bal.ncols;
        Mloc.locfile = matrix_get_derived_cache_filename_stem(Mloc.mname, mmt.pi, Mloc.bal.checksum);

        mmt.matrices.push_back(Mloc);
    }

    mmt.n[0] = mmt.matrices[0].n[0];
    mmt.n0[0] = mmt.matrices[0].n0[0];
    mmt.n[1] = mmt.matrices[mmt.matrices.size()-1].n[1];
    mmt.n0[1] = mmt.matrices[mmt.matrices.size()-1].n0[1];

    /* in the second loop below, get_local_permutations_ranges uses
     * mmt.n[], so we do it in a second pass.
     * Now, given that double matrices only barely work at the moment,
     * I'm not absolutely sure that it's really needed.
     */
    for(size_t i = 0 ; i < mmt.matrices.size() ; i++) {
        matmul_top_matrix & Mloc = mmt.matrices[i];

        matmul_top_init_prepare_local_permutations(mmt, i);

        if (!mmt.pi->interleaved) {
            matmul_top_read_submatrix(mmt, i, pl, optimized_direction );
        } else {
            /* Interleaved threads will share their matrix data. The first
             * thread to arrive will do the initialization, and the second
             * one will just grab the pointer. The trick is to be able to
             * pick the pointer in the right location ! We add a generic
             * pointer dictionary feature in the parallelizing_info
             * interface for this purpose.
             */

#define MMT_MM_MAGIC_KEY        0xaa000000UL

            if (mmt.pi->interleaved->idx == 0) {
                matmul_top_read_submatrix(mmt, i, pl, optimized_direction);
                pi_store_generic(mmt.pi, MMT_MM_MAGIC_KEY + i, mmt.pi->m->trank, Mloc.mm);
            } else {
                Mloc.mm = (matmul_ptr) pi_load_generic(mmt.pi, MMT_MM_MAGIC_KEY + i, mmt.pi->m->trank);
            }
        }
    }
}

unsigned int matmul_top_rank_upper_bound(matmul_top_data & mmt)
{
    unsigned int r = MAX(mmt.n0[0], mmt.n0[1]);
    for(auto const & Mloc : mmt.matrices) {
        r = MAX(r, Mloc.bal.nrows - Mloc.bal.nzrows);
        r = MAX(r, Mloc.bal.ncols - Mloc.bal.nzcols);
    }
    return r;
}


static int export_cache_list_if_requested(matmul_top_matrix & Mloc, parallelizing_info_ptr pi, param_list_ptr pl)
{
    const char * cachelist = param_list_lookup_string(pl, "export_cachelist");
    if (!cachelist) return 0;

    char * myline = NULL;
    int rc;
    rc = asprintf(&myline, "%s %s", pi->nodename, Mloc.mm->cachefile_name);
    ASSERT_ALWAYS(rc >= 0);
    ASSERT_ALWAYS(myline != NULL);
    char const ** tlines = NULL;
    tlines = (char const **) shared_malloc_set_zero(pi->m, pi->m->ncores * sizeof(const char *));
    tlines[pi->m->trank] = myline;
    serialize_threads(pi->m);

    /* Also, just out of curiosity, try to see what we have currently */
    struct stat st[1];
    int * has_cache = (int *) shared_malloc_set_zero(pi->m, pi->m->totalsize * sizeof(int));
    rc = stat(Mloc.mm->cachefile_name, st);
    unsigned int mynode = pi->m->ncores * pi->m->jrank;
    has_cache[mynode + pi->m->trank] = rc == 0;
    serialize_threads(pi->m);

    int len = 0;
    if (pi->m->trank == 0) {
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                has_cache, pi->m->ncores, MPI_INT, pi->m->pals);

        for(unsigned int j = 0 ; j < pi->m->ncores ; j++) {
            int s = strlen(tlines[j]);
            if (s >= len) len = s;
        }
        MPI_Allreduce(MPI_IN_PLACE, &len, 1, MPI_INT, MPI_MAX, pi->m->pals);
        char * info = (char *) malloc(pi->m->totalsize * (len + 1));
        char * mybuf = (char *) malloc(pi->m->ncores * (len+1));
        memset(mybuf, 0, pi->m->ncores * (len+1));
        for(unsigned int j = 0 ; j < pi->m->ncores ; j++) {
            memcpy(mybuf + j * (len+1), tlines[j], strlen(tlines[j])+1);
        }
        MPI_Allgather(mybuf, pi->m->ncores * (len+1), MPI_BYTE,
                info, pi->m->ncores * (len+1), MPI_BYTE,
                pi->m->pals);
        if (pi->m->jrank == 0) {
            FILE * f = fopen(cachelist, "wb");
            DIE_ERRNO_DIAG(f == NULL, "fopen(%s)", cachelist);
            for(unsigned int j = 0 ; j < pi->m->njobs ; j++) {
                unsigned int j0 = j * pi->m->ncores;
                fprintf(f, "get-cache ");
                for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            for(unsigned int j = 0 ; j < pi->m->njobs ; j++) {
                unsigned int j0 = j * pi->m->ncores;
                fprintf(f, "has-cache ");
                for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    if (!has_cache[pi->m->ncores * j + k]) continue;
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        free(info);
        free(mybuf);
    }
    serialize_threads(pi->m);
    shared_free(pi->m, tlines);
    shared_free(pi->m, has_cache);
    serialize_threads(pi->m);
    free(myline);
    serialize(pi->m);

    return 1;
}

static unsigned int local_fraction(unsigned int normal, pi_comm_ptr wr)
{
    unsigned int i = wr->jrank * wr->ncores + wr->trank;
    return normal / wr->totalsize + (i < (normal % wr->totalsize));
}


static void matmul_top_read_submatrix(matmul_top_data & mmt, int midx, param_list_ptr pl, int optimized_direction)
{
    int rebuild = 0;
    param_list_parse_int(pl, "rebuild_cache", &rebuild);
    int can_print = (mmt.pi->m->jrank == 0 && mmt.pi->m->trank == 0);

    matmul_top_matrix & Mloc = mmt.matrices[midx];

    Mloc.mm = matmul_init(mmt.abase,
            Mloc.n[0] / mmt.pi->wr[1]->totalsize,
            Mloc.n[1] / mmt.pi->wr[0]->totalsize,
            Mloc.locfile.empty() ? NULL : Mloc.locfile.c_str(),
            NULL  /* means: choose mm_impl from pl */,
            pl, optimized_direction);

    // *IF* we need to do a collective read of the matrix, we need to
    // provide the pointer *now*.
    unsigned int sqread = 0;
    param_list_parse_uint(pl, "sequential_cache_read", &sqread);

    int cache_loaded = 0;

    if (export_cache_list_if_requested(Mloc, mmt.pi, pl)) {
        /* If we are being called from dispatch, once all submatrices
         * have had their list of required files printed, the program
         * will exit. */
        return;
    }

    if (!rebuild) {
        if (can_print) {
            printf("Now trying to load matrix cache files\n");
        }
        if (sqread) {
            for(unsigned int j = 0 ; j < mmt.pi->m->ncores ; j += sqread) {
                serialize_threads(mmt.pi->m);
                double t_read = -wct_seconds();
                if (j / sqread == mmt.pi->m->trank / sqread)
                    cache_loaded = matmul_reload_cache(Mloc.mm);
                serialize(mmt.pi->m);
                t_read += wct_seconds();
                if (mmt.pi->m->jrank == 0 && mmt.pi->m->trank == j && cache_loaded) {
                    printf("[%s] J%uT%u-%u: read cache %s (and others) in %.2fs (round %u/%u)\n",
                    mmt.pi->nodenumber_s,
                    mmt.pi->m->jrank,
                    mmt.pi->m->trank,
                    MIN(mmt.pi->m->ncores, mmt.pi->m->trank + sqread) - 1,
                    Mloc.mm->cachefile_name,
                    t_read,
                    j / sqread, iceildiv(mmt.pi->m->ncores, sqread)
                    );
                }
            }
        } else {
            double t_read = -wct_seconds();
            cache_loaded = matmul_reload_cache(Mloc.mm);
            serialize(mmt.pi->m);
            t_read += wct_seconds();
            if (mmt.pi->m->jrank == 0 && mmt.pi->m->trank == 0 && cache_loaded) {
                printf("[%s] J%u: read cache %s (and others) in %.2fs\n",
                        mmt.pi->nodenumber_s,
                        mmt.pi->m->jrank,
                        Mloc.mm->cachefile_name,
                        t_read
                      );
            }
        }
        if (!mmt.pi->m->trank) {
            printf("J%u %s done reading (result=%d)\n", mmt.pi->m->jrank, mmt.pi->nodename, cache_loaded);
        }
    }

    if (!pi_data_eq(&cache_loaded, 1, BWC_PI_INT, mmt.pi->m)) {
        if (can_print) {
            fprintf(stderr, "Fatal error: cache files not present at expected locations\n");
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt.pi->m) {
            fprintf(stderr, "[%s] J%uT%u: cache %s: %s\n",
                    mmt.pi->nodenumber_s,
                    mmt.pi->m->jrank,
                    mmt.pi->m->trank,
                    Mloc.mm->cachefile_name,
                    cache_loaded ? "ok" : "not ok");
        }
        SEVERAL_THREADS_PLAY_MPI_END();
        serialize(mmt.pi->m);
        abort();
    }

    unsigned int sqb = 0;
    param_list_parse_uint(pl, "sequential_cache_build", &sqb);

    matrix_u32 m;
    memset(m, 0, sizeof(matrix_u32));
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (!cache_loaded) {
        // the mm layer is informed of the higher priority computations
        // that will take place. Depending on the implementation, this
        // may cause the direct or transposed ordering to be preferred.
        // Thus we have to read this back from the mm structure.
        m->mfile = Mloc.mname.empty() ? NULL : Mloc.mname.c_str();
        m->bfile = Mloc.bname.empty() ? NULL : Mloc.bname.c_str();
        m->transpose = Mloc.mm->store_transposed;
        m->withcoeffs = !mmt.abase->is_characteristic_two();
        if (Mloc.mname.empty()) {
            if (can_print) {
                printf("Begin creation of fake matrix data in parallel\n");
            }
            /* Mloc.mm->dim[0,1] contains the dimensions of the padded
             * matrix. This is absolutely fine in the normal case. But in
             * the case of staged matrices, it's a bit different. We must
             * make sure that we generate matrices which have zeroes in
             * the padding area.
             */

            unsigned int data_nrows = local_fraction(Mloc.n0[0], mmt.pi->wr[1]);
            unsigned int data_ncols = local_fraction(Mloc.n0[1], mmt.pi->wr[0]);
            unsigned int padded_nrows = Mloc.n[0] / mmt.pi->wr[1]->totalsize;
            unsigned int padded_ncols = Mloc.n[1] / mmt.pi->wr[0]->totalsize;

            random_matrix_get_u32(mmt.pi, pl, m,
                    data_nrows, data_ncols, padded_nrows, padded_ncols);
        } else {
            if (can_print) {
                printf("Matrix dispatching starts\n");
            }

            /* It might be that the leader and the other nodes do _not_
             * share a common filesystem, in which case we must do this
             * also here.
             */
            if (mmt.pi->m->trank == 0)
                matrix_create_derived_cache_subdir(Mloc.mname, mmt.pi);
            serialize_threads(mmt.pi->m);

            balancing_get_matrix_u32(mmt.pi, pl, m);

            int ssm = 0;
            param_list_parse_int(pl, "save_submatrices", &ssm);
            if (ssm) {
                std::string submat = matrix_get_derived_submatrix_filename(Mloc.mname, mmt.pi);
                fprintf(stderr, "DEBUG: creating %s\n", submat.c_str());
                FILE * f = fopen(submat.c_str(), "wb");
                fwrite(m->p, sizeof(uint32_t), m->size, f);
                fclose(f);
            }
        }
    }

    if (can_print && !Mloc.bname.empty()) {
        balancing bal;
        balancing_init(bal);
        balancing_read_header(bal, Mloc.bname);
        balancing_set_row_col_count(bal);
        printf("Matrix: total %" PRIu32 " rows %" PRIu32 " cols "
                "%" PRIu64 " coeffs\n",
                bal.nrows, bal.ncols, bal.ncoeffs);
        balancing_clear(bal);
    }

    if (!sqb) {
        if (!cache_loaded) {
            // everybody does it in parallel
            if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                printf("[%s] J%uT%u building cache for %s\n",
                        mmt.pi->nodenumber_s,
                        mmt.pi->m->jrank,
                        mmt.pi->m->trank,
                        Mloc.locfile.c_str());
            matmul_build_cache(Mloc.mm, m);
            matmul_save_cache(Mloc.mm);
        }
    } else {
        if (can_print)
            printf("Building local caches %d at a time\n", sqb);
        for(unsigned int j = 0 ; j < mmt.pi->m->ncores + sqb ; j += sqb) {
            serialize_threads(mmt.pi->m);
            if (cache_loaded) continue;
            if (j / sqb == mmt.pi->m->trank / sqb) {
                if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                    printf("[%s] J%uT%u building cache for %s\n",
                            mmt.pi->nodenumber_s,
                            mmt.pi->m->jrank,
                            mmt.pi->m->trank,
                            Mloc.locfile.c_str());
                matmul_build_cache(Mloc.mm, m);
            } else if (j / sqb == mmt.pi->m->trank / sqb + 1) {
                matmul_save_cache(Mloc.mm);
            }
        }
    }

    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (Mloc.mm->cachefile_name && verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO)) {
        my_pthread_mutex_lock(mmt.pi->m->th->m);
        printf("[%s] J%uT%u uses cache file %s\n",
                mmt.pi->nodenumber_s,
                mmt.pi->m->jrank, mmt.pi->m->trank,
                /* cache for mmt.locfile, */
                Mloc.mm->cachefile_name);
        my_pthread_mutex_unlock(mmt.pi->m->th->m);
    }
}

void matmul_top_report(matmul_top_data & mmt, double scale, int full)
{
    for(auto const & Mloc : mmt.matrices) {
        matmul_report(Mloc.mm, scale);
        size_t max_report_size = 0;
        pi_allreduce(&Mloc.mm->report_string_size, &max_report_size, 1, BWC_PI_SIZE_T, BWC_PI_MAX, mmt.pi->m);
        char * all_reports = (char *) malloc(mmt.pi->m->totalsize * max_report_size);
        memset(all_reports, 0, mmt.pi->m->totalsize * max_report_size);
        memcpy(all_reports + max_report_size * (mmt.pi->m->jrank * mmt.pi->m->ncores + mmt.pi->m->trank), Mloc.mm->report_string, Mloc.mm->report_string_size);
        pi_allgather(NULL, 0, 0,
                all_reports, max_report_size, BWC_PI_BYTE, mmt.pi->m);

        if (max_report_size > 1 && mmt.pi->m->jrank == 0 && mmt.pi->m->trank == 0) {
            for(unsigned int j = 0 ; j < mmt.pi->m->njobs ; j++) {
                for(unsigned int t = 0 ; t < mmt.pi->m->ncores ; t++) {
                    char * locreport = all_reports + max_report_size * (j * mmt.pi->m->ncores + t);
                    if (full || (j == 0 && t == 0))
                        printf("##### J%uT%u timing report:\n%s", j, t, locreport);
                }
            }
        }
        serialize(mmt.pi->m);
        free(all_reports);
    }
}

matmul_top_data::~matmul_top_data()
{
    matmul_top_data & mmt = *this;
    pi_free_arith_datatype(mmt.pi, mmt.pitype);
    serialize_threads(mmt.pi->m);
    serialize(mmt.pi->m);
    serialize_threads(mmt.pi->m);

    for(auto const & Mloc : mmt.matrices) {
        serialize(mmt.pi->m);
        if (!mmt.pi->interleaved) {
            matmul_clear(Mloc.mm);
        } else if (mmt.pi->interleaved->idx == 1) {
            matmul_clear(Mloc.mm);
            /* group 0 is the first to leave, thus it doesn't to freeing.
            */
        }
    }
    serialize(mmt.pi->m);
}
