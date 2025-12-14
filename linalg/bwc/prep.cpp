#include "cado.h" // IWYU pragma: keep

#include <cstdint> // for uint32_t
#include <cstdio>
#include <cstdlib>
#include <cstring> // for memset
#include <climits>
#include <ctime>   // for time

#include <algorithm>
#include <ios>
#include <istream>
#include <fstream>
#include <memory>
#include <set>
#include <utility>
#include <vector>
#include <string>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "abase_proxy.hpp"
#include "arith-generic.hpp"
#include "bw-common.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "matmul_top_vec.hpp"
#include "mmt_vector_pair.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "select_mpi.h"
#include "xdotprod.hpp"
#include "xvectors.hpp"
#include "bwc_filenames.hpp"

static void bw_rank_check(matmul_top_data & mmt, cxx_param_list & pl)
{
    int const tcan_print = bw->can_print && mmt.pi->m->trank == 0;
    unsigned int const r = matmul_top_rank_upper_bound(mmt);
    if (tcan_print) {
        fmt::print("Matrix rank is at most {} (based on zero columns and rows "
               "encountered)\n",
               r);
    }
    int skip = 0;
    param_list_parse_int(pl, "skip_bw_early_rank_check", &skip);
    if (bw->m + r < mmt.n0[0]) {
        fmt::print(
            stderr,
            "Based on the parameter m (={}) and the rank defect of the matrix "
            "(>={}), we can't expect to compute solutions reliably.\n",
            bw->m, mmt.n0[0] - r);
        if (skip) {
            fmt::print(stderr,
                    "Proceeding anyway as per skip_bw_early_rank_check=1\n");
        } else {
            fmt::print(stderr, "Aborting. Use skip_bw_early_rank_check=1 to "
                            "proceed nevertheless.\n");
            exit(EXIT_FAILURE);
        }
    }
}

static unsigned int matrix_rank(arith_generic * A, arith_generic::elt * M,
                         unsigned int nrows, unsigned int ncols)
{
    /* This is an abstract Gauss based on the arith_generic layer. It's
     * certainly not meant to be fast. If performance is an issue, then
     * we should specialize and bring it down to the arith_hard level.
     *
     * This _destructively_ computes the rank of the matrix pointed by M.
     * The matrix is expected to be in row-major order, with ncols
     * grouped by A->simd_groupsize() items
     */
    ASSERT_ALWAYS(ncols % A->simd_groupsize() == 0);
    unsigned int const row_stride = ncols / A->simd_groupsize();
    unsigned int rank = 0;
    auto elt = A->alloc_vector(1);
    auto buf = A->alloc_vector(row_stride);
    for (unsigned int r = 0; r < nrows; r++) {
        /* find the row with least significant first set bit among the
         * rows M[r:]
         */
        int jpiv = INT_MAX;
        unsigned int rpiv = UINT_MAX;
        for (unsigned int r1 = r; r1 < nrows; r1++) {
            auto * row1 = A->vec_subvec(M, r1 * row_stride);
            size_t const j1 = A->vec_simd_find_first_set(*elt, row1, row_stride);
            if (j1 == SIZE_MAX)
                continue;
            if (int(j1) < jpiv) {
                jpiv = int(j1);
                rpiv = r1;
            }
        }
        if (jpiv == INT_MAX)
            return rank;
        /* swap rows r and rpiv */
        if (r != rpiv) {
            auto * row = A->vec_subvec(M, r * row_stride);
            auto * rowpiv = A->vec_subvec(M, rpiv * row_stride);
            A->vec_set(buf.get(), row, row_stride);
            A->vec_set(row, rowpiv, row_stride);
            A->vec_set(rowpiv, buf.get(), row_stride);
        }
        int const j = jpiv;
        auto * row = A->vec_subvec(M, r * row_stride);
        /* find the smallest element again */
        A->vec_simd_find_first_set(*elt, row, row_stride);
        A->inverse(*elt, *elt);
        A->vec_set_zero(buf.get(), row_stride);
        A->vec_addmul_and_reduce(buf.get(), row, *elt, row_stride);
        A->vec_neg(buf.get(), buf.get(), row_stride);
        for (unsigned int r1 = r + 1; r1 < nrows; r1++) {
            auto * row1 = A->vec_subvec(M, r1 * row_stride);
            size_t const j1 = A->vec_simd_find_first_set(*elt, row1, row_stride);
            ASSERT_ALWAYS(j1 == SIZE_MAX || int(j1) >= j);
            if (int(j1) == j) // we have a coefficient to cancel
                A->vec_addmul_and_reduce(row1, buf.get(), *elt, row_stride);
        }
        rank++;
    }

    return rank;
}

struct rhs_header {
    uint32_t nrows = 0;
    unsigned int ncols = 0;
    cxx_mpz p = 0;
};

static std::istream & operator>>(std::istream & is, rhs_header & hdr)
{
    return is >> hdr.nrows >> hdr.ncols >> hdr.p;
}

/* We have several options to read the rhs coefficients: read the full
 * vector from every node/thread, or read from only one place and
 * broadcast the result. We prefer the latter, as it's probably better in
 * terms of filesystem pressure for large jobs.
 */

static void read_rhs_from_file(std::vector<mmt_vec> & rhs_vecs, std::istream * is,
                        size_t itemsondisk)
{
    if (rhs_vecs.empty())
        return;

    /* {std::istream * is} is actually nullptr if pi->m->jrank != 0 */

    arith_generic * A = rhs_vecs[0].abase;
    parallelizing_info_ptr pi = rhs_vecs[0].pi;
    int const d = rhs_vecs[0].d;
    auto i1 = rhs_vecs[0].i1;
    auto i0 = rhs_vecs[0].i0;
    auto eitems = i1 - i0;

    pi_comm_ptr m = pi->m;
    pi_comm_ptr wr = pi->wr[d];
    pi_comm_ptr xwr = wr->xwr;

    size_t const nrhs = rhs_vecs.size() * A->simd_groupsize();

    for (auto & r: rhs_vecs)
        mmt_full_vec_set_zero(r);

    /* similar to {std::istream * is}, create a reference leader_vecs to
     * a vector that has size 0 except at rank 0
     */
    std::vector<arith_generic::owned_vector> local_vecs;
    if (m->trank == 0 && m->jrank == 0) {
        for (auto const & r MAYBE_UNUSED: rhs_vecs)
            local_vecs.emplace_back(A->alloc_vector(eitems));
    }
    /* make leader_vecs accessible to all threads at pi->m->jrank == 0
     * and all thread ranks along xwr */
    auto * p_local_vecs = &local_vecs;
    pi_bcast(&p_local_vecs, sizeof(void *), BWC_PI_BYTE, 0, 0, xwr);
    auto & leader_vecs(*p_local_vecs);

    cxx_mpz c;

    for (unsigned int jpeer = 0; jpeer < xwr->njobs; jpeer++) {
        if (wr->trank || wr->jrank) {
            /* what happens along wr will be dealt with by a
             * broadcast along this communicator
             */
            continue;
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN (xwr) {
            /* we'll use this as an MPI tag. It's important that we
             * include a thread id in there, since there will be multiple
             * transfers with same (src, dst) job ids!
             */
            unsigned int const tpeer = xwr->trank;
            unsigned int const round = jpeer * xwr->ncores + tpeer;

            if (xwr->jrank == 0) {
                /* read eitems from the file. These will go to the own
                 * subvec of (jpeer, tpeer), and correspond to indices
                 * that start at (jpeer * wr->ncores + tpeer) * eitems.
                 *
                 * This can later be made consistent with all
                 * jobs/threads along wr by an allreduce operation.
                 */

                unsigned int const offset = round * eitems;
                unsigned int const bound = MIN(offset + eitems, itemsondisk);

                for (auto & r: leader_vecs)
                    A->vec_set_zero(r.get(), eitems);

                for (unsigned int i = offset ; i < bound ; i++) {
                    for (auto & r: leader_vecs) {
                        /* TODO: a binary version would only have to
                         * update these three lines
                         */
                        (*is) >> c;
                        ASSERT_ALWAYS(bool(*is));
                        A->set(A->vec_item(r.get(), i - offset), c);
                    }
                }
            }

            /* send to the right recipient */
            for (size_t j = 0; j * A->simd_groupsize() < nrhs; j++) {
                if (jpeer != 0) {
                    if (xwr->jrank == 0)
                        MPI_Send(leader_vecs[j].get(),
                                 (int)A->vec_elt_stride(eitems), MPI_BYTE,
                                 (int)jpeer, (int)round, xwr->pals);
                    else if (xwr->jrank == jpeer)
                        MPI_Recv(rhs_vecs[j].v, (int)A->vec_elt_stride(eitems),
                                 MPI_BYTE, 0, (int)round, xwr->pals,
                                 MPI_STATUS_IGNORE);
                } else if (xwr->jrank == 0) {
                    A->vec_set(rhs_vecs[j].v, leader_vecs[j].get(), eitems);
                }
            }
        }
        SEVERAL_THREADS_PLAY_MPI_END();
    }

    serialize_threads(xwr);

    /* it's an easy enough way to get everything consistent. */
    for (auto & r: rhs_vecs) {
        r.consistency = 1;
        mmt_vec_allreduce(r);
    }
}

struct prep_object {
    parallelizing_info_ptr pi;
    int const tcan_print;
    int const char2;
    int const splitwidth;
    unsigned int const A_multiplex;
    std::unique_ptr<arith_generic> A;
    matmul_top_data mmt;
    cxx_gmp_randstate rstate;
    mmt_vector_pair ymy;
    mmt_vec & y;
    unsigned int const unpadded;
    unsigned int nrhs = 0;
    std::vector<uint32_t> xvecs;
    std::vector<unsigned int> Z;

    prep_object(parallelizing_info_ptr pi, cxx_param_list & pl)
        : pi(pi)
        , tcan_print(bw->can_print && pi->m->trank == 0)
        , char2(mpz_cmp_ui(bw->p, 2) == 0)
        , splitwidth(char2 ? 64 : 1)
        , A_multiplex(bw->n / splitwidth)
        , A(arith_generic::instance(bw->p, splitwidth))
        , mmt(A.get(), pi, pl, bw->dir)
        , ymy(mmt, bw->dir)
        , y(ymy[0])
        , unpadded(MAX(mmt.n0[0], mmt.n0[1]))
    {
        /* Interleaving does not make sense for this program. So the second
         * block of threads just leave immediately */
        ASSERT_ALWAYS(!pi->interleaved);

        // Doing the ``hello world'' test is a very good way of testing the
        // global mpi/pthreads setup. So despite its apparent irrelevance, I
        // suggest leaving it here as a cheap sanity check.
        pi_hello(pi);

        // I don't think multi-matrix was ever tested beyond the case p==2.
        ASSERT_ALWAYS(char2 || mmt.matrices.size() == 1);

        bw_rank_check(mmt, pl);

        set_common_seed();

        nrhs = load_and_prepare_rhs_vectors(pl);

        Z = zero_columns();
    }

    void set_common_seed()
    { // {{{
        if (pi->m->trank == 0 && !bw->seed) {
            /* note that bw is shared between threads, thus only thread 0 should
             * test and update it here.
             * at pi->m->jrank > 0, we don't care about the seed anyway
             */
            bw->seed = int(time(nullptr));
            MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
        }
        serialize_threads(pi->m);
        gmp_randseed_ui(rstate, bw->seed);

        int const tcan_print = bw->can_print && pi->m->trank == 0;

        if (tcan_print)
            fmt::print("// Random generator seeded with {}\n", bw->seed);
    } // }}}

    /* this function reads the text file that gives the rhs vector, and
     * transforms it into an mmt vector in the correct vector space (or at
     * least we hope so).
     *
     * TODO: read binary data as well in order to do binary systems too.
     *
     * TODO: we must probably add a layer of twisting / untwisting here.
     *
     */
    unsigned int load_and_prepare_rhs_vectors(cxx_param_list & pl)
    { // {{{
        parallelizing_info_ptr pi = mmt.pi;
        arith_generic * A = mmt.abase;

        /* First create all RHS vectors -- these are just splits of the big
         * RHS block. Those files get created together. */
        char const * rhs_name = param_list_lookup_string(pl, "rhs");

        if (!rhs_name)
            return 0;

        rhs_header hdr;

        // open rhs only at job 0, thread 0
        std::unique_ptr<std::ifstream> rhs;
        if (pi->m->jrank == 0 && pi->m->trank == 0) {
            rhs = std::make_unique<std::ifstream>(rhs_name);
            (*rhs) >> hdr;
            ASSERT_ALWAYS(bool(*rhs));
            if (char2) {
                /* as it turns out, it doesn't take a whole lot to use
                 * read_rhs_from_file in both cases (I think).
                 */
                (*rhs) >> std::hex;
            }
        }
        // all threads at job 0 get access to it (even though only
        // threads along one communicator will deal with it)
        // other mpi ranks will just be sharing a null pointer, here.
        auto * p_rhs = rhs.get();
        pi_thread_bcast(&p_rhs, sizeof(void *), BWC_PI_BYTE, 0, pi->m);
        pi_bcast(&hdr.ncols, sizeof(hdr.ncols), BWC_PI_BYTE, 0, 0, pi->m);
        pi_bcast(&hdr.nrows, sizeof(hdr.nrows), BWC_PI_BYTE, 0, 0, pi->m);
        /* don't share the cxx_mpz hdr.p */

        unsigned int const nrhs = hdr.ncols;

        ASSERT_ALWAYS(nrhs);
        ASSERT_ALWAYS(nrhs <= mmt.n[!bw->dir]);

        abase_proxy const natural = abase_proxy::most_natural(pi);
        arith_generic * Av = natural.A.get();
        pi_datatype_ptr Av_pi = natural.A_pi;

        /* XXX I'm pretty sure that this assert holds, in which case Av can
         * be A == mmt.abase, and Av_pi can be mmt.pitype
         */
        ASSERT_ALWAYS(Av->simd_groupsize() == A->simd_groupsize());

        std::vector<mmt_vec> rhs_vecs;

        /* the rhs vectors have mmt.n0[!bw->dir] coordinates, and they
         * must be considered as vectors in the direction !bw->dir if we
         * want to get a chance to apply the T (decorrelating)
         * permutation on them.
         *
         * This being said, for the purposes of serving as init vectors
         * for the iteration, we obviously want to store them as vectors
         * with respect to the larger dimension unpadded =
         * max(mmt->n0[*]).
         */

        for (unsigned int j = 0; j * splitwidth < nrhs; j++) {
            rhs_vecs.emplace_back(mmt, Av, Av_pi, !bw->dir, /* shared ! */ 1,
                                  mmt.n[!bw->dir]);

            mmt_full_vec_set_zero(rhs_vecs.back());

            ASSERT_ALWAYS(mmt.n0[bw->dir] <= mmt.n[bw->dir]);
        }

        ASSERT_ALWAYS(hdr.nrows == mmt.n0[!bw->dir]);

        read_rhs_from_file(rhs_vecs, p_rhs, mmt.n0[!bw->dir]);

        for (unsigned int j = 0; j * splitwidth < nrhs; j++) {
            /* Twist the rhs vectors! It's important for inhomogenous systems
             * with nullspace=left
             * XXX true ? twist ? untwist ? Always ?
             */
            mmt_vec_apply_T(mmt, rhs_vecs[j]);

            auto pat = bwc_V_file::pattern(0);

            std::string const name = fmt::format(fmt::runtime(pat),
                    j * splitwidth, (j + 1) * splitwidth);

            if (tcan_print)
                fmt::print("// Creating {} (extraction from {})\n", name,
                           rhs_name);

            /* although the RHS vectors have mmt->n0[!bw->dir], we save
             * then wrt the unpadded dimension, which is max(mmt->n0[*])
             */
            mmt_vec_save(rhs_vecs[j], pat, unpadded, j);
        }

        return nrhs;
    } // }}}

    void complete_with_random_y_vectors()
    { // {{{
        /* Create purely random vectors for V */
        size_t const splitwidth = y.abase->simd_groupsize();
        parallelizing_info_ptr pi = y.pi;
        int const tcan_print = bw->can_print && pi->m->trank == 0;

        for (unsigned int j = nrhs; j < (unsigned int)bw->n; j += splitwidth) {
            auto pat = bwc_V_file::pattern(0);
            mmt_vec_set_random_through_file(y, pat, unpadded, rstate, j);
            if (tcan_print)
                fmt::print("// generated {}\n",
                        fmt::format(fmt::runtime(pat), j, j + splitwidth));
        }
    } // }}}

    void find_init_vectors()
    { // {{{
        /*
         * XXX beyond the zero_columns fix that is coded above, there's
         * also another way in which we may fail. We don't want any of
         * the x vectors to be in a small degree characteristic subspace
         * of M. How we want to guard against that is not entirely clear.
         *
         * The q&d approach (which, by the way, was the one that was in
         * place before we reimplemented prep) is to set my_nx to a
         * not-too-small value, and hope for the best.
         *
         * A better way would probably be to verify the absence of this
         * phenomenon up to a certain bounded degree.
         */
        unsigned int my_nx = std::max((int)iceildiv(Z.size(), bw->m), 4);
        if (my_nx > 4 && tcan_print) {
            fmt::print("Because of zero columns, we need"
                       " at least {} coordinates per x vector\n",
                       my_nx);
        }

        for (unsigned ntri = 0;; ntri++) {
            if (nrhs == (unsigned int)bw->n) {
                if (ntri)
                    ++my_nx;
                if (ntri >= 4) {
                    fmt::print(stderr,
                               "Cannot find a satisfactory initialization, "
                               "and your RHS vectors leave with no leeway for "
                               "randomness "
                               "(nrhs={}, n={}). "
                               "Maybe your RHS vectors are bad ?\n",
                               nrhs, bw->n);
                    exit(EXIT_FAILURE);
                }
            } else if (ntri >= my_nx * 10) {
                ++my_nx;
                if (tcan_print)
                    fmt::print("// Getting bored. Trying {} x vectors\n", my_nx);
            }
            serialize_threads(pi->m);

            if (tcan_print)
                fmt::print(
                    "// Generating new vectors for x and y[{}:] (trial # {})\n",
                    nrhs, ntri);

            unsigned int const rk = do_one_trial(my_nx);

            if (tcan_print)
                fmt::print("// Dimension of kernel: {}\n", bw->m - rk);

            if (rk == (unsigned int)bw->m) {
                if (tcan_print)
                    fmt::print("// Found good x,y vector pair after {} trials\n",
                           ntri + 1);
                save_x(xvecs, bw->m, my_nx, pi);
                return;
            }
        }
    } // }}}

    arith_generic::owned_vector
    compute_first_few_matrices(unsigned int nx, unsigned int iterations)
    { // {{{
        auto xymats = A->alloc_vector(bw->m * iterations * A_multiplex,
                                      ALIGNMENT_ON_ALL_BWC_VECTORS);

        // we have indices mmt.wr[1]->i0..i1 available.
        A->vec_set_zero(xymats.get(), bw->m * iterations * A_multiplex);

        // compute x^T M y, x^T M^2 y, and so on. Since we do that
        // piecewise for the different vectors, we first collect
        // everything in the xymats array, and compute the rank later on.

        for (size_t j = 0; j < A_multiplex; j++) {
            auto pat = bwc_V_file::pattern(0);
            int const ok = mmt_vec_load(y, pat, unpadded, j * splitwidth);
            ASSERT_ALWAYS(ok);

            /*
            // XXX Note that x^Ty does not count here, because it does not
            // take part to the sequence computed by lingen !
            //
            // XXX -> it's no longer true if we have a RHS.
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy.vectors(), nullptr);
            mmt_vec_untwist(mmt, y);
            */

            for (unsigned int k = 0; k < iterations; k++) {
                x_dotprod(
                    A->vec_subvec(xymats.get(), (k * A_multiplex + j) * bw->m),
                    xvecs, 0, bw->m, nx, y, 1);
                mmt_vec_twist(mmt, y);
                matmul_top_mul(mmt, ymy.vectors(), nullptr);
                mmt_vec_untwist(mmt, y);
            }
        }

        return xymats;
    } // }}}

    std::vector<unsigned int> zero_columns()
    {
        /* XXX There's the following catch: if nullspace==right, we're
         * going to build a vector that is orthogonal to x, x*M, x*M^2,
         * and so on. If our matrix has zero columns, which may well
         * happen in tests, the only vector among these that may have
         * non-zero stuff at these coordinates is x.
         *
         * This use case does fit in a possible scenario, after all: say
         * 1000 relation-sets, 997 ideals, and 3 schirokauer maps computed
         * for each relation-set. We _should_ get a solution. And yet, we
         * don't always find it, because if x does not touch the last 3
         * coordinates, and the vector of logs and logs of SMs might give
         * us a non-zero result at the last 3 relation-sets.
         *
         * So the fix that we want to make is that if we can identify
         * zero coordinates in things that are essentially the same
         * vectors as the ones that get computed by bwc/secure, then we
         * want to start over.
         */
        mmt_vector_pair myy(mmt, !bw->dir);
        mmt_vec & my = myy[0];


        // mmt_vec_set_x_indices(my, gxvecs.get(), MIN(nchecks, bw->m), nx);
        mmt_vec_set_random_inconsistent(my, rstate);
        mmt_vec_twist(mmt, my);
        matmul_top_mul(mmt, myy.vectors(), nullptr);
        mmt_vec_untwist(mmt, my);

        /* because the vector is fully consistent, there's nothing to
         * fear: everyone will compute the same data.
         */
        std::set<unsigned int> ret0;
        for (unsigned int i = mmt.n0[bw->dir]; i < unpadded; i++)
            ret0.insert(i);
        for (unsigned int i = my.i0; i < my.i1 && i < mmt.n0[bw->dir]; i++) {
            if (A->is_zero(A->vec_item(my.v, i - my.i0)))
                ret0.insert(i);
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN(my.pi->wr[my.d]) {
            parallelizing_info_experimental::allgather(ret0, my.pi->wr[my.d]->xwr);
        }
        SEVERAL_THREADS_PLAY_MPI_END();

        /* do it a second time: this gives us a chance to remove some
         * zero coordinates from our list */
        for ( ; ; ) {
            mmt_vec_twist(mmt, my);
            matmul_top_mul(mmt, myy.vectors(), nullptr);
            mmt_vec_untwist(mmt, my);

            std::set<unsigned int> ret1;
            for (unsigned int i = mmt.n0[bw->dir]; i < unpadded; i++)
                ret1.insert(i);
            for (unsigned int i = my.i0; i < my.i1 && i < mmt.n0[bw->dir]; i++) {
                if (A->is_zero(A->vec_item(my.v, i - my.i0)))
                    ret1.insert(i);
            }
            SEVERAL_THREADS_PLAY_MPI_BEGIN(my.pi->wr[my.d]) {
                parallelizing_info_experimental::allgather(ret1, my.pi->wr[my.d]->xwr);
            }
            SEVERAL_THREADS_PLAY_MPI_END();

            std::set<unsigned int> early_nz;
            for(auto const i : ret1) {
                if (ret0.find(i) == ret0.end())
                    /* i wasn't zero in an earlier iteration */
                    early_nz.insert(i);
            }

            for(auto const i : early_nz)
                ret1.erase(ret1.find(i));
            /* at this point, ret1 is a subset of ret0. We'll keep
             * working on it if it is a struct subset */
            if (ret1.size() == ret0.size())
                break;
            std::swap(ret1, ret0);
        }

        std::vector<unsigned int> ret(ret0.begin(), ret0.end());

        if (tcan_print && !ret.empty()) {
            std::string s;
            for (auto j: ret)
                s += fmt::format(" {}", j);
            fmt::print("found zero {}(s):{}\n", bw->dir ? "column" : "row", s);
        }

        return ret;
    }

    unsigned int do_one_trial(unsigned int nx)
    { // {{{
        if (tcan_print)
            fmt::print("// Choosing new projection vectors with"
                    " {} non-zero coordinates per x vector\n", nx);
        // generate indices w.r.t *unpadded* dimensions !
        xvecs = setup_x_random(bw->m, nx, unpadded, pi, rstate, Z);

        complete_with_random_y_vectors();

        /* we now want to make sure that this choice of x and y makes the
         * lingen initialization feasible
         */

        /* Minimal number of m by n matrices to use in order to obtain a
         * matrix of rank m.
         *
         * Note that it must be at least m/n, otherwise we stand no chance !
         *
         * The +1 is here if we have a RHS.
         */
        unsigned int const prep_iterations = iceildiv(bw->m, bw->n);

        auto xymats = compute_first_few_matrices(nx, prep_iterations + 1);

        /* now compute the rank of xymats. It's a bit tricky because we
         * laid it out in striped format. So the first thing we have to
         * do is copy it.
         */
        auto xymats_not_striped =
            A->alloc_vector(bw->m * prep_iterations * A_multiplex,
                            ALIGNMENT_ON_ALL_BWC_VECTORS);

        for (unsigned int k = 0; k < (unsigned int)iceildiv(bw->m, bw->n);
             k++) {
            ASSERT_ALWAYS(A->simd_groupsize() == (size_t)splitwidth);
            for (size_t j = 0; j < A_multiplex; j++) {
                unsigned int k1 = k;
                k1 += j * splitwidth >= nrhs;
                for (int r = 0; r < bw->m; r++) {
                    memcpy(A->vec_subvec(
                               xymats_not_striped.get(),
                               (r * prep_iterations + k) * A_multiplex + j),
                           A->vec_subvec(xymats.get(),
                                         (k1 * A_multiplex + j) * bw->m + r),
                           A->elt_stride());
                }
            }
        }

        /* Make sure computation is over for everyone ! */
        serialize_threads(pi->m);

        /* Now all threads and jobs must collectively reduce the zone
         * pointed to by xymats */
        pi_allreduce(nullptr, xymats_not_striped.get(),
                     bw->m * prep_iterations * A_multiplex, mmt.pitype,
                     BWC_PI_SUM, pi->m);

        /* OK -- now everybody has the same data */

        return matrix_rank(A.get(), xymats_not_striped.get(), bw->m,
                           prep_iterations * bw->n);
    } // }}}
};

static void * prep_prog(parallelizing_info_ptr pi, cxx_param_list & pl,
                 void * arg MAYBE_UNUSED)
{
    prep_object P(pi, pl);
    P.find_init_vectors();
    return nullptr;
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches */
    param_list_decl_usage(pl, "rhs",
                          "file with the right-hand side vectors for "
                          "inhomogeneous systems mod p");

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    /* This program does not support interleaving */
    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters: none here (so far). */
    if (mpz_cmp_ui(bw->p, 2) != 0)
        param_list_lookup_string(pl, "rhs");

    ASSERT_ALWAYS(!param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(!param_list_lookup_string(pl, "solutions"));

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank)
            param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(prep_prog, pl, nullptr);

    parallelizing_info_finish();
    bw_common_clear(bw);

    return 0;
}
