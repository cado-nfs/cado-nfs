#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <memory>
#include <string>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "async.hpp"
#include "bw-common.h"
#include "bwc_filenames.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"
#include "mmt_vector_pair.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "rolling.hpp"
#include "select_mpi.h"
#include "utils_cxx.hpp"
#include "xdotprod.hpp"
#include "xvectors.hpp"

using namespace fmt::literals;

struct check_data {
    matmul_top_data & mmt;
    parallelizing_info_ptr pi;
    int nchecks;
    arith_generic * A;
    std::unique_ptr<arith_generic> Ac;
    pi_datatype_ptr Ac_pi;
    std::unique_ptr<arith_cross_generic> AxAc;
    mmt_vec check_vector;
    arith_generic::elt * Tdata = nullptr;
    arith_generic::elt * ahead = nullptr;

    int tcan_print = 0;

    bool leader() const {
        return pi->m->trank == 0 && pi->m->jrank == 0;
    }

    check_data(matmul_top_data & mmt, arith_generic * A)
        : mmt(mmt)
        , pi(mmt.pi)
        , nchecks(mpz_cmp_ui(bw->p, 2) > 0 ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2)
        , A(A)
        , Ac(arith_generic::instance(bw->p, nchecks))
        , Ac_pi(pi_alloc_arith_datatype(pi, Ac.get()))
        , AxAc(arith_cross_generic::instance(A, Ac.get()))
        , check_vector(mmt, Ac.get(), Ac_pi,
                  bw->dir, THREAD_SHARED_VECTOR, mmt.n[bw->dir])
        , tcan_print(bw->can_print && pi->m->trank == 0)
      {
      }

    void load() {
        /* We do the dot product by working on the local vector chunks.
         * Therefore, we must really understand the check vector as
         * playing a role in the very same direction of the y vector!
         */
        std::string const Cv_filename_pattern = bwc_Cv_file::pattern(bw->interval);
        const int ok = mmt_vec_load(check_vector, Cv_filename_pattern, mmt.n0[bw->dir], 0);
        if (!ok) {
            if (tcan_print)
                fmt::print(stderr, "check file {} not found", Cv_filename_pattern);
            pi_abort(EXIT_FAILURE, pi->m);
        }
        std::string const Ct_filename = fmt::format("Ct0-{}.0-{}", nchecks, bw->m);
        Tdata = Ac->alloc(bw->m, ALIGNMENT_ON_ALL_BWC_VECTORS);
        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            auto Tfile = fopen_helper(Ct_filename, "rb");
            size_t const rc = fread(Tdata, Ac->vec_elt_stride(bw->m), 1, Tfile.get());
            ASSERT_ALWAYS(rc == 1);
        }
        if (tcan_print) fmt::print("loaded {}\n", Ct_filename);
        pi_bcast(Tdata, bw->m, Ac_pi, 0, 0, pi->m);

        ahead = A->alloc(nchecks, ALIGNMENT_ON_ALL_BWC_VECTORS);
    }
    ~check_data() {
        A->free(ahead);
        Ac->free(Tdata);
        pi_free_arith_datatype(pi, Ac_pi);
    }

    void plan_ahead(mmt_vec const & y) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.

        /* Note that the check vector is always stored untwisted in
         * memory */

        /* create a matrix of size nchecks * nbys with the dot
         * product with Cv -- in the case where nbys != nchecks,
         * dealing with that data will require some care.
         *
         */
        A->vec_set_zero(ahead, nchecks);
        /* The syntax of ->dotprod is a bit weird. We compute
         * transpose(data-operand-0)*data-operand1, but data-operand0
         * (check_vector here) actually refers to field-operand1 (Ac
         * here).
         */
        AxAc->add_dotprod(ahead,
                mmt_my_own_subvec(y),
                mmt_my_own_subvec(check_vector),
                mmt_my_own_size_in_items(y));
    }


    bool verify(mmt_vec const & y, std::vector<uint32_t> const & gxvecs, int nx) const
    {
        /* Last dot product. This must cancel ! */
            arith_generic::elt * tmp1 = nullptr;
            tmp1 = A->alloc(nchecks, ALIGNMENT_ON_ALL_BWC_VECTORS);
            for(int c = 0 ; c < bw->m ; c += nchecks) {
                /* First zero out the matrix of size nchecks * nbys.  */
                A->vec_set_zero(tmp1, nchecks);
                x_dotprod(tmp1, gxvecs, c, c + nchecks, nx, y, -1);
                /* And now compute the product transpose(part of
                 * T)*ahead_tmp, and subtract that from our check value
                 */
                AxAc->add_dotprod(
                        ahead,
                        tmp1,
                        Ac->vec_subvec(Tdata, c),
                        nchecks);
            }
            A->free(tmp1);

        pi_allreduce(nullptr, ahead, nchecks, mmt.pitype, BWC_PI_SUM, pi->m);
        return A->vec_is_zero(ahead, nchecks);
    }
};

static void * krylov_prog(parallelizing_info_ptr pi, cxx_param_list & pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != nullptr;
    fake = fake || param_list_lookup_string(pl, "static_random_matrix") != nullptr;
    if (fake) bw->skip_online_checks = 1;
    int const tcan_print = bw->can_print && pi->m->trank == 0;
    struct timing_data timing[1];

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->ys[1]-bw->ys[0]) % 2 == 0);
        ys[0] = bw->ys[0] + pi->interleaved->idx * (bw->ys[1]-bw->ys[0])/2;
        ys[1] = ys[0] + (bw->ys[1]-bw->ys[0])/2;
    }

    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, ys[1]-ys[0]));
    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    std::shared_ptr<check_data> C;
    if (!bw->skip_online_checks)
        C = std::make_shared<check_data>(mmt, A.get());

    mmt_vector_pair ymy(mmt, bw->dir);

    unsigned int const unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    serialize(pi->m);
    
    std::vector<uint32_t> gxvecs;
    unsigned int nx = 0;
    if (!fake) {
        gxvecs = load_x(bw->m, nx, pi);
    } else {
        gxvecs = set_x_fake(bw->m, nx, pi);
    }
    indices_twist(mmt, gxvecs, bw->dir);

    /* let's be generous with interleaving protection. I don't want to be
     * bothered, really */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);
    matmul_top_comm_bench(mmt, bw->dir);
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    /* I have absolutely no idea why, but the two --apparently useless--
     * serializing calls around the next block seem to have a beneficial
     * impact on the SEGv's we see every now and then with --mca
     * mpi_leave_pinned 1
     */
    serialize(pi->m);
    if (!fake) {
        int const ok = mmt_vec_load(ymy[0], bwc_V_file::pattern(bw->start), unpadded, ys[0]);
        ASSERT_ALWAYS(ok);
    } else {
        cxx_gmp_randstate rstate;
#if 0
        /* This is for setting the source vector to something consistent
         * across mappings, so that given a fixed (fake, here) matrix, any
         * splitting will give the same source vector. This is just a
         * copy of the mechanism which exists in prep for doing exactly
         * this. Alas, in what we denote as a fake situation, there is
         * no chance of course that two different splittings lesad to
         * identical matrices ! Hence, we'd rather not bother with
         * generating something consistent.
         */
        if (pi->m->trank == 0 && !bw->seed) {
            bw->seed = time(nullptr);
            MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
        }
        serialize_threads(pi->m);
        gmp_randseed_ui(rstate, bw->seed);
        if (tcan_print) {
            printf("// Random generator seeded with %d\n", bw->seed);
        }
        if (tcan_print) { printf("Creating fake %s...", v_name); fflush(stdout); }
        mmt_vec_set_random_through_file(mmt, nullptr, v_name, bw->dir, bw->start, unpadded, rstate);
        if (tcan_print) { printf("done\n"); }
#else
        unsigned long const g = pi->m->jrank * pi->m->ncores + pi->m->trank;
        gmp_randseed_ui(rstate, bw->seed + g);
        mmt_vec_set_random_inconsistent(ymy[0], rstate);
        mmt_vec_truncate(mmt, ymy[0]);
#endif
    }
    serialize(pi->m);

    serialize_threads(pi->m);
    if (pi->m->trank == 0) {
        /* the bw object is global ! */
        bw_set_length_and_interval_krylov(bw, mmt.n0);
    }
    serialize_threads(pi->m);
    if (tcan_print) {
        fmt::print("Target iteration is {}\n", bw->end);
    }
    ASSERT_ALWAYS(bw->end % bw->interval == 0);


    if (C) C->load();

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    arith_generic::elt * xymats;

    if (tcan_print) {
        fmt::print("Each thread allocates {} kb for the A matrices\n",
                A->vec_elt_stride(bw->m*bw->interval) >> 10);
    }
    xymats = A->alloc(bw->m*bw->interval, ALIGNMENT_ON_ALL_BWC_VECTORS);
   
#if 0
    /* FIXME -- that's temporary ! only for debugging */
    pi_log_init(pi->m);
    pi_log_init(pi->wr[0]);
    pi_log_init(pi->wr[1]);
#endif

    timing_init(timing, 4 * mmt.matrices.size(), bw->start, bw->end);
    auto clean_timing = call_dtor([&]() { timing_clear(timing); });

    for(size_t i = 0 ; i < mmt.matrices.size(); i++) {
        timing_set_timer_name(timing, 4*i, "CPU%zu", i);
        timing_set_timer_items(timing, 4*i, mmt.matrices[i].mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%zu", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%zu", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%zu", i);
    }

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {

        if (C) C->plan_ahead(ymy[0]);

        /* Create an empty slot in program execution, so that we don't
         * impose strong constraints on twist/untwist_vector being free of
         * MPI calls (well, it's *not* free of MPI calls, to start
         * with...).
         */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);
        mmt_vec_twist(mmt, ymy[0]);

        A->vec_set_zero(xymats, bw->m*bw->interval);
        serialize(pi->m);
        pi_interleaving_flip(pi);
        for(int i = 0 ; i < bw->interval ; i++) {
            /* Compute the product by x */
            x_dotprod(A->vec_subvec(xymats, i * bw->m),
                    gxvecs, 0, bw->m, nx, ymy[0], 1);

            matmul_top_mul(mmt, ymy.vectors(), timing);

            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        /* See remark above. */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);

        if (C && !C->verify(ymy[0], gxvecs, nx)) {
            fmt::print("Failed check at iteration {}\n",
                    s + bw->interval);
            exit(1);
        }


        mmt_vec_untwist(mmt, ymy[0]);

        /* Now (and only now) collect the xy matrices */
        pi_allreduce(nullptr, xymats,
                bw->m * bw->interval,
                mmt.pitype, BWC_PI_SUM, pi->m);

        if (pi->m->trank == 0 && pi->m->jrank == 0 && !fake) {
            std::string const tmp = fmt::format("A{}-{}.{}-{}", ys[0], ys[1], s, s+bw->interval);
            std::string const tmptmp = tmp + ".tmp";
            FILE * f = fopen(tmptmp.c_str(), "wb");
            int rc = fwrite(xymats, A->elt_stride(), bw->m*bw->interval, f);
            fclose(f);
            if (rc != bw->m*bw->interval) {
                fmt::print(stderr, "Ayee -- short write\n");
                // make sure our input data won't be deleted -- this
                // chunk will have to be redone later, maybe the disk
                // failure is temporary (?)
                // 
                // anyway, better not rename the file immediately in this
                // case.
            } else {
                rc = rename(tmptmp.c_str(), tmp.c_str());
                ASSERT_ALWAYS(rc == 0);
            }
        }

        if (!fake) {
            mmt_vec_save(ymy[0], bwc_V_file::pattern(s + bw->interval), unpadded, ys[0]);
        }

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            keep_rolling_checkpoints(fmt::format("V{}-{}", ys[0], ys[1]), s + bw->interval);
        }

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "krylov");
    }

    timing_final_tally(pi, timing, tcan_print, "krylov");

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (tcan_print)
        fmt::print("Done krylov.\n");
    serialize(pi->m);

    A->free(xymats);

    int want_full_report = 0;
    param_list_parse_int(pl, "full_report", &want_full_report);
    matmul_top_report(mmt, 1.0, want_full_report);

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
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */
    if (bw->ys[0] < 0) { fmt::print(stderr, "no ys value set\n"); exit(1); }

    ASSERT_ALWAYS(param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(!param_list_lookup_string(pl, "solutions"));

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(krylov_prog, pl, 0);

    parallelizing_info_finish();

    bw_common_clear(bw);

    return 0;
}

