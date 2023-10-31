#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <ctime>                // for time
#include <string>                // for string
#include <gmp.h>                 // for gmp_randclear, gmp_randinit_default
#include "matmul.hpp"              // for matmul_public_s
#include "parallelizing_info.hpp"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "select_mpi.h"
#include "params.h"
#include "misc.h"
#include "bw-common.h"
#include "async.hpp"
#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "fmt/format.h"
#include "macros.h"
#include "mmt_vector_pair.hpp"
using namespace fmt::literals;

void * mksol_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;
    if (fake) bw->skip_online_checks = 1;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    struct timing_data timing[1];

    unsigned int solutions[2] = { bw->solutions[0], bw->solutions[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->solutions[1]-bw->solutions[0]) % 2 == 0);
        solutions[0] = bw->solutions[0] + pi->interleaved->idx * (bw->solutions[1]-bw->solutions[0])/2;
        solutions[1] = solutions[0] + (bw->solutions[1]-bw->solutions[0])/2;
    }

    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    int splitwidth = char2 ? 64 : 1;

    /* Define and initialize our arithmetic back-ends. Because simd group
     * size differs, we have two distinct backends to create. One for the
     * vectors iterates we read from disk (width being splitwidth),
     * and one for the matrix times vector multiplications we perform
     * (width being solutions[1]-solutions[0]).
     */

    /* {{{ First: only relative to the vectors we read */
    unsigned int Av_width = splitwidth;
    unsigned int Av_multiplex = bw->n / Av_width;
    std::unique_ptr<arith_generic> Av(arith_generic::instance(bw->p, Av_width));
    pi_datatype_ptr Av_pi = pi_alloc_arith_datatype(pi, Av.get());
    /* }}} */

    /* {{{ Second: We intend to perform only a single spmv per iteration,
     * which constrains solutions[1]-solutions[0] to being a type width we can
     * handle profitably.
     */
    unsigned int As_multiplex = 1;
    unsigned int As_width = solutions[1]-solutions[0];
    if ((char2 && (As_width != 64 && As_width != 128 && As_width != 256))
            || (!char2 && As_width > 1))
    {
        fprintf(stderr,
                "We cannot support computing %u solutions at a time "
                "with one single Spmv operation, given the currently "
                "implemented code\n",
                As_width);
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<arith_generic> As(arith_generic::instance(bw->p, As_width));
    /* How many F files do we need to read simultaneously to form
     * solutions[1]-solutions[0] columns ? */
    unsigned int Af_multiplex = As_width / Av_width;
    /* }}} */

    /* {{{ ... and the combined operations */
    std::unique_ptr<arith_cross_generic> AvxAs(arith_cross_generic::instance(Av.get(), As.get()));
    /* }}} */

    block_control_signals();

    /* Now that we do this in Horner fashion, we multiply on vectors
     * whose width is the number of solutions we compute. */
    matmul_top_data mmt(As.get(), pi, pl, bw->dir);
    pi_datatype_ptr As_pi = mmt.pitype;

    /* allocate vectors (two batches): */

    mmt_vector_pair ymy(mmt, bw->dir);

    /* {{{ For the vectors which we read from disk and which participate in
     *   the coefficients which get added to the computation at each
     *   iteration: we need n vectors -- or n/64 for the binary case.
     */
    mmt_vec * vi = new mmt_vec[bw->n / splitwidth];
    matmul_top_matrix * mptr = &mmt.matrices[bw->dir ? (mmt.matrices.size() - 1) : 0];
    for(int i = 0 ; i < bw->n / splitwidth ; i++) {
        mmt_vec_setup(vi[i], mmt, Av.get(), Av_pi, bw->dir, 1, mptr->n[bw->dir]);
        mmt_full_vec_set_zero(vi[i]);
    }
    /* }}} */


    unsigned int unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    serialize(pi->m);
    
    /* {{{ i/o stats.
     * let's be generous with interleaving protection. I don't want to be
     * bothered, really */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);
    matmul_top_comm_bench(mmt, bw->dir);
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);
    /* }}} */

    /* {{{ Read all vi's */
    gmp_randstate_t rstate;
    if (fake) {
        gmp_randinit_default(rstate);
        if (pi->m->trank == 0 && !bw->seed) {
            bw->seed = time(NULL);
            MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
        }
        serialize_threads(pi->m);
        gmp_randseed_ui(rstate, bw->seed);
        if (tcan_print) 
            printf("// Random generator seeded with %d\n", bw->seed);
    }
    /* }}} */

    unsigned int expected_last_iteration;

    serialize_threads(pi->m);
    if (pi->m->trank == 0) {
        /* the bw object is global ! */
        expected_last_iteration = bw_set_length_and_interval_mksol(bw, mmt.n0);
    }
    pi_thread_bcast(&expected_last_iteration, 1, BWC_PI_UNSIGNED, 0, pi->m);
    serialize_threads(pi->m);
    if (bw->end == INT_MAX) {
        if (tcan_print)
            printf ("Target iteration is unspecified ;"
                    " going to end of F file\n");
    } else {
        if (tcan_print)
            printf ("Target iteration is %u\n", bw->end);
        expected_last_iteration = bw->end;
    }
    ASSERT_ALWAYS(bw->end == INT_MAX || bw->end % bw->interval == 0);

    pi_interleaving_flip(pi);
    if (bw->checkpoint_precious) {
        if (tcan_print) {
            printf("As per interval=%d checkpoint_precious=%d, we'll load vectors every %d iterations, and print timings every %d iterations\n",
                    bw->interval,
                    bw->checkpoint_precious,
                    bw->checkpoint_precious,
                    bw->interval);
        }
    } else {
        bw->checkpoint_precious = bw->interval;
    }
    pi_interleaving_flip(pi);

    /* {{{ Prepare temp space for F coefficients */
    /* F plays the role of a right-hand-side in mksol. Vector iterates
     * have to be multiplied on the right by a matrix corresponding to
     * some coefficient of F. Because F has been written to disk in the
     * "proper" order (namely: reversed order, leading coeff first) by
     * lingen, iterate i is to be multiplied by the i-th coefficient of F.
     *
     * The F-coefficients, from the overall description of the algorithm,
     * are square matrices having bw->n rows and columns.  Columns
     * correspond to candidate solutions. In some cases we are interested
     * in only a few solutions, maybe even one for the prime field case.
     * More precisely, we are intereseted in columns whose indices are
     * given by the interval solutions[0]..solutions[1].
     *
     * The bw->n rows correspond to the fact that we have
     * bw->n/splitwidth vectors: at each iteration, we will perform that
     * many multiplications by coefficients of F.
     */

    // XXX remove ?
    // ASSERT(Av->vec_elt_stride(Av,As_width) == As->vec_elt_stride(As, Av->simd_groupsize(Av)));

    /* We'll load all the F coefficient matrices before the main loop
     * (and for a full interval).
     * It's mostly dual to the treatment of the xy matrices in krylov.
     * Same considerations apply w.r.t the memory footprint.
     *
     * We have Av_multiplex "sets of rows", and As_multiplex "sets of columns".
     *
     * This is complicated further by the facts that the "sets of
     * columns" above (on which we do spmv) are actually divided into
     * several smaller sets. The number of these subsets is Af_multiplex,
     * and each holds splitwidth columns.
     */
    size_t one_fcoeff = As->vec_elt_stride(Av->simd_groupsize());
    if (tcan_print) {
        char buf[20];
        printf("Each thread allocates %d*%u*%u*%zu*%d=%s for the F matrices\n",
                Af_multiplex > 1 ? 2 : 1,
                Av_multiplex,
                As_multiplex,
                one_fcoeff,
                bw->interval,
                size_disp(
                    (Af_multiplex > 1 ? 2 : 1) *
                    Av_multiplex *
                    As_multiplex *
                    one_fcoeff *
                    bw->interval, buf));
    }
    arith_generic::elt ** fcoeffs = new arith_generic::elt*[Av_multiplex * As_multiplex];
    for(unsigned int k = 0 ; k < Av_multiplex * As_multiplex ; k++) {
        (fcoeffs[k]) = As->alloc(one_fcoeff * bw->interval, ALIGNMENT_ON_ALL_BWC_VECTORS);
        As->vec_set_zero(fcoeffs[k], one_fcoeff * bw->interval);
    }
    ASSERT_ALWAYS(Av_width * Af_multiplex == (unsigned int) As->simd_groupsize());
    arith_generic::elt * fcoeff_tmp = NULL;
    if (Af_multiplex > 1) {
            (fcoeff_tmp) = As->alloc(one_fcoeff / Af_multiplex * bw->interval, ALIGNMENT_ON_ALL_BWC_VECTORS);
            As->vec_set_zero(fcoeff_tmp, one_fcoeff / Af_multiplex * bw->interval);
    }
    /* }}} */
    
    /* {{{ bless our timers */
    timing_init(timing, 4 * mmt.matrices.size(), bw->start, expected_last_iteration);
    for(size_t i = 0 ; i < mmt.matrices.size(); i++) {
        timing_set_timer_name(timing, 4*i, "CPU%zu", i);
        timing_set_timer_items(timing, 4*i, mmt.matrices[i].mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%zu", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%zu", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%zu", i);
    }
    /* }}} */

    int bw_end_copy = bw->end; /* avoid race conditions w/ interleaving */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw_end_copy ; s += bw->checkpoint_precious ) {
        serialize(pi->m);
        for(int i = 0 ; i < bw->n / splitwidth ; i++) {
            int ys[2] = { i * splitwidth, (i + 1) * splitwidth };
            std::string v_name = fmt::format(FMT_STRING("V%u-%u.{}"), s);
            if (fake) {
                mmt_vec_set_random_through_file(vi[i], v_name, unpadded, rstate, ys[0]);
            } else {
                int ok = mmt_vec_load(vi[i], v_name, unpadded, ys[0]);
                ASSERT_ALWAYS(ok);
            }
            mmt_vec_twist(mmt, vi[i]);
        }

        serialize(pi->m);

        mmt_full_vec_set_zero(ymy[0]);

        int sx = MIN(bw_end_copy, s + bw->checkpoint_precious);
        if (tcan_print) {
            /*
            printf("// bw->start=%d bw_end_copy=%d sx=%d s=%d\n",
                    bw->start,
                    bw_end_copy,
                    sx, s);
                    */
            printf("about to do %d iterations starting from vectors at iteration %d to handle the coefficients of degree [%d..%d] in F\n", sx-s-1, s, s, sx-1);
        }

        /* read coefficients of F by windows */
        int n_windows = bw->checkpoint_precious / bw->interval;

        for(int i_window = 0 ; i_window < n_windows ; i_window++) {
            /* We'll read from coefficient s0 */
            int s0 = s + (n_windows - 1 - i_window) * bw->interval;
            int s1 = s + (n_windows     - i_window) * bw->interval;
            if (s0 >= bw_end_copy)
                continue;

            /* This is used only on the leader node */
            bool short_read = false;

            serialize(pi->m);
            if (pi->m->trank == 0 && pi->m->jrank == 0) {
                int rc0 = 0, rc = 0;
                for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                    for(unsigned int j = 0 ; j < As_multiplex ; j++) {
                        arith_generic::elt * ff = fcoeffs[j * Av_multiplex + i];
                        /* points to bw->interval * one_fcoeff */
                        /* Zero out first, then read; when we reach EOF while
                         * reading F, it is crucially important that we have
                         * a zero area past the end of the file ! */
                        As->vec_set_zero(ff, one_fcoeff * bw->interval);

                        /* Now read piece by piece */
                        for(unsigned int k = 0 ; k < Af_multiplex ; k++) {
                            arith_generic::elt * buffer;
                            if (Af_multiplex == 1) {
                                buffer = ff;
                            } else {
                                buffer = fcoeff_tmp;
                            }

                            unsigned int sol0, sol1;
                            sol0 = (j * Af_multiplex + k) * Av_width;
                            sol0 += solutions[0];
                            sol1 = sol0 + Av_width;
                            std::string f_name = fmt::format(FMT_STRING("F.sols{}-{}.{}-{}"), 
                                    sol0, sol1,
                                    i * Av_width, (i + 1) * Av_width);

                            // printf("[%d] reading from %s\n", pi->interleaved ? pi->interleaved->idx : -1, tmp);
                            FILE * f = fopen(f_name.c_str(), "rb");
                            DIE_ERRNO_DIAG(f == NULL, "fopen(%s)", f_name.c_str());
                            rc = fseek(f, one_fcoeff / Af_multiplex * s0, SEEK_SET);
                            if (rc >= 0) {
                                /* Read everything in one go. We might want to
                                 * reconsider how the coefficients inside the
                                 * individual F files are written -- transposed or
                                 * not. Of course that does not matter for the prime
                                 * case where individual F files are polynomials, but
                                 * surely it does in the binary case where individual
                                 * F files are a priori made of 64*64 matrices.
                                 */
                                rc = fread(buffer, one_fcoeff / Af_multiplex, bw->interval, f);
                                ASSERT_ALWAYS(rc >= 0 && rc <= bw->interval);
                                if (Af_multiplex > 1) {
                                    /* spread to fcoeff */
                                    const arith_generic::elt * src = buffer;
                                    arith_generic::elt * dst = Av->vec_subvec(ff, k);
                                    for(unsigned int row = 0 ; row < Av_width *  rc; row++) {
                                        memcpy(dst, src, Av->elt_stride());
                                        src = Av->vec_subvec(src, 1);
                                        dst = Av->vec_subvec(dst, Af_multiplex);
                                    }
                                }
                            } else {
                                /* Otherwise we're off bounds already, don't try
                                 * to read anything...
                                 */
                                rc = 0;
                            }

                            if (i == 0 && j == 0 && k == 0) rc0 = rc;

                            if (rc != rc0) {
                                fprintf(stderr, "Inconsistency in number of "
                                        "coefficients for F files\n");
                                exit(EXIT_FAILURE);
                            }
                            /* TODO: *maybe* transpose the F coefficients
                             * at this point */
                            fclose(f);
                        }
                    }
                }

                if (rc0 < bw->interval) {
                    short_read = true;
                    if (s1 != sx) {
                        fprintf(stderr, "Problem while reading coefficients of f for degrees [%d..%d[ ; we should not have a short read given that bw->end=%d\n",
                                s0, s1, bw_end_copy);
                        exit(EXIT_FAILURE);
                    }
                    sx = bw_end_copy = s0+rc;
                }
            }
            serialize(pi->m);
            pi_bcast(&s0, 1, BWC_PI_INT, 0, 0, pi->m);
            pi_bcast(&s1, 1, BWC_PI_INT, 0, 0, pi->m);
            pi_bcast(&sx, 1, BWC_PI_INT, 0, 0, pi->m);
            pi_bcast(&bw_end_copy, 1, BWC_PI_INT, 0, 0, pi->m);

            serialize_threads(mmt.pi->m);

            if (s0 == bw_end_copy)
                continue;

            if (bw_end_copy < s1)
                s1 = bw_end_copy;

            if (tcan_print && short_read)
                printf("We read %d coefficients from F in total\n", bw_end_copy);

            /* broadcast f */
            for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                for(unsigned int j = 0 ; j < As_multiplex ; j++) {
                    arith_generic::elt * ff = fcoeffs[i * As_multiplex + j];
                    pi_bcast(ff, Av->simd_groupsize() * (s1 - s0), As_pi, 0, 0, pi->m);
                }
            }

            serialize(pi->m);

                /* Despite the fact that the bw->end value might lead us
                 * to stop earlier, we want to stop at multiples of the checking
                 * interval. That's important, otherwise the last bunch of
                 * computations won't be checked.
                 */

                pi_interleaving_flip(pi);

                size_t eblock = mmt_my_own_size_in_items(ymy[0]);

                for(int k = 0 ; k < s1 - s0 ; k++) {

                    serialize_threads(pi->m);

                    /*
                    if (tcan_print) {
                        printf("// %d\n", bw->start + sx - s1 + k);
                        printf("v:=v+f[%d];\n", s1 - 1 - k);
                    }
                    */

                    for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                        arith_generic::elt * ff = fcoeffs[i * As_multiplex /* + j */];
                        /* or maybe transpose here instead ?? */
                        AvxAs->addmul_tiny(
                                mmt_my_own_subvec(ymy[0]),
                                mmt_my_own_subvec(vi[i]),
                                As->vec_subvec(ff,
                                    (s1 - s0 - 1 - k) * Av->simd_groupsize()),
                                eblock);
                    }
                    /* addmul_tiny degrades consistency ! */
                    ymy[0].consistency = 1;
                    mmt_vec_broadcast(ymy[0]);
                    /* we're doing something which we normally avoid: write
                     * on the next input vector. This means a race condition
                     * with the forthcoming spmv, so we have to serialize.
                     *
                     * I believe we could probably get away with this.
                     */

                    serialize_threads(pi->m);

                    /* It's equivalent to doing the multiply at the
                     * beginning of the loop with the condition
                     *  if (i_window || k)
                     * except that by doing it here, we'll avoid nasty
                     * surprises with displayed timings.
                     *
                     * (and also it might well be more correct like this
                     * at the end of the computation, when we begin
                     * halway through a range of coefficients).
                     */
                    if (i_window < n_windows-1 || k < s1 - s0 - 1) {
                        // if (tcan_print) printf("v:=M*v;\n");
                        matmul_top_mul(mmt, ymy.vectors(), timing);

                        timing_check(pi, timing, s + sx - s1 + k + 1, tcan_print);
                    }
                }

                serialize(pi->m);
                /* See remark above. */
                pi_interleaving_flip(pi);
                pi_interleaving_flip(pi);

            // reached s + bw->interval. Count our time on cpu, and compute the sum.
            timing_disp_collective_oneline(pi, timing, s + sx - s0, tcan_print, "mksol");
        }

        mmt_vec_untwist(mmt, ymy[0]);

        if (!fake) {
            /* We have only one (block of) vectors at a time, so j=0,
             * really (and As_multiplex == 1)
             */
            int j = 0;
            std::string s_name = fmt::format(FMT_STRING("S.sols%u-%u.{}-{}"), s, s + bw->checkpoint_precious);
            ASSERT_ALWAYS(ymy[0].abase->simd_groupsize() == As_width);
            mmt_vec_save(ymy[0], s_name, unpadded,
                    solutions[0] + j * As_width);
        }
    }

    timing->end_mark = bw->start + bw->interval * iceildiv(bw_end_copy - bw->start, bw->interval);

    timing_final_tally(pi, timing, tcan_print, "mksol");

    if (tcan_print) {
        printf("Done mksol.\n");
    }
    serialize(pi->m);
    if (fake) gmp_randclear(rstate);
    delete[] vi;

    for(unsigned int k = 0 ; k < Av_multiplex * As_multiplex ; k++) {
        As->free((fcoeffs[k]));
    }
    delete[] fcoeffs;

    if (Af_multiplex > 1) {
        As->free((fcoeff_tmp));
    }

    pi_free_arith_datatype(pi, Av_pi);

    timing_clear(timing);

    return NULL;
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters: none here (so far). */

    ASSERT_ALWAYS(!param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(param_list_lookup_string(pl, "solutions"));

    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(mksol_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}


