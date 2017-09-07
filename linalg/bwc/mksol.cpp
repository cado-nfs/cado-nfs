#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "async.h"
#include "xdotprod.h"
#include "rolling.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"

void * mksol_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;
    if (fake) bw->skip_online_checks = 1;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
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
    mpfq_vbase Av;
    unsigned int Av_width = splitwidth;
    unsigned int Av_multiplex = bw->n / Av_width;
    mpfq_vbase_oo_field_init_byfeatures(Av, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, Av_width,
            MPFQ_DONE);
    pi_datatype_ptr Av_pi = pi_alloc_mpfq_datatype(pi, Av);
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
    mpfq_vbase As;
    mpfq_vbase_oo_field_init_byfeatures(As,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, As_width,
            MPFQ_DONE);
    /* }}} */

    /* {{{ ... and the combined operations */
    mpfq_vbase_tmpl AvxAs;
    mpfq_vbase_oo_init_templates(AvxAs, Av, As);
    /* }}} */

    block_control_signals();

    /* Now that we do this in Horner fashion, we multiply on vectors
     * whose width is the number of solutions we compute. */
    matmul_top_init(mmt, As, pi, pl, bw->dir);
    pi_datatype_ptr As_pi = mmt->pitype;

    /* allocate vectors (two batches): */

    /* {{{ For the vectors which get multiplied by matrices:
     *   we need nmatrices vectors, rounded to the next even integer
     *   (well, in truth two would always suffice, but the code isn't
     *   ther yet).
     *
     * XXX note that As_multiplex is equal to 1 above
     */
    int nmats_odd = mmt->nmatrices & 1;
    mmt_vec * ymy = new mmt_vec[mmt->nmatrices + nmats_odd];
    matmul_top_matrix_ptr mptr;
    mptr = (matmul_top_matrix_ptr) mmt->matrices + (bw->dir ? (mmt->nmatrices - 1) : 0);
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        int shared = (i == 0) & nmats_odd;
        mmt_vec_init(mmt,0,0, ymy[i], bw->dir ^ (i&1), shared, mptr->n[bw->dir]);
        mmt_full_vec_set_zero(ymy[i]);

        mptr += bw->dir ? -1 : 1;
    }
    if (nmats_odd) {
        mmt_vec_init(mmt,0,0, ymy[mmt->nmatrices], !bw->dir, 0, mmt->matrices[0]->n[bw->dir]);
        mmt_full_vec_set_zero(ymy[mmt->nmatrices]);
    }
    /* }}} */

    /* {{{ For the vectors which we read from disk and which participate to
     *   the coefficients which get added to the computation at each
     *   iteration: we need n vectors -- or n/64 for the binary case.
     */
    mmt_vec * vi = new mmt_vec[bw->n / splitwidth];
    mptr = (matmul_top_matrix_ptr) mmt->matrices + (bw->dir ? (mmt->nmatrices - 1) : 0);
    for(int i = 0 ; i < bw->n / splitwidth ; i++) {
        mmt_vec_init(mmt, Av, Av_pi, vi[i], bw->dir, 1, mptr->n[bw->dir]);
        mmt_full_vec_set_zero(vi[i]);
    }
    /* }}} */


    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

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
    // ASSERT(Av->vec_elt_stride(Av,As_width) == As->vec_elt_stride(As, Av->groupsize(Av)));

    /* We'll load all the F coefficient matrices before the main loop
     * (and for a full interval).
     * It's mostly dual to the treatment of the xy matrices in krylov.
     * Same considerations apply w.r.t the memory footprint.
     *
     * We have Av_multiplex "sets of rows", and As_multiplex "sets of columns".
     */
    void ** fcoeffs = new void*[Av_multiplex * As_multiplex];
    size_t one_fcoeff = As->vec_elt_stride(As, Av->groupsize(Av));
    if (tcan_print) {
        printf("Each thread allocates %zd kb for the F matrices\n",
                (Av_multiplex *
                 As_multiplex *
                 one_fcoeff *
                 bw->interval) >> 10);
    }
    for(unsigned int k = 0 ; k < Av_multiplex * As_multiplex ; k++) {
        cheating_vec_init(As, &(fcoeffs[k]), one_fcoeff * bw->interval);
        As->vec_set_zero(As, fcoeffs[k], one_fcoeff * bw->interval);
    }
    /* }}} */
    
    /* {{{ tell something about the end value. Both relative to the commabd
     * line, and to what we expect as a total duration */
    if (bw->end == 0) {
        // for mksol we use EOF as an ending indication.
        serialize_threads(pi->m);
        if (pi->m->trank == 0)
            bw->end = INT_MAX;
        serialize_threads(pi->m);
        if (tcan_print) {
            printf ("Target iteration is unspecified ;"
                    " going to end of F file\n");
        }
    } else if (tcan_print) {
        printf ("Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }
    
    /* We can either look up the F files to get an idea of the number of
     * coefficients to be considered, or make our guess based on the
     * expected degree of the generator. The latter is obviously less
     * accurate, but not by much, and anyway for the purpose of making up
     * an ETA, it's good enough.
     */

    int exp_end = iceildiv(MAX(mmt->n[0], mmt->n[1]), bw->n);
    /* }}} */

    /* {{{ bless our timers */
    timing_init(timing, 4 * mmt->nmatrices, bw->start, bw->interval * iceildiv(exp_end, bw->interval));
    for(int i = 0 ; i < mmt->nmatrices; i++) {
        timing_set_timer_name(timing, 4*i, "CPU%d", i);
        timing_set_timer_items(timing, 4*i, mmt->matrices[i]->mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%d", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%d", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%d", i);
    }
    /* }}} */

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        serialize(pi->m);
        for(int i = 0 ; i < bw->n / splitwidth ; i++) {
            int ys[2] = { i * splitwidth, (i + 1) * splitwidth };
            char * v_name = NULL;
            int rc = asprintf(&v_name, "V%u-%u.%u", ys[0], ys[1], s);
            ASSERT_ALWAYS(rc >= 0);
            if (tcan_print) {
                printf(fake ? "Creating fake %s..." : "Loading %s ...",
                        v_name);
                fflush(stdout);
            }
            if (fake) {
                mmt_vec_set_random_through_file(vi[i], v_name, unpadded, rstate);
            } else {
                mmt_vec_load(vi[i], v_name, unpadded);
                mmt_vec_reduce_mod_p(vi[i]);
            }
            mmt_vec_twist(mmt, vi[i]);
            if (tcan_print) { printf("done\n"); }
            free(v_name);
        }
        serialize(pi->m);
        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            int rc0 = 0, rc = 0;
            for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                for(unsigned int j = 0 ; j < As_multiplex ; j++) {
                    void * ff = fcoeffs[i * As_multiplex + j];

                    /* Zero out first, then read; when we reach EOF while
                     * reading F, it is crucially important that we have
                     * a zero area past the end of the file ! */
                    As->vec_set_zero(As, ff, one_fcoeff * bw->interval);

                    char * tmp;
                    rc = asprintf(&tmp, "F.sols%u-%u.%u-%u",
                            solutions[0] + j * As_width,
                            solutions[0] + (j + 1) * As_width,
                            i * Av_width, (i + 1) * Av_width);

                    FILE * f = fopen(tmp, "rb");
                    DIE_ERRNO_DIAG(f == NULL, "fopen", tmp);
                    rc = fseek(f, one_fcoeff * s, SEEK_SET);
                    if (rc >= 0) {
                        /* Read everything in one go. We might want to
                         * reconsider how the coefficients inside the
                         * individual F files are written -- transposed or
                         * not. Of course that does not matter for the prime
                         * case where individual F files are polynomials, but
                         * surely it does in the binary case where individual
                         * F files are a priori made of 64*64 matrices.
                         */
                        rc = fread(ff, one_fcoeff, bw->interval, f);
                        ASSERT_ALWAYS(rc >= 0 && rc <= bw->interval);
                    } else {
                        /* Otherwise we're off bounds already, don't try
                         * to read anything...
                         */
                        rc = 0;
                    }

                    if (i == 0 && j == 0) rc0 = rc;

                    if (rc != rc0) {
                        fprintf(stderr, "Inconsistency in number of "
                                "coefficients for F files "
                                "F.0-%u.sols%u-%u and %s\n",
                                Av_width,
                                solutions[0], solutions[0] + As_width,
                                tmp);
                        exit(EXIT_FAILURE);
                    }
                    /* TODO: *maybe* transpose the F coefficients at this
                     * point */
                    fclose(f);
                    free(tmp);
                }
            }

            if (rc0 < bw->interval) {
                printf("We read %d coefficients from F in total\n", s+rc);
                bw->end = s+rc;
            }
        }
        if (pi->m->trank == 0) {
            /* we're good friends with the global leader. Try to grasp
             * the information about the end value -- it might have
             * changed because of EOF. We're not so disturbed
             * by serialization points here. */
            int err = MPI_Bcast(&(bw->end), 1, MPI_INT, 0, pi->m->pals);
            ASSERT_ALWAYS(!err);
        }
        serialize(pi->m);

        /* broadcast f */
        for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
            for(unsigned int j = 0 ; j < As_multiplex ; j++) {
                void * ff = fcoeffs[i * As_multiplex + j];
                pi_bcast(ff, Av->groupsize(Av) * bw->interval, As_pi, 0, 0, pi->m);
            }
        }

        serialize(pi->m);

        /* this loop on j is somewhat ridiculous, because we know that we
         * force As_multiplex to be equal to 1 above. But well,
         * conceivably, we *could* do that. It's a bit crazy because it
         * has no advantage over multiple mksol runs (hmm, except maybe
         * with interleaving on ?).
         *
         * (also, given the way it is presently, As_multiplex > 1 will
         * wreak havoc with timings).
         */
        for(unsigned int j = 0 ; j < As_multiplex ; j++) {
            /* Create an empty slot in program execution, so that we don't
             * impose strong constraints on twist/untwist_vector being free of
             * MPI calls.
            pi_interleaving_flip(pi);
            pi_interleaving_flip(pi);
            mmt_vec_twist(mmt, ymy[0]);
             */

            mmt_full_vec_set_zero(ymy[0]);

            serialize(pi->m);

            /* Despite the fact that the bw->end value might lead us
             * to stop earlier, we want to stop at multiples of the checking
             * interval. That's important, otherwise the last bunch of
             * computations won't be checked.
             */

            pi_interleaving_flip(pi);

            size_t eblock = mmt_my_own_size_in_items(ymy[0]);

            for(int k = 0 ; k < bw->interval ; k++) {
                if (k)
                    matmul_top_mul(mmt, ymy, timing);

                serialize_threads(pi->m);

                for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                        void * ff = fcoeffs[i * As_multiplex + j];
                        /* or maybe transpose here instead ?? */
                        AvxAs->addmul_tiny(Av, As,
                                mmt_my_own_subvec(ymy[0]),
                                mmt_my_own_subvec(vi[i]),
                                As->vec_subvec(As, ff,
                                    (bw->interval - 1 - k) * Av->groupsize(Av)),
                                eblock);
                }
                ymy[0]->consistency = 1;
                mmt_vec_broadcast(ymy[0]);
                /* we're doing something which we normally avoid: write
                 * on the next input vector. This means a race condition
                 * with the forthcoming spmv, so we have to serialize.
                 *
                 * I believe we could probably get away with this.
                 */

                serialize_threads(pi->m);

                timing_check(pi, timing, s+k+1, tcan_print);
            }

            serialize(pi->m);
            /* See remark above. */
            pi_interleaving_flip(pi);
            pi_interleaving_flip(pi);

            mmt_vec_untwist(mmt, ymy[0]);

            // FIXME. serialize here ??
            // serialize(pi->m);

            if (!fake) {
                char * s_name;
                int rc = asprintf(&s_name, "S.sols%u-%u.%u-%u",
                        solutions[0] + j * As_width,
                        solutions[0] + (j + 1) * As_width,
                        s, s + bw->interval);
                ASSERT_ALWAYS(rc >= 0);
                mmt_vec_save(ymy[0], s_name, unpadded);
                free(s_name);
            }

        }
        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "mksol");
    }

    timing->end_mark = bw->start + bw->interval * iceildiv(bw->end - bw->start, bw->interval);

    timing_final_tally(pi, timing, tcan_print, "mksol");

    if (tcan_print) {
        printf("Done mksol.\n");
    }
    serialize(pi->m);
    if (fake) gmp_randclear(rstate);
    for(int i = 0 ; i < bw->n / splitwidth ; i++) {
        mmt_vec_clear(mmt, vi[i]);
    }
    delete[] vi;

    for(unsigned int k = 0 ; k < Av_multiplex * As_multiplex ; k++) {
        cheating_vec_clear(As, &(fcoeffs[k]), one_fcoeff * bw->interval);
    }
    delete[] fcoeffs;

    for(int i = 0 ; i < mmt->nmatrices + nmats_odd ; i++) {
        mmt_vec_clear(mmt, ymy[i]);
    }
    delete[] ymy;

    matmul_top_clear(mmt);
    pi_free_mpfq_datatype(pi, Av_pi);

    As->oo_field_clear(As);
    Av->oo_field_clear(Av);

    timing_clear(timing);

    return NULL;
}

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
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(mksol_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}


