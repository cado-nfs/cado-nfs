#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "debug.h"
#include "params.h"
#include "xvectors.h"
#include "xymats.h"
#include "misc.h"
#include "bw-common-mpi.h"
#include "async.h"
#include "filenames.h"
#include "xdotprod.h"
#include "rolling.h"
#include "mpfq/mpfq.h"
#include "mpfq/abase_vbase.h"

void * mksol_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    struct timing_data timing[1];

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->ys[1]-bw->ys[0]) % 2 == 0);
        ys[0] = bw->ys[0] + pi->interleaved->idx * (bw->ys[1]-bw->ys[0])/2;
        ys[1] = ys[0] + (bw->ys[1]-bw->ys[0])/2;
    }

    int withcoeffs = param_list_lookup_string(pl, "prime") != NULL;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpz_t p;
    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "prime", p);
    abase_vbase A;
    abase_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, p,
            MPFQ_GROUPSIZE, ys[1]-ys[0],
            MPFQ_DONE);
    /* Hmmm. This would deserve better thought. Surely we don't need 64
     * in the prime case. Anything which makes checks relevant will do.
     * For the binary case, we used to work with 64 as a constant, but
     * for the prime case we want to make this tunable (or maybe 1 ?)
     */
    abase_vbase Ac;
    abase_vbase_oo_field_init_byfeatures(Ac,
            MPFQ_PRIME_MPZ, p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    abase_vbase Ar;
    abase_vbase_oo_field_init_byfeatures(Ar,
            MPFQ_PRIME_MPZ, p,
            MPFQ_GROUPSIZE, bw->n,
            MPFQ_DONE);
    mpz_clear(p);

    if (pi->m->trank == 0) Ar->mpi_ops_init(Ar);

    block_control_signals();

    matmul_top_init(mmt, A, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];
    pi_wiring_ptr picol = mmt->pi->wr[bw->dir];

    serialize(pi->m);
    
    A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);

    uint32_t * gxvecs = NULL;
    unsigned int nx = 0;
    
    if (!bw->skip_online_checks) {
        load_x(&gxvecs, bw->m, &nx, pi);
        /*
        indices_apply_S(mmt, gxvecs, nx * bw->m, bw->dir);
        if (bw->dir == 0)
            indices_apply_P(mmt, gxvecs, nx * bw->m, bw->dir);
            */
    }

    int rc;

    char * v_name;
    rc = asprintf(&v_name, V_FILE_BASE_PATTERN, ys[0], ys[1]);
    char * s_name;
    rc = asprintf(&s_name, S_FILE_BASE_PATTERN, ys[0], ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, bw->dir, bw->start, unpadded);
    if (tcan_print) { printf("done\n"); }

    mmt_vec check_vector;
    /* Our ``ahead zone'' will also be used for storing the data prior to
     * transposing */
    mmt_vec ahead;
    vec_init_generic(pi->m, A, ahead, 0, bw->n);

    abase_vbase_tmpl AxAc;
    abase_vbase_oo_init_templates(AxAc, A, Ac);

    abase_vbase_tmpl ArxA;
    abase_vbase_oo_init_templates(ArxA, Ar, A);

    abase_vbase_tmpl AxAr;
    abase_vbase_oo_init_templates(AxAr, A, Ar);

    if (!bw->skip_online_checks) {
        matmul_top_vec_init_generic(mmt, Ac,
                check_vector, !bw->dir, THREAD_SHARED_VECTOR);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        matmul_top_load_vector_generic(mmt,
                check_vector, CHECK_FILE_BASE, !bw->dir, bw->interval, unpadded);
        if (tcan_print) { printf("done\n"); }
    }

    /* F plays the role of a right-hand-side in mksol. Vector iterates
     * have to be multiplied on the right by a matrix corresponding to
     * some coefficient of F. Because F has been written to disk in the
     * proper order (=reversed order, leading coeff first) by lingen,
     * iterate i is to be multiplied by the i-th coefficient of F.
     *
     * The F-coefficients, from the overall description of the
     * algorithm, are square matrices having bw->n rows and columns.
     * Columns correspond to candidate solutions. Rows correspond to the
     * fact that contributions from all vector iterates have to be added
     * in order to form a candidate solution. However, for
     * the purpose of splitting the work across sites, it's stored on
     * disk in the transposed order, so that the already existing split
     * program is able to do the split. This means that the data we read
     * here consists of matrices having unconditionally bw->n rows, and
     * ys[1]-ys[0] columns. Because our multiply-by-smallish-matrix code
     * does not like this ordering, we have to transpose again in order
     * to obtain a matrix with ys[1]-ys[0] rows, and bw->n columns. This
     * matrix is called rhs in the text below.
     */

    ASSERT(A->vec_elt_stride(A,bw->n) == Ar->vec_elt_stride(Ar, A->groupsize(A)));

    /* We'll load all the F coefficient matrices before the main loop.
     * It's dual to the treatment of the xy matrices in krylov. Same
     * considerations apply w.r.t the memory footprint.
     */
    mmt_vec fcoeffs;
    if (tcan_print) {
        printf("Each thread allocates %zd kb for the F matrices\n",
                Ar->vec_elt_stride(Ar, A->groupsize(A)*bw->interval) >> 10);
    }
    vec_init_generic(pi->m, Ar, fcoeffs, 0, A->groupsize(A)*bw->interval);
    
    /* Our sum vector is only part of the story of course. All jobs, all
     * threads, will participate to the computation of the products by
     * tiny matrices. When dir==1, within a column, all threads have
     * access to the mcol->i0..mcol->i1 interval (which is roughly a
     * fraction equal to 1 / pirow->totalsize of the total
     * vector size. Now this interval will be split in picol->totalsize
     * parts.
     */

    unsigned int ii0, ii1;
    unsigned int di = mcol->i1 - mcol->i0;

    ii0 = mcol->i0 + di * (picol->jrank * picol->ncores + picol->trank) /
        picol->totalsize;
    ii1 = mcol->i0 + di * (picol->jrank * picol->ncores + picol->trank + 1) /
        picol->totalsize;

    mmt_vec sum;
    vec_init_generic(mmt->pi->m, Ar, sum, 0, ii1-ii0);

    if (bw->end == 0) {
        // for mksol we use EOF as an ending indication.
        serialize_threads(pi->m);
        bw->end = INT_MAX;
        if (tcan_print) {
            fprintf(stderr, "Target iteration is unspecified ;"
                    " going to end of F file\n");
        }
    } else if (tcan_print) {
        fprintf(stderr, "Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }
    
    /* We can either look up the F files to get an idea of the number of
     * coefficients to be considered, or make our guess based on the
     * expected degree of the generator. The latter is obiviously less
     * accurate, but not by much, and anyway for the purpose of making up
     * an ETA, it's good enough.
     */

    int exp_end;
    exp_end = MAX(mmt->n[0], mmt->n[1]);
    exp_end = iceildiv(exp_end, bw->n);

    timing_init(timing, bw->start, exp_end);

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        serialize(pi->m);
        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            Ar->vec_set_zero(Ar, fcoeffs->v, A->groupsize(A)*bw->interval);
            char * tmp;
            rc = asprintf(&tmp, F_FILE_SLICE_PATTERN, ys[0], ys[1]);

            FILE * f = fopen(tmp, "r");
            DIE_ERRNO_DIAG(f == NULL, "fopen", tmp);
            rc = fseek(f, A->vec_elt_stride(A, s * bw->n), SEEK_SET);
            if (rc < 0) {
                bw->end = s;
            }

            for(int i = 0 ; i < bw->interval ; i++) {
                rc = fread(ahead->v, A->vec_elt_stride(A,1), bw->n, f);
                ASSERT_ALWAYS(rc == 0 || rc == bw->n);
                if (rc == 0) {
                    printf("Exhausted F data after %d coefficients\n", s+i);
                    bw->end = s+i;
                    break;
                }

#if 0
                {
                    void * test;
                    test = electric_alloc(rstride * abnbits(abase));
                    abvtranspose(abase_rhs, abase, test, ahead->v);
                    electric_free(test,rstride * abnbits(abase));
                }
#endif
                ArxA->transpose(Ar->obj, A->obj,
                        SUBVEC(fcoeffs, v, i * A->groupsize(A)),
                        ahead->v);
            }
            fclose(f);
            free(tmp);
        }
        if (pi->m->trank == 0) {
            /* we're good friends with the global leader. Try to grasp
             * the information about the end value -- it might have
             * changed because of EOF. We're not so disturbed
             * by serialization points here. */
            int err = MPI_Bcast(&(bw->end), 1, MPI_INT, 0, pi->m->pals);
            ASSERT_ALWAYS(!err);
        }
        serialize_threads(pi->m);

        /* broadcast f */
        broadcast_generic(fcoeffs, pi->m, A->groupsize(A) * bw->interval, 0, 0);

        Ar->vec_set_zero(Ar, sum->v, ii1 - ii0);

        // Plan ahead. The check vector is here to predict the final A matrix.
        // Our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        if (!bw->skip_online_checks) {
            A->vec_set_zero(A, ahead->v, nchecks);
            unsigned int how_many;
            unsigned int offset_c;
            unsigned int offset_v;
            how_many = intersect_two_intervals(&offset_c, &offset_v,
                    mrow->i0, mrow->i1,
                    mcol->i0, mcol->i1);

            if (how_many) {
                AxAc->dotprod(A->obj, Ac->obj, ahead->v,
                        SUBVEC(check_vector, v, offset_c),
                        SUBVEC(mcol->v, v, offset_v),
                        how_many);
            }
        }

        serialize(pi->m);
        /* Create an empty slot in program execution, so that we don't
         * impose strong constraints on twist/untwist_vector being free of
         * MPI calls.
         */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);
        matmul_top_twist_vector(mmt, bw->dir);
        serialize(pi->m);

        /* Despite the fact that the bw->end value might lead us
         * to stop earlier, we want to stop at multiples of the checking
         * interval. That's important, otherwise the last bunch of
         * computations won't be checked.
         */
        for(int i = 0 ; i < bw->interval ; i++) {
            pi_interleaving_flip(pi);
            timing_flip_timer(timing);
            /* This segment must be guaranteed to be free of any mpi
             * calls */

            AxAr->addmul_tiny(A->obj, Ar->obj,
                    sum->v,
                    SUBVEC(mcol->v, v, ii0 - mcol->i0),
                    SUBVEC(fcoeffs, v, i * A->groupsize(A)),
                    ii1 - ii0);

            matmul_top_mul_cpu(mmt, bw->dir);
            pi_interleaving_flip(pi);
            timing_flip_timer(timing);
            /* Now we can resume MPI communications. */
            matmul_top_mul_comm(mmt, bw->dir);
            timing_check(pi, timing, s+i+1, tcan_print);
        }

        /* See remark above. */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);

        matmul_top_untwist_vector(mmt, bw->dir);

        // FIXME. serialize here ??
        // serialize(pi->m);

        if (!bw->skip_online_checks) {
            /* Last dot product. This must cancel ! Recall that x_dotprod
             * adds to the result. */
            x_dotprod(mmt, gxvecs, nx, ahead, 0, nchecks, -1);

            allreduce_generic(ahead, pi->m, nchecks);
            if (!A->vec_is_zero(A, ahead->v, nchecks)) {
                printf("Failed check at iteration %d\n", s + bw->interval);
                exit(1);
            }
        }

        matmul_top_save_vector(mmt, v_name, bw->dir, s + bw->interval, unpadded);

        /* We need to write our S files, untwisted. There are several
         * ways to do this untwisting. In order to keep the memory
         * footprint to a minimum, we reuse the vector area mcol->v. 
         */
        /* Since S is possibly wider than one single vector, we need
         * several passes. Each corresponds to one batch of potential
         * solutions in the end.
         *
         * We first make sure that an integral number of passes is
         * needed.
         */
        serialize(pi->m);
        size_t  stride =  A->vec_elt_stride(A, 1);
        size_t rstride =  Ar->vec_elt_stride(Ar, 1);
        ASSERT_ALWAYS(rstride % stride == 0);
        int npasses = rstride / stride;
        for(int i = 0 ; i < npasses ; i++) {
            serialize(pi->m);
            A->vec_set_zero(A, mcol->v->v, mcol->i1 - mcol->i0);
            serialize(pi->m);
            /* Each job/thread copies its data share to mcol->v */
            void * dst;
            const void * src;
            src = pointer_arith_const(sum->v, i * stride);
            dst = pointer_arith(mcol->v->v, (ii0 - mcol->i0) * stride);
            for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
                /* XXX If we dereference a function pointer for each such
                 * iteration here, we're going to pay a lot. And the mpfq
                 * API does not contain this ``strided copy''. So for the
                 * moment, since all types considered are flat anyway,
                 * this looks good enough.
                 */
                memcpy(dst, src, stride);
                dst = pointer_arith(dst, stride);
                src = pointer_arith_const(src, rstride);
            }
            /* Now we collect the data in mcol->v. Note that each
             * location is non-zero only at one location, so unreduced
             * addition is unnecessary.
             */
            allreduce_across(mmt, bw->dir);
            /* untwist the result */
            matmul_top_untwist_vector(mmt, bw->dir);

            serialize_threads(pi->m);
            dst = pointer_arith(sum->v, i * stride);
            src = pointer_arith_const(mcol->v->v, (ii0 - mcol->i0) * stride);
            for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
                memcpy(dst, src, stride);
                dst = pointer_arith(dst, rstride);
                src = pointer_arith_const(src, stride);
            }
        }
        pi_save_file_2d(pi, bw->dir, s_name, s+bw->interval, sum->v, Ar->vec_elt_stride(Ar, ii1 - ii0), Ar->vec_elt_stride(Ar, unpadded));
        serialize(pi->m);

        /* recover the vector */
        matmul_top_load_vector(mmt, v_name, bw->dir, s + bw->interval, unpadded);

        if (pi->m->trank == 0 && pi->m->jrank == 0)
            keep_rolling_checkpoints(v_name, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, mmt->mm->ncoeffs, tcan_print, 1);
    }

    if (tcan_print) {
        printf("Done.\n");
    }
    serialize(pi->m);

    matmul_top_vec_clear_generic(mmt, sum, bw->dir);
    if (!bw->skip_online_checks) {
        matmul_top_vec_clear_generic(mmt, check_vector, !bw->dir);
        vec_clear_generic(pi->m, ahead, nchecks);
        free(gxvecs);
    }
    vec_clear_generic(pi->m, fcoeffs, A->groupsize(A)*bw->interval);
    free(v_name);
    free(s_name);

    matmul_top_clear(mmt);

    if (pi->m->trank == 0) Ar->mpi_ops_clear(Ar);
    Ar->oo_field_clear(Ar);
    Ac->oo_field_clear(Ac);
    A->oo_field_clear(A);

    timing_clear(timing);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./mksol <options>\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix interval start ys\n");
    fprintf(stderr, "Note: data files must be found in wdir !\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    bw_common_init_mpi(bw, pl, &argc, &argv);
    if (param_list_warn_unused(pl)) usage();

    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    catch_control_signals();
    pi_go(mksol_prog, pl, 0);
    param_list_clear(pl);
    bw_common_clear_mpi(bw);

    return 0;
}

