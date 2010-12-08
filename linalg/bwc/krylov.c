#include "cado.h"

#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "debug.h"
#include "params.h"
#include "xvectors.h"
#include "xymats.h"
#include "bw-common-mpi.h"
#include "async.h"
// #include "rusage.h"
#include "filenames.h"
#include "xdotprod.h"
#include "rolling.h"

void * krylov_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    struct timing_data timing[1];
    abobj_t abase;

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->ys[1]-bw->ys[0]) % 2 == 0);
        ys[0] = bw->ys[0] + pi->interleaved->idx * (bw->ys[1]-bw->ys[0])/2;
        ys[1] = ys[0] + (bw->ys[1]-bw->ys[0])/2;
    }

    abobj_init(abase);
    abobj_set_nbys(abase, ys[1]-ys[0]);

    block_control_signals();

    matmul_top_init(mmt, abase, pi, flags, pl, bw->dir);

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    serialize(pi->m);
    
    abzero(abase, mrow->v->v, mrow->i1 - mrow->i0);

    uint32_t * gxvecs = malloc(bw->nx * bw->m * sizeof(uint32_t));

    load_x(gxvecs, bw->m, bw->nx, pi, mmt->bal);

    char * v_name;

    int rc = asprintf(&v_name, V_FILE_BASE_PATTERN, ys[0], ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, bw->dir, bw->start);
    if (tcan_print) { printf("done\n"); }

    size_t stride =  abbytes(abase, 1);

    // We must also fetch the check vector C.
    abvobj_t abase_check;
    mmt_generic_vec check_vector;
    mmt_vec ahead;

    if (!bw->skip_online_checks) {
        abvobj_set_nbys(abase_check, NCHECKS_CHECK_VECTOR);
        size_t vstride = abvbytes(abase_check, 1);
        matmul_top_vec_init_generic(mmt, abvbytes(abase_check, 1),
                check_vector, !bw->dir, THREAD_SHARED_VECTOR);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        matmul_top_load_vector_generic(mmt, vstride,
                check_vector, CHECK_FILE_BASE, !bw->dir, bw->interval);
        if (tcan_print) { printf("done\n"); }
    }

    if (!bw->skip_online_checks) {
        vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, 0, NCHECKS_CHECK_VECTOR);
    }

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    mmt_vec xymats;
    if (tcan_print) {
        printf("Each thread allocates %zu kb for the A matrices\n",
                abbytes(abase, bw->m*bw->interval) >> 10);
    }
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, bw->m*bw->interval);
    
    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        length = MAX(mmt->n[0], mmt->n[1]);
        length = iceildiv(length, bw->m) + iceildiv(length, bw->n);
        length += 2 * iceildiv(bw->m, bw->n);
        length += 2 * iceildiv(bw->n, bw->m);
        length += 10;
        /* Because bw is a global variable, we protect its use */
        if (serialize_threads(pi->m)) {
            bw->end = length;
        }
        serialize(pi->m);
    }

    if (tcan_print) {
        fprintf(stderr, "Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }

#if 0
    /* FIXME -- that's temporary ! only for debugging */
    pi_log_init(pi->m);
    pi_log_init(pi->wr[0]);
    pi_log_init(pi->wr[1]);
#endif

    timing_init(timing, bw->start, bw->interval * iceildiv(bw->end, bw->interval));

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                mrow->i0, mrow->i1,
                mcol->i0, mcol->i1);

        if (!bw->skip_online_checks) {
            abzero(abase, ahead->v, NCHECKS_CHECK_VECTOR);
            if (how_many) {
                size_t bytes_c =  abvbytes(abase_check, offset_c);
                abvt * c = abase_generic_ptr_add(check_vector->v, bytes_c);
                abt * v = mcol->v->v + aboffset(abase, offset_v);
                abvdotprod(abase, abase_check, ahead->v, c, v, how_many);
            }
        }

        abzero(abase, xymats->v, bw->m*bw->interval);
        serialize(pi->m);
        for(int i = 0 ; i < bw->interval ; i++) {
            pi_interleaving_flip(pi);
            timing_flip_timer(timing);
            /* This segment must be guaranteed to be free of any mpi
             * calls */
            /* Compute the product by x */
            x_dotprod(mmt, gxvecs,
                    xymats->v + aboffset(abase, i * bw->m), bw->m);

            
            matmul_top_mul_cpu(mmt, bw->dir);
            pi_interleaving_flip(pi);
            timing_flip_timer(timing);
            /* Now we can resume MPI communications. */
            matmul_top_mul_comm(mmt, bw->dir);
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        if (!bw->skip_online_checks) {
            /* Last dot product. This must cancel ! */
            x_dotprod(mmt, gxvecs, ahead->v, NCHECKS_CHECK_VECTOR);

            allreduce_generic(abase, ahead, pi->m, NCHECKS_CHECK_VECTOR);
            if (!abis_zero(abase, ahead->v, NCHECKS_CHECK_VECTOR)) {
                printf("Failed check at iteration %d\n", s + bw->interval);
                exit(1);
            }
        }

        /* Now (and only now) collect the xy matrices */
        allreduce_generic(abase, xymats, pi->m, bw->m * bw->interval);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            rc = asprintf(&tmp, A_FILE_PATTERN, ys[0], ys[1], s, s+bw->interval);
            FILE * f = fopen(tmp, "w");
            rc = fwrite(xymats->v, stride, bw->m*bw->interval, f);
            if (rc != bw->m*bw->interval) {
                fprintf(stderr, "Ayee -- short write\n");
                // make sure our input data won't be deleted -- this
                // chunk will have to be redone later, maybe the disk
                // failure is temporary (?)
            }
            fclose(f);
            free(tmp);
        }

        matmul_top_save_vector(mmt, v_name, bw->dir, s + bw->interval);
        if (pi->m->trank == 0 && pi->m->jrank == 0)
            keep_rolling_checkpoints(mmt->bal, v_name, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, mmt->mm->ncoeffs, tcan_print, 0);
    }

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (tcan_print) {
        printf("Done.\n");
    }
    serialize(pi->m);

    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, bw->m*bw->interval);

    if (!bw->skip_online_checks) {
        matmul_top_vec_clear_generic(mmt, abvbytes(abase_check, 1), check_vector, !bw->dir);
        vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, NCHECKS_CHECK_VECTOR);
    }

    free(gxvecs);
    free(v_name);

    matmul_top_clear(mmt, abase);

    timing_clear(timing);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./krylov <options>\n");
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

    if (bw->nx == 0) { fprintf(stderr, "no nx value set\n"); exit(1); } 
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    catch_control_signals();
    pi_go(krylov_prog, pl, 0);
    param_list_clear(pl);
    bw_common_clear_mpi(bw);

    return 0;
}

