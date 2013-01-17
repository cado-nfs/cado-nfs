#include "cado.h"
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "gauss.h"
#include "gauss.h"
#include "params.h"
#include "xvectors.h"
#include "xymats.h"
#include "bw-common-mpi.h"
#include "filenames.h"
#include "mpfq/mpfq.h"
#include "mpfq/abase_vbase.h"

void * prep_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    if (pi->interleaved && pi->interleaved->idx)
        return NULL;

    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    hello(pi);

    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    mpz_t p;
    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "prime", p);
    abase_vbase A;
    abase_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, p,
            MPFQ_GROUPSIZE, bw->n,
            MPFQ_DONE);
    mpz_clear(p);


    matmul_top_init(mmt, A, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    /* Number of copies of m by n matrices to use for trying to obtain a
     * full-rank matrix (rank m).
     *
     * Note that it must be at least m/n, otherwise we stand no chance !
     */
    unsigned int prep_lookahead_iterations;
    prep_lookahead_iterations = iceildiv(bw->m, bw->m) + 1;

    unsigned int my_nx = 1;

    mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    uint32_t * xvecs = malloc(my_nx * bw->m * sizeof(uint32_t));

    mmt_vec xymats;

    /* We're cheating on the generic init routines */
    vec_init_generic(pi->m,
            A,
            xymats,
            0,
            bw->m * prep_lookahead_iterations);

    for (unsigned ntri = 0;; ntri++) {
        serialize_threads(pi->m);

        if (tcan_print) {
            printf("// Generating new x,y vector pair (trial # %u -- seed %lu)\n", ntri, (unsigned long) rand());
        }
        if (ntri >= my_nx * 10) {
            ++my_nx;
            if (tcan_print) {
                printf("// Getting bored. Trying %u x vectors\n", my_nx);
            }
            xvecs = realloc(xvecs, my_nx * bw->m * sizeof(uint32_t));
            ASSERT_ALWAYS(xvecs != NULL);
        }

        // if we're looking for the right nullspace, then x is on the left.
        // Otherwise, it's on the right.

        // generate indices w.r.t *unpadded* dimensions !
        setup_x_random(xvecs, bw->m, my_nx, mmt->n0[bw->dir], pi);

        /* Random generation + save is better done as writing random data
         * to a file followed by reading it: this way, seeding works
         * better.
         */

        if (pi->m->trank == 0) {
            const char * filename = Y_FILE_BASE ".0";

            void * y;
            A->vec_init(A, &y, mmt->n[bw->dir]);
            A->vec_set_zero(A, y, mmt->n[bw->dir]);
            /* Again, important. Generate zero coordinates for padding !
             * This provides reproducibility of random choices.
             */
            if (pi->m->jrank == 0)
                A->vec_random(A, y, mmt->n0[bw->dir]);
            int err = MPI_Bcast(y,
                    mmt->n[bw->dir],
                    A->mpi_datatype(A),
                    0, pi->m->pals);
            ASSERT_ALWAYS(!err);
            FILE * f = fopen(filename, "w");
            ASSERT_ALWAYS(f);
            int rc = fwrite(y, A->vec_elt_stride(A,1), unpadded, f);
            ASSERT_ALWAYS(rc == (int) unpadded);
            fclose(f);
            A->vec_clear(A, &y, mmt->n[bw->dir]);
        }
        matmul_top_load_vector(mmt, Y_FILE_BASE, bw->dir, 0, unpadded);

        if (tcan_print) {
            printf("// vector generated and dispatched (trial # %u)\n", ntri);
        }

#if 0
        // Compute y.
        matmul_top_fill_random_source(mmt, bw->dir);

        // we need to save this starting vector for later use if it turns out
        // that we need to save it for real.
        matmul_top_save_vector(mmt, Y_FILE_BASE, bw->dir, 0, unpadded);
#endif

        // We must compute x^T M y, x^T M^2 y, and so on.
        // XXX Note that x^Ty does not count here, because it does not
        // take part to the sequence computed by lingen !
        matmul_top_twist_vector(mmt, bw->dir);
        matmul_top_mul(mmt, bw->dir);
        matmul_top_untwist_vector(mmt, bw->dir);
        
        // we have indices mmt->wr[1]->i0..i1 available.
        A->vec_set_zero(A, xymats->v, bw->m * prep_lookahead_iterations);

        for(unsigned int k = 0 ; k < prep_lookahead_iterations ; k++) {
            for(int j = 0 ; j < bw->m ; j++) {
                void * where = SUBVEC(xymats, v, j * prep_lookahead_iterations + k);
                for(unsigned int t = 0 ; t < my_nx ; t++) {
                    uint32_t i = xvecs[j*my_nx+t];
                    if (i < mcol->i0 || i >= mcol->i1)
                        continue;
                    /* We want the data to match our bw->interval on both
                     * directions, because otherwise we'll end up
                     * computing rubbish -- recall that no broadcast_down
                     * has occurred yet.
                     */
                    if (i < mrow->i0 || i >= mrow->i1)
                        continue;
                    A->add(A, where, where, SUBVEC(mcol->v, v, i - mcol->i0));
                }
            }
            matmul_top_twist_vector(mmt, bw->dir);
            matmul_top_mul(mmt, bw->dir);
            matmul_top_untwist_vector(mmt, bw->dir);
        }

        /* Make sure computation is over for everyone ! */
        serialize_threads(pi->m);

        /* Now all threads and jobs must collectively reduce the zone
         * pointed to by xymats */
        allreduce_generic(xymats, pi->m,
                bw->m * prep_lookahead_iterations);

        /* OK -- now everybody has the same data */

        int dimk;
        int * pdimk;
        
        /* the kernel() call is not reentrant */
        if (pi->m->trank == 0) {
            dimk = kernel((mp_limb_t *) xymats->v, NULL,
                    bw->m, prep_lookahead_iterations * A->groupsize(A),
                    A->vec_elt_stride(A, prep_lookahead_iterations)/sizeof(mp_limb_t),
                    0);
            pdimk = &dimk;
        }
        thread_broadcast(pi->m, (void **) &pdimk, 0);
        dimk = * pdimk;

        if (tcan_print)
            printf("// Dimension of kernel: %d\n", dimk);

        if (dimk == 0) {
            if (tcan_print)
                printf("// Found good x,y vector pair after %u trials\n",
                        ntri+1);
            break;
        }
    }

    save_x(xvecs, bw->m, my_nx, pi);

    matmul_top_clear(mmt);

    /* clean up xy mats stuff */
    vec_clear_generic(pi->m, xymats, bw->m * prep_lookahead_iterations);

    A->oo_field_clear(A);


    free(xvecs);
    return NULL;
}

void usage()
{
    fprintf(stderr, "Usage: ./prep <options>\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);
    bw_common_init_mpi(bw, pl, &argc, &argv);
    if (param_list_warn_unused(pl)) usage();

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);
    
    srand(bw->seed ? bw->seed : time(NULL));

    pi_go(prep_prog, pl, 0);

    param_list_remove_key(pl, "sequential_cache_build");
    param_list_remove_key(pl, "rebuild_cache");

    param_list_clear(pl);

    bw_common_clear_mpi(bw);
    return 0;
}

