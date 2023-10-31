#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cinttypes>
#include <cstdint>              // for uint32_t
#include <cstring>              // for memset
#include <ctime>                // for time
#include <cstdlib>
#include <memory>
#include <gmp.h>
#include "balancing.hpp"           // for balancing_pre_shuffle
#include "parallelizing_info.hpp"
#include "matmul_top.hpp"
#include "select_mpi.h"
#include "bblas_gauss.h"
#include "params.h"
#include "xvectors.hpp"
#include "bw-common.h"
#include "arith-generic.hpp"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "macros.h"
#include "cxx_mpz.hpp"
#include "mmt_vector_pair.hpp"
#include "utils_cxx.hpp"


void bw_rank_check(matmul_top_data & mmt, param_list_ptr pl)
{
    int tcan_print = bw->can_print && mmt.pi->m->trank == 0;
    unsigned int r = matmul_top_rank_upper_bound(mmt);
    if (tcan_print) {
        printf("Matrix rank is at most %u (based on zero columns and rows encountered)\n", r);
    }
    int skip=0;
    param_list_parse_int(pl, "skip_bw_early_rank_check", &skip);
    if (bw->m + r < mmt.n0[0]) {
        fprintf(stderr, "Based on the parameter m (=%u) and the rank defect of the matrix (>=%u), we can't expect to compute solutions reliably.\n",
                bw->m, mmt.n0[0]-r);
        if (skip) {
            fprintf(stderr, "Proceeding anyway as per skip_bw_early_rank_check=1\n");
        } else {
            fprintf(stderr, "Aborting. Use skip_bw_early_rank_check=1 to proceed nevertheless.\n");
            exit(EXIT_FAILURE);
        }
    }
}

void * prep_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    ASSERT_ALWAYS(!pi->interleaved);

    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    pi_hello(pi);

    int tcan_print = bw->can_print && pi->m->trank == 0;

    unsigned int nrhs = 0;
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    int splitwidth = char2 ? 64 : 1;

    unsigned int A_width = splitwidth;
    unsigned int A_multiplex = bw->n / A_width;
    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, A_width));
    ASSERT_ALWAYS(A->simd_groupsize() * A_multiplex == (unsigned int) bw->n);

    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    bw_rank_check(mmt, pl);


    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = time(NULL);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }
    serialize_threads(pi->m);

    gmp_randseed_ui(rstate, bw->seed);
    if (tcan_print) {
        printf("// Random generator seeded with %d\n", bw->seed);
    }

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    FILE * rhs = NULL;
    if (rhs_name) {
        rhs = fopen(rhs_name, "r");
        get_rhs_file_header_stream(rhs, NULL, &nrhs, NULL);
        ASSERT_ALWAYS(rhs != NULL);
        ASSERT_ALWAYS(nrhs <= mmt.n[!bw->dir]);
    }

    mmt_vector_pair ymy(mmt, bw->dir);

    mmt_vec & y = ymy[0];

    unsigned int unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    /* Number of copies of m by n matrices to use for trying to obtain a
     * matrix of rank m.
     *
     * Note that it must be at least m/n, otherwise we stand no chance !
     */
    unsigned int prep_lookahead_iterations = iceildiv(bw->m, bw->n) + 1;

    unsigned int my_nx = 1;
    uint32_t * xvecs = (uint32_t*) malloc(my_nx * bw->m * sizeof(uint32_t));

    arith_generic::elt * xymats = A->alloc(bw->m * prep_lookahead_iterations * A_multiplex, ALIGNMENT_ON_ALL_BWC_VECTORS);

    for (unsigned ntri = 0;; ntri++) {
        if (nrhs == A_multiplex) {
            if (ntri) ++my_nx;
            if (ntri >= 4) {
                fprintf(stderr, "Cannot find a satisfactory initialization. "
                        "Maybe your RHS vectors are bad ?\n");
                exit(EXIT_FAILURE);
            }
        } else if (ntri >= my_nx * 10) {
            ++my_nx;
            if (tcan_print) {
                printf("// Getting bored. Trying %u x vectors\n", my_nx);
            }
            xvecs = (uint32_t*) realloc(xvecs, my_nx * bw->m * sizeof(uint32_t));
            ASSERT_ALWAYS(xvecs != NULL);
        }
        serialize_threads(pi->m);

        if (tcan_print) {
            printf("// Generating new x,y vector pair (trial # %u)\n", ntri);
        }

        // if we're looking for the right nullspace, then x is on the left.
        // Otherwise, it's on the right.

        // generate indices w.r.t *unpadded* dimensions !
        setup_x_random(xvecs, bw->m, my_nx, mmt.n0[bw->dir], pi, rstate);

        // we have indices mmt.wr[1]->i0..i1 available.
        A->vec_set_zero(xymats, bw->m * prep_lookahead_iterations * A_multiplex);

        ASSERT_ALWAYS(nrhs <= A_multiplex);

        for(unsigned int j = 0 ; j < A_multiplex ; j++) {
            /* Random generation + save is better done as writing random data
             * to a file followed by reading it: this way, seeding works
             * better.
             */
            if (j < nrhs) {
                /* create it as an extraction from the rhs file */
                ASSERT_ALWAYS(0);       /* implement me... See GF(p) below. */
            } else {
                mmt_vec_set_random_through_file(y, "V%u-%u.0", unpadded, rstate, j * A_width);
            }
            if (tcan_print) {
                printf("// generated V%u-%u.0 (trial # %u)\n", 
                        j * A_width, (j + 1) * A_width, ntri);
            }

            // compute x^T M y, x^T M^2 y, and so on. Since we do that
            // piecewise for the different vectors, we first collect
            // everything in the xymats array, and compute the rank later on.
            
            // XXX Note that x^Ty does not count here, because it does not
            // take part to the sequence computed by lingen !
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy.vectors(), NULL);
            mmt_vec_untwist(mmt, y);
            

            /* XXX it's really like x_dotprod, except that we're filling
             * the coefficients in another order (but why ?) */
            for(unsigned int k = 0 ; k < prep_lookahead_iterations ; k++) {
                for(int r = 0 ; r < bw->m ; r++) {
                    arith_generic::elt & where = A->vec_item(xymats, (r * prep_lookahead_iterations + k) * A_multiplex + j);
                    for(unsigned int t = 0 ; t < my_nx ; t++) {
                        uint32_t row = xvecs[r*my_nx+t];
                        unsigned int vi0 = y.i0 + mmt_my_own_offset_in_items(y);
                        unsigned int vi1 = vi0 + mmt_my_own_size_in_items(y);
                        if (row < vi0 || row >= vi1)
                            continue;

                        arith_generic::elt const & coeff = y.abase->vec_item(y.v, row - y.i0);
                        A->add_and_reduce(where, coeff);
                    }
                }
                mmt_vec_twist(mmt, y);
                matmul_top_mul(mmt, ymy.vectors(), NULL);
                mmt_vec_untwist(mmt, y);
            }
        }

        /* Make sure computation is over for everyone ! */
        serialize_threads(pi->m);

        /* Now all threads and jobs must collectively reduce the zone
         * pointed to by xymats */
        pi_allreduce(NULL, xymats,
                bw->m * prep_lookahead_iterations * A_multiplex,
                mmt.pitype, BWC_PI_SUM, pi->m);

        /* OK -- now everybody has the same data */

        int dimk;
        
        /* the kernel() call is not reentrant */
        if (pi->m->trank == 0) {
            dimk = kernel((mp_limb_t *) xymats, NULL,
                    bw->m, prep_lookahead_iterations * A->simd_groupsize() * A_multiplex,
                    A->vec_elt_stride(prep_lookahead_iterations * A_multiplex)/sizeof(mp_limb_t),
                    0);
        }
        pi_thread_bcast((void *) &dimk, 1, BWC_PI_INT, 0, pi->m);

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

    if (rhs_name)
        free(rhs);

    gmp_randclear(rstate);

    /* clean up xy mats stuff */
    A->free(xymats);

    free(xvecs);
    return NULL;
}

/* The GF(p) case is significantly different:
 *  - this is *NOT* a parallel program, although we're happy to use the
 *    same tools as everywhere else for some common operations, namely
 *    the matmul_top layer.
 *  - we do *NOT* perform the * consistency check.
 *  - and we possibly inject the RHS in the data produced.
 */
void * prep_prog_gfp(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    ASSERT_ALWAYS(!pi->interleaved);

    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    pi_hello(pi);

    int tcan_print = bw->can_print && pi->m->trank == 0;

    unsigned int nrhs = 0;
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    int splitwidth = char2 ? 64 : 1;

    const char * rhs_name;
    FILE * rhs;

    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, splitwidth));

    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    // I don't think this was ever tested.
    ASSERT_ALWAYS(mmt.matrices.size() == 1);

    bw_rank_check(mmt, pl);

    if (pi->m->trank || pi->m->jrank) {
        /* as said above, this is *NOT* a parallel program.  */
        return NULL;
    }

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = time(NULL);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }

    gmp_randseed_ui(rstate, bw->seed);

    if (tcan_print) {
        printf("// Random generator seeded with %d\n", bw->seed);
    }

    rhs_name = param_list_lookup_string(pl, "rhs");
    rhs = NULL;
    if (rhs_name) {
        rhs = fopen(rhs_name, "r");
        get_rhs_file_header_stream(rhs, NULL, &nrhs, NULL);
        ASSERT_ALWAYS(rhs != NULL);
        ASSERT_ALWAYS(nrhs <= mmt.n[!bw->dir]);
    }

    /* First create all RHS vectors -- these are just splits of the big
     * RHS block. Those files get created together. */
    if (nrhs) {
        char ** vec_names = (char**) malloc(nrhs * sizeof(char *));
        FILE ** vec_files = (FILE**) malloc(nrhs * sizeof(FILE *));
        for(unsigned int j = 0 ; j < nrhs ; j++) {
            int rc = asprintf(&vec_names[j], "V%d-%d.0", j, j+1);
            ASSERT_ALWAYS(rc >= 0);
            vec_files[j] = fopen(vec_names[j], "wb");
            ASSERT_ALWAYS(vec_files[j] != NULL);
            printf("// Creating %s (extraction from %s)\n", vec_names[j], rhs_name);
        }
        arith_generic::elt * coeff = A->alloc(1);
        cxx_mpz c;
        for(unsigned int i = 0 ; i < mmt.n0[!bw->dir] ; i++) {
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                int rc;
                memset(coeff, 0, A->elt_stride());
                rc = gmp_fscanf(rhs, "%Zd", (mpz_ptr) c);
                ASSERT_ALWAYS(rc == 1);
                A->set(A->vec_item(coeff, 0), c);
                rc = fwrite(coeff, A->elt_stride(), 1, vec_files[j]);
                ASSERT_ALWAYS(rc == 1);
            }
        }
        A->free(coeff);

        for(unsigned int j = 0 ; j < nrhs ; j++) {
            fclose(vec_files[j]);
            free(vec_names[j]);
        }
        free(vec_files);
        free(vec_names);
    }
    /* Now create purely random vectors */
    for(int j = (int) nrhs ; j < bw->n ; j++) {
        char * vec_name;
        FILE * vec_file;
        int rc = asprintf(&vec_name, "V%d-%d.0", j, j+1);
        ASSERT_ALWAYS(rc >= 0);
        vec_file = fopen(vec_name, "wb");
        ASSERT_ALWAYS(vec_file != NULL);
        printf("// Creating %s\n", vec_name);
        unsigned int unpadded = MAX(mmt.n0[0], mmt.n0[1]);
        auto vec = A->alloc(unpadded);
        A->vec_set_random(vec, unpadded, rstate);
        A->vec_set_zero(A->vec_subvec(vec, mmt.n0[bw->dir]), unpadded - mmt.n0[bw->dir]);
        rc = fwrite(vec, A->elt_stride(), unpadded, vec_file);
        ASSERT_ALWAYS(rc >= 0 && ((unsigned int) rc) == unpadded);
        A->free(vec);
        fclose(vec_file);
        free(vec_name);
    }


    gmp_randclear(rstate);

    {
        /* initialize x -- make it completely deterministic. */
        unsigned int my_nx = 4;
        uint32_t * xvecs = (uint32_t*) malloc(my_nx * bw->m * sizeof(uint32_t));
        /* with rhs, consider the strategy where the matrix is kept with
         * its full size, but the SM block is replaced with zeros. Here
         * we just force the x vectors to have data there, so as to avoid
         * the possibility that a solution vector found by BW still has
         * non-zero coordinates at these locations.
         */
        if (bw->m < (int) nrhs) {
            fprintf(stderr, "m < nrhs is not supported\n");
            exit(EXIT_FAILURE);
        }
        /* I am not sure that using balancing_pre_shuffle is right both
         * for bw->dir == 0 and bw->dir == 1. Let's make sure we're in
         * the case where this has been tested and seems to work
         * correctly.
         */
        ASSERT_ALWAYS(bw->dir == 1);
        ASSERT_ALWAYS(mmt.matrices.size() == 1);
        for(unsigned int i = 0 ; i < nrhs ; i++) {
            xvecs[i * my_nx] = balancing_pre_shuffle(mmt.matrices[0].bal, mmt.n0[!bw->dir]-nrhs+i);
            printf("Forced %d-th x vector to be the %" PRIu32"-th canonical basis vector\n", i, xvecs[i * my_nx]);
            ASSERT_ALWAYS(xvecs[i * my_nx] >= (uint32_t) (bw->m - nrhs));
            for(unsigned int j = 1 ; j < my_nx ; j++) {
                xvecs[i * my_nx + j] = (1009 * (i * my_nx + j)) % mmt.n0[!bw->dir];
            }
        }
        for(int i = (int) nrhs ; i < bw->m ; i++) {
            xvecs[i * my_nx] = i - nrhs;
            for(unsigned int j = 1 ; j < my_nx ; j++) {
                xvecs[i * my_nx + j] = (1009 * (i * my_nx + j)) % mmt.n0[!bw->dir];
            }
        }
        /* save_x operates only on the leader thread */
        save_x(xvecs, bw->m, my_nx, pi);
        free(xvecs);
    }

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
    /* declare local parameters and switches */
    param_list_decl_usage(pl, "rhs",
            "file with the right-hand side vectors for inhomogeneous systems mod p");

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
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (mpz_cmp_ui(bw->p, 2) == 0) {
        pi_go(prep_prog, pl, 0);
    } else {
        pi_go(prep_prog_gfp, pl, 0);
    }

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}

