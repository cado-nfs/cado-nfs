#include "cado.h"
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "bw-common.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "async.h"
#include "portability.h"

/* We merely have to build up a check vector. We've got a fairly easy
 * candidate for that: the x vector.  */

#if 0
void save_untwisted_transposed_vector(matmul_top_data_ptr mmt, const char * name, int d, int iter)
{
    mmt_comm_ptr mcol = mmt->wr[!d];
    mmt_comm_ptr mrow = mmt->wr[d];
    pi_comm_ptr picol = mmt->pi->wr[!d];
    pi_comm_ptr pirow = mmt->pi->wr[d];

    // v ---> Pr^-1*Sr^-1*v
    matmul_top_twist_vector(mmt, d);

    serialize_threads(mmt->pi->m);

    // v ---> Pr*v          [[ Sr^-1 * v ]]
    if (pirow->trank == 0 || !(mrow->v->flags & THREAD_SHARED_VECTOR)) {
        if (pirow->jrank == 0 && pirow->trank == 0) {
            /* ok, we keep the data */
        } else {
            abzero(mmt->abase, mrow->v->v, mrow->i1 - mrow->i0);
        }
    }
    matmul_top_mul_comm(mmt, !d);
    serialize_threads(mmt->pi->m);

    mmt_vec_save(mmt, NULL, name, !d, iter);
    // v --> Pr^-1 * v      [[ Pr^-1*Sr^-1*v ]]
    if (picol->trank == 0 || !(mcol->v->flags & THREAD_SHARED_VECTOR)) {
        if (picol->jrank == 0 && picol->trank == 0) {
            /* ok, we keep the data */
        } else {
            abzero(mmt->abase, mcol->v->v, mcol->i1 - mcol->i0);
        }
    }
    matmul_top_mul_comm(mmt, d);
    serialize_threads(mmt->pi->m);

    // v --> Sr*Pr*v        [[ v ]]
    matmul_top_untwist_vector(mmt, d);
    serialize_threads(mmt->pi->m);
}
#endif

void * sec_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;

    ASSERT_ALWAYS(!pi->interleaved);

    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    matmul_top_init(mmt, A, pi, pl, bw->dir);

    mmt_vec myy[2];
    mmt_vec_ptr my = myy[0];
    mmt_vec_ptr y = myy[1];

    /* we work in the opposite direction compared to other programs */
    mmt_vec_init(mmt,0,0, y,   bw->dir,                0, mmt->n[bw->dir]);
    mmt_vec_init(mmt,0,0, my, !bw->dir, /* shared ! */ 1, mmt->n[!bw->dir]);

    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    /* Because we're a special case, we _expect_ to work opposite to
     * optimized direction. So we pass bw->dir even though _we_ are going
     * to call mmt_mul with !bw->dir.
     */

    serialize(pi->m);

    if (bw->start == 0) {
        if (tcan_print)
            printf("We have start=0: creating C.0 as an expanded copy of X\n");
        uint32_t * gxvecs = NULL;
        unsigned int nx = 0;

        if (!fake) {
            load_x(&gxvecs, bw->m, &nx, pi);
        } else {
            set_x_fake(&gxvecs, bw->m, &nx, pi);
        }

        mmt_vec_set_x_indices(my, gxvecs, MIN(nchecks, bw->m), nx);
        mmt_vec_save(my, "C.0", unpadded);

        free(gxvecs);
    } else {
        char * tmp;
        int rc = asprintf(&tmp, "C.%d", bw->start);
        ASSERT_ALWAYS(rc >= 0);
        mmt_vec_load(my, tmp, unpadded);
        if (tcan_print) {
            printf("loaded %s\n", tmp);
        }
        free(tmp);
    }

    /*
    for(int j = 0 ; j < nchecks ; j++) {
        for(unsigned int k = 0 ; k < nx ; k++) {
            uint32_t i = gxvecs[j*nx+k];
            // set bit j of entry i to 1.
            if (i < mrow->i0 || i >= mrow->i1)
                continue;
            abt * where;
            where = mrow->v->v + aboffset(abase, i-mrow->i0);
            abset_ui(abase, where, j, 1);
        }
    }
    */
    // mmt_vec_allreduce(mmt, !bw->dir);

    // matmul_top_twist_vector(mmt, !bw->dir);
    // mmt_vec_save(mmt, NULL, "ux", !bw->dir, 0);
    // save_untwisted_transposed_vector(mmt, "tx", !bw->dir, 0);


    if (tcan_print) {
        printf("Computing trsp(x)*M^k for check stops k=");
        for(int s = 0 ; s < bw->number_of_check_stops ; s++) {
            int next = bw->check_stops[s];
            if (s) printf(",");
            printf("%d", next);
        }
        printf("\n");
    }

    serialize(pi->m);

    // kill the warning.
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        mmt->matrices[i]->mm->iteration[!bw->dir] = INT_MIN;
    }

    int k = bw->start;
    for(int s = 0 ; s < bw->number_of_check_stops ; s++) {
        int next = bw->check_stops[s];
        mmt_vec_twist(mmt, my);
        for( ; k < next ; k++) {
            pi_log_op(mmt->pi->m, "iteration %d", k);
            matmul_top_mul(mmt, myy, NULL);

            if (tcan_print) {
                putchar('.');
                fflush(stdout);
            }
        }
        serialize(pi->m);
        serialize_threads(mmt->pi->m);
        mmt_vec_untwist(mmt, my);
        {
            char * tmp;
            int rc = asprintf(&tmp, "C.%d", k);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(my, tmp, unpadded);
            if (tcan_print) {
                printf("saved %s\n", tmp);
            }
            free(tmp);
        }
    }

    mmt_vec_clear(mmt, y);
    mmt_vec_clear(mmt, my);
    matmul_top_clear(mmt);

    A->oo_field_clear(A);

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

    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(sec_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);
    return 0;
}

