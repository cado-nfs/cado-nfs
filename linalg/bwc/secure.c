#include "cado.h"
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "bw-common-mpi.h"
#include "filenames.h"

/* We merely have to build up a check vector. We've got a fairly easy
 * candidate for that: the x vector.  */

void * sec_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    /* Interleaving does not make sense for this program. So the second
     * block of threads just leave immediately */
    if (pi->interleaved && pi->interleaved->idx)
        return NULL;

    abobj_t abase;
    abobj_init(abase);
    abobj_set_nbys(abase, NCHECKS_CHECK_VECTOR);

    matmul_top_data mmt;

    /* XXX Here we're roking in the opposite direction compared to
     * prep/krylov/mksol ! */
    int flags[2];
    flags[!bw->dir] = THREAD_SHARED_VECTOR;
    flags[bw->dir] = 0;

    int tcan_print = bw->can_print && pi->m->trank == 0;

    /* Because we're a special case, we _expect_ to work opposite to
     * optimized direction. So we pass bw->dir even though _we_ are going
     * to call mmt_mul with !bw->dir.
     */
    matmul_top_init(mmt, abase, pi, flags, pl, bw->dir);

    // mmt_wiring_ptr mcol = mmt->wr[bw->dir];
    mmt_wiring_ptr mrow = mmt->wr[!bw->dir];

    abzero(abase, mrow->v->v, mrow->i1 - mrow->i0);

    uint32_t * gxvecs = malloc(bw->nx * bw->m * sizeof(uint32_t));

    load_x(gxvecs, bw->m, bw->nx, pi, mmt->bal);

    ASSERT_ALWAYS(bw->m >= NCHECKS_CHECK_VECTOR);
    for(int j = 0 ; j < NCHECKS_CHECK_VECTOR ; j++) {
        for(unsigned int k = 0 ; k < bw->nx ; k++) {
            // set bit j of entry gxvecs[j*bw->nx+k] to 1.
            uint32_t i = gxvecs[j*bw->nx+k];
            if (i < mrow->i0 || i >= mrow->i1)
                continue;
            abt * where;
            where = mrow->v->v + aboffset(abase, i-mrow->i0);
            abset_ui(abase, where, j, 1);
        }
    }

    if (tcan_print) {
        printf("Computing trsp(x)*M^%d\n",bw->interval);
    }

    serialize(pi->m);

#if 0
    /* FIXME -- that's temporary ! only for debugging */
    pi_log_init(pi->m);
    pi_log_init(pi->wr[0]);
    pi_log_init(pi->wr[1]);
#endif

    // kill the warning.
    mmt->mm->iteration[!bw->dir] = INT_MIN;

    for(int k = 0 ; k < bw->interval ; k++) {
        pi_log_op(mmt->pi->m, "iteration %d", k);
        matmul_top_mul(mmt, !bw->dir);
        if (tcan_print) {
            putchar('.');
            fflush(stdout);
        }
    }
    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_save_vector(mmt, CHECK_FILE_BASE, !bw->dir, bw->interval);

    matmul_top_clear(mmt, abase);

    free(gxvecs);

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./secure <options>\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr matrix interval\n");
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

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);


    pi_go(sec_prog, pl, 0);

    param_list_clear(pl);
    bw_common_clear_mpi(bw);
    return 0;
}

