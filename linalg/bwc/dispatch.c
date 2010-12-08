#include "cado.h"

/* This is a very silly program which merely reads the matrix and then
 * exits. It must come before the prep program, and coalescing this one
 * and prep into one would be difficult and/or artificial, because they
 * handle different data widths.
 */
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "debug.h"
#include "params.h"
#include "bw-common-mpi.h"
#include "async.h"
// #include "rusage.h"
#include "filenames.h"

void * dispatch_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
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

    matmul_top_clear(mmt, abase);

    timing_clear(timing);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./dispatch <options>\n");
    fprintf(stderr, "%s", bw_common_usage_string());
    fprintf(stderr, "Relevant options here: wdir cfg m n mpi thr ys\n");
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
    pi_go(dispatch_prog, pl, 0);
    param_list_clear(pl);
    bw_common_clear_mpi(bw);

    return 0;
}

