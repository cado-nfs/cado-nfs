#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include "portability.h"
#include "utils.h"
#include "mf.h"
#include "balancing.h"
#include "rowset_heap.h"
#include "cheating_vec_init.h"
#include "mf_bal.h"

/* This program computes how a matrix would have to be balanced for
 * fitting a job grid of given size. This does _not_ read the matrix,
 * only the row-weight and col-weight files are read.
 */

/* TODO: Now that all files input and output by bwc are permutation
 * independent, there is room for computing this permutation on the fly
 * if the balancing argument is not provided (and only in this case).
 * Essentially we would hook onto the mmt_fill_fields_from_balancing()
 * function called from mmt_finish_init().
 */
int main(int argc, char * argv[])
{
    param_list pl;
    struct mf_bal_args mba[1];
    memset(mba, 0, sizeof(struct mf_bal_args));
    mba->do_perm[0] = MF_BAL_PERM_AUTO;
    mba->do_perm[1] = MF_BAL_PERM_AUTO;
    // int display_correlation = 0;

    param_list_init(pl);

    mf_bal_decl_usage(pl);
    mf_bal_configure_switches(pl, mba);
    mf_bal_parse_cmdline(mba, pl, &argc, &argv);
    mf_bal_interpret_parameters(mba, pl);

    mf_bal(mba);

    /* mba holds pointers to inside the pl struct, so beware ! */
    param_list_clear(pl);

    return 0;
}
