#include "cado.h" // IWYU pragma: keep

#include <cstdio>   // for fflush, stdout
#include <cstring>  // for memset

#include "mf_bal.hpp"  // for mf_bal and friends.
#include "params.h"  // for param_list_clear, param_list_init, param_list_pr...

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
int main(int argc, char const * argv[])
{
    mf_bal_args mba;
    mba.do_perm[0] = mf_bal_args::MF_BAL_PERM_AUTO;
    mba.do_perm[1] = mf_bal_args::MF_BAL_PERM_AUTO;
    // int display_correlation = 0;

    cxx_param_list pl;

    mf_bal_decl_usage(pl);
    mf_bal_configure_switches(pl, &mba);
    mf_bal_parse_cmdline(&mba, pl, &argc, &argv);
    mf_bal_interpret_parameters(&mba, pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    mf_bal(&mba);

    return 0;
}
