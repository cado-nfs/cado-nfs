#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "finding_good_strategy.hpp"
#include "params.h"
#include "strategy.hpp"
#include "tab_strategy.hpp"
#include "macros.h"

/************************************************************************/
/*             USAGE                                                    */
/************************************************************************/

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage_header("This binary allow to choose the good strategy for each pair\n"
			    "\t\t of cofactor according to the distribution of these pairs\n"
			    "\t\t and the time required to establish our cofactors.\n");
    pl.declare_usage("st",
			  "the pathname of our directory which contains our strategies.");
    pl.declare_usage("dist",
			  "the pathname of our file which contains the distribution\n"
			  "\t\t of our pairs of cofactors.");
    pl.declare_usage("t",
			  "specify the time (seconds) to optain cofactors in the file\n"
			  "\t\t given by the option 'dist'.");
    pl.declare_usage("mfb0", "set the first cofactor bound to 2^mfb0");
    pl.declare_usage("mfb1",
			  "set the second cofactor bound to 2^mfb1");
    pl.declare_usage("out",
			  "the output file which contain final strategies.");

}

/************************************************************************/
/*             MAIN                                                     */
/************************************************************************/

int main(int argc, char const * argv[])
{
    cxx_param_list pl;
    declare_usage(pl);

    param_list_process_command_line_and_extra_parameter_files(pl, &argc, &argv);


    int mfb0 = -1, mfb1 = -1;
    double time_C = -1;
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_double(pl, "t", &time_C);
    if (mfb0 < 0 || mfb1 < 0 || time_C < 0)
	pl.fail("The following parameters are mandatory:\n"
	      "\t\t -mfb0 -mfb1 -time_C\n");

    printf("time_C = %lf seconds!\n", time_C);
    //convert the time in micro-s. because all previous binaries
    //compute their times in micro-s.
    time_C *= 1000000;		//s-->ms 

    const char *pathname_C;
    if ((pathname_C = param_list_lookup_string(pl, "dist")) == NULL)
	pl.fail("missing option -dist"
	      " followed by the pathname of the file which stores the"
	      " distribution of our cofactors.\n");

    const char *pathname_st;
    if ((pathname_st = param_list_lookup_string(pl, "st")) == NULL)
	pl.fail("missing option -st"
	      " followed by the pathname of the directory which"
	      " contains our strategies.\n");

    const char *pathname_output;
    if ((pathname_output = param_list_lookup_string(pl, "out")) == NULL)
	pl.fail("missing option -out"
	      " to specify where the result will be stored!\n");

    printf("EXTRACT DATA STRAT\n");
    tabular_strategy_t ***matrix_strat =
	extract_matrix_strat(pathname_st, mfb0 + 1, mfb1 + 1);
    if (matrix_strat == NULL) {
	exit(EXIT_FAILURE);
    }

    printf("EXTRACT DATA C\n");
    FILE *file_C = fopen(pathname_C, "r");
    DIE_ERRNO_DIAG(!file_C, "fopen(%s)", pathname_C);
    unsigned long **matrix_C = extract_matrix_C(file_C, mfb0 + 1, mfb1 + 1);
    if (matrix_C == NULL) {
	fprintf(stderr, "Error while reading file %s\n", pathname_C);
	exit(EXIT_FAILURE);
    }
    fclose(file_C);

    printf("GENERATE FINAL STRATEGIES\n");
    strategy_t ***matrix_strat_res =
	compute_best_strategy(matrix_strat, matrix_C, mfb0 + 1, mfb1 + 1,
			      time_C);

    FILE *file_output = fopen(pathname_output, "w");
    DIE_ERRNO_DIAG(!file_output, "fopen(%s)", pathname_output);
    int err = fprint_final_strategy(file_output, matrix_strat_res,
				    mfb0 + 1, mfb1 + 1);
    if (err == -1) {
	fprintf(stderr, "Error when i want to write in '%s'\n",
		pathname_output);
	exit(EXIT_FAILURE);
    }
    fclose(file_output);

    //free
    for (int r0 = 0; r0 <= mfb0; r0++) {
	for (int r1 = 0; r1 <= mfb1; r1++) {
	    tabular_strategy_free(matrix_strat[r0][r1]);
	    strategy_free(matrix_strat_res[r0][r1]);
	}
	free(matrix_strat[r0]);
	free(matrix_strat_res[r0]);
	free(matrix_C[r0]);
    }
    free(matrix_strat);
    free(matrix_strat_res);
    free(matrix_C);

    return EXIT_SUCCESS;
}
