#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>

#include "finding_good_strategy.h"
#include "params.h"     // param_list
#include "strategy.h"
#include "tab_strategy.h"

/************************************************************************/
/*             USAGE                                                    */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allow to choose the good strategy for each pair\n"
			    "\t\t of cofactor according to the distribution of these pairs\n"
			    "\t\t and the time required to establish our cofactors.\n");
    param_list_decl_usage(pl, "st",
			  "the pathname of our directory which contains our strategies.");
    param_list_decl_usage(pl, "dist",
			  "the pathname of our file which contains the distribution\n"
			  "\t\t of our pairs of cofactors.");
    param_list_decl_usage(pl, "t",
			  "specify the time (seconds) to optain cofactors in the file\n"
			  "\t\t given by the option 'dist'.");
    param_list_decl_usage(pl, "mfb0", "set the first cofactor bound to 2^mfb0");
    param_list_decl_usage(pl, "mfb1",
			  "set the second cofactor bound to 2^mfb1");
    param_list_decl_usage(pl, "out",
			  "the output file which contain final strategies.");

}

/************************************************************************/
/*             MAIN                                                     */
/************************************************************************/

int main(int argc, char *argv[])
{
    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    /* 
       Passing NULL is allowed here. Find value with
       param_list_parse_switch later on 
     */
    if (argc <= 1) {
	param_list_print_usage(pl, argv[0], stderr);
	exit(EXIT_FAILURE);
    }

    argv++, argc--;
    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv)) {
	    continue;
	}
	/* Could also be a file */
	FILE *f;
	if ((f = fopen(argv[0], "r")) != NULL) {
	    param_list_read_stream(pl, f, 0);
	    fclose(f);
	    argv++, argc--;
	    continue;
	}

	fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
	param_list_print_usage(pl, argv[0], stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    int mfb0 = -1, mfb1 = -1;
    double time_C = -1;
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_double(pl, "t", &time_C);
    if (mfb0 < 0 || mfb1 < 0 || time_C < 0) {
	fputs("The following parameters are mandatory:\n"
	      "\t\t -mfb0 -mfb1 -time_C\n", stderr);

	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    printf("time_C = %lf seconds!\n", time_C);
    //convert the time in micro-s. because all previous binaries
    //compute their times in micro-s.
    time_C *= 1000000;		//s-->ms 

    const char *pathname_C;
    if ((pathname_C = param_list_lookup_string(pl, "dist")) == NULL) {
	fputs("Parser error: Please re-run with the option -dist"
	      "followed by the pathname of the file which stores the "
	      "distribution of our cofactors.\n", stderr);
	exit(EXIT_FAILURE);
    }

    const char *pathname_st;
    if ((pathname_st = param_list_lookup_string(pl, "st")) == NULL) {
	fputs("Parser error: Please re-run with the option -st"
	      "followed by the pathname of the directory which"
	      " contains our strategies.\n", stderr);
	exit(EXIT_FAILURE);
    }

    const char *pathname_output;
    if ((pathname_output = param_list_lookup_string(pl, "out")) == NULL) {
	fputs("Parser error: Please re-run with the option -out "
	      "to specify where the result will be stored!\n", stderr);
	exit(EXIT_FAILURE);
    }

    printf("EXTRACT DATA STRAT\n");
    tabular_strategy_t ***matrix_strat =
	extract_matrix_strat(pathname_st, mfb0 + 1, mfb1 + 1);
    if (matrix_strat == NULL) {
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    printf("EXTRACT DATA C\n");
    FILE *file_C = fopen(pathname_C, "r");
    unsigned long **matrix_C = extract_matrix_C(file_C, mfb0 + 1, mfb1 + 1);
    if (matrix_C == NULL) {
	fprintf(stderr, "Error while reading file %s\n", pathname_C);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    fclose(file_C);

    printf("GENERATE FINAL STRATEGIES\n");
    strategy_t ***matrix_strat_res =
	compute_best_strategy(matrix_strat, matrix_C, mfb0 + 1, mfb1 + 1,
			      time_C);

    FILE *file_output = fopen(pathname_output, "w");
    int err = fprint_final_strategy(file_output, matrix_strat_res,
				    mfb0 + 1, mfb1 + 1);
    if (err == -1) {
	fprintf(stderr, "Error when i want to write in '%s'\n",
		pathname_output);
	param_list_clear(pl);
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

    param_list_clear(pl);

    return EXIT_SUCCESS;
}
