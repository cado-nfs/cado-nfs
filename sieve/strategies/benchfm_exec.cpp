#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <gmp.h>

#include "tab_fm.h"     // tabular_fm_free
#include "generate_factoring_method.hpp"
#include "params.h"     // param_list

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to do a bench on the probabilities and/or the times\n"
			    "for each size of prime number from 'lb'.\n");

    param_list_decl_usage(pl, "lb",
			  "to begin the benchmark with prime numbers of 'lb' bits.");
    param_list_decl_usage(pl, "p", "to bench the probabilities.");
    param_list_decl_usage(pl, "t", "to bench the times.");
    param_list_decl_usage(pl, "N", "number of tests to run for each bench");
    param_list_decl_usage(pl, "in",
			  "to locate the file which contains our factoring methods.");
    param_list_decl_usage(pl, "f",
			  "to keep only a number of factoring methods.");
    param_list_decl_usage(pl, "out",
			  "to locate the file which contains our benchmark.");

}

/************************************************************************/
/*                      MAIN */
/************************************************************************/

// coverity[root_function]
int main(int argc, char *argv[])
{
    int nb_test = 0;
    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    /* 
       Passing NULL is allowed here. Find value with
       param_list_parse_switch later on 
     */
    param_list_configure_switch(pl, "p", NULL);
    param_list_configure_switch(pl, "t", NULL);

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

    //default values
    int len_p_min = -1; //default_value
    int final_nb_fm = -1; //default value
    
    int opt_proba = param_list_parse_switch(pl, "-p");
    int opt_time = param_list_parse_switch(pl, "-t");
    param_list_parse_int(pl, "lb", &len_p_min);
    param_list_parse_int(pl, "f", &final_nb_fm);
    param_list_parse_int(pl, "N", &nb_test);

    const char *pathname_in;
    const char *pathname_out;
    if ((pathname_in = param_list_lookup_string(pl, "in")) == NULL ||
	(pathname_out = param_list_lookup_string(pl, "out")) == NULL)
	{
	    fputs("Parse error: Please re-run with options "
		  "-in and -out with each one a valid file name.\n", stderr);
	    param_list_clear(pl);
	    exit(EXIT_FAILURE);
	}

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    FILE *file_in = fopen(pathname_in, "r");
    tabular_fm_t *c = tabular_fm_fscan(file_in);
    if (c == NULL) {
	fprintf(stderr, "impossible to read %s\n", pathname_in);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    fclose(file_in);

    if (opt_proba)
	{
	    if (len_p_min == -1)
		{
		    fprintf(stderr, "error: -lb is missing.\n");
		    param_list_clear(pl);
		    exit(EXIT_FAILURE);
		}
	    bench_proba(state, c, len_p_min, 0, nb_test);
	}
    if (opt_time)
	bench_time(state, c, nb_test);

    if (final_nb_fm != -1)
	{
	    tabular_fm_t* res = filtering (c, final_nb_fm);
	    tabular_fm_free (c);
	    c = res;
	}
    
    FILE *file_out = fopen(pathname_out, "w");
    int err = tabular_fm_fprint(file_out, c);
    if (err < 0) {
	fprintf(stderr, "error:: try to write in the file %s.\n", pathname_out);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    fclose(file_out);

    //free
    param_list_clear(pl);
    tabular_fm_free(c);

    gmp_randclear(state);

    return EXIT_SUCCESS;
}
