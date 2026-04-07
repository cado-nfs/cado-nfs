#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <gmp.h>

#include "tab_fm.hpp"
#include "generate_factoring_method.hpp"
#include "params.h"

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage_header("This binary allows to do a bench on the probabilities and/or the times\n"
			    "for each size of prime number from 'lb'.\n");

    pl.declare_usage("lb",
			  "to begin the benchmark with prime numbers of 'lb' bits.");
    pl.declare_usage("p", "to bench the probabilities.");
    pl.declare_usage("t", "to bench the times.");
    pl.declare_usage("N", "number of tests to run for each bench");
    pl.declare_usage("in",
			  "to locate the file which contains our factoring methods.");
    pl.declare_usage("f",
			  "to keep only a number of factoring methods.");
    pl.declare_usage("out",
			  "to locate the file which contains our benchmark.");

}

/************************************************************************/
/*                      MAIN */
/************************************************************************/

// coverity[root_function]
int main(int argc, char const * argv[])
{
    int nb_test = 0;
    cxx_param_list pl;
    declare_usage(pl);
    /* 
       Passing NULL is allowed here. Find value with
       param_list_parse_switch later on 
     */
    pl.configure_switch("p");
    pl.configure_switch("t");

    param_list_process_command_line_and_extra_parameter_files(pl, &argc, &argv);

    //default values
    int len_p_min = -1; //default_value
    int final_nb_fm = -1; //default value
    
    int const opt_proba = param_list_parse_switch(pl, "-p");
    int const opt_time = param_list_parse_switch(pl, "-t");
    param_list_parse_int(pl, "lb", &len_p_min);
    param_list_parse_int(pl, "f", &final_nb_fm);
    param_list_parse_int(pl, "N", &nb_test);

    const char *pathname_in;
    const char *pathname_out;
    if ((pathname_in = param_list_lookup_string(pl, "in")) == NULL)
        pl.fail("missing argument -in");
    if ((pathname_out = param_list_lookup_string(pl, "out")) == NULL)
        pl.fail("missing argument -out");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    FILE *file_in = fopen(pathname_in, "r");
    tabular_fm_t *c = tabular_fm_fscan(file_in);
    if (c == NULL) {
	fprintf(stderr, "impossible to read %s\n", pathname_in);
	exit(EXIT_FAILURE);
    }
    fclose(file_in);

    if (opt_proba)
	{
	    if (len_p_min == -1)
                pl.fail("missing argument -lb");
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
    int const err = tabular_fm_fprint(file_out, c);
    if (err < 0) {
	fprintf(stderr, "error:: try to write in the file %s.\n", pathname_out);
	exit(EXIT_FAILURE);
    }
    fclose(file_out);

    //free
    tabular_fm_free(c);

    gmp_randclear(state);

    return EXIT_SUCCESS;
}
