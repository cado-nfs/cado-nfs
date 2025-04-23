#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "params.h"     // param_list
#include "facul_ecm.h"
#include "facul_method.hpp" // EC_METHOD ...
#include "generate_factoring_method.hpp"
#include "macros.h"
#include "tab_fm.hpp"

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to choose the best parameters (B1, B2 (B2=c*B1))\n"
			    "for each factoring methods: PM1, PP1-27, PP1-65, ECM-M12, ECM-M16, ECM-B12,\n"
			    "to find a prime number in [2^lb, 2^ub]. Note that you can : \n"
			    "-make your study only with one method (-m).\n"
			    "-specify the sieving region of B1 and c instead of using default values.\n"
			    "Then, You can apply the convex hull to keep only the best methods\n "
			    "(-fch from a file or just -ch if you want to apply it after the bench\n).");
    param_list_decl_usage(pl, "ub", "to set large prime bound to 2^ub.");
    param_list_decl_usage(pl, "lb", "to set factor base bound to 2^lb.");
    param_list_decl_usage(pl, "n", "to give the lenght of n.");
    param_list_decl_usage(pl, "out",
			  "to specify the file which contain our factoring methods.");
    param_list_decl_usage(pl, "m",
			  "to specify the method : PM1, PP1-27, PP1-65, ECM-M12, ECM-M16,\n ECM-B12. By default, we use all methods one after the other.");
    param_list_decl_usage(pl, "ch", "to apply the convex hull.");

    param_list_decl_usage(pl, "b1min", "to set b1_min (sieve region).");
    param_list_decl_usage(pl, "b1max", "to set b1_max (sieve region).");
    param_list_decl_usage(pl, "b1step", "to set b1_step (sieve region).");
    param_list_decl_usage(pl, "cmin", "to set c_min (sieve region).");
    param_list_decl_usage(pl, "cmax", "to set c_max (sieve region).");
    param_list_decl_usage(pl, "cstep", "to set c_step (sieve region).");

    param_list_decl_usage(pl, "fch", "to apply the convex hull");
    param_list_decl_usage(pl, "fch_in",
			  "to specify the input file which contains \n the factoring methods (default name :'default_fch_in').");
    param_list_decl_usage(pl, "fch_out",
			  "to specify the output file which contains \n the convex hull (default name :'default_fch_out').");

}

/************************************************************************/
/*                      MAIN */
/************************************************************************/

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    declare_usage(pl);
    /* 
       Passing nullptr is allowed here. Find value with
       param_list_parse_switch later on 
     */
    param_list_configure_switch(pl, "ch", nullptr);
    param_list_configure_switch(pl, "fch", nullptr);

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
	if ((f = fopen(argv[0], "r")) != nullptr) {
	    param_list_read_stream(pl, f, 0);
	    fclose(f);
	    argv++, argc--;
	    continue;
	}

	fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
	param_list_print_usage(pl, argv[0], stderr);
	exit(EXIT_FAILURE);
    }

    int const opt_fch = param_list_parse_switch(pl, "-fch");
    if (opt_fch) {
	const char *pathname_fch_in;
	const char *pathname_fch_out;
	if ((pathname_fch_in = param_list_lookup_string(pl, "fch_in")) == nullptr
	    || (pathname_fch_out =
		param_list_lookup_string(pl, "fch_out")) == nullptr) {
	    fputs("Parse error: Please re-run with the options "
		  "-in and -out with each one a valid file name.\n", stderr);
	    exit(EXIT_FAILURE);
	}

        auto file_in = fopen_helper(pathname_fch_in, "r");
        auto file_out = fopen_helper(pathname_fch_out, "w");

	tabular_fm_t *res_ch = convex_hull_from_file(file_in.get(), file_out.get());
	if (res_ch == nullptr) {
	    fprintf(stderr, "impossible to read %s\n"
		    "impossible to write in the file %s\n",
		    pathname_fch_in, pathname_fch_out);
	    exit(EXIT_FAILURE);
	}
	tabular_fm_free(res_ch);
    } else {
	//default values
	int lb = -1, ub = -1, len_n = -1;
        ec_parameterization_t curve = MONTY12;
	// {b1min, b1max, b1step, cmin, cmax, cstep}
        std::vector<int> param(6, 0);

	int const opt_ch = param_list_parse_switch(pl, "-ch");
	param_list_parse_int(pl, "ub", &ub);
	param_list_parse_int(pl, "lb", &lb);
	param_list_parse_int(pl, "n", &len_n);
	param_list_parse_int(pl, "b1min", &param[0]);
	param_list_parse_int(pl, "b1max", &param[1]);
	param_list_parse_int(pl, "b1step", &param[2]);
	param_list_parse_int(pl, "cmin", &param[3]);
	param_list_parse_int(pl, "cmax", &param[4]);
	param_list_parse_int(pl, "cstep", &param[5]);

	if (lb == -1 || ub == -1 || ub <= lb) {
	    fprintf(stderr, "Error: options -lb, -ub are mandatory here,\n"
		    "and lb must be less than ub.\n\n");
	    param_list_print_usage(pl, argv[0], stderr);
	    exit(EXIT_FAILURE);
	}

	/*
	   Default value for len_n:
	   So n will be a composite such that n=pq, with :
	   size of p equal to len_p bits and for q 3*len_p bits.
	   Thus, q will not be a problem when we are going to p!
	 */
	if (len_n == -1)
	    len_n = 60 + ub;

	/* Extract the method and the curve\n" */
	const char *name_fm = param_list_lookup_string(pl, "m");

        facul_method_code method = NO_METHOD;

	if (name_fm != nullptr) {
	    if (strcmp(name_fm, "PM1") == 0)
		method = PM1_METHOD;
	    else if (strcmp(name_fm, "PP1-27") == 0)
		method = PP1_27_METHOD;
	    else if (strcmp(name_fm, "PP1-65") == 0)
		method = PP1_65_METHOD;
	    else if (strcmp(name_fm, "ECM-M12") == 0) {
		method = EC_METHOD;
		curve = MONTY12;
	    } else if (strcmp(name_fm, "ECM-B12") == 0) {
		method = EC_METHOD;
		curve = BRENT12;
	    } else if (strcmp(name_fm, "ECM-M16") == 0) {
		method = EC_METHOD;
		curve = MONTY16;
	    } else {
		fprintf(stderr,
			"Your factoring method '%s' doesn't exist.\n"
			"Choose an other method from the below file:\n"
			"PM1, PP1-27, PP1-65, ECM-M12, ECM-M16, ECM-B12\n",
			name_fm);
		exit(EXIT_FAILURE);
	    }
	}

	cxx_gmp_randstate state;
	/* Initializing radom generator */
	gmp_randseed_ui(state, time(nullptr));

	/*
	   To generate our factoring methods!
	 */

	tabular_fm_t *res;

	int *param_sieve = nullptr;
	if (param[0] && param[1] && param[2] &&
	    param[3] && param[4] && param[5])
        {
	    param_sieve = param.data();
            /*
	    printf("param_sieve: b1min= %d, b1max=%d, b1step=%d"
		   "       cmin= %d, cmax=%d, cstep=%d\n",
		   param_sieve[0], param_sieve[1], param_sieve[2],
		   param_sieve[3], param_sieve[4], param_sieve[5]);
                   */
	}

	if (method == NO_METHOD)
	    res = generate_factoring_methods
		(state, lb, ub, len_n, opt_ch, param_sieve);
	else {
	    res = generate_factoring_methods_mc
		(state, lb, ub, len_n, method, curve, opt_ch, param_sieve);
	}

	//print the result in file!
	//to do: change the file output!
	const char *name_file_out = param_list_lookup_string(pl, "out");


        /* if -out is not given, print to stdout */
	FILE *file_out = (name_file_out == nullptr)
          ? stdout : fopen(name_file_out, "w");
        DIE_ERRNO_DIAG(file_out == nullptr, "fopen(%s)", name_file_out);
	int const err = tabular_fm_fprint(file_out, res);
	if (err < 0) {
	    fprintf(stderr,
		    "error:: try to write in the file '%s'.\n", name_file_out);
	    exit(EXIT_FAILURE);
	}
	fclose(file_out);
	tabular_fm_free(res);
    }

    return EXIT_SUCCESS;
}
