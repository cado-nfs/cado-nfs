#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>

#include "facul_ecm.h"
#include "facul_method.hpp"
#include "fm.hpp"
#include "gen_decomp.hpp"
#include "generate_strategies.hpp"
#include "macros.h"
#include "params.h" // param_list
#include "tab_decomp.hpp"
#include "tab_fm.hpp"
#include "tab_strategy.hpp"

/************************************************************************/
/*                            USAGE                                     */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(
        pl, "This binary allows to build the best strategies for each couple "
            "(r0, r1)\n"
            "where (r0,r1) are the bits size for our couple of cofactors.\n");

    param_list_decl_usage(pl, "gdc",
                          "to precompute all decompositions of cofactors of "
                          "mfb bits given that \n "
                          "\t \t it has no prime divisors less than lim. So, "
                          "you must specify these options:\n"
                          "\t \t -lim0, -mfb0\n");
    param_list_decl_usage(pl, "gst_r",
                          "to precompute the best strategies \n"
                          "\t \t for one bit size cofactor.\n "
                          "\t \t You must specify these options:\n"
                          "\t \t -lim0, -lpb0, -r0, -ncurves, -decomp\n");
    param_list_decl_usage(
        pl, "gst",
        "to merge two (or all) precomputing did by the option 'gst_r',\n "
        "\t \t and thus find the best strategie(s) for one (or each) couple "
        "(r0,r1).\n"
        "\t \t So, you must specify these options:\n"
        "\t \t -lim0, -lim1 ,-in, and ((-r0 -r1) or (-mfb0, -mfb1))\n");

    param_list_decl_usage(pl, "ncurves", "controls number of curves.\n");
    param_list_decl_usage(pl, "r0",
                          "set the bit size of the studied cofactor to r0.\n");
    param_list_decl_usage(pl, "r1",
                          "set the bit size of the studied cofactor to r1.\n");
    param_list_decl_usage(pl, "lim0",
                          "set rationnal factor base bound to lim0\n");
    param_list_decl_usage(pl, "lim1",
                          "set algebraic factor base bound to lim1\n");
    param_list_decl_usage(pl, "lpb0",
                          "set rational large prime bound to 2^lpb0");
    param_list_decl_usage(pl, "lpb1",
                          "set algebraic large prime bound to 2^lpb1");
    param_list_decl_usage(pl, "mfb0", "set the first cofactor bound to 2^mfb0");
    param_list_decl_usage(pl, "mfb1",
                          "set the second cofactor bound to 2^mfb1");
    param_list_decl_usage(
        pl, "in",
        "to locate the file which contains\n "
        "\t \t our factoring methods, or locate the directory \n"
        "\t \t where the precomputed files for option 'gst' are stored");
    param_list_decl_usage(pl, "out",
                          "to locate the directory where the "
                          "file(s) will be stored.");
    param_list_decl_usage(
        pl, "decomp",
        "to locate the file or the directory , according to\n"
        "\t \t if you need one or several files,\n"
        "\t \t which contain(s) the file(s) of cofactors decompositions.");
}

/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;
    declare_usage(pl);
    /*
       Passing NULL is allowed here. Find value with
       param_list_parse_switch later on
     */
    param_list_configure_switch(pl, "gdc", nullptr);
    param_list_configure_switch(pl, "gst", nullptr);
    param_list_configure_switch(pl, "gst_r", nullptr);

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
        FILE * f;
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

    /*default values */
    unsigned long lim0 = 0;
    unsigned long lim1 = 0;
    int mfb1 = -1;
    int mfb0 = -1;
    int lpb0 = -1;
    int lpb1 = -1;
    int ncurves = -1;

    param_list_parse_ulong(pl, "lim0", &lim0);
    param_list_parse_ulong(pl, "lim1", &lim1);
    param_list_parse_int(pl, "lpb0", &lpb0);
    param_list_parse_int(pl, "lpb1", &lpb1);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_int(pl, "ncurves", &ncurves);

    int const gdc = param_list_parse_switch(pl, "-gdc");
    int const gst_r = param_list_parse_switch(pl, "-gst_r");
    int const gst = param_list_parse_switch(pl, "-gst");

    // and precompute just for one side!
    if (gdc) { // precompute all decompositions!

        if (lim0 == 0) {
            fprintf(stderr, "Error: parameter -lim0 is mandatory\n");
            param_list_print_usage(pl, argv[0], stderr);
            exit(EXIT_FAILURE);
        }
        if (mfb0 == -1) {
            fprintf(stderr, "Error: parameter -mfb0 is mandatory\n");
            param_list_print_usage(pl, argv[0], stderr);
            exit(EXIT_FAILURE);
        }

        char const * file_out = param_list_lookup_string(pl, "out");

        FILE * file = stdout;
        if (file_out) {
            file = fopen(file_out, "w");
            DIE_ERRNO_DIAG(!file, "fopen(%s)", file_out);
        }

        tabular_decomp res = generate_all_decomp(mfb0, lim0);
        fmt::print(file, "{}", res);
        if (file_out)
            fclose(file);

    } else {

        // store data
        char const * directory_out;
        if ((directory_out = param_list_lookup_string(pl, "out")) == nullptr) {
            directory_out = "./";
        }

        if (gst) {
            // check parameters!
            char const * directory_in;
            if ((directory_in = param_list_lookup_string(pl, "in")) ==
                nullptr) {
                fprintf(stderr, "Error: parameter -in is mandatory\n");
                param_list_print_usage(pl, argv[0], stderr);
                exit(EXIT_FAILURE);
            }
            if (lim0 == 0) {
                fprintf(stderr, "Error: parameter -lim0 is mandatory\n");
                param_list_print_usage(pl, argv[0], stderr);
                exit(EXIT_FAILURE);
            }
            if (lim1 == 0) {
                fprintf(stderr, "Error: parameter -lim1 is mandatory\n");
                param_list_print_usage(pl, argv[0], stderr);
                exit(EXIT_FAILURE);
            }

            int r0 = -1, r1 = -1;
            param_list_parse_int(pl, "r0", &r0);
            param_list_parse_int(pl, "r1", &r1);

            if (r0 != -1 && r1 != -1) {
                // just one couple will be studied!
                char name_file_in[strlen(directory_in) + 64];
                FILE * file_in;

                // get back the best strategies for r0!
                snprintf(name_file_in, sizeof(name_file_in),
                         "%s/strategies%lu_%d", directory_in, lim0, r0);
                file_in = fopen(name_file_in, "r");
                DIE_ERRNO_DIAG(!file_in, "fopen(%s)", name_file_in);
                tabular_strategy_t * strat_r0 = tabular_strategy_fscan(file_in);
                ASSERT_ALWAYS(strat_r0);
                fclose(file_in);

                // get back the best strategies for r1!
                snprintf(name_file_in, sizeof(name_file_in),
                         "%s/strategies%lu_%d", directory_in, lim1, r1);
                file_in = fopen(name_file_in, "r");
                DIE_ERRNO_DIAG(!file_in, "fopen(%s)", name_file_in);
                tabular_strategy_t * strat_r1 = tabular_strategy_fscan(file_in);
                ASSERT_ALWAYS(strat_r1);
                fclose(file_in);

                // compute the best strategies for (r0, r1)!
                char res_file[200];
                tabular_strategy_t * res =
                    generate_strategy_r0_r1(strat_r0, strat_r1);
                // print
                snprintf(res_file, sizeof(res_file), "%s/strategies_%d_%d",
                         directory_out, r0, r1);
                FILE * file = fopen(res_file, "w");
                DIE_ERRNO_DIAG(!file, "fopen(%s)", res_file);
                tabular_strategy_fprint(file, res);
                fclose(file);
                // free
                tabular_strategy_free(strat_r0);
                tabular_strategy_free(strat_r1);
                tabular_strategy_free(res);

            } else {
                /* For the both side and each size of cofactor the best
                   strategies had been choosed. Now, it's time to merge these
                   datas to compute the best strategies for each couple (r0,
                   r1)!
                */

                if (mfb0 == 0) {
                    fprintf(stderr, "Error: parameter -mfb0 is mandatory\n");
                    param_list_print_usage(pl, argv[0], stderr);
                    exit(EXIT_FAILURE);
                }
                if (mfb1 == 0) {
                    fprintf(stderr, "Error: parameter -mfb1 is mandatory\n");
                    param_list_print_usage(pl, argv[0], stderr);
                    exit(EXIT_FAILURE);
                }
                auto * data_rat = new tabular_strategy_t *[mfb0 + 1];

                char name_file_in[strlen(directory_in) + 64];
                for (int r0 = 0; r0 <= mfb0; r0++) {
                    snprintf(name_file_in, sizeof(name_file_in),
                             "%s/strategies%lu_%d", directory_in, lim0, r0);
                    FILE * file_in = fopen(name_file_in, "r");
                    DIE_ERRNO_DIAG(!file_in, "fopen(%s)", name_file_in);
                    data_rat[r0] = tabular_strategy_fscan(file_in);
                    ASSERT_ALWAYS(data_rat[r0]);
                    fclose(file_in);
                }

                for (int r1 = 0; r1 <= mfb1; r1++) {
                    snprintf(name_file_in, sizeof(name_file_in),
                             "%s/strategies%lu_%d", directory_in, lim1, r1);
                    FILE * file_in = fopen(name_file_in, "r");
                    DIE_ERRNO_DIAG(!file_in, "fopen(%s)", name_file_in);
                    tabular_strategy_t * strat_r1 =
                        tabular_strategy_fscan(file_in);
                    ASSERT_ALWAYS(strat_r1);
                    fclose(file_in);

                    for (int r0 = 0; r0 <= mfb0; r0++) {
                        char res_file[200];
                        tabular_strategy_t * res =
                            generate_strategy_r0_r1(data_rat[r0], strat_r1);
                        // print
                        snprintf(res_file, sizeof(res_file),
                                 "%s/strategies_%d_%d", directory_out, r0, r1);
                        FILE * file = fopen(res_file, "w");
                        DIE_ERRNO_DIAG(!file, "fopen(%s)", res_file);
                        tabular_strategy_fprint(file, res);
                        fclose(file);

                        tabular_strategy_free(res);
                    }
                    tabular_strategy_free(strat_r1);
                }
                for (int r0 = 0; r0 <= mfb0; r0++)
                    tabular_strategy_free(data_rat[r0]);
                delete[] data_rat;
            }
        } else {
            char const * name_file_in;
            if ((name_file_in = param_list_lookup_string(pl, "in")) ==
                nullptr) {
                fprintf(stderr, "Error: parameter -in is mandatory\n");
                param_list_print_usage(pl, argv[0], stderr);
                exit(EXIT_FAILURE);
            }

            FILE * file_in = fopen(name_file_in, "r");
            DIE_ERRNO_DIAG(!file_in, "fopen(%s)", name_file_in);
            tabular_fm_t * c = tabular_fm_fscan(file_in);
            ASSERT_ALWAYS(c);
            fclose(file_in);

            tabular_fm_t * data_pp1_27 = extract_fm_method(c, PP1_27_METHOD, 0);
            tabular_fm_t * data_pp1_65 = extract_fm_method(c, PP1_65_METHOD, 0);

            tabular_fm_t * data_pm1 = extract_fm_method(c, PM1_METHOD, 0);

            tabular_fm_t * data_ecm_m16 =
                extract_fm_method(c, EC_METHOD, MONTY16);

            tabular_fm_t * data_ecm_m12 =
                extract_fm_method(c, EC_METHOD, MONTY12);

            tabular_fm_t * data_ecm_b12 =
                extract_fm_method(c, EC_METHOD, BRENT12);

            tabular_fm_t * data_pp1 = tabular_fm_create();
            tabular_fm_concat(data_pp1, data_pp1_27);
            tabular_fm_concat(data_pp1, data_pp1_65);

            // add zero method if doesn't exist!
            fm_t * zero = fm_create();
            unsigned long method_zero[4] = {0, 0, 0, 0};
            fm_set_method(zero, method_zero, 4);

            if (data_pm1->tab[0]->method[2] != 0) // B1=0
                tabular_fm_add_fm(data_pm1, zero);
            if (data_pp1->tab[0]->method[2] != 0) // B1=0
                tabular_fm_add_fm(data_pp1, zero);
            if (data_ecm_b12->tab[0]->method[2] != 0) // B1=0
                tabular_fm_add_fm(data_ecm_b12, zero);
            if (data_ecm_m12->tab[0]->method[2] != 0) // B1=0
                tabular_fm_add_fm(data_ecm_m12, zero);
            if (data_ecm_m16->tab[0]->method[2] != 0) // B1=0
                tabular_fm_add_fm(data_ecm_m16, zero);

            tabular_fm_t * data_ecm = tabular_fm_create();
            tabular_fm_concat(data_ecm, data_ecm_m16);
            tabular_fm_concat(data_ecm, data_ecm_b12);
            tabular_fm_concat(data_ecm, data_ecm_m12);

            tabular_fm_sort(data_pm1);
            tabular_fm_sort(data_pp1);
            tabular_fm_sort(data_ecm);

            if (gst_r) {

                int r0 = -1;
                param_list_parse_int(pl, "r0", &r0);
                if (r0 == -1 || lpb0 == -1 || lim0 == 0 || ncurves == -1) {
                    fprintf(stderr, "Error: parameters -r0 -lim0 "
                                    "-lpb0 -ncurves are mandatories.\n");
                    param_list_print_usage(pl, argv[0], stderr);
                    exit(EXIT_FAILURE);
                }

                // precompute the convex hull for one bit size of cofactor
                int const fbb0 = ceil(log2((double)(lim0 + 1)));
                int const lim_is_prime = 2 * fbb0 - 1;
                tabular_decomp tab_decomp;
                if (r0 >= lim_is_prime) {
                    char const * name_file_decomp;
                    if ((name_file_decomp = param_list_lookup_string(
                             pl, "decomp")) == nullptr) {

                        fprintf(stderr,
                                "Error: parameter -decomp is mandatory\n");
                        param_list_print_usage(pl, argv[0], stderr);
                        exit(EXIT_FAILURE);
                    }
                    std::ifstream is(name_file_decomp);
                    is >> tab_decomp;
                    if (is.fail())
                        ASSERT_ALWAYS(0);
                }

                tabular_strategy_t * res = generate_strategies_oneside(
                    tab_decomp, zero, data_pm1, data_pp1, data_ecm, ncurves,
                    lim0, lpb0, r0);

                char name_file[strlen(directory_out) + 64];
                snprintf(name_file, sizeof(name_file), "%s/strategies%lu_%d",
                         directory_out, lim0, r0);
                FILE * file_out = fopen(name_file, "w");
                DIE_ERRNO_DIAG(!file_out, "fopen(%s)", name_file);
                tabular_strategy_fprint(file_out, res);
                fclose(file_out);
                tabular_strategy_free(res);
                fm_free(zero);

            } else {
                /*
                  default parameter: this function generates our best
                  strategies for each couple (r0, r1) without to have
                  precomputed the best strategies for each side and bit
                  size r(0/1).
                 */
                // Check parameters

                if (lim0 == 0 || lpb0 == -1 || mfb0 == -1 || lim1 == 0 ||
                    lpb1 == -1 || mfb1 == -1 || ncurves) {
                    fprintf(stderr, "Error: parameters -lim0 -lpb0 -mfb0"
                                    " -lim1 -lpb1 -mfb1 are mandatories\n");
                    param_list_print_usage(pl, argv[0], stderr);
                    exit(EXIT_FAILURE);
                }
                /*
                  computes the matrix of strategies where for each couple
                  (r0, r1) we will optain the best strategies to find a
                  relation.
                */
                char const * name_directory_decomp;
                if ((name_directory_decomp =
                         param_list_lookup_string(pl, "decomp")) == nullptr) {
                    fprintf(stderr, "Error: parameter -decomp is mandatory\n");
                    param_list_print_usage(pl, argv[0], stderr);
                    exit(EXIT_FAILURE);
                }

                tabular_strategy_t *** matrix = generate_matrix(
                    name_directory_decomp, data_pm1, data_pp1, data_ecm,
                    ncurves, lim0, lpb0, mfb0, lim1, lpb1, mfb1);

                char res_file[strlen(directory_out) + 64];
                for (int r0 = 0; r0 <= mfb0; r0++)
                    for (int r1 = 0; r1 <= mfb1; r1++) {
                        snprintf(res_file, sizeof(res_file),
                                 "%s/strategies_%d_%d", directory_out, r0, r1);
                        FILE * file = fopen(res_file, "w");
                        DIE_ERRNO_DIAG(!file, "fopen(%s)", res_file);
                        tabular_strategy_fprint(file, matrix[r0][r1]);
                        fclose(file);
                    }

                // free
                for (int r0 = 0; r0 <= mfb0; r0++) {
                    for (int r1 = 0; r1 <= mfb1; r1++)
                        tabular_strategy_free(matrix[r0][r1]);
                    delete[] matrix[r0];
                }
                delete[] matrix;
            }
            // free
            tabular_fm_free(data_pp1_65);
            tabular_fm_free(data_pp1_27);
            tabular_fm_free(data_pp1);
            tabular_fm_free(data_pm1);
            tabular_fm_free(data_ecm_m16);
            tabular_fm_free(data_ecm_m12);
            tabular_fm_free(data_ecm_b12);
            tabular_fm_free(data_ecm);
            tabular_fm_free(c);
        }
    }
    return EXIT_SUCCESS;
}
