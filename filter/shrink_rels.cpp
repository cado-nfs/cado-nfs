#include "cado.h" // IWYU pragma: keep
#include <sstream>
#include <istream>
#include <ostream>
#include <iostream>
#include <memory>
#include <type_traits>
#include "fmt/core.h"
#include "fmt/format.h"
#include "macros.h"
#include "misc.h"
#include "params.h"
#include "relation.hpp"
#include "relation-tools.h"
#include "indexed_relation.hpp"
#include "gzip.h"

/*
 * The goal of this binary is to read relations that come out of
 * dup1/dup2, and shrink them in a way similar to what fake_rels does
 * with the "shrink" option.
 *
 * Two distinct shrink parameters are provided:
 *      - an integer sigma that denotes the shrink factor with the same meaning
 *        as for fake_rels: it divides column indices by sigma.
 *      - optionally, a second parameter is passed (a floating point
 *        between 0 and 1) that gives the fraction of the input rows that
 *        are kept. If unspecified, this parameter defaults to 1/sigma.
 * Note that if dup1/dup2 split the input into several slices (say 2),
 * the caller has the following options:
 *      - apply a shrink factor sigma to columns, and leave the second
 *      parameter unspecified (i.e. keep one row every sigma), and do
 *      that to the two input files.
 *      - assuming sigma>1, read only the first file, apply a shrink
 *      factor sigma to columns, and 2/sigma to rows.
 * Assuming the split that is done by dup1 introduces no statistical
 * bias, the second option is faster since it reads only part of the
 * input.
 *
 *
 * Note that relations that are read are post-dup2, and hence these are
 * renumbered relations!
 *
 * A random seed can be provided (this has an impact on which rows
 * are kept).
 */


struct shrink_action {
    double row_fraction = 0;
    double shrink_factor = 0;
    index_t shrink_threshold = 0;
    int dl = 0;
    gmp_randstate_t rstate;

    shrink_action() {
        gmp_randinit_default(rstate);
    }
    shrink_action(shrink_action const &) = delete;
    shrink_action& operator=(shrink_action const &) = delete;
    ~shrink_action() {
        gmp_randclear(rstate);
    }

    void process(std::ostream & os, std::istream& is) 
    {
        for(std::string line ; std::getline(is, line) ; ) {
            if (line[0] == '#')
                continue;
            indexed_relation rel;
            if (!(std::istringstream(line) >> rel))
                throw std::runtime_error(fmt::format(FMT_STRING("parse error while reading {}"), line));

            double rnd = double(u64_random(rstate)) / double(UINT64_MAX);
            if (rnd >= row_fraction)
                continue;

            rel.shrink(shrink_factor, shrink_threshold);
            rel.sort();
            rel.compress(dl);

            os << rel << std::endl;
        }
    }
};


static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "out", "output file (defaults to stdout)");
    param_list_decl_usage(pl, "in", "input file (defaults to stdin)");
    param_list_decl_usage(pl, "shrink-factor", "divide all column indices by n");
    param_list_decl_usage(pl, "shrink-threshold", "colums below threshold are not shrunk");
    param_list_decl_usage(pl, "row-fraction", "ratio of rows to keep");
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "dl", "DL mode (do not reduce valuations mod 2)");
}

// coverity[root_function]
int
main (int argc, char *argv[])
{
    char * argv0 = argv[0];
    cxx_param_list pl;

    declare_usage(pl);

    argv++, argc--;

    shrink_action A;

    param_list_configure_switch(pl, "-dl", &A.dl);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

        /* Could also be a file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++,argc--;
            continue;
        }

        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit (EXIT_FAILURE);
    }
    // param_list_print_command_line(stdout, pl);
    //

    param_list_parse_double(pl, "shrink-factor", &A.shrink_factor);
    if (A.shrink_factor < 1) {
        fprintf(stderr, "Error: shrink factor must be an integer >= 1\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    
    {
        unsigned int thresh = 0;
        param_list_parse_uint(pl, "shrink-threshold", &thresh);
        A.shrink_threshold = thresh;
    }
    
    unsigned long seed;

    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(A.rstate, seed);
    }

    if (param_list_parse_double(pl, "row-fraction", &A.row_fraction)) {
        if (A.row_fraction < 0 || A.row_fraction > 1) {
            fprintf(stderr, "Error: row-fraction must be an real number in [0,1]\n");
            param_list_print_usage(pl, argv0, stderr);
            exit(EXIT_FAILURE);
        }
    } else {
        A.row_fraction = 1 / (double) A.shrink_factor;
    }

    std::istream * ptr_in = &std::cin;
    const char * in = param_list_lookup_string(pl, "in");
    std::unique_ptr<std::istream> p_in;
    if (in) {
        p_in = std::unique_ptr<std::istream>(new ifstream_maybe_compressed(in));
        ptr_in = p_in.get();
    }

    std::ostream * ptr_out = &std::cout;
    const char * out = param_list_lookup_string(pl, "out");
    std::unique_ptr<std::ostream> p_out;
    if (out) {
        p_out = std::unique_ptr<std::ostream>(new ofstream_maybe_compressed(out));
        ptr_out = p_out.get();
    }

    A.process(*ptr_out, *ptr_in);

    return 0;
}
