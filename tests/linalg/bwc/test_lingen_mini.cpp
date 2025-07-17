#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "fmt/format.h"

#include "bw-common.h"                    // for bw, bw_common_clear, bw_com...
#include "misc.h"
#include "lingen_checkpoints.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_matpoly_select.hpp"
#include "lingen_io_matpoly.hpp"
#include "lingen_tuning.hpp"
#include "params.h"
#include "lingen.hpp"

static int verbose;

static void declare_usage(cxx_param_list & pl)
{
    param_list_usage_header(pl,
            "This program takes an E checkpoint as computed by lingen,"
            " and computes the pi matrix that can be obtained from it."
            " This matrix should obey the following constraints."
            " We assume that E is defined mod X^d and that the delta values"
            " are delta_0 to delta_{n-1}. We must have E*pi = O(X^d),"
            " the columns of ((X^{delta_i})_i) * pi must be bounded"
            " by the (new) delta_i's, and the total increase from the delta_i"
            " to the new delta_i's is at most m*d");
    param_list_decl_usage(pl, "v", "be verbose");
}

template<bool is_binary>
static void do_one_lingen(std::string const & filename, cxx_param_list &pl)
{
    typename lingen_checkpoint<is_binary>::header_info h;
    if (!(std::ifstream(filename) >> h))
        throw std::runtime_error(
                fmt::format("bad aux file for {}", filename));

    bmstatus<is_binary> bm(h.m, h.n, h.p);
    bm.set_t0(h.t0);
    bm.delta = h.delta;
    bm.lucky = h.lucky;
    bm.done = h.done;
    // bm.depth = h.level;

    /* because we're reading E directly, the length that we must pass to
     * the tuner is simply the number of coefficients. And it order to
     * make sure that the configuration that is picked is exactly as we
     * expect, we must pass a depth of zero, irrespective of the depth
     * at which this E file was created originally.
     */
    bm.depth = 0;
    bm.hints = lingen_tuning(bm.d, h.ncoeffs, bm.com[0], pl);

    lingen_checkpoint<is_binary> cp(bm, h.t0, h.t1, 0, filename);

    matpoly<is_binary> E(&bm.d.ab, bm.d.m, bm.d.m + bm.d.n, (int) h.ncoeffs);
    E.set_size(h.ncoeffs);

    // lingen checkpoints are definitely a mess, because they include
    // snapshots of the whole tree.
    /*
    size_t Xsize;
    cp.load_aux_file(Xsize);
    E.set_size(Xsize);
    bool ok = cp.load_data_file(E);
    ASSERT_ALWAYS(ok);
    */

    std::unique_ptr<FILE, delete_FILE> data(fopen(cp.datafile.c_str(), "rb"));
    if (!data) {
        fmt::print(stderr, "Warning: cannot open {}\n", cp.datafile);
        return;
    }
    int rc = lingen_io_matpoly<is_binary>::read(&bm.d.ab, data.get(), E, 0, E.get_size(), 0, 0);
    if (rc != (int) E.get_size()) {
        throw std::runtime_error(fmt::format("{}: short read", filename));
        return;
    }

    E.clear_high_word();
    // do one debug print
    // matpoly_write(&bm.d.ab, std::cout, E, 0, E.get_size(), 1, 0);

    if (!lingen_checkpoint<is_binary>::default_directory.empty()) {
        lingen_checkpoint<is_binary>::threshold = 0;
        save_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E, h.t0, h.t1);
    }

    matpoly pi = bw_lingen_single(bm, E);
    // do another debug print
    std::cout << "\n\n\n";
    // matpoly_write(&bm.d.ab, std::cout, pi, 0, pi.get_size(), 1, 0);
    std::ofstream fpi(filename + ".pi");
    // matpoly_write(&bm.d.ab, std::cout, pi, 0, pi.get_size(), 1, 0);
    lingen_io_matpoly<is_binary>::write(&bm.d.ab, fpi, pi, 0, pi.get_size(), 0, 0);

    if (!lingen_checkpoint<is_binary>::default_directory.empty()) {
        lingen_checkpoint<is_binary>::threshold = 0;
        save_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, h.t0, h.t1);
    }
}


int main(int argc, char const * argv[])
{
#ifdef LINGEN_BINARY
    /* well, we don't compile lingen_mini in the binary case, it seems */
    constexpr bool is_binary = true;
#else
    constexpr bool is_binary = false;
#endif

    bw_common_init(bw, &argc, &argv);

    cxx_param_list pl;
    cxx_param_list pl2;

    char const * argv0 = argv[0];

    declare_usage(pl);
    lingen_tuning_decl_usage(pl);
    lingen_checkpoint<is_binary>::decl_usage(pl);

    param_list_configure_switch(pl, "-v", &verbose);

    std::vector<std::string> E_files;

    for (argc--, argv++; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }

        if (has_suffix(argv[0], ".aux")) {
            E_files.emplace_back(argv[0]);
            argc--, argv++;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    lingen_tuning_lookup_parameters(pl);
    lingen_checkpoint<is_binary>::interpret_parameters(pl);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    for(auto const & E : E_files) {
        do_one_lingen<is_binary>(E, pl);
    }

    bw_common_clear(bw);
}
