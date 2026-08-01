#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "parallelizing_info.hpp"
#include "select_mpi.h"
#include "params.hpp"
#include "macros.h"

static int verbose=0;

static void * program(parallelizing_info & pi, cxx_param_list & pl MAYBE_UNUSED, void * arg MAYBE_UNUSED)
{
    if (verbose) {
        pi.m.log_init();
        pi.wr[0].log_init();
        pi.wr[1].log_init();
    }

    // it is here as a cheap sanity check.
    pi.hello();

    if (verbose) {
        pi.m.log_op("serialize");
        pi.m.serialize(__FILE__, __LINE__);

        /* note that in order to do pi.wr[0].serialize(__FILE__, __LINE__), we need to make
         * sure that only one thread in the intersecting communicator
         * executes.
         */
        if (pi.wr[1].trank == 0) {
            pi.wr[0].log_op("serialize(2nd)");
            pi.wr[0].serialize(__FILE__, __LINE__);
        }
        pi.wr[1].serialize_threads(__FILE__, __LINE__);

        if (pi.wr[0].trank == 0) {
            pi.wr[1].log_op("serialize(3rd)");
            pi.wr[1].serialize(__FILE__, __LINE__);
        }
        pi.wr[0].serialize_threads(__FILE__, __LINE__);

        pi.log_print_all();

        pi.m.log_clear();
        pi.wr[0].log_clear();
        pi.wr[1].log_clear();
    }

    return nullptr;
}

int main(int argc, char const * argv[])
{
    int rank;
    int size;
    cxx_param_list pl;

    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    MPI_Init(&argc, (char ***) &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    parallelizing_info::init_attribute_things();

    parallelizing_info::declare_usage(pl);
    pl.declare_usage("v", "turn on some demo logging");

    pl.configure_switch("v");

    pl.process_command_line(argc, argv, false);

    pl.parse("v", verbose);
    parallelizing_info::lookup_parameters(pl);

    if (verbose)
        pl.display_debug(stderr);

    if (pl.warn_unused())
        pl.fail("Unused parameters are given");

    pi_go(program, pl, nullptr);

    parallelizing_info::clear_attribute_things();

    MPI_Finalize();

    return 0;
}

