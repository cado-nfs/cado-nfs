#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "parallelizing_info.hpp"
#include "select_mpi.h"
#include "params.h"
#include "macros.h"

static int verbose=0;

static void * program(parallelizing_info_ptr pi, cxx_param_list & pl MAYBE_UNUSED, void * arg MAYBE_UNUSED)
{
    if (verbose) {
        pi_log_init(pi->m);
        pi_log_init(pi->wr[0]);
        pi_log_init(pi->wr[1]);
    }

    // it is here as a cheap sanity check.
    pi_hello(pi);

    if (verbose) {
        pi_log_op(pi->m, "serialize");
        serialize(pi->m);

        /* note that in order to do serialize(pi->wr[0]), we need to make
         * sure that only one thread in the intersecting communicator
         * executes.
         */
        if (pi->wr[1]->trank == 0) {
            pi_log_op(pi->wr[0], "serialize(2nd)");
            serialize(pi->wr[0]);
        }
        serialize_threads(pi->wr[1]);

        if (pi->wr[0]->trank == 0) {
            pi_log_op(pi->wr[1], "serialize(3rd)");
            serialize(pi->wr[1]);
        }
        serialize_threads(pi->wr[0]);

        pi_log_print_all(pi);

        pi_log_clear(pi->m);
        pi_log_clear(pi->wr[0]);
        pi_log_clear(pi->wr[1]);
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


    parallelizing_info_init();

    parallelizing_info_decl_usage(pl);
    param_list_decl_usage(pl, "v", "turn on some demo logging");

    param_list_configure_switch(pl, "v", &verbose);

    param_list_process_command_line(pl, &argc, &argv, false);

    parallelizing_info_lookup_parameters(pl);

    if (verbose)
        param_list_display (pl, stderr);

    if (param_list_warn_unused(pl))
        pl.fail("Unused parameters are given");

    pi_go(program, pl, nullptr);

    parallelizing_info_finish();

    MPI_Finalize();

    return 0;
}

