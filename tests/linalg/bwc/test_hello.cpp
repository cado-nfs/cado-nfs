#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "parallelizing_info.hpp"
#include "select_mpi.h"
#include "params.h"
#include "macros.h"

int verbose=0;

void * program(parallelizing_info_ptr pi, param_list pl MAYBE_UNUSED, void * arg MAYBE_UNUSED)
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

    return NULL;
}

int main(int argc, char * argv[])
{
    int rank;
    int size;
    param_list pl;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    param_list_init (pl);
    parallelizing_info_init();

    parallelizing_info_decl_usage(pl);
    param_list_decl_usage(pl, "v", "turn on some demo logging");

    const char * programname = argv[0];

    argv++, argc--;
    param_list_configure_switch(pl, "v", &verbose);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, programname, stderr);
        exit(EXIT_FAILURE);
    }

    parallelizing_info_lookup_parameters(pl);

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, programname, stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(program, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);

    MPI_Finalize();

    return 0;
}

