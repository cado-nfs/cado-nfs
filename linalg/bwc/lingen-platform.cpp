#include "cado.h"
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif
#include "lingen-platform.hpp"

void lingen_platform::lookup_parameters(cxx_param_list & pl) {
    param_list_lookup_string(pl, "max_ram");
    param_list_lookup_string(pl, "tuning_thr");
    param_list_lookup_string(pl, "tuning_mpi");
    param_list_lookup_string(pl, "mpi");
    param_list_lookup_string(pl, "thr");
}

void lingen_platform::declare_usage(cxx_param_list & pl) {
    /* TODO: this shall supersede mpi= and thr= that are currently
     * parsed from within plingen.cpp */
    param_list_decl_usage(pl, "max_ram",
            "Maximum local memory to be used for transforms and matrices, in GB");
    param_list_decl_usage(pl, "tuning_thr",
            "Number of threads to be used for tuning only (if different from real thr=)");
    param_list_decl_usage(pl, "tuning_mpi",
            "Number of mpi jobs to be used for tuning only (if different from real mpi=)");
}

lingen_platform::lingen_platform(MPI_Comm comm, cxx_param_list & pl) : comm(comm) {

    int mpi[2] = { 1, 1 };
    int thr[2] = { 1, 1 };

    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "tuning_mpi", mpi);

    param_list_parse_intxint(pl, "thr", thr);
    param_list_parse_intxint(pl, "tuning_thr", thr);

    T = thr[0] * thr[1];
    r = mpi[0];

    int rank;
    MPI_Comm_rank(comm, &rank);

    if (mpi[0] != mpi[1]) {
        if (!rank)
            fprintf(stderr, "The current plingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
                    mpi[0], mpi[1]);
        abort();
    }

    double dtmp = 1;
    param_list_parse_double(pl, "max_ram", &dtmp);
    available_ram = dtmp * (1 << 30);

    openmp_threads = 1;
#ifdef HAVE_OPENMP
    openmp_threads = omp_get_max_threads();
#endif
}
