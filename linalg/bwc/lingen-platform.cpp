#include "cado.h"
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif
#include "lingen-platform.hpp"

void lingen_platform::lookup_parameters(cxx_param_list & pl) {/*{{{*/
    param_list_lookup_string(pl, "max_ram");
    param_list_lookup_string(pl, "tuning_T");
    param_list_lookup_string(pl, "tuning_r");
    param_list_lookup_string(pl, "mpi");
    param_list_lookup_string(pl, "thr");
}/*}}}*/

void lingen_platform::declare_usage(cxx_param_list & pl) {/*{{{*/
    /* TODO: this shall supersede mpi= and thr= that are currently
     * parsed from within plingen.cpp */
    param_list_decl_usage(pl, "max_ram",
            "Maximum local memory to be used for transforms and matrices, in GB");
    param_list_decl_usage(pl, "tuning_T",
            "For --tune only: target number of threads (if for different platform)");
    param_list_decl_usage(pl, "tuning_r",
            "For --tune only: size of the mpi grid (the grid would be r times r)");
}/*}}}*/

lingen_platform::lingen_platform(MPI_Comm comm, cxx_param_list & pl) : comm(comm) {/*{{{*/

    int mpi[2];
    int thr[2];

    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);

    T = thr[0] * thr[1];
    r = mpi[0];

    double dtmp = 1;
    param_list_parse_double(pl, "max_ram", &dtmp);
    available_ram = dtmp * (1 << 30);

    param_list_parse_uint(pl, "tuning_T", &T);
    param_list_parse_uint(pl, "tuning_r", &r);

    openmp_threads = 1;
#ifdef HAVE_OPENMP
    openmp_threads = omp_get_max_threads();
#endif
}/*}}}*/
