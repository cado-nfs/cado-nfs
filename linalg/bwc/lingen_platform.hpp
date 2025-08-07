#ifndef CADO_LINGEN_PLATFORM_HPP
#define CADO_LINGEN_PLATFORM_HPP

#include <cstddef>
#include "select_mpi.h"
struct cxx_param_list;

struct lingen_platform {
    /* input characteristics -- the ones we have to live with */

    /* We give timings for a run on r*r nodes, with T threads per node */
    /* Note that all of this can also be auto-detected, a priori */
    MPI_Comm comm;
    unsigned int r;
    unsigned int T;     /* **PHYSICAL** cores, or we say rubbish  */
    int openmp_threads;

    size_t available_ram;       // zero means unlimited.

    /* Assume we output something like one gigabyte per second. This is
     * rather conservative for HPC networks */
    /*
     * Note that is the reason for one of the inaccuracies in the timing
     * estimations for lingen.
     */
    double mpi_xput = 1e9;

    static void lookup_parameters(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);
    lingen_platform(MPI_Comm comm, cxx_param_list & pl);

    lingen_platform single() const {
        lingen_platform res = *this;
        res.r = 1;
        return res;
    }
};


#endif	/* LINGEN_PLATFORM_HPP_ */
