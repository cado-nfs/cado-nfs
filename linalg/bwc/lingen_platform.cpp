#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#ifdef HAVE_OPENMP
#include "omp_proxy.h"
#endif
#include "lingen_platform.hpp"
#include "params.hpp"

/* TODO: fetch data from hwloc */

/* TODO tuning_XXX here are mostly leftovers, they should go. OTOH, this
 * platform class is where we should declare, parse, and host the mpi and
 * thr parameters.
 *
 * As to tuning specifics, the question is if/how we retain the
 * posibility to compute the optimal schedules on a machine given a cache
 * that was computed elsewhere. A priori, this all remains possible, but
 * hasn't been used much. At least we should bear this in mind when using
 * hwloc here, because this makes a difference.
 */
void lingen_platform::lookup_parameters(cxx_param_list & pl) {
    pl.lookup("max_ram");
    pl.lookup("tuning_thr");
    pl.lookup("tuning_mpi");
    pl.lookup("mpi");
    pl.lookup("thr");
}

void lingen_platform::declare_usage(cxx_param_list & pl) {
    /* TODO: this shall supersede mpi= and thr= that are currently
     * parsed from within lingen.cpp */
    pl.declare_usage("max_ram",
            "Maximum local memory to be used for transforms and matrices, in GB");
    pl.declare_usage("tuning_thr",
            "Number of threads to be used for tuning only (if different from real thr=)");
    pl.declare_usage("tuning_mpi",
            "Number of mpi jobs to be used for tuning only (if different from real mpi=)");
}

lingen_platform::lingen_platform(MPI_Comm comm, cxx_param_list & pl) : comm(comm) {

    std::array<int, 2> mpi = { 1, 1 };
    std::array<int, 2> thr = { 1, 1 };

    pl.parse("mpi", mpi, "x");
    pl.parse("tuning_mpi", mpi, "x");

    pl.parse("thr", thr, "x");
    pl.parse("tuning_thr", thr, "x");

    T = thr[0] * thr[1];
    r = mpi[0];

    int rank;
    MPI_Comm_rank(comm, &rank);

    if (mpi[0] != mpi[1]) {
        if (!rank)
            fprintf(stderr, "The current lingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
                    mpi[0], mpi[1]);
        abort();
    }

    /* default for the RAM amount is zero, which means: take everything.
     * (see #30022)
     * note that if hwloc is available, we could conceivably do a quick
     * peek at the real RAM amount.
     */
    double dtmp = 0;
    pl.parse("max_ram", dtmp);
    available_ram = dtmp * (1 << 30);

    openmp_threads = 1;
#ifdef HAVE_OPENMP
    openmp_threads = omp_get_max_threads();
#endif
}
