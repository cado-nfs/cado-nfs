#ifndef CADO_SELECT_MPI_H
#define CADO_SELECT_MPI_H

// scan-headers: no prototypes
// scan-headers: stop here

#include "cado_mpi_config.h"    // IWYU pragma: keep

#ifdef WITH_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>        // IWYU pragma: export
#else
#include "fakempi.h"    // IWYU pragma: export
#endif

#endif	/* CADO_SELECT_MPI_H */
