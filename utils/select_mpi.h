#ifndef SELECT_MPI_H_
#define SELECT_MPI_H_

// scan-headers: no prototypes
// scan-headers: stop here

#include "cado_mpi_config.h"
#include "macros.h"

#ifdef WITH_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#else
#include "fakempi.h"
#endif

#endif	/* SELECT_MPI_H_ */
