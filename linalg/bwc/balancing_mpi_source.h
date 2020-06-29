#ifndef BALANCING_MPI_SOURCE_H_
#define BALANCING_MPI_SOURCE_H_

#include <stddef.h>                 // for size_t
#include <stdint.h>                 // for uint32_t
#include "balancing_data_source.h"
#include "select_mpi.h"

struct mpi_source_s {
    data_source b;
    int tag;
    int peer;
    int nparallel;

    int over;
    size_t size;
    size_t tailsize;
    uint32_t * buf;
    MPI_Comm comm;
};
typedef struct mpi_source_s mpi_source[1];
typedef struct mpi_source_s *mpi_source_ptr;

#ifdef __cplusplus
extern "C" {
#endif

extern size_t mpi_source_get(mpi_source_ptr s, uint32_t ** p, size_t avail);
extern data_source_ptr mpi_source_alloc(MPI_Comm comm, int peer, size_t queue_size);
extern void mpi_source_free(data_source_ptr q);

#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_MPI_SOURCE_H_ */
