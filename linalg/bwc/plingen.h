#ifndef PLINGEN_H_
#define PLINGEN_H_

/* This contains some definitions for lingen mod p.
 *
 * This code is to be understood as scotch tape around an old and quirky
 * implementation.
 */

#include <stdlib.h>

#include "mpfq_layer.h"
#include "select_mpi.h"

/* TODO: Rename ! */
typedef struct {
    unsigned int m,n,nrhs;
    abfield ab;
} dims;

struct bmstatus_s {
    dims d[1];
    unsigned int t;
    int * lucky;

    double t_basecase;
    double t_mp;
    double t_mul;
    double t_cp_io;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm com[3]; /* [0]: MPI_COMM_WORLD, reordered.
                        [1]: row-wise
                        [2]: column-wise */
};
typedef struct bmstatus_s bmstatus[1];
typedef struct bmstatus_s *bmstatus_ptr;

extern void bmstatus_init(bmstatus_ptr bm, unsigned int m, unsigned int n);
extern void bmstatus_clear(bmstatus_ptr bm);


#endif	/* PLINGEN_H_ */
