#ifndef CADO_ROPT_SUBLATTICE_CRT_H
#define CADO_ROPT_SUBLATTICE_CRT_H

#include "ropt_single_sublattice_priority_queue.h"
#include "ropt_sublattice_priority_queue.h"
#include "ropt_str.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned int ropt_sublattice_combine_all_crt(unsigned int nprimes, const unsigned int * primes, single_sublattice_priority_queue const * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_ROPT_SUBLATTICE_CRT_H */
