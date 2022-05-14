#ifndef ROPT_SUBLATTICE_CRT_H_
#define ROPT_SUBLATTICE_CRT_H_

#include "ropt_single_sublattice_priority_queue.h"
#include "ropt_sublattice_priority_queue.h"
#include "ropt_str.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned int ropt_sublattice_combine_all_crt(unsigned int nprimes, unsigned int * primes, single_sublattice_priority_queue * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue);

#ifdef __cplusplus
}
#endif

#endif	/* ROPT_SUBLATTICE_CRT_H_ */
