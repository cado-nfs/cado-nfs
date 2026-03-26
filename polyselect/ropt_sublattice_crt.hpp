#ifndef CADO_ROPT_SUBLATTICE_CRT_HPP
#define CADO_ROPT_SUBLATTICE_CRT_HPP

#include "ropt_single_sublattice_priority_queue.h"
#include "ropt_sublattice_priority_queue.hpp"
#include "ropt_str.hpp"

unsigned int ropt_sublattice_combine_all_crt(unsigned int nprimes, const unsigned int * primes, single_sublattice_priority_queue const * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue);

#endif	/* CADO_ROPT_SUBLATTICE_CRT_HPP */
