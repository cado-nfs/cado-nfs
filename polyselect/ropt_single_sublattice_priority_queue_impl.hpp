#ifndef CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_IMPL_HPP
#define CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_IMPL_HPP

#include "ropt_single_sublattice_priority_queue.h"
#include "min_max_heap.hpp"

/**
 * Priority queue for single_sublattices over a product of p^e.
 */

struct single_sublattice_priority_queue_impl : min_max_heap<single_sublattice_info> {
    size_t max_count;
    static single_sublattice_priority_queue_impl * cast(single_sublattice_priority_queue_ptr s) { return (single_sublattice_priority_queue_impl *) s->impl; }
    static single_sublattice_priority_queue_impl const * cast(single_sublattice_priority_queue_srcptr s) { return (single_sublattice_priority_queue_impl const *) s->impl; }
};

#endif	/* CADO_ROPT_SINGLE_SUBLATTICE_PRIORITY_QUEUE_IMPL_HPP */
