#include "cado.h" // IWYU pragma: keep

#include "ropt_single_sublattice_priority_queue.h"
#include "ropt_single_sublattice_priority_queue_impl.hpp"

void single_sublattice_priority_queue_init(single_sublattice_priority_queue_ptr q,
                         size_t max_len)
{
    auto * qi = new single_sublattice_priority_queue_impl;
    q->impl = (void*) qi;
    if (max_len == 0)
        max_len = 1;
    qi->max_count = max_len;
    q->impl = (void *) qi;

}

void single_sublattice_priority_queue_push(single_sublattice_priority_queue_ptr q, const struct single_sublattice_info * z)
{
    auto * qi = single_sublattice_priority_queue_impl::cast(q);
    qi->push(*z);
    if (qi->size() > qi->max_count)
        qi->popMin();
}

void single_sublattice_priority_queue_clear (single_sublattice_priority_queue_ptr q)
{
    delete single_sublattice_priority_queue_impl::cast(q);
}

size_t single_sublattice_priority_queue_size (single_sublattice_priority_queue_srcptr q)
{
    return single_sublattice_priority_queue_impl::cast(q)->size();
}

int single_sublattice_priority_queue_empty (single_sublattice_priority_queue_srcptr q)
{
    return single_sublattice_priority_queue_impl::cast(q)->empty();
}

/* This calls the fucntion f on all the queue elements, in order. We
 * could as well pass each element's rank and score, but presently we do
 * not.
 */
void single_sublattice_priority_queue_pop_all(single_sublattice_priority_queue_srcptr q, struct single_sublattice_info * out)
{
    auto const * qi = single_sublattice_priority_queue_impl::cast(q);
    /* copy, and consume the copy */
    single_sublattice_priority_queue_impl cpi = *qi;
    for( ; !cpi.empty() ; )
        *out++ = cpi.popMax();
}

