#ifndef CHEATING_VEC_INIT_HPP_
#define CHEATING_VEC_INIT_HPP_

#include <stdlib.h>
#include "macros.h"
#include "memory.h"
#include "arith-generic.hpp"

#define FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES      32
#define MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES     4

/* The mpfq routines for doing vec_init rely on simple malloc() to do their
 * job.
 *
 * Unfortunately our alignment constraints are stricter.
 *
 * We could consider modifying mpfq to allow an alignment constraint in
 * vec_init, but that would be clutter.
 *
 * Instead, we cheat on vec_init, and leverage the fact that cado-nfs
 * already has all sorts of aligned_malloc's and friends
 */

static inline void cheating_vec_init(arith_generic * A, arith_generic::elt ** p, size_t nitems)
{
    *p = (arith_generic::elt *) malloc_aligned(A->vec_elt_stride(nitems), FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES);
}

static inline void cheating_vec_clear(arith_generic * A MAYBE_UNUSED, arith_generic::elt ** p, size_t nitems MAYBE_UNUSED)
{
    free_aligned(*p);
    *p=NULL;
}

#endif	/* CHEATING_VEC_INIT_HPP_ */
