#include <stdlib.h>
#ifndef CHEATING_VEC_INIT_H_
#define CHEATING_VEC_INIT_H_

#include "macros.h"
#include "memory.h"
#include "mpfq/mpfq_vbase.h"


/* If you change stuff here, make sure that you also change the identical
 * definitions in
 tests/linalg/bwc/bwc-ptrace.sh
 tests/linalg/bwc/convert_magma.pl
 */
#define FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES      64
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
#ifdef __cplusplus
extern "C" {
#endif

static inline void cheating_vec_init(mpfq_vbase_ptr A, void ** p, size_t nitems)
{
    *p = malloc_aligned(A->vec_elt_stride(A, nitems), FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES);
}

static inline void cheating_vec_clear(mpfq_vbase_ptr A MAYBE_UNUSED, void ** p, size_t nitems MAYBE_UNUSED)
{
    free_aligned(*p);
    *p=NULL;
}


#ifdef __cplusplus
}
#endif

#endif	/* CHEATING_VEC_INIT_H_ */
