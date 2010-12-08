#ifndef ABASE_U64_H_
#define ABASE_U64_H_

/* This is used in order to enable GF(2)-only code */
#define ABASE_BINARY

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* I know it's ugly */
#include "bwc_config.h"
#include "macros.h"
#include "random_generation.h"
#include "pad.h"
// #include "electric_alloc.h"

#include "abase-generic.h"

/* provides an interface for the arithmetic base of computations.
 * Typically, we have plain old datatypes here, e.g. uint64_t's
 *
 * The rather contrived interface must be regarded while having in mind
 * how all this could carry over to data types stored within mpns, even
 * possibly when the type width is going to be determined at runtime.
 */

typedef void * abase_u64_obj_t;
typedef void * abase_u64_obj_ptr;
typedef const void * abase_u64_obj_srcptr;

#define abase_u64_obj_init(x) (x=NULL)
#define abase_u64_obj_clear(x)
#define abase_u64_obj_init_set(y,x) (y=x)

typedef uint64_t abase_u64_base_type;

/* Our `best-friend' variable-width type is u64n. We ask not to bind it,
 * but our call to abase-api.h will do vbinding though */
#define ABASE_DONTBIND_u64n
#include "abase-u64n.h"

// repeating is just the fact of putting several base type elements next
// to each other. in the case of this header, it's an easy constant. But
// the api has provision for it because in hard cases it could be handy
// (albeit never needed in critical paths -- only handy as a convenient
// way of doing bookkeeping tasks in e.g. prep).
#define abase_u64_repeat(x)  1

// nbits takes into account the repeat count. So it's repeat(x) times the
// number of bits in abase_u64_base_type.
#define abase_u64_nbits(x)  64

#define abase_u64_max_accumulate(x) UINT_MAX
#define abase_u64_max_accumulate_wide(x) UINT_MAX

#define P(X)    PAD(abase_u64,X)
#define PV(X)    PAD(abase_u64n,X)
#define ABASE_F(t,n,a) static inline t P(n) a
#include "abase-binary-generic.h"

/* This one is in the C file. */
extern void
P(dotprod)(P(obj_srcptr) x,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n);
extern void
P(vdotprod)(P(obj_srcptr) x,
        PV(obj_srcptr) y,
        P(base_type) * w,
        const PV(base_type) * u,
        const P(base_type) * v,
        unsigned int n);
extern void
P(vaddmul_tiny)(P(obj_srcptr) x,
        PV(obj_srcptr) y,
        PV(base_type) * w,
        const P(base_type) * u,
        const PV(base_type) * v,
        unsigned int n);
extern void
P(vtranspose)(PV(obj_srcptr) y,
        P(obj_srcptr) x,
        PV(base_type) * w,
        const P(base_type) * u);

#ifndef ABASE_DONTBIND_u64
/* Bind our interface as the default one */
#define ABASE_BIND(X)    PAD(abase_u64,X)
#define ABASE_VBIND(X)    PAD(abase_u64n,X)
#include "abase-api.h"
#endif

#undef  P
#undef  ABASE_F

#if !(defined(ABASE_UNIVERSAL_READAHEAD_ITEMS) && ABASE_UNIVERSAL_READAHEAD_ITEMS % 2 == 0)
#error "The universal readahead count must be a multiple of 2"
#endif

/* proxy for the generic functions */
#define abase_u64_init(x,n) abase_generic_init(abase_u64_bytes(x,1),n)
#define abase_u64_clear(x,p,n) abase_generic_clear(abase_u64_bytes(x,1),p,n)
#define abase_u64_initf(x,n) abase_generic_initf(abase_u64_bytes(x,1),n)
#define abase_u64_clearf(x,p,n) abase_generic_clearf(abase_u64_bytes(x,1),p,n)
#define abase_u64_zero(x,p,n) abase_generic_zero(abase_u64_bytes(x,1),p,n)
#define abase_u64_copy(x,q,p,n) abase_generic_copy(abase_u64_bytes(x,1),q,p,n)
#define abase_u64_write(x,p,n) abase_generic_write(abase_u64_bytes(x,1),p,n)
#define abase_u64_read(x,p,n) abase_generic_read(abase_u64_bytes(x,1),p,n)
#define abase_u64_random(x,p,n) abase_generic_random(abase_u64_bytes(x,1),p,n)

#endif	/* ABASE_U64_H_ */
