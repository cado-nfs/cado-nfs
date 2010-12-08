#ifndef ABASE_U64N_H_
#define ABASE_U64N_H_

/* This provides a version of abase_u64, with a runtime-determined repeat
 * value. The runtime value is the only thing to be seen in the abobj
 * data. Note that this structure is bound to be pretty slow, as the
 * runtime-specified repeat count is checked in all tight loops !!
 *
 * For another version which does the same for a constant repeat count,
 * hence not suffering from this speed penalty, see abase-u64k.h
 */

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

typedef unsigned int abase_u64n_obj_t[1];
typedef unsigned int * abase_u64n_obj_ptr;
typedef const unsigned int * abase_u64n_obj_srcptr;

#define abase_u64n_obj_init(x) (x[0] = 2)
#define abase_u64n_obj_clear(x)
#define abase_u64n_obj_init_set(y,x) (y[0]=x[0])

typedef uint64_t abase_u64n_base_type;

/* Our `best-friend' variable-width type is u64n. Funny, that's
 * ourselves !! */

#define abase_u64n_repeat(x)  (x[0])

// nbits takes into account the repeat count. So it's repeat(x) times the
// number of bits in abase_u64n_base_type.
#define abase_u64n_nbits(x)  (64*x[0])

#define abase_u64n_max_accumulate(x) UINT_MAX
#define abase_u64n_max_accumulate_wide(x) UINT_MAX

#define P(X)    PAD(abase_u64n,X)
#define PV(X)    PAD(abase_u64n,X)
#define ABASE_F(t,n,a) static inline t P(n) a

#define AVOID_GENERIC_obj_set_nbys
/* This generic version assumes compile-time-constant stride */
ABASE_F(void,obj_set_nbys,(P(obj_ptr) x MAYBE_UNUSED, unsigned int nbys))
{
    unsigned int q = nbys / 64;
    unsigned int r = nbys % 64;
    ASSERT_ALWAYS(r == 0);
    x[0] = q;
}
#include "abase-binary-generic.h"
#undef AVOID_GENERIC_obj_set_nbys


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


#ifndef ABASE_DONTBIND_u64n
/* Bind our interface as the default one */
#define ABASE_BIND(X)    PAD(abase_u64n,X)
#define ABASE_VBIND(X)    PAD(abase_u64n,X)
#include "abase-api.h"
#endif

#undef  P
#undef  ABASE_F

#if !(defined(ABASE_UNIVERSAL_READAHEAD_ITEMS) && ABASE_UNIVERSAL_READAHEAD_ITEMS % 2 == 0)
#error "The universal readahead count must be a multiple of 2"
#endif

/* proxy for the generic functions */
#define abase_u64n_init(x,n) abase_generic_init(abase_u64n_bytes(x,1),n)
#define abase_u64n_clear(x,p,n) abase_generic_clear(abase_u64n_bytes(x,1),p,n)
#define abase_u64n_initf(x,n) abase_generic_initf(abase_u64n_bytes(x,1),n)
#define abase_u64n_clearf(x,p,n) abase_generic_clearf(abase_u64n_bytes(x,1),p,n)
#define abase_u64n_zero(x,p,n) abase_generic_zero(abase_u64n_bytes(x,1),p,n)
#define abase_u64n_copy(x,q,p,n) abase_generic_copy(abase_u64n_bytes(x,1),q,p,n)
#define abase_u64n_write(x,p,n) abase_generic_write(abase_u64n_bytes(x,1),p,n)
#define abase_u64n_read(x,p,n) abase_generic_read(abase_u64n_bytes(x,1),p,n)
#define abase_u64n_random(x,p,n) abase_generic_random(abase_u64n_bytes(x,1),p,n)

#endif	/* ABASE_U64N_H_ */
