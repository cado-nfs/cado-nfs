#ifndef ABASE_U128_H_
#define ABASE_U128_H_

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
 * Typically, we have plain old datatypes here, e.g. uint128_t's
 *
 * The rather contrived interface must be regarded while having in mind
 * how all this could carry over to data types stored within mpns, even
 * possibly when the type width is going to be determined at runtime.
 */

typedef void * abase_u128_obj_t;
typedef void * abase_u128_obj_ptr;
typedef const void * abase_u128_obj_srcptr;

#define abase_u128_obj_init(x) (x=NULL)
#define abase_u128_obj_clear(x)
#define abase_u128_obj_init_set(y,x) (y=x)

#include <emmintrin.h>
typedef __v2di abase_u128_base_type;

/* Our `best-friend' variable-width type is u64n. We ask not to bind it,
 * but our call to abase-api.h will do vbinding though */
#define ABASE_DONTBIND_u64n
#include "abase-u64n.h"

#define abase_u128_nbits(x)  128

// repeating is just the fact of putting several base type elements next
// to each other.
#define abase_u128_repeat(x)  1

#define abase_u128_max_accumulate(x) UINT_MAX
#define abase_u128_max_accumulate_wide(x) UINT_MAX

#define P(X)    PAD(abase_u128,X)
#define PV(X)    PAD(abase_u64n,X)
#define ABASE_F(t,n,a) static inline t P(n) a

#define AVOID_GENERIC_is_zero
#include "abase-binary-generic.h"
#undef  AVOID_GENERIC_is_zero

ABASE_F(int, is_zero, (P(obj_srcptr) x MAYBE_UNUSED,
            const P(base_type) * p,
            unsigned int n MAYBE_UNUSED))
{
    for(unsigned int i = 0 ; i < P(bytes)(x,n) ; i++) {
        if (((unsigned char *)p)[i])
            return 0;
    }
    return 1;
}

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


#ifndef ABASE_DONTBIND_u128
/* Bind our interface as the default one */
#define ABASE_BIND(X)    PAD(abase_u128,X)
#define ABASE_VBIND(X)    PAD(abase_u64n,X)
#include "abase-api.h"
#endif

#undef  P
#undef  ABASE_F

/* Contrary to most u64-based abases, this one does not require readahead
 * padding */

/* proxy for the generic functions */
#define abase_u128_init(x,n) abase_generic_init(abase_u128_bytes(x,1),n)
#define abase_u128_clear(x,p,n) abase_generic_clear(abase_u128_bytes(x,1),p,n)
#define abase_u128_initf(x,n) abase_generic_initf(abase_u128_bytes(x,1),n)
#define abase_u128_clearf(x,p,n) abase_generic_clearf(abase_u128_bytes(x,1),p,n)
#define abase_u128_zero(x,p,n) abase_generic_zero(abase_u128_bytes(x,1),p,n)
#define abase_u128_copy(x,q,p,n) abase_generic_copy(abase_u128_bytes(x,1),q,p,n)
#define abase_u128_write(x,p,n) abase_generic_write(abase_u128_bytes(x,1),p,n)
#define abase_u128_read(x,p,n) abase_generic_read(abase_u128_bytes(x,1),p,n)
#define abase_u128_random(x,p,n) abase_generic_random(abase_u128_bytes(x,1),p,n)

#endif	/* ABASE_U128_H_ */
