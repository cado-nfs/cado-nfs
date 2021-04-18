#ifndef GMP_HACKS_H_
#define GMP_HACKS_H_

#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"

/* TODO: remove this. It's only here by lack of something better for the
 * needs of crtalgsqrt. */

/* GMP field access macros from gmp-impl.h */
/* Starting with gmp-6.0.0 there is a better interface that we can use.
 * Actually we _have_ to use it from gmp-6.2.0 onwards.
 *
 * At some point we'll get rid of these ugly macros for good (i.e. when
 * we safely assume that we have gmp-6+)
 */
#ifndef PTR
#define PTR(x) ((x)->_mp_d)
#endif
#ifndef SIZ
#define SIZ(x) ((x)->_mp_size)
#endif
#ifndef ALLOC
#define ALLOC(x) ((x)->_mp_alloc)
#endif

#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {                                                                  \
    while ((NLIMBS) > 0)                                                \
      {                                                                 \
        if ((DST)[(NLIMBS) - 1] != 0)                                   \
          break;                                                        \
        (NLIMBS)--;                                                     \
      }                                                                 \
  } while (0)
#endif

#ifndef	MPN_COPY
#define	MPN_COPY(dst, src, n)	memcpy((dst), (src), (n) * sizeof(mp_limb_t))
#endif

#ifndef	MPN_ZERO
#define	MPN_ZERO(dst, n)	memset((dst), 0, (n) * sizeof(mp_limb_t))
#endif

/* Useful for the lazy boyz */

#ifdef __cplusplus
extern "C" {
#endif

static inline void MPZ_INIT_SET_MPN(mpz_ptr DST, const mp_limb_t * SRC, size_t NLIMBS)
{
#if GMP_VERSION_ATLEAST(6, 0, 0)
    mpz_t foo;
    mpz_init_set(DST, mpz_roinit_n(foo, SRC, NLIMBS));
#else
    ALLOC(DST) = (NLIMBS);
    SIZ(DST) = (NLIMBS);
    PTR(DST) = (mp_limb_t*) malloc((NLIMBS) * sizeof(mp_limb_t));
    memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));
    MPN_NORMALIZE(PTR(DST),SIZ(DST));
#endif
}

static inline void MPZ_SET_MPN(mpz_ptr DST, const mp_limb_t * SRC, size_t NLIMBS)
{
#if GMP_VERSION_ATLEAST(6, 0, 0)
    mpn_copyi(mpz_limbs_write(DST, NLIMBS),(SRC),(NLIMBS));
    mpz_limbs_finish(DST, NLIMBS);
#else
    /* MPZ_GROW_ALLOC(DST, NLIMBS); */
    {
        if (ALLOC(DST) < (int) (NLIMBS)) {
            ALLOC(DST) = (NLIMBS);
            PTR(DST)=(mp_limb_t *) realloc(PTR(DST),
                    (NLIMBS) * sizeof(mp_limb_t));
        }
    }
    SIZ(DST) = (NLIMBS);
    memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));
    MPN_NORMALIZE(PTR(DST),SIZ(DST));
#endif
}

static inline void MPN_SET_MPZ(mp_limb_t * DST, size_t NLIMBS, mpz_srcptr SRC)
{
#if GMP_VERSION_ATLEAST(6, 0, 0)
    mp_size_t r = MIN((size_t) mpz_size(SRC), NLIMBS);
    mpn_copyi(DST, mpz_limbs_read(SRC), r);
    mpn_zero(DST+mpz_size(SRC), NLIMBS-r);
#else
    mp_size_t r = MIN((size_t) ABS(SIZ(SRC)), NLIMBS);
    memcpy((DST),PTR(SRC),r * sizeof(mp_limb_t));
    memset((DST)+ABS(SIZ(SRC)),0,((NLIMBS)-r) * sizeof(mp_limb_t));
#endif
}

#ifdef __cplusplus
}
#endif


#endif /* GMP_HACKS_H_ */	
