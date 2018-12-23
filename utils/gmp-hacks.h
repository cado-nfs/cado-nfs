#ifndef GMP_HACKS_H_
#define GMP_HACKS_H_

#include <string.h>
#include <gmp.h>
#include "macros.h"

/* TODO: remove this. It's only here by lack of something better for the
 * needs of crtalgsqrt. */

/* GMP field access macros from gmp-impl.h */
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

static inline void MPZ_INIT_SET_MPN(mpz_ptr DST, const mp_limb_t * SRC, size_t NLIMBS)
{
    ALLOC(DST) = (NLIMBS);
    SIZ(DST) = (NLIMBS);
    PTR(DST) = (mp_limb_t*) malloc((NLIMBS) * sizeof(mp_limb_t));
    memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));
    MPN_NORMALIZE(PTR(DST),SIZ(DST));
}

static inline void MPZ_SET_MPN(mpz_ptr DST, const mp_limb_t * SRC, size_t NLIMBS)
{
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
}

static inline void MPN_SET_MPZ(mp_limb_t * DST, size_t NLIMBS, mpz_srcptr SRC)
{
    mp_size_t r = MIN((size_t) ABS(SIZ(SRC)), NLIMBS);
    memcpy((DST),PTR(SRC),r * sizeof(mp_limb_t));
    memset((DST)+ABS(SIZ(SRC)),0,((NLIMBS)-r) * sizeof(mp_limb_t));
}

#endif /* GMP_HACKS_H_ */	
