#ifndef CADO_RELATION_TOOLS_H
#define CADO_RELATION_TOOLS_H

#include <stdint.h>    // for int64_t, uint64_t
#include "gmp_aux.h"
#include "typedefs.h"  // for p_r_values_t

#ifdef __cplusplus
extern "C" {
#endif

p_r_values_t relation_compute_r (int64_t a, uint64_t b, p_r_values_t p);
extern char * u64toa16 (char *p, uint64_t m);
extern char * u64toa10 (char *p, uint64_t m);
extern char * d64toa10 (char *p, int64_t m);
extern char * d64toa16 (char *p, int64_t m);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* Only used in dup2 which is a C++ binary */
inline p_r_values_t
relation_compute_r(mpz_srcptr a, mpz_srcptr b, p_r_values_t p)
{
    /* assumes p will never be > 63 bits */
    int64_t sa = (int64_t) mpz_tdiv_uint64(a, p);
    if (mpz_sgn(a) < 0) {
        sa = -sa;
    }
    return relation_compute_r(sa, mpz_tdiv_uint64(b, p), p);
}

inline char *
u64toa16 (char *p, mpz_srcptr m)
{
    ASSERT(mpz_sgn(m) >= 0); /* m should not be negative */
    mpz_get_str(p, 16, m);
    return p + mpz_sizeinbase(m, 16);
}

inline char *
d64toa16 (char *p, mpz_srcptr m)
{
    mpz_get_str(p, 16, m);
    return p + mpz_sizeinbase(m, 16) + (mpz_sgn(m) < 0 ? 1 : 0);
}

#endif

#endif	/* CADO_RELATION_TOOLS_H */
