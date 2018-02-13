#ifndef CADO_UTILS_GMP_AUX_H_
#define CADO_UTILS_GMP_AUX_H_

#include <gmp.h>
#include <stdint.h>
#include <macros.h>
#include <stdbool.h>
#include "getprime.h"

/* the following function are missing in GMP */
#ifndef mpz_addmul_si
#define mpz_addmul_si(a, b, c)                  \
  do {                                          \
    if (c >= 0)                                 \
      mpz_addmul_ui (a, b, c);                  \
    else                                        \
      mpz_submul_ui (a, b, -(c));               \
  }                                             \
  while (0)
#endif

#ifndef mpz_add_si
#define mpz_add_si(a,b,c)                       \
  if (c >= 0) mpz_add_ui (a, b, c);             \
  else mpz_sub_ui (a, b, -(c))
#endif

#ifndef mpz_submul_si
#define mpz_submul_si(a,b,c)                    \
  if (c >= 0) mpz_submul_ui (a, b, c);          \
  else mpz_addmul_ui (a, b, -(c))
#endif
  
#ifdef __cplusplus
extern "C" {
#endif

/* gmp_aux */
extern void mpz_set_uint64 (mpz_t, uint64_t);
extern void mpz_set_int64 (mpz_t, int64_t);
extern uint64_t mpz_get_uint64 (mpz_srcptr);
extern int64_t mpz_get_int64 (mpz_srcptr);
extern void mpz_mul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c);
extern void mpz_addmul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c);
extern void mpz_submul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c);
extern void mpz_submul_int64 (mpz_t a, mpz_srcptr b, int64_t c);
extern void mpz_divexact_uint64 (mpz_t a, mpz_srcptr b, uint64_t c);
extern void mpz_mul_int64 (mpz_t a, mpz_srcptr b, int64_t c);
extern void mpz_addmul_int64 (mpz_t a, mpz_srcptr b, int64_t c);
extern int mpz_fits_uint64_p(mpz_srcptr);
extern int mpz_fits_int64_p(mpz_srcptr);
extern unsigned long ulong_nextprime (unsigned long);
extern uint64_t uint64_nextprime (uint64_t);
extern int ulong_isprime (unsigned long);
extern unsigned long ulong_nextcomposite (unsigned long, unsigned long);
extern void mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, const mpz_t d);
extern void mpz_ndiv_qr_ui (mpz_t q, mpz_t r, mpz_t n, unsigned long int d);
extern void mpz_ndiv_q (mpz_t q, mpz_t n, const mpz_t d);
extern void mpz_ndiv_q_ui (mpz_t q, mpz_t n, unsigned long int d);
extern int mpz_divisible_uint64_p (mpz_t a, uint64_t c);
extern int mpz_coprime_p (mpz_t a, mpz_t b);

/* Put in r the smallest legitimate value that it at least s + diff (note
   that if s+diff is already legitimate, then r = s+diff will result.

   Here, legitimate means prime or squarefree composite, with the constraint
   that all the prime factors must be in [pmin, pmax[ .

   The prime factors of r are put in factor_r, and the number of them is
   returned. The caller must have allocated factor_r with enough space.
   */
int next_mpz_with_factor_constraints(mpz_t r, unsigned long factor_r[],
        const mpz_t s, const unsigned long diff, unsigned long pmin,
        unsigned long pmax);

/* return the number of bits of p, counting from the least significant end */
extern int nbits (uintmax_t p);
extern long double mpz_get_ld (mpz_t z);

extern int mpz_p_valuation(mpz_srcptr a, mpz_srcptr p);
extern int mpz_p_valuation_ui(mpz_srcptr a, unsigned long p);

#if !GMP_VERSION_ATLEAST(5,0,0)
mp_limb_t mpn_neg (mp_limb_t *rp, const mp_limb_t *sp, mp_size_t n);
void mpn_xor_n (mp_limb_t *rp, const mp_limb_t *s1p, const mp_limb_t *s2p,
		mp_size_t n);
void mpn_zero(mp_limb_t *rp, mp_size_t n);
void mpn_copyi(mp_limb_t *rp, const mp_limb_t * up, mp_size_t n);
void mpn_copyd(mp_limb_t *rp, const mp_limb_t * up, mp_size_t n);
#endif

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GMP_AUX_H_ */
