#ifndef CADO_UTILS_GMP_AUX_H_
#define CADO_UTILS_GMP_AUX_H_

#include <gmp.h>
#include <stdint.h>

/* the following function is missing in GMP */
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

#ifdef __cplusplus
extern "C" {
#endif

/* gmp_aux */
extern void mpz_set_uint64 (mpz_t, uint64_t);
extern void mpz_set_int64 (mpz_t, int64_t);
extern uint64_t mpz_get_uint64 (mpz_t);
extern int64_t mpz_get_int64 (mpz_t);
extern void mpz_mul_uint64 (mpz_t a, mpz_t b, uint64_t c);
extern void mpz_mul_int64 (mpz_t a, mpz_t b, int64_t c);
extern uint64_t uint64_nextprime (uint64_t);
extern unsigned long ulong_nextprime (unsigned long);
extern int isprime (unsigned long);
extern void mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d);
extern void mpz_ndiv_q (mpz_t q, mpz_t n, mpz_t d);

/* return the number of bits of p, counting from the least significant end */
extern int nbits (uintmax_t p);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GMP_AUX_H_ */
