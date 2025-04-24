#ifndef CADO_ROOTS_MOD_H
#define CADO_ROOTS_MOD_H

#include <stdint.h>
#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Structure to hold factorization of an unsigned long, and to iterate through
   its proper divisors. An integer < 2^64 has at most 15 prime factors. The 
   smallest integer with 16 prime factors is 53# =~ 2^64.8. */
typedef struct {
  unsigned long p[15];
  unsigned char e[15];
  unsigned char c[15];
  unsigned int ndiv;
} enumeratediv_t;

unsigned int roots_mod_uint64 (uint64_t * r, uint64_t a, int d, uint64_t p, gmp_randstate_ptr);
unsigned char factor_ul (unsigned long *, unsigned char *, unsigned long n);
void enumeratediv_init (enumeratediv_t *, unsigned long n);
unsigned long enumeratediv (enumeratediv_t *);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_ROOTS_MOD_H */
