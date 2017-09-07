#ifndef TRIALDIV_H_
#define TRIALDIV_H_

#include <gmp.h>

/* The maximum number of words in the numbers to be trial-divided.
   The $l$ value from the thesis text is equal to TRIALDIV_MAXLEN - 1 */
#ifndef TRIALDIV_MAXLEN
#define TRIALDIV_MAXLEN 8
#endif

typedef struct {
  unsigned long p;
  unsigned long w[TRIALDIV_MAXLEN - 1]; /* w[i] = w^{i+1} mod p */
  unsigned long pinv; /* pinv == 1/p (mod w) */
  unsigned long plim; /* plim = (w-1)/p */
} trialdiv_divisor_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Divides primes in d out of N and stores them (with multiplicity) in f */
size_t trialdiv (unsigned long * f, mpz_t N, const trialdiv_divisor_t *d,
		 const size_t max_div);

/* Get the largest candidate factor that can be trial divided */
unsigned long trialdiv_get_max_p();

/* Initialise a trialdiv_divisor_t array with the nr primes stored in an 
   unsigned long array. Allocates memory. */
trialdiv_divisor_t *trialdiv_init (const unsigned long *f, const unsigned int nr);

/* Frees the memory allocated by trialdiv_init() */
void trialdiv_clear (trialdiv_divisor_t *d);

#ifdef __cplusplus
}
#endif

#endif	/* TRIALDIV_H_ */
