#ifndef CADO_UTILS_GETPRIME_H_
#define CADO_UTILS_GETPRIME_H_
#ifndef MAIN
#include "macros.h"
#else
#define ATTRIBUTE_DEPRECATED
#endif

struct prime_info_s {
  unsigned long offset;  /* offset for current primes */
  long current;          /* index of previous prime */
  unsigned int *primes;  /* small primes up to sqrt(p) */
  unsigned long nprimes; /* length of primes[] */
  unsigned char *sieve;  /* sieving table */
  long len;              /* length of sieving table */
  unsigned int *moduli;  /* offset for small primes */
};
typedef struct prime_info_s prime_info[1];

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime function returns successive odd primes, starting with 3. */
extern void prime_info_init (prime_info);
extern void prime_info_clear (prime_info);
extern unsigned long getprime_mt (prime_info);

extern unsigned long getprime (unsigned long) ATTRIBUTE_DEPRECATED;


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GETPRIME_H_ */
