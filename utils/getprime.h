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
typedef struct prime_info_s * prime_info_ptr;
typedef const struct prime_info_s * prime_info_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime function returns successive odd primes, starting with 3. */
extern void prime_info_init (prime_info_ptr);

/* this second intialization function makes it possible to specify a
 * lower bound that is different from 3
 */
extern void prime_info_init_seek(prime_info_ptr pi, unsigned long lower_bound);
extern void prime_info_clear (prime_info_ptr);
extern void prime_info_seek (prime_info_ptr, unsigned long lower_bound);
extern unsigned long getprime_mt (prime_info_ptr);

extern unsigned long getprime (unsigned long) ATTRIBUTE_DEPRECATED;


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GETPRIME_H_ */
