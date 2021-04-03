/* Dynamic Eratosthenes sieve.

  Copyright 2001, 2002, 2003, 2005 Paul Zimmermann and Alexander Kruppa.
  (Modified wrt GMP-ECM to use 'unsigned long' instead of 'double'.)

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

/* compile with -DMAIN to use as a standalone program */

#ifndef MAIN
#include "cado.h"		// IWYU pragma: keep
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "getprime.h"
#include "macros.h"
#ifndef MAIN
#endif

/* provided for in cado.h, but we want getprime.c to be standalone */
#ifndef ASSERT
#define ASSERT(x)
#endif

/* This function returns successive odd primes, starting with 3.
   To perform a loop over all primes <= B1, do the following
   (compile this file with -DMAIN to count primes):

      prime_info pi;
      prime_info_init (pi);
      for (p = 2; p <= B1; p = getprime_mt (pi))
         {
            ...
         }

      prime_info_clear (pi);

  If you prefer to have 2 returned as well, you may use the
  prime_info_init_seek call as follows:

  prime_info pi;
  prime_info_init_seek(pi, 2);  // depart from default lower bound, 3
  for (unsigned long p; (p = getprime_mt (pi)) <= B1; ) {
    ...
  }
  prime_info_clear (pi);

*/

void prime_info_init(prime_info_ptr pi)
{
    /* pi->offset must always be odd */
    pi->offset = 3;
    pi->current = 0;
    pi->primes = NULL;
    pi->nprimes = 0;
    pi->sieve = NULL;
    pi->len = 0;
    pi->moduli = NULL;
}

void prime_info_init_seek(prime_info_ptr pi, unsigned long lower_bound)
{
    prime_info_init(pi);
    prime_info_seek(pi, lower_bound);
}

void prime_info_clear(prime_info_ptr pi)
{
    free(pi->primes);
    free(pi->sieve);
    free(pi->moduli);
}

/* this function is not thread-safe, as it uses a global state */
unsigned long getprime(unsigned long p)
{
    static prime_info pi;
    static int initialized = 0;

    if (p == 0) {
	prime_info_clear(pi);
	initialized = 0;
	return p;
    }

    if (initialized == 0) {
	prime_info_init(pi);
	initialized = 1;
    }

    return getprime_mt(pi);
}


/*
  unsigned long offset;  // offset for current primes
  long current;          // index of previous prime
  unsigned int *primes;  // small primes up to sqrt(p)
  unsigned long nprimes; // length of primes[]
  unsigned char *sieve;  // sieving table
  long len;              // length of sieving table
  unsigned int *moduli;  // offset for small primes
 */

void prime_info_seek(prime_info_ptr pi, unsigned long lower_bound)
{
    /* note that pi->offset must be odd
     */
    if (lower_bound >= 3 && !(lower_bound & 1))
	lower_bound++;
    pi->offset = lower_bound;

    /* Setting pi->len and pi->primes to zero will trigger
     * reinitialization of the sieve array and prime list
     */
    pi->len = 0;
    pi->nprimes = 0;
}

/* This initializes only the list of small primes to something minimal.
 *
 * The list of moduli is initialized according to pi->offset.
 *
 * This does not touch the sieve table
 */
void prime_info_small_primes_table_bootstrap(prime_info_ptr pi)
{
    pi->nprimes = 1;

    pi->primes =
	(unsigned int *) realloc(pi->primes,
				 pi->nprimes * sizeof(unsigned int));
    /* assume this "small" malloc will not fail in normal usage */
    ASSERT(pi->primes != NULL);

    pi->moduli =
	(unsigned int *) realloc(pi->moduli,
				 pi->nprimes * sizeof(unsigned int));
    /* assume this "small" malloc will not fail in normal usage */
    ASSERT(pi->moduli != NULL);

    if (pi->offset <= 3)
	pi->offset = 3;

    pi->primes[0] = 3;
    unsigned int p = 3;
    unsigned int j = pi->offset % p;
    j = (j == 0) ? j : p - j;	/* -offset mod p */
    if ((j % 2) != 0)
	j += p;			/* ensure j is even */
    pi->moduli[0] = j / 2;
    if (pi->offset == 3)
	pi->moduli[0] += 3;

}

unsigned long prime_info_sieving_max_reach(prime_info_srcptr pi)
{
    unsigned long p = pi->primes[pi->nprimes - 1];
    return p * p;
}

/* This doubles (possibly several times, if we used prime_info_seek) the
 * size of the small primes, and adjust the moduli according to
 * pi->offset
 */
void prime_info_small_primes_table_expand(prime_info_ptr pi)
{
    for (; prime_info_sieving_max_reach(pi) < pi->offset + 2 * pi->len;) {
	unsigned int k, p, j, ok;
	k = pi->nprimes;
	pi->nprimes *= 2;
	pi->primes = (unsigned int *) realloc(pi->primes, pi->nprimes *
					      sizeof(unsigned int));
	pi->moduli = (unsigned int *) realloc(pi->moduli, pi->nprimes *
					      sizeof(unsigned int));
	/* assume those "small" realloc's will not fail in normal usage */
	ASSERT(pi->primes != NULL && pi->moduli != NULL);
	for (p = pi->primes[k - 1]; k < pi->nprimes; k++) {
	    /* find next (odd) prime > p ; trial divide by all smaller
	     * primes, which is sub-optimal, but shouldn't hurt too much. */
	    do {
		for (p += 2, ok = 1, j = 0; (ok != 0) && (j < k); j++)
		    ok = p % pi->primes[j];
	    }
	    while (ok == 0);
	    pi->primes[k] = p;
	    /* moduli[k] is the smallest m such that
	     * offset + 2*m = 0 mod primes[k] */
	    j = pi->offset % p;
	    j = (j == 0) ? j : p - j;	/* -offset mod p */
	    if ((j % 2) != 0)
		j += p;		/* ensure j is even */
	    pi->moduli[k] = j / 2;
	}
    }
}

/* we want to enforce:
 *    len >= 1,
 *    len^2 >= offset,
 *    pi->sieve allocated of size len+1, with an end mark
 *
 * note that this is a no-op if post-conditions are already met.
 */
void prime_info_sieving_table_expand(prime_info_ptr pi)
{
    /* enlarge sieving table if too small */
    if (pi->len && (unsigned long) pi->len * pi->len >= pi->offset)
	return;

    pi->len += pi->len == 0;
    for (; (unsigned long) pi->len * pi->len < pi->offset; pi->len *= 2);

    pi->sieve = (unsigned char *) realloc(pi->sieve, pi->len + 1);
    /* assume this "small" malloc will not fail in normal usage */
    ASSERT(pi->sieve != NULL);

    pi->sieve[pi->len] = 1;	/* End mark */
}

void prime_info_sieve_more(prime_info_ptr pi)
{
    /* now sieve for new primes */
    long k;
    unsigned long j, p;

    memset(pi->sieve, 1, sizeof(unsigned char) * (pi->len + 1));
    for (j = 0; j < pi->nprimes; j++) {
	p = pi->primes[j];
	for (k = pi->moduli[j]; k < pi->len; k += p)
	    pi->sieve[k] = 0;
	pi->moduli[j] = k - pi->len;	/* for next sieving array */
    }
}

unsigned long prime_info_sieve_find(prime_info_ptr pi)
{
    for (; pi->current < pi->len; pi->current++) {
	if (pi->sieve[pi->current]) {
	    unsigned long p = pi->offset + 2 * pi->current;
	    /* for next call */
	    pi->current++;
	    /* most calls will end here */
	    return p;
	}
    }
    return 0;
}

/* this function is thread-safe, provided of course that no two threads
 * tinker with the same prime_info data. */
unsigned long getprime_mt(prime_info_ptr pi)
{
    unsigned long p;
    p = prime_info_sieve_find(pi);
    if (p)
	return p;

    if (pi->offset <= 2) {
	pi->offset = 3;
	return 2;
    }

    /* otherwise we have to sieve, and maybe find new primes first. */

    pi->offset += 2 * pi->len;

    /* we want to enforce:
     *    len >= 1,
     *    len^2 >= offset,
     *    last_prime^2 >= offset+2*len
     */

    if (pi->nprimes == 0)
	prime_info_small_primes_table_bootstrap(pi);

    prime_info_sieving_table_expand(pi);

    prime_info_small_primes_table_expand(pi);

    prime_info_sieve_more(pi);

    pi->current = 0;

    p = prime_info_sieve_find(pi);
    if (p)
	return p;

    /* otherwise we found a prime gap >= sqrt(x) around x */
    ASSERT_ALWAYS(0);

    return 0;
}

#ifdef MAIN
int main(int argc, char *argv[])
{
    unsigned long p, B;
    unsigned long ii = 0;
    int diff = 0;
    prime_info pi;

    if (argc != 2) {
	fprintf(stderr, "Usage: getprime <bound>\n");
	exit(EXIT_FAILURE);
    }

    B = strtoul(argv[1], NULL, 0);

    prime_info_init(pi);

    for (ii = 0, p = 2; p <= B; ii++) {
	unsigned long newp = getprime_mt(pi);
	if (newp - p > diff)
	    printf("firstdiff(%d)=%lu\n", diff = newp - p, p);
	p = newp;
    }
    printf("pi(%lu)=%lu, maxdiff=%d\n", B, ii, diff);

    prime_info_clear(pi);	/* free the tables */

    return 0;
}
#endif
