/*****************************************************************
 *         Additional functions for the factor base              * 
 *****************************************************************/
/* The functions here used to be used by the linesiever, which is the
 * predecessor of the lattice siever las. Since las has its own data
 * types for most tasks, the typedefs here are outdated.
 */

#ifndef LINESIEVE_FB_H
#define LINESIEVE_FB_H

#include <stdio.h>
#include "cado_poly.h"
#include "fb.h"

/*********************************************************************/

#define SIEVE_BLOCKING_MAX 5

/* Factorbase entry with prime < 2^24 and exactly 1 root */
typedef struct {
  fbprime_t p;            /* A prime or a prime power < 2^24 */
  fbroot_t root_and_log;  /* The root and approximate log (upper 8 bits) */
} factorbase_small_t;

/* Factorbase entry with prime < 2^24 and index to next sieve array
   location where this prime divides  */
typedef struct {
  fbprime_t p;            /* A prime or a prime power < 2^24 */
  fbprime_t loc_and_log;  /* Next location in sieve array where this prime
			     will divide (low 24 bits) and approximate log
			     (upper 8 bits) */
} factorbase_small_inited_t;


typedef struct {
  factorbase_degn_t *fullfb; /* The complete factor base */
  factorbase_degn_t *fblarge; /* Pointer to the entries in fullfb with primes
  				 > fbL2bound */
  factorbase_degn_t *firstlog[256]; /* Pointers to the first entry in 
                                       fullfb with plog == i */
  factorbase_small_t *fbsmall[SIEVE_BLOCKING_MAX];
  factorbase_small_inited_t *fbinit[SIEVE_BLOCKING_MAX];
  unsigned int fbsmallsize[SIEVE_BLOCKING_MAX];  /* Number of entries in small
                                                    fb, incl. stop marker */
  fbprime_t fbsmallbound[SIEVE_BLOCKING_MAX]; /* Upper bound on primes in 
                                                 small fb */
} factorbase_t[1];

void            fb_init_firstlog (factorbase_t);
void		fb_disable_roots (factorbase_degn_t *, const unsigned long, 
                                  const int);
void		fb_restore_roots (factorbase_degn_t *, const unsigned long, 
                                  const int);
void 		fb_extract_small (factorbase_t, const unsigned int, const int,
                                  const int);
int             fb_check (factorbase_t, cado_poly, int);
void            fb_clear (factorbase_t);

/* Number of blocking levels for small factor base primes, should
   correspond to cache levels. Sieving will be done in SIEVE_BLOCKING + 1
   passes: SIEVE_BLOCKING passes updating directly, and one pass with
   bucket sorting. Bucket sorting not implemented atm. */
#define SIEVE_BLOCKING 2
static const unsigned long CACHESIZES[SIEVE_BLOCKING] = {32768, 1048576};

#endif  /* LINESIEVE_FB_H */
