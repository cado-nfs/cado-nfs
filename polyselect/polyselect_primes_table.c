#include "cado.h"
#include <stdlib.h>
#include <inttypes.h>
#include "getprime.h"
#include "misc.h"
#include "polyselect_primes_table.h"
#include "timing.h"


/* init prime array */
/* initialize primes in [P,2*P] */
static unsigned long
initPrimes ( unsigned long P,
             uint32_t **primes )
{
  unsigned long p, nprimes = 0;
  unsigned long Pmax = 2*P;
#ifdef LESS_P // if impatient for root finding
  Pmax = P + P/2;
#endif
  unsigned long maxprimes = nprimes_interval(P, Pmax);

  *primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  prime_info pi;
  prime_info_init (pi);

  /* It's now fairly trivial to parallelize this prime search if we want
   * to, but I think that it's a trivial computation anyway, and most
   * probably not worth the work.
   */
  prime_info_seek(pi, P);

  for (p = P, nprimes = 0; (p = getprime_mt (pi)) <= Pmax; nprimes++) {
    if (nprimes + 1 >= maxprimes) {
      maxprimes += maxprimes / 10;
      *primes = (uint32_t*) realloc (*primes, maxprimes * sizeof (uint32_t));
      if ( (*primes) == NULL) {
        fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
        exit (1);
      }
    }
    (*primes)[nprimes] = p;
  }

  prime_info_clear (pi);

  *primes = (uint32_t*) realloc (*primes, (nprimes) * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  return nprimes;
}


/* clear prime array */
void
polyselect_primes_table_print (polyselect_primes_table_srcptr pt)
{
    uint32_t *primes = pt->Primes;
    unsigned long size = pt->lenPrimes;
    unsigned long i;
    for (i = 0; i < size; i++) {
        fprintf (stderr, "(%lu, %" PRIu32 ") ", i, primes[i]);
        if ((i+1) % 5 == 0)
            fprintf (stderr, "\n");
    }
    fprintf (stderr, "\n");
}

void polyselect_primes_table_init(polyselect_primes_table_ptr pt, unsigned long P)
{
    unsigned long st = milliseconds();

    pt->lenPrimes = initPrimes(P, &pt->Primes);

    printf("# Info: initializing %lu P primes took %lums\n",
            pt->lenPrimes, milliseconds() - st);
}

void polyselect_primes_table_clear(polyselect_primes_table_ptr pt)
{
    free(pt->Primes);
}
