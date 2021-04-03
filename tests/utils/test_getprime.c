#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "getprime.h"
#include "macros.h"

  int verbose = 0;

void
test_getprime (unsigned long count, unsigned long exp_max)
{
  unsigned long p;
  unsigned long i;
  prime_info pi;

  /* there is no really something random we can test here, since
     getprime can start only at 2 */
  prime_info_init (pi);
  for (p = 2, i = 0; i < count; p = getprime_mt (pi), i++) {
      if (verbose)
          printf("%lu\n", p);
  }
  if (p != exp_max) {
      fprintf(stderr, "Reached %lu-th prime = %lu, while we expected %lu instead\n",
              count, p, exp_max);
      ASSERT_ALWAYS (p == exp_max);
  }
  prime_info_clear (pi);
}

void
test_getprime_range (unsigned long lower_bound, unsigned long upper_bound, unsigned long exp_count)
{
  unsigned long p;
  unsigned long i;
  prime_info pi;

  /* This returns the number of primes such that lower_bound <= p <
   * upper_bound.
   */
  prime_info_init (pi);

  prime_info_seek(pi, lower_bound);

  for (p = lower_bound, i = 0; (p = getprime_mt (pi)) < upper_bound; i++) {
      if (verbose)
          printf("%lu\n", p);
  }
  if (i != exp_count) {
      fprintf(stderr, "Found %lu primes in interval (%lu,%lu) while we expected %lu instead\n",
              i, lower_bound, upper_bound, exp_count);
      ASSERT_ALWAYS (i == exp_count);
  }
  prime_info_clear (pi);
}

/* either use:
 *      - count and exp-max (default)
 *      - bound and count (and, optionally, seek)
 */
int
main (int argc, char * argv[])
{
  unsigned long count = 1000000;
  unsigned long exp_max = 15485867;
  unsigned long bound = 0;
  unsigned long seek = 0;

    for(argc--,argv++ ; argc ; argc--,argv++) {
        if (strcmp(argv[0], "-v") == 0) {
            verbose = 1;
            continue;
        }
        if (argc > 1 && strcmp(argv[0], "-seek") == 0) {
            seek = atol(argv[1]);
            argc--,argv++;
            continue;
        }
        if (argc > 1 && strcmp(argv[0], "-bound") == 0) {
            bound = atol(argv[1]);
            argc--,argv++;
            continue;
        }
        if (argc > 1 && strcmp(argv[0], "-count") == 0) {
            count = atol(argv[1]);
            argc--,argv++;
            continue;
        }
        if (argc > 1 && strcmp(argv[0], "-exp-max") == 0) {
            exp_max = atol(argv[1]);
            argc--,argv++;
            continue;
        }
        fprintf(stderr, "Unexpected argument: %s\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (bound && count) {
        test_getprime_range(seek, bound, count);
    } else if (count && exp_max) {
        test_getprime (count, exp_max);

        if (count <= 1000000) {
            /* Do it a second time for tiny tests. We use this to check that
             * prim_info_{init,clear} behave as they should */
            test_getprime (count, exp_max);
        }
    }

  exit (EXIT_SUCCESS);
}
