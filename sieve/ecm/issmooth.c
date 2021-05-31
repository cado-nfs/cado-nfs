/* Reads integers "N", or "p N", or "p O N" from stdin and checks N for 
     smoothness, where smooth is defined as
     For given B1, B2, k, d,
     N is smooth if c=N/gcd(N,E) is 1 or is prime and B1 < c <= B2,
     or c is coprime to d and <d/2, or c divides d,
     where E = k * lcm(1, 2, 3, ..., B1).
   By default, k=1 and d=1, so that they have no effect.
   The purpose of the value d is including those N where ECM would find 
     the factor during stage 2 initialisation. 
   With -pmin and -pmax, test only those input where pmin <= p <= pmax,
     (or pmin <= N <= pmax if input has no p value).
   Prints the number of smooth numbers, and with -v option each "p O N" 
     where N is smooth. 
   Several B2 can be specified on the command line (but only one B1 value).
     For each B2, the number of smooth N are counted.
   If an "m" value can be specified, the number of input numbers and smooth 
     (w.r.t. the largest B2) numbers with p==r (mod m) for the possible r
     are also counted and printed.
   Also prints the average exponent of primes up to 19 in N and in O/N 
     (which requires N|O).
*/

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>
#include "macros.h"
#include "gcd.h"

#define NR_EXPONENTS 8
#define EXP_PRIMES {2,3,5,7,11,13,17,19}

static const unsigned long exp_primes[NR_EXPONENTS] = EXP_PRIMES;

static int 
isprime (const unsigned long N, mpz_t t)
{
  mpz_set_ui (t, N);
  return mpz_probab_prime_p (t, 1);
}


static unsigned int
valuation (const unsigned long n, const unsigned long p)
{
  unsigned long c = n;
  unsigned int i = 0;
  while (c % p == 0)
    {
      c /= p;
      i++;
    }
  return i;
}

typedef struct {
  unsigned long input, *smooth;
  unsigned int exp_order[NR_EXPONENTS];
  unsigned int exp_index[NR_EXPONENTS];
} stats_t;


int 
main (int argc, char **argv)
{
  mpz_t E, m_c;
  unsigned long N, O, c, c_gcdiv_d, p, i, j, imin = 0, imax = 0;
  /* Parameters */
  unsigned long B1, *B2, maxB2, k = 1, d = 1, pmin = 0, pmax = ULONG_MAX, 
                m = 1;
  int quiet = 0, verbose = 0;
  unsigned int nr_B2 = 0;
  stats_t *stats;
  int scanf_code;
  char *d_coprime;
  char buf[256];
  
  /* Parse arguments */

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argv[1][1] == 'v')
        {
          verbose = 1;
          argc--;
          argv++;
        }
      else if (argv[1][1] == 'q')
        {
          quiet = 1;
          argc--;
          argv++;
        }
      else if (argc > 2 && argv[1][1] == 'k')
        {
          k = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && argv[1][1] == 'd')
        {
          d = strtoul (argv[2], NULL, 10);
          if (d == 0) {
              fprintf (stderr, "d == 0 does not make sense\n");
              exit (EXIT_FAILURE);
          }
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && argv[1][1] == 'm')
        {
          m = strtoul (argv[2], NULL, 10);
          if (m == 0) {
              fprintf (stderr, "m == 0 does not make sense\n");
              exit (EXIT_FAILURE);
          }
#ifdef __COVERITY__
          /* any non-zero value will do */
          __coverity_mark_pointee_as_sanitized__(&m, GENERIC); // allocation, divisor, loop bound
#endif
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-pmin") == 0)
        {
          pmin = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-pmax") == 0)
        {
          pmax = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (EXIT_FAILURE);
        }
    }

  if (argc < 3)
    {
      printf ("Please provide B1 and B2 values\n");
      exit (EXIT_FAILURE);
    }

  /* Get B1 and all the B2 values */
  B1 = strtoul (argv[1], NULL, 10);
  nr_B2 = argc - 2;
  B2 = malloc(nr_B2 * sizeof(unsigned long));
  maxB2 = 0;
  for (i = 0; i < nr_B2; i++)
    {
      B2[i] = strtoul (argv[2 + i], NULL, 10);
      maxB2 = MAX(B2[i], maxB2);
    }

  d_coprime = malloc (d * sizeof(char));
  for (i = 1; i < d; i++)
    d_coprime[i] = (gcd_ul(i,d) == 1) ? 1 : 0;
  if (d > 1)
    {
      imin = (B1 + d / 2) / d;
      imax = (maxB2 - d / 2) / d;
    }
  
  /* Allocate a stat_t for each residue class (mod m) */
  stats = malloc (m * sizeof (stats_t));
  for (i = 0; i < m; i++)
    {
      stats[i].input = 0;
      for (j = 0; j < NR_EXPONENTS; j++)
        stats[i].exp_order[j] = stats[i].exp_index[j] = 0;

      stats[i].smooth = malloc (nr_B2 * sizeof(unsigned long));
      for (j = 0; j < nr_B2; j++)
        stats[i].smooth[j] = 0;
    }


  /* Compute E */
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  for (i = 2UL; i <= B1; i = (i + 1UL) | 1UL)
    if (mpz_gcd_ui (NULL, E, i) == 1UL) /* A prime? */
      {
        for (j = 1UL; j <= B1 / i; j *= i);
        mpz_mul_ui (E, E, j);
      }
  mpz_mul_ui (E, E, k);

  if (!quiet)
    {
      fprintf (stderr, "B1 = %lu, ", B1);
      for (i = 0; i < nr_B2; i++)
        fprintf (stderr, "B2[%lu] = %lu, ", i, B2[i]);
      gmp_fprintf (stderr, "k = %lu, pmin = %lu, pmax = %lu, E = %Zd\n", 
                   k, pmin, pmax, E);
    }
  

  /* Test input numbers */
  mpz_init (m_c);
  while (!feof(stdin))
    {
      if (fgets (buf, 256, stdin) == NULL)
        break;
      scanf_code = sscanf (buf, "%lu %lu %lu\n", &p, &O, &N);
      if (pmin <= p && p <= pmax)
        {
          unsigned long r;
          int was_smooth = 0, ispr = 0;
          /* input line may contain only "N", or "p N", or "p O N" */
          if (scanf_code == 1)
            N = p;
          else if (scanf_code == 2)
            N = O;

          r = p % m;
          stats[r].input++;

          c = N / mpz_gcd_ui (NULL, E, N);
          c_gcdiv_d = c / gcd_ul (c, d);
          ispr = (B1 < c && c <= maxB2 && isprime(c, m_c));

          for (i = 0; i < nr_B2; i++)
            if (c_gcdiv_d == 1UL || /* c == 1 (factor found after stage 1) || c | d */
                (imin <= c_gcdiv_d && c_gcdiv_d <= imax) || /* factor appears in list of idP */
                (c < d / 2 && d_coprime[c] == 1) || /* factor appears in list of jP */
                (ispr && c <= B2[i])) /* factor appears normally in stage 2 */
              {
                stats[r].smooth[i]++;
                was_smooth = 1;
              }

          /* If order was smooth for any B2 and verbose is set, print it */
          if (was_smooth && verbose)
            printf ("%lu %lu\n", p, N);

          if (scanf_code == 3)
            {
              ASSERT (O % N == 0);
              for (i = 0; i < NR_EXPONENTS; i++)
                {
                  stats[r].exp_order[i] += valuation (O, exp_primes[i]);
                  stats[r].exp_index[i] += valuation (O/N, exp_primes[i]);
                }
            }
        }
    }
  mpz_clear (E);
  mpz_clear (m_c);

  if (!quiet)
    {
      unsigned long r;
      for (r = 0; r < m; r++)
        {
          if (stats[r].input > 0)
            {
              printf ("%lu (mod %lu): %lu input numbers\n", 
                r, m, stats[r].input);
              for (i = 0; i < nr_B2; i++)
                printf ("%lu (mod %lu): %lu %lu,%lu-smooth numbers, ratio %f\n", 
                  r, m, stats[r].smooth[i], B1, B2[i], 
                  (double) stats[r].smooth[i] / (double) stats[r].input);
            }
              
          if (stats[r].exp_order[0] > 0)
            {
              printf ("Sum of valuation of q in O:\n");
              for (i = 0; i < NR_EXPONENTS; i++)
                if (stats[r].exp_order[i])
                  printf ("q = %lu: %u  ", exp_primes[i], stats[r].exp_order[i]);
              printf ("\n");
            }
          if (stats[r].exp_index[0] > 0)
            {
              printf ("Sum of valuation of q in N/O (index of starting element):\n");
              for (i = 0; i < NR_EXPONENTS; i++)
                if (stats[r].exp_index[i])
                  printf ("q = %lu: %u  ", exp_primes[i], stats[r].exp_index[i]);
              printf ("\n");
            }
        }
    }
  
  /* Free memory */
  free (B2);
  B2 = NULL;
  free (d_coprime);
  d_coprime = NULL;
  for (i = 0; i < m; i++)
    {
      free(stats[i].smooth);
      stats[i].smooth = NULL;
    }
  free (stats);
  stats = NULL;

  return 0;
}
