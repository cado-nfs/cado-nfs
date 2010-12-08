#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include "utils.h"

/* #define DEBUG */

#ifdef DEBUG /* compute exponent of given prime */
#define DEBUG_PRIME 101921
#endif

#define FAST /* uses product tree instead of naive accumulation */

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 + 
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
mpz_t*
accumulate_fast (mpz_t *prd, mpz_t a, unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  mpz_mul (prd[0], prd[0], a);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (mpz_t*) realloc (prd, *lprd * sizeof (mpz_t));
          mpz_init_set_ui (prd[i + 1], 1);
        }
      mpz_mul (prd[i + 1], prd[i + 1], prd[i]);
      mpz_set_ui (prd[i], 1);
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
void
accumulate_fast_end (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    mpz_mul (prd[0], prd[0], prd[i]);
}

int
main (int argc, char **argv)
{
  FILE * matfile, *kerfile = NULL;
  mpz_t m, a, b, c, *prd;
  unsigned long lprd; /* number of elements in prd[] */
  unsigned long nprd; /* number of accumulated products in prd[] */
  int ret;
  unsigned long w;
  int i, j, nlimbs, nab;
  char str[1024];
  int depnum = 0;
  cado_poly pol;
#ifdef DEBUG
  unsigned long debug_exponent = 0;
#endif

  fprintf (stderr, "%s revision %s\n", argv[0], CADO_REV);

  if (argc > 2 && strcmp (argv[1], "-depnum") == 0)
    {
      depnum = atoi (argv[2]);
      argc -= 2;
      argv += 2;
    }

  if (argc != 3 && argc != 4)
    {
      fprintf (stderr, "usage: %s [-depnum nnn] matfile kerfile polyfile\n", argv[0]);
      /* The following is a way to perform the square root after msieve's
         linear algebra. It assumes ab_file contains all (a,b) pairs of the
         dependency, one per line, with a and b separated by a space,
         and a header "nrows ncols" where nrows is the number of (a,b) pairs,
         and ncols is any integer.
         To produce ab_file, do: "msieve -nc3 k,k N > ab_file", where k
         is the dependency number (1 <= k <= 64).
      */
      fprintf (stderr, "usage: %s ab_file polyfile\n", argv[0]);
      exit (1);
    }

  matfile = fopen(argv[1], "r");
  ASSERT (matfile != NULL);
  if (argc == 4)
    {
      kerfile = fopen (argv[2], "r");
      ASSERT (kerfile != NULL);
    }
  /* otherwise kerfile = NULL */
  cado_poly_init (pol);
  ret = cado_poly_read (pol, argv[argc - 1]);
  ASSERT (ret);

  mpz_init (m);
  mpz_neg (m, pol->g[0]);

  gmp_fprintf (stderr, "m = %Zd\n", m);

  lprd = 1;
  nprd = 0;
  prd = (mpz_t*) malloc (lprd * sizeof (mpz_t));
  mpz_init_set_ui(prd[0], 1);
  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  if (kerfile != NULL)
    {
      int nrows, ncols;
      ret = fscanf (matfile, "%d %d", &nrows, &ncols);
      ASSERT (ret == 2);
      fgets (str, 1024, matfile); // read end of first line
      nlimbs = (nrows / GMP_NUMB_BITS) + 1;

      /* go to dependency depnum */
      while (depnum > 0)
        {
          int c;
          /* read one line */
          while ((c = fgetc (kerfile)) != '\n')
            if (c == EOF)
              break;
          depnum --;
        }

      if (depnum > 0)
        {
          fprintf (stderr, "Error, not enough dependencies\n");
          exit (1);
        }
    }
  else
    nlimbs = INT_MAX;

  for (nab = 0, i = 0; i < nlimbs; ++i) {
    if (kerfile != NULL)
      {
        ret = fscanf (kerfile, "%lx", &w);
        ASSERT (ret == 1);
      }
    else
      w = ~0UL; /* trick so that all relations are considered */
    for (j = 0; j < GMP_NUMB_BITS; ++j) {
      if (fgets(str, 1024, matfile) != NULL) {
	if (w & 1UL) {
	  ret = gmp_sscanf(str, "%Zd %Zd", a, b);
          nab ++;
          if(!(nab % 100000))
            fprintf (stderr, "# Read ab pair #%d at %2.2lf\n", nab, seconds ());
	  ASSERT (ret == 2);
          /* FIXME: instead of accumulating a-b*m, where m = -g0/g1 (mod N),
             and accumulate g1^nab on the algebraic side, we could accumulate
             g1*a+g0*b on the rational side. */
	  mpz_mul (c, b, m);
	  mpz_sub (c, a, c);
#ifndef FAST
          mpz_mul (prd[0], prd[0], c);
#else
          prd = accumulate_fast (prd, c, &lprd, nprd++);
#endif
#ifdef DEBUG
          if (mpz_divisible_ui_p (c, DEBUG_PRIME))
            gmp_printf ("prime %lu appears in relation (%Zd,%Zd)\n",
                        DEBUG_PRIME, a, b);
          while (mpz_divisible_ui_p (c, DEBUG_PRIME))
            {
              mpz_divexact_ui (c, c, DEBUG_PRIME);
              debug_exponent ++;
            }
#endif          
	}
      }
      else
        {
          i = nlimbs - 1;
          break; /* end of file */
        }
      w >>= 1;
    }
  }
  fprintf (stderr, "# Read %d relations\n", nab);

  /* check the number of relations is even */
  if (nab & 1)
    {
      fprintf (stderr, "Error, odd number of relations\n");
      exit (1);
    }

#ifdef DEBUG
  printf ("exponent of %lu is %lu\n", DEBUG_PRIME, debug_exponent);
#endif

#ifdef FAST
  accumulate_fast_end (prd, lprd);
#endif

  fprintf(stderr, "Size of product = %lu bits\n", mpz_sizeinbase (prd[0], 2));

  if (mpz_sgn (prd[0]) < 0)
    {
      fprintf (stderr, "Error, product is negative: try another dependency\n");
      exit (1);
    }

#ifdef DEBUG2 /* print all primes with odd exponent */
  {
#define MAX_ACCU 1048576
    unsigned long p, P[MAX_ACCU], E[MAX_ACCU], f, F[MAX_ACCU];
    int i, j;
    mpz_t t, q, u;
    double d;
    mpz_init (t);
    mpz_init (q);
    mpz_init (u);
    p = 2;
    /* since we know we have a square, take the square root */
    mpz_sqrtrem (prd[0], t, prd[0]);
    ASSERT_ALWAYS(mpz_cmp_ui (t, 0) == 0);
    while (mpz_cmp_ui (prd[0], 1) > 0)
      {
        fprintf (stderr, "# Remains %lu bits\n", mpz_sizeinbase (prd[0], 2));
        /* accumulate small primes in q */
        i = 0;
        d = 0.0;
        mpz_set_ui (q, 1);
        while (i < MAX_ACCU && d < 1.0)
          {
            mpz_mul_ui (q, q, p);
            P[i] = p;
            E[i] = 0;
            i ++;
            d += log((double) p) * log((double) p) / (double) p;
            do
              {
                p = p + 1 + (p > 2);
              }
            while (isprime (p) == 0);
          }
        f = 1; /* invariant: q = (P[0] * P[1] * ... * P[i-1])^f */
        while (mpz_divisible_p (prd[0], q))
          {
            mpz_mul (q, q, q);
            f *= 2;
          }
        if (f > 1)
          {
            mpz_sqrt (q, q);
            f /= 2; /* now f >= 1 */
          }
        for (j = 0; j < i; j++)
          F[j] = f;
        /* invariant: q = P[0]^F[0] * P[1]^F[1] * ... * P[i-1]^F[i-1] */
        while (mpz_cmp_ui (q, 1) > 0)
          {
            mpz_fdiv_qr (t, u, prd[0], q);
            if (mpz_cmp_ui (u, 0) == 0) /* division is exact */
              {
                for (j = 0; j < i; j++)
                  E[j] += F[j];
                mpz_set (prd[0], t);
              }
            else /* some P[j]^F[j] do not divide fully: divide F[j] by 2 */
              {
                for (j = 0; j < i; j++)
                  {
                    if (F[j] == 0)
                      continue;
                    mpz_ui_pow_ui (t, P[j], F[j]);
                    while (mpz_divisible_p (u, t) == 0 && F[j] > 0)
                      {
                        if (F[j] >= 2)
                          mpz_sqrt (t, t);
                        /* else F[j]=1: t = P[j] */
                        mpz_divexact (q, q, t);
                        F[j] /= 2;
                      }
                  }
              }
          }
        for (j = 0; j < i; j++)
          printf ("%lu %lu\n", P[j], 2 * E[j]);
        fflush (stdout);
      }
    mpz_clear (t);
    mpz_clear (q);
    mpz_clear (u);
  }
#endif
  mpz_sqrtrem (prd[0], a, prd[0]);
  {
    size_t la = mpz_sizeinbase (a, 2);
    if (la <= 100)
      gmp_fprintf (stderr, "remainder is %Zd\n", a);
    else
      fprintf (stderr, "remainder has %lu bits\n", la);
  }

  mpz_mod(prd[0], prd[0], pol->n);
  gmp_fprintf(stderr, "rational square root is %Zd\n", prd[0]);
  gmp_printf("%Zd\n", prd[0]);

  for (i = 0; i < (int) lprd; i++)
    mpz_clear (prd[i]);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (c);
  free (prd);

  return 0;
}
