#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <gmp.h>
#include "getprime.h"
#include "gmp_aux.h"
#include "rootfinder.h"

int
roots_mod_uint64 (uint64_t *r, uint64_t a, int d, uint64_t p);

/* sort the roots r[0], ..., r[n-1] in increasing order */
static void
sort_roots (uint64_t *r, int n)
{
  int i, j;
  uint64_t t;

  for (i = 1; i < n; i++)
    {
      t = r[i];
      for (j = i; j > 0 && r[j-1] > t; j--)
        r[j] = r[j-1];
      r[j] = t;
    }
}

int 
roots_mod_mpz(uint64_t *r, uint64_t a, int d, uint64_t p)
{
  mpz_t *f;
  int n, i;
  
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  for (i = 0; i <= d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[d], 1);
  mpz_set_uint64 (f[0], p - a);
  n = poly_roots_uint64 (r, f, d, p);
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
  sort_roots (r, n);
  return n;
}


int main(int argc, char **argv) {
  uint64_t *r1, *r2;
  unsigned long a, p, d;
  int n1, n2, i;
  unsigned long minp = 100, maxp=10000, mina=1, maxa=100, mind=1, maxd=10;
  int check = 1;
  
  if (argc > 1 && strcmp(argv[1], "-nc") == 0) {
    printf ("Checking results disabled.\n");
    check = 0;
    argc--;
    argv++;
  }
  if (argc > 1)
    minp = strtoul (argv[1], NULL, 10);
  if (argc > 2)
    maxp = strtoul (argv[2], NULL, 10);
  if (argc > 3)
    mina = strtoul (argv[3], NULL, 10);
  if (argc > 4)
    maxa = strtoul (argv[4], NULL, 10);
  if (argc > 5)
    mind = strtoul (argv[5], NULL, 10);
  if (argc > 6)
    maxd = strtoul (argv[6], NULL, 10);
  
  printf ("minp = %lu, maxp = %lu, mina = %lu, maxa = %lu, mind = %lu, maxd = %lu\n",
          minp, maxp, mina, maxa, mind, maxd);

  r1 = malloc (sizeof(uint64_t) * maxd);
  r2 = malloc (sizeof(unsigned long) * maxd);
  
  for (p = 2; p < minp; p = getprime(p));

  while (p <= maxp) {
    for (a = mina; a <= maxa && a < p; a++) {
      for (d = mind; d <= maxd; d++) {
        n1 = roots_mod_uint64 (r1, a, d, p);
        if (check) {
          n2 = roots_mod_mpz (r2, a, d, p);
          if (n1 != n2) {
            fprintf (stderr, "Error: for a=%lu, d=%lu, p=%lu, roots_mod_uint64()"
                     " reports %d roots, roots_mod_mpz() reports %d\n", 
                     a, d, p, n1, n2);
            exit (EXIT_FAILURE);
          }
          for (i = 0; i < n1 && r1[i] == r2[i]; i++);
          if (i != n1) {
            fprintf (stderr, "Error: for a=%lu, d=%lu, p=%lu, roots_mod_uint64()"
                     " reports roots: ", a, d, p);
            for (i = 0; i < n1; i++)
              fprintf (stderr, "%" PRIu64 " ", r1[i]);
            fprintf (stderr, ", roots_mod_mpz() reports ");
            for (i = 0; i < n2; i++)
              fprintf (stderr, "%" PRIu64 " ", r2[i]);
            fprintf (stderr, "\n");
            exit (EXIT_FAILURE);
          }
        }
      }
    }
    p = getprime(p);
  }
  free (r1);
  free (r2);
  exit (EXIT_SUCCESS);
}
