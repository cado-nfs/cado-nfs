#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int
main (int argc, char *argv[])
{
  if (argc != 2) {
      fprintf(stderr, "Usage: nextprime p\n");
      exit(1);
  }

  mpz_t P;
  mpz_init(P);
  mpz_set_str(P, argv[1], 0);
  mpz_nextprime (P, P);
  gmp_printf ("%Zd\n", P);
  mpz_clear (P);
  return 0;
}
