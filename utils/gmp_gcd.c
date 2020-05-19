#include "cado.h"
#include <gmp.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "usage: %s integer integer\n", argv[0]);
		return 1;
	}

	mpz_t a, b, gcd;
	mpz_init_set_str (a,argv[1],10);
	mpz_init_set_str (b,argv[2],10);
	mpz_init (gcd);

	mpz_gcd( gcd, a, b);
	gmp_printf ("%Zd\n", gcd);

    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (gcd);
	return 0;
}

