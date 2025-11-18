#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "rootfinder.h"
#include "tests_common.h"
#include "mpz_poly.h"
#include "macros.h"
#include "gmp_aux.h"
#include "misc.h"
#include "portability.h" // IWYU pragma: keep

int cmp(mpz_t * a, mpz_t * b)
{
    return mpz_cmp(*a, *b);
}

struct test_rootfinder_example {
    const char * p;
    const char * poly;
    int expected_nroots;
};

void test(mpz_srcptr p, mpz_poly_srcptr F, int expected_nroots)
{
    mpz_t * r = NULL;
    unsigned long n = mpz_poly_roots_gen (&r, F, p, state);

    if (tests_common_get_verbose()) {
        char * a;
        mpz_poly_asprintf(&a, F);
        gmp_printf("%s mod %Zd\n", a, p);
        free(a);
    }

    if (mpz_probab_prime_p (p, 5) && mpz_sizeinbase (p, 2) <= 64) {
        unsigned long n1 = mpz_poly_roots_uint64 (NULL, F, mpz_get_uint64 (p), state);
        ASSERT_ALWAYS(n1 == n);
    }

    if (tests_common_get_verbose()) {
        for (unsigned int i = 0; i < n; i++)
            gmp_printf("\t%Zd\n", r[i]);
    }

    if (expected_nroots != -1) {
        mpz_t v;
        mpz_init(v);
        ASSERT_ALWAYS((int) n == expected_nroots);
        for (unsigned int i = 0; i < n; i++) {
            mpz_poly_eval (v, F, r[i]);
            ASSERT_ALWAYS(mpz_divisible_p (v, p));
        }
        mpz_clear (v);
    }

    for (unsigned int i = 0; i < n; i++)
        mpz_clear (r[i]);
    free (r);
}

/* Check roots of polynomial ff[d]*x^d + ... + ff[1]*x + ff[0] mod pp.
   If pp is the empty string, generate a random integer (and then generate
   random coefficients for ff[]).
   If nroots <> -1, it is the expected number of roots. */
void
test_one_example (const struct test_rootfinder_example * ex)
{
    mpz_t p;
    mpz_poly F;

    mpz_init (p);
    mpz_poly_init(F, -1);

    ASSERT_ALWAYS(ex->p);
    ASSERT_ALWAYS(ex->poly);

    mpz_set_from_expression(p, ex->p);
    mpz_poly_set_from_expression(F, ex->poly);

    test(p, F, ex->expected_nroots);

    mpz_clear(p);
    mpz_poly_clear(F);
}

void test_random_p_and_poly(int d, int max_bits, gmp_randstate_t state)
{
    mpz_t p;
    mpz_poly F;

    mpz_init(p);
    mpz_poly_init(F, -1);

    for(;;) {
        mpz_urandomb (p, state, max_bits);
        if (mpz_even_p(p)) continue;
        mpz_poly_set_randomb(F, d, state, max_bits,
                MPZ_POLY_URANDOM |
                MPZ_POLY_DEGREE_EXACT
                );
        if (F->deg == -1) continue;
        if (mpz_divisible_p(mpz_poly_lc(F), p)) continue;
        break;
    }
    test(p, F, -1);
    mpz_clear(p);
    mpz_poly_clear(F);
}

void usage()
{
    fprintf(stderr, "Usage: test_rootfinder [-seed nnn] [-iter nnn] [-v] [[-degree <d>] | <p> <polynomial in x>]\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char const * argv[])
{
    const struct test_rootfinder_example examples[] = {
	{ "4294967291", "-3*x^2 + 1", 2 },
	{ "18446744073709551557", "5*x^3 + 3*x^2 + 2*x + 1", 1 },
	{ "18446744073709551629", "-x^3 + 7*x^2 - x + 1", 0 },
	{ "246089", "4*x^3 + 3*x^2 + 2*x + 1", 27 },
	{ "9", "x^3 + 1", 3 },
	{ NULL, NULL, -1 },
    };
    unsigned long iter = 100;
    int max_bits = 40;

    tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);
    tests_common_get_iter (&iter);

    const char * pstr = NULL;
    const char * polystr = NULL;
    int d = -1;

    for(int i = 1 ; i < argc ; i++) {
        if (strcmp(argv[i], "-degree") == 0 && i+1 < argc) {
            char * t;
            d = (int) strtol(argv[i+1], &t, 0);
            if (*t || d <= 0) usage();
            i++;
        } else if (strcmp(argv[i], "-bits") == 0 && i+1 < argc) {
            char * t;
            max_bits = (int) strtol(argv[i+1], &t, 0);
            if (*t || max_bits <= 0) usage();
            i++;
        } else if (!pstr) {
            pstr = argv[i];
        } else if (!polystr) {
            polystr = argv[i];
        } else
            usage();
    }

    if ((pstr == NULL) != (polystr == NULL))
        usage();

    if (d > -1 && pstr)
        usage();

    if (pstr && polystr) {
        const struct test_rootfinder_example ex = { pstr, polystr, -1 };
        test_one_example(&ex);
    } else if (d >= 1) {
        for( ; iter-- ; )
            test_random_p_and_poly(d, max_bits, state);
    } else {
        for(size_t i = 0 ; examples[i].p ; i++)
            test_one_example(&examples[i]);

        while (iter--) {
            d = 1 + (int) gmp_urandomm_ui(state, 15);
            test_random_p_and_poly(d, max_bits, state);
        }

    }

    tests_common_clear ();
    return EXIT_SUCCESS;
}

#if 0
// magma code for producing test cases.
s:=1.2;            
p:=10;
while p lt 2^200 do
    for i in [1..100] do
        p:=NextPrime(p);
        d:=Random([2..7]);
        coeffs:=[Random(GF(p)):i in [0..d]];
        F:=PolynomialRing(GF(p))!coeffs;
        printf "in %o", p;
        for c in Reverse(coeffs) do printf " %o", c; end for;
        printf "\n";
        r:=Sort([x[1]: x in Roots(F) ]);
        printf "out";
        for c in r do printf " %o", c; end for;
        printf "\n";
    end for;
    p := Ceiling(p*s);
end while;


#endif

