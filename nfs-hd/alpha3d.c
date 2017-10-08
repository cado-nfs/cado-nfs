#include "cado.h"
#include "utils.h"
#include <math.h>
#include <time.h>
#include "alpha3d.h"

static double expect_val_p(mpz_poly_srcptr f, uint64_t p, gmp_randstate_t state)
{
  mpz_t p_Z;
  mpz_init(p_Z);
  mpz_set_ui(p_Z, p);
  mpz_poly_factor_list lf;
  mpz_poly_factor_list_init(lf);

  mpz_poly_factor(lf, f, p_Z, state);
  unsigned int nr1 = 0;
  unsigned int nr2 = 0;

  for (int i = 0; i < lf->size; i++) {
    if (lf->factors[i]->f->deg == 1) {
      nr1++;
    } else if (lf->factors[i]->f->deg == 2) {
      nr2++;
    }
  }

  double val_1 = (double)nr1 * ((double) (p * p + p) / (double)(p * p * p - 1));
  double val_2 = 2.0 * (double)nr2 * (double)(p * p) / (double)((p * p - 1) * (p
        * p + p + 1));

  mpz_clear(p_Z);
  mpz_poly_factor_list_clear(lf);
  return val_1 + val_2;
}

/* Return the largest integer i such that p^i divides n.
   Warning: this function modifies n (powers of p are removed).
   FIXME: if/when a function mpz_remove_ui exists in GMP
   (analogous to mpz_remove for an unsigned long p) we should use it. */
static unsigned int
mpz_valuation (mpz_ptr n, unsigned long p)
{
  unsigned int i = 0;

  while (mpz_divisible_ui_p (n, p)) {
    mpz_divexact_ui (n, n, p);
    i++;
  }

  return i;
}

/* a is a degree-2 polynomial */
static unsigned int is_irreducible(mpz_poly_srcptr a)
{
  mpz_t delta;
  int res;

  ASSERT(a->deg == 2);
  mpz_init(delta);
  mpz_mul (delta, a->coeff[2], a->coeff[0]);
  mpz_mul_2exp (delta, delta, 2);
  mpz_neg (delta, delta);
  mpz_addmul(delta, a->coeff[1], a->coeff[1]);

  res = mpz_perfect_square_p(delta) == 0;

  mpz_clear(delta);
  return res;
}

//Generate a random number between [offset, offset + length[
static inline void sub_rand_mpz(mpz_ptr rand_Z, mpz_srcptr length,
    mpz_srcptr offset, gmp_randstate_t state)
{
  mpz_urandomm(rand_Z, state, length);
  mpz_add(rand_Z, rand_Z, offset);
}

static void rand_mpz(mpz_ptr rand_Z, mpz_srcptr min, mpz_srcptr max,
    gmp_randstate_t state)
{
  mpz_t length;
  mpz_init(length);
  mpz_sub(length, max, min);
  sub_rand_mpz(rand_Z, length, min, state);
  mpz_clear(length);
}

static void random_mpz_poly(mpz_poly_ptr g, int degree, mpz_srcptr min,
    mpz_srcptr max, gmp_randstate_t state)
{
  mpz_t rand_Z;
  mpz_init(rand_Z);

  for (int i = 0; i < degree; i++) {
    rand_mpz(rand_Z, min, max, state);
    mpz_poly_setcoeff(g, i, rand_Z);
  }
  rand_mpz(rand_Z, min, max, state);
  while(mpz_cmp_ui(rand_Z, 0) == 0) {
    rand_mpz(rand_Z, min, max, state);
  }
  mpz_poly_setcoeff(g, degree, rand_Z);

  mpz_clear(rand_Z);
}


static void mpz_poly_irred(mpz_poly_ptr a, mpz_srcptr one, mpz_srcptr B,
    gmp_randstate_t state)
{
  mpz_t content;
  mpz_init(content);

  do {
    random_mpz_poly(a, 2, one, B, state);
    mpz_poly_content(content, a);
  } while (mpz_cmp_ui(content, 1) != 0 || !is_irreducible(a));

#ifndef NDEBUG
  mpz_poly_content(content, a);
  ASSERT(mpz_cmp_ui(content, 1) == 0);
  ASSERT(is_irreducible(a));
#endif // NDEBUG

  mpz_clear(content);
}

static void monte_carlo_average_value(double * V, mpz_poly_srcptr f,
    unsigned long * bad_p, unsigned int length, unsigned int N,
    gmp_randstate_t state)
{
  for (unsigned int i = 0; i < length; i++) {
    V[i] = 0.0;
  }
  mpz_t B;
  mpz_init(B);
  mpz_set_ui(B, N);
  mpz_mul(B, B, B);
  mpz_add_ui(B, B, 1);

  mpz_t one;
  mpz_init(one);
  mpz_set_ui(one, 1);

  mpz_poly a;
  mpz_poly_init(a, 2);

  mpz_t resultant;
  mpz_init(resultant);

  mpz_t prod; /* product of all bad primes */
  mpz_t g;
  mpz_init_set_ui (prod, 1);
  mpz_init (g);
  for (unsigned int j = 0; j < length; j++)
    mpz_mul_ui (prod, prod, bad_p[j]);

  for (unsigned int i = 0; i < N; i++) {
    mpz_poly_irred(a, one, B, state);
    mpz_poly_resultant(resultant, f, a);
    mpz_gcd (g, resultant, prod);
    /* if g=1, no bad prime can divide the resultant */
    for (unsigned int j = 0; j < length && mpz_cmp_ui (g, 1) > 0; j++) {
      unsigned int v = mpz_valuation(resultant, bad_p[j]);
      V[j] += (double)v;
      if (v)
        mpz_divexact_ui (g, g, bad_p[j]);
      /* now g is the gcd of 'resultant' and the remaining bad primes */
    }
  }

  mpz_clear (g);
  mpz_clear (prod);
  mpz_clear(resultant);
  mpz_poly_clear(a);
  mpz_clear(one);
  mpz_clear(B);

  for (unsigned int i = 0; i < length; i++) {
    V[i] = V[i] / N;
  }
}

/* p_end is the bound on primes, N is the number of iterations in
   monte_carlo_average_value */
double alpha3d(mpz_poly_srcptr f, unsigned long p_end, gmp_randstate_t rstate, unsigned int N)
{
  mpz_t discriminant;
  mpz_init(discriminant);
  mpz_t lc;
  mpz_init(lc);

  mpz_poly_discriminant(discriminant, f);
  mpz_set(lc, mpz_poly_lc_const(f));

  prime_info pi;
  prime_info_init(pi);
  unsigned long p = 2;

  double alpha = 0.0;

  //bad_p is very long, but why not?
  unsigned long * bad_p = (unsigned long *) malloc(sizeof(unsigned long) *
      p_end);
  unsigned int index = 0;

  for ( ; p < p_end; p = getprime_mt(pi)) {
    if (mpz_divisible_ui_p(discriminant, p) ||
        mpz_divisible_ui_p(lc, p)) {
      bad_p[index] = p;
      index++;
    } else {
      alpha += log((double)p) * (1.0 / (double)(p - 1) - expect_val_p(f, p,
            rstate));
    }
  }

  double * V = (double *) malloc(sizeof(double) * index);
  monte_carlo_average_value(V, f, bad_p, index, N, rstate);

  for (unsigned int i = 0; i < index; i++) {
    alpha += log((double)bad_p[i]) * (1 / (double)(bad_p[i] - 1) - V[i]);
  }

  free (V);
  free(bad_p);
  prime_info_clear(pi);
  mpz_clear(lc);
  mpz_clear(discriminant);

  return alpha;
}

#ifdef MAIN_ALPHA3D
void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "f", "a polynomial");
  param_list_decl_usage(pl, "p", "bound of the primes");
  param_list_decl_usage(pl, "N", "bound for Monte Carlo computation");
  param_list_decl_usage(pl, "seed", "seed for Monte Carlo computation");
}

void initialise_parameters(int argc, char * argv[], mpz_poly_ptr f,
    unsigned long * p, gmp_randstate_t rstate, unsigned int * N)
{
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  FILE * fpl;
  char * argv0 = argv[0];

  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

    /* Could also be a file */
    if ((fpl = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, fpl, 0);
      fclose(fpl);
      argv++,argc--;
      continue;
    }

    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit (EXIT_FAILURE);
  }

  param_list_parse_mpz_poly(pl, "f", f);
  ASSERT(f->deg > 0);
  
  * p = 2000;
  param_list_parse_ulong(pl, "p", p);

  * N = 10000;
  param_list_parse_uint(pl, "N", N);

  unsigned long seed;
  if (param_list_parse_ulong(pl, "seed", &seed))
      gmp_randseed_ui(rstate, seed);

  param_list_clear(pl);
}
int main(int argc, char ** argv)
{
  mpz_poly f;
  unsigned long p;
  unsigned int N;
  gmp_randstate_t rstate;

  mpz_poly_init(f, -1);
  gmp_randinit_default(rstate);

  initialise_parameters(argc, argv, f, &p, rstate, &N);

  printf("%f\n", alpha3d(f, p, rstate, N));

  gmp_randclear(rstate);
  mpz_poly_clear(f);

  return 0;
}
#endif // MAIN_ALPHA3D
