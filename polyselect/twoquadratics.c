#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "auxiliary.h"
#include "area.h"
#include "cado_poly.h"
#include "murphyE.h"
#include "mpz_poly.h"
#include "mpz_vector.h"
#include "macros.h"
#include "gmp_aux.h"
#include "params.h"
#include "verbose.h"
#include "polyselect_norms.h"

typedef struct {
  cado_poly cpoly;
  mpz_t m;
  mpz_t p;  // common root is m/p mod N
  mpz_t skew;
  double E;
} cado_poly_extended_s;

typedef cado_poly_extended_s cado_poly_extended[1];
typedef cado_poly_extended_s *cado_poly_extended_ptr;

void
cado_poly_extended_init (cado_poly_extended poly)
{
  cado_poly_init (poly->cpoly);
  cado_poly_provision_new_poly (poly->cpoly);
  cado_poly_provision_new_poly (poly->cpoly);
  mpz_init (poly->m);
  mpz_init (poly->p);
  mpz_init (poly->skew);
  poly->E = 0.0;
}

void
cado_poly_extended_clear (cado_poly_extended poly)
{
  cado_poly_clear (poly->cpoly);
  mpz_clear (poly->m);
  mpz_clear (poly->p);
  mpz_clear (poly->skew);
  poly->E = 0.0;
}

void
cado_poly_set2 (cado_poly_ptr cpoly, mpz_poly f, mpz_poly g, mpz_t N, 
                mpz_t skew)
{
  for( ; cpoly->nb_polys < 2 ; )
      cado_poly_provision_new_poly(cpoly);
  mpz_poly_set (cpoly->pols[0], f);
  mpz_poly_set (cpoly->pols[1], g);
  mpz_set (cpoly->n, N);
  cpoly->skew = mpz_get_d (skew);
}

void
cado_poly_extended_set (cado_poly_extended poly, mpz_poly f, mpz_poly g,
                        mpz_t N, mpz_t p, mpz_t skew, double E, mpz_t m)
{
  mpz_poly_set (poly->cpoly->pols[0], f);
  mpz_poly_set (poly->cpoly->pols[1], g);
  mpz_set (poly->cpoly->n, N);
  poly->cpoly->skew = mpz_get_d (skew);

  mpz_set (poly->p, p);
  mpz_set (poly->skew, skew);
  poly->E = E;
  mpz_set(poly->m, m);
}

double
cado_poly_extended_get_E (cado_poly_extended poly)
{
  return poly->E;
}

void cado_poly_extended_print_cado_format(FILE *out, cado_poly_extended poly,
    mpz_srcptr N) {
  mpz_poly_ptr f0 = poly->cpoly->pols[0];
  mpz_poly_ptr f1 = poly->cpoly->pols[1];

  gmp_fprintf(out, "n: %Zd\n", N);
  for (int i = 0; i < 3; i++)
    gmp_fprintf(out, "c%d: %Zd\n", i, f0->_coeff[i]);
  for (int i = 0; i < 3; i++)
    gmp_fprintf(out, "Y%d: %Zd\n", i, f1->_coeff[i]);
  fprintf(out, "skew: %1.2f\n", L2_combined_skewness2(f0, f1));
}

void
cado_poly_extended_print (FILE *out, cado_poly_extended poly, char *pre)
{
  mpz_t tmp;
  mpz_init (tmp);
  double s = poly->cpoly->skew;

  mpz_poly_ptr f0 = poly->cpoly->pols[0];
  mpz_poly_ptr f1 = poly->cpoly->pols[1];

  gmp_fprintf (out, "%sN = %Zd\n", pre, poly->cpoly->n);

  fprintf(out, "%sf0 = ", pre);
  mpz_poly_fprintf (out, f0);
  fprintf(out, "%sf1 = ", pre);
  mpz_poly_fprintf (out, f1);

  if (f0->deg != -1 && f1->deg != -1)
  {
    fprintf (out, "%sE = %e\n", pre, poly->E);
    fprintf (out, "%salpha_f0 = %.2f\n", pre, get_alpha (f0, get_alpha_bound ()));
    fprintf (out, "%salpha_f1 = %.2f\n", pre, get_alpha (f1, get_alpha_bound ()));
    fprintf (out, "%sskew_f0 = %.2f\n", pre, L2_skewness (f0));
    fprintf (out, "%sskew_f1 = %.2f\n", pre, L2_skewness (f1));
    gmp_fprintf (out, "%sskewness = %Zd\n", pre, poly->skew);
    fprintf (out, "%sL2_skew_norm_f0 = %.2f\n", pre, L2_lognorm (f0, s));
    fprintf (out, "%sL2_skew_norm_f1 = %.2f\n", pre, L2_lognorm (f1, s));
    gmp_fprintf (out, "%sp = %Zd\n", pre, poly->p);
    gmp_fprintf (out, "%sm = %Zd\n", pre, poly->m);
    mpz_mod (tmp, poly->m, poly->p);
    gmp_fprintf (out, "%sr = %Zd\n", pre, tmp);
  }

  mpz_clear (tmp);
}



/* Set m to the nearest integer to m0, such that
   m = r mod P, i.e. find the integer in [m-p/2,m+p/2] which is congruent to
   r modulo P.
 */
void
compute_m (mpz_t m, mpz_t m0, mpz_t r, mpz_t P)
{
  mpz_t t;
  mpz_init (t);
  mpz_tdiv_q_2exp (m, P, 1); // m = floor(p/2)
  mpz_add (m, m ,m0);        // m = m0 + floor(p/2)
  mpz_sub (t, m ,r);         // t = m0 + floor(p/2) - r
  mpz_mod (t, t, P);         // t = (m0 + floor(p/2) -r) mod P
  mpz_sub (m, m, t);         // m = m0 + floor(p/2) - t
                             // By definition, (m0 + floor(p/2) - t) = r mod P
  mpz_clear (t);
}

/* Compute maximun skewness, which in floor(N^(1/d^2)) */
void
compute_default_max_skew (mpz_t skew, mpz_t N, int d)
{
  mpz_root (skew, N, (unsigned long) d*d);
}

/* Returns the square root of N modulo P, using Tonelli-Shanks' algorithm.
   Assume P is prime.
   Returns 0 in case of error (e.g. N is not a square), else non-zero.
*/
int
mpz_msqrt (mpz_t r, mpz_t N, mpz_t P)
{
  if (mpz_cmp_ui (P, 2) == 0)
  {
    mpz_mod (r, N, P);
    return mpz_sgn(r); // Return 0 if r == 0
  }
  else if (mpz_legendre (N, P) != 1)
    return 0;
  else
  {
    mpz_t Q, n, z, c, R, t, t0, b;

    mpz_init(n);
    mpz_init_set_ui(z, 2);
    mpz_init(c);
    mpz_init(R);
    mpz_init(t);
    mpz_init(t0);
    mpz_init(b);
    mpz_init(Q);

    mpz_mod(n, N, P);

    mpz_sub_ui (Q, P, 1); // Q = P-1
    int M, S;
    for (S = 0; mpz_even_p(Q); S++)
      mpz_tdiv_q_2exp (Q, Q, 1);
    M = S;

    while (mpz_legendre(z, P) != -1)
        mpz_add_ui(z, z, 1);

    // We want c = z^Q mod P , t = n^Q mod P and R = n^((Q+1)/2) mod P
    mpz_powm (c, z, Q, P); // c = z^Q mod P
    mpz_sub_ui (Q, Q, 1);
    mpz_tdiv_q_2exp (Q, Q, 1);
    mpz_powm (R, n, Q, P); // R = n ^ ((Q-1)/2) mod P
    mpz_mul (t, R, R);
    mpz_mod (t, t, P);       // t = n^(Q-1) mod P
    mpz_mul (R, R, n);
    mpz_mod (R, R, P);       // R = n ^ ((Q+1)/2) mod P
    mpz_mul (t, t, n);
    mpz_mod (t, t, P);       // t = n^Q mod P

    while (mpz_cmp_ui (t, 1) != 0)
    {
      int i = 0;
      mpz_set(t0, t);
      while(mpz_cmp_ui(t, 1) != 0)
      {
        mpz_pow_ui(t, t, 2);
        mpz_mod(t, t, P);
        i++;
      }
      mpz_set(b, c);
      int j;
      for (j = 0; j < M - i - 1; j++)
      {
        mpz_pow_ui(b, b, 2);
        mpz_mod(b, b, P);
      }
      mpz_mul(R, R, b);
      mpz_mod(R, R, P);
      mpz_mul(c, b, b);
      mpz_mod(c, c, P);
      mpz_mul(t, t0, c);
      mpz_mod(t, t, P);
      M = i;
    }

  mpz_set (r, R);

  mpz_clear(Q);
  mpz_clear(n);
  mpz_clear(z);
  mpz_clear(c);
  mpz_clear(t);
  mpz_clear(t0);
  mpz_clear(b);
  mpz_clear(R);
  return 1;
  }
}

/* Compute f and g the two polynomials found via Montgomery's Two
   Quadratics Method, given N, P and m.
   P and N should be coprime and m^2-N should be divisible by P.
   The geometric progression is [c0, c1, c2] = [ p, m, (m^2-N)/p ]
   We start with a = [m, -p, 0] and b = [(m*t-c2)/P, -t, 1]
      where t = c2/m mod P
   We then compute a reduced basis of the lattice spanned by a and b with the
   Lagrange algorithm using the maximun possible skewness (bound by maxS)
 */
void MontgomeryTwoQuadratics (mpz_poly f, mpz_poly g, mpz_t skew, mpz_t N,
                              mpz_t P, mpz_t m, mpz_t maxS)
{
  ASSERT_ALWAYS (mpz_coprime_p(P, N));

  mpz_vector_t a, b, reduced_a, reduced_b;
  mpz_vector_init (a, 3);
  mpz_vector_init (b, 3);
  mpz_vector_init (reduced_a, 3);
  mpz_vector_init (reduced_b, 3);

  mpz_t tmp, c2, t;
  mpz_init(tmp);
  mpz_init(c2);
  mpz_init(t);

  // compute c2
  mpz_mul (tmp, m, m);
  mpz_sub (tmp, tmp, N);
  ASSERT_ALWAYS (mpz_divisible_p(tmp, P));
  mpz_divexact (c2, tmp, P);    // c2 = (m^2 - N) / P

  //compute t
  int ret = mpz_invert (tmp, m, P);
  ASSERT_ALWAYS (ret != 0);
  mpz_mul (t, tmp, c2);          // t = c2 / m mod P

  // Set vector a
  mpz_vector_setcoordinate (a, 0, m); // a[0] = m
  mpz_neg (tmp, P);
  mpz_vector_setcoordinate (a, 1, tmp); // a[1] = -P
  mpz_vector_setcoordinate_ui (a, 2, 0); // a[2] = 0

  // Set vector b
  mpz_mul (tmp, m, t);
  mpz_sub (tmp, tmp, c2);
  ASSERT_ALWAYS (mpz_divisible_p(tmp, P));
  mpz_divexact (tmp, tmp, P);
  mpz_vector_setcoordinate (b, 0, tmp); // b[0] = (m*t - c2) / P
  mpz_neg (tmp, t);
  mpz_vector_setcoordinate (b, 1, tmp); // b[1] = -t
  mpz_vector_setcoordinate_ui (b, 2, 1); // b[2] = 1

  mpz_vector_reduce_with_max_skew (reduced_a, reduced_b, skew, a, b, maxS, 2);

  mpz_vector_get_mpz_poly(f, reduced_a);
  mpz_vector_get_mpz_poly(g, reduced_b);

  mpz_clear(tmp);
  mpz_clear(c2);
  mpz_clear(t);
  mpz_vector_clear (a);
  mpz_vector_clear (b);
  mpz_vector_clear (reduced_a);
  mpz_vector_clear (reduced_b);
}


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "N", "input number (default c59)");
  param_list_decl_usage(pl, "minP", "Use P > minP (default 2)");
  param_list_decl_usage(pl, "maxP", "Use P <= maxP (default nextprime(minP))");
  param_list_decl_usage(pl, "skewness", "maximun skewness possible "
                                        "(default floor(N^(1/4))");
  param_list_decl_usage(pl, "v", "verbose output (print all polynomials)");
  param_list_decl_usage(pl, "q", "quiet output (print only best polynomials)");
  char str[200];
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
  param_list_decl_usage(pl, "print", "print in cado_poly format");
  verbose_decl_usage(pl);
}

static void usage (const char *argv, param_list pl)
{
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

int main(int argc, char const * argv[])
{
  const char * argv0 = argv[0];
  int verbose = 0, quiet = 0, print = 0;

  mpz_t N, minP, maxP, max_skewness;
  mpz_init (N);
  mpz_init (minP);
  mpz_init (maxP);
  mpz_init (max_skewness);

  param_list pl;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-q", &quiet);
  param_list_configure_switch (pl, "-print", &print);

  if (argc == 1)
    usage (argv[0], pl);

  argc--,argv++;
  for ( ; argc ; )
  {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, pl);
  }

  if (!param_list_parse_mpz(pl, "skewness", max_skewness))
    mpz_set_ui (max_skewness, 0);
  else if (mpz_cmp_ui (max_skewness, 1) < 0)
  {
    gmp_fprintf(stderr, "Error, skewness (%Zd) should be positive\n",
                         max_skewness);
    abort();
  }

  if (param_list_parse_double (pl, "area", &area) == 0) /* no -area */
    area = AREA;
  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;
  if (!param_list_parse_mpz(pl, "minP", minP))
    mpz_set_ui (minP, 2);
  if (!param_list_parse_mpz(pl, "maxP", maxP))
    mpz_nextprime(maxP, minP);
  if (!param_list_parse_mpz(pl, "N", N))
    mpz_set_str (N,
        "71641520761751435455133616475667090434063332228247871795429", 10);

  if (mpz_cmp_ui (minP, 0) < 0)
  {
    gmp_fprintf(stderr, "Error, minP (%Zd) should be greater or equal to 0\n",
                        minP);
    abort();
  }
  if (mpz_cmp (maxP, minP) < 0)
  {
    gmp_fprintf(stderr, "Error, maxP (%Zd) should be greater or equal to minP "
                        "(%Zd)\n", maxP, minP);
    abort();
  }
  if (mpz_cmp_ui (N, 3) < 0)
  {
    gmp_fprintf(stderr, "Error, N (%Zd) should be greater or equal to 3\n", N);
    abort();
  }
  if (mpz_divisible_ui_p (N, 2))
  {
    gmp_fprintf(stderr, "Error, N (%Zd) should not be divisible by 2\n", N);
    abort();
  }

  if (quiet)
    verbose = -1;

  verbose_interpret_parameters(pl);
  if (param_list_warn_unused(pl))
    usage (argv0, pl);
  param_list_print_command_line (stdout, pl);

  if (verbose >= 0)
  {
    gmp_printf("### N = %Zd\n", N);
    printf("### Bf = %.1e , Bg = %.1e , area = %.1e\n", bound_f, bound_g, area);
  }

  /* Compute max_skewness: use max_skewness if max_skewness argument is greater
     than 0 and lesser than default value */
  mpz_t tmp;
  mpz_init (tmp);
  compute_default_max_skew (tmp, N, 2);
  if (mpz_cmp_ui(max_skewness, 0) == 0 || mpz_cmp(max_skewness, tmp) > 0)
    mpz_set (max_skewness, tmp);
  mpz_clear (tmp);

/* Given an integer N and bounds on P, tests all pairs of polynomials
   with P such that minP < P <= maxP and print the pair having the best
   rating, along with its skewness, norm, alpha(f), alpha(g), and the
   prime number P which gave it
 */
  cado_poly_extended best_poly, poly;
  cado_poly cur_poly;
  mpz_t P, sqrtN, r, m, skew_used;
  mpz_poly f, g;

  cado_poly_extended_init (best_poly);
  cado_poly_extended_init (poly);
  cado_poly_init (cur_poly);
  mpz_init(P);
  mpz_init(sqrtN);
  mpz_init(r);
  mpz_init(m);
  mpz_init(skew_used);
  mpz_poly_init (f, 2);
  mpz_poly_init (g, 2);

  // sqrtN = floor(sqrt(N))
  mpz_sqrt(sqrtN, N);

  mpz_nextprime (P, minP);
  while(mpz_cmp (P, maxP) <= 0)
  {
    if (mpz_kronecker(N, P) == 1)
    {
      for (unsigned int k = 0; k < 2; k++)
      {
        if (k == 0) // For the first time compute r = sqrt of N modulo P
        {
          int ret = mpz_msqrt (r, N, P);
          ASSERT_ALWAYS (ret != 0);
        }
        else // The second root mod P is P-r
          mpz_sub(r, P, r);

        if (verbose >= 0)
          gmp_printf("### P = %Zd: compute polynomials for r = %Zd\n", P, r);
        // m is the first integer >= sqrtN such that m = r (mod P)
        compute_m (m, sqrtN, r, P);
        MontgomeryTwoQuadratics (f, g, skew_used, N, P, m, max_skewness);
        cado_poly_set2 (cur_poly, f, g, N, skew_used);
        double E = MurphyE (cur_poly, bound_f, bound_g, area, MURPHY_K,
                            get_alpha_bound ());
        if(E > cado_poly_extended_get_E(best_poly))
          cado_poly_extended_set (best_poly, f, g, N, P, skew_used, E, m);
        if (verbose >= 1)
        {
          cado_poly_extended_set (poly, f, g, N, P, skew_used, E, m);
          cado_poly_extended_print (stdout, poly, "# ");
        }
      }
    }
    else if (verbose >= 0)
      gmp_printf("### P = %Zd: N is not a square modulo P, skipping it.\n", P);
    mpz_nextprime(P, P);
  }

  if (print) {
    cado_poly_extended_print_cado_format(stdout, best_poly, N);
    cado_poly_extended_print (stdout, best_poly, "# ");
  } else {
    printf("### Best polynomials found:\n");
    cado_poly_extended_print (stdout, best_poly, "");
  }

  mpz_clear(P);
  mpz_clear(sqrtN);
  mpz_clear(r);
  mpz_clear(m);

  mpz_poly_clear (f);
  mpz_poly_clear (g);

  cado_poly_clear (cur_poly);
  cado_poly_extended_clear (poly);
  cado_poly_extended_clear (best_poly);
  mpz_clear (N);
  mpz_clear (skew_used);
  mpz_clear (max_skewness);
  mpz_clear (minP);
  mpz_clear (maxP);
  param_list_clear(pl);
  return EXIT_SUCCESS;
}
