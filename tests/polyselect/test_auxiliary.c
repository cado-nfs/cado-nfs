#include "cado.h"
#include "utils.h"
#include "polyselect/auxiliary.h"
#include "polyselect/size_optimization.h"
#include "tests_common.h"

static int
check_num (double x, double y, double emax)
{
  double e = fabs (x - y) / fabs (y);

  if (e > emax)
    {
      printf ("expected %.16e, got %.16e (rel. error %e)\n", y, x, e);
      return 0;
    }
  return 1;
}

static void
test_L2_lognorm (void)
{
  mpz_poly p;
  double n;

  mpz_poly_init (p, MAX_DEGREE);

  /* degree 1 */
  mpz_set_ui (p->coeff[0], 1);
  mpz_set_ui (p->coeff[1], 2);
  p->deg = 1;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.68393671871075337595, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 0.94925084419690070780, 2e-10);

  /* degree 2 */
  mpz_set_ui (p->coeff[2], 3);
  p->deg = 2;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.82777775480769542870, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.3718482492081025724, 2e-10);

  /* degree 3 */
  mpz_set_si (p->coeff[3], -4);
  p->deg = 3;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.73159180848396739500, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 1.7170713330509004847, 2e-10);

  /* degree 4 */
  mpz_set_si (p->coeff[4], 5);
  p->deg = 4;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.88625243226159843905, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.1476529791715452362, 2e-10);

  /* degree 5 */
  mpz_set_si (p->coeff[5], -6);
  p->deg = 5;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.90490889517894624795, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.5284494790416146421, 2e-10);

  /* degree 6 */
  mpz_set_si (p->coeff[6], 7);
  p->deg = 6;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.94130996783641860020, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 2.9064780318501317870, 2e-10);

  /* degree 7 */
  mpz_set_si (p->coeff[7], -8);
  p->deg = 7;
  n = L2_lognorm (p, 1.0);
  check_num (n, 0.95405307631727864405, 2e-10);
  n = L2_lognorm (p, 2.0);
  check_num (n, 3.2796366844822268800, 2e-10);

  mpz_poly_clear (p);
}

/* t=0: generate polynomial with skewness < 1
   t=1: generate polynomial with skewness > 1
   t=2: generate random polynomial */
static void
test_L2_skewness (int t)
{
  mpz_poly p;
  int d, i;
  double s, n, sl, nl, sh, nh, eps;
  int prec = 10;

  mpz_poly_init (p, MAX_DEGREE);
  if (t == 0)
    mpz_set_ui (p->coeff[0], 1);
  else if (t == 1)
    mpz_set_ui (p->coeff[0], 4294967295UL);
  else
    mpz_set_ui (p->coeff[0], mrand48 ());
  for (d = 1; d <= 7; d++)
    {
      if ((t == 0 || t == 1) && d > 1)
        mpz_set_ui (p->coeff[d-1], 1);
      if (t == 0)
        mpz_set_ui (p->coeff[d], 4294967295UL);
      else if (t == 1)
        mpz_set_ui (p->coeff[d], 1);
      else
        {
          do
            mpz_set_ui (p->coeff[d], mrand48 ());
          while (mpz_cmp_ui (p->coeff[d], 0) == 0);
        }
      p->deg = d;
      s = L2_skewness (p, prec);
      n = L2_lognorm (p, s);
      eps = ldexp (fabs (n), -prec);
      /* check that skewness to the left and to the right is worse */
      for (i = 0; i < 53; i++)
        {
          sl = s - ldexp (s, i - 53);
          nl = L2_lognorm (p, sl);
          if (nl < n - eps)
            {
              printf ("Non-optimal skewness for polynomial\n");
              mpz_poly_fprintf (stdout, p);
              printf ("For skewness %.16e, norm is %.16e\n", s, n);
              printf ("For skewness %.16e, norm is %.16e\n", sl, nl);
              abort ();
            }
          sh = s + ldexp (s, i - 53);
          nh = L2_lognorm (p, sh);
          if (nh < n - eps)
            {
              printf ("Non-optimal skewness for polynomial\n");
              mpz_poly_fprintf (stdout, p);
              printf ("For skewness %.16e, norm is %.16e\n", s, n);
              printf ("For skewness %.16e, norm is %.16e\n", sh, nh);
              abort ();
            }
        }
    }

  /* non-regression test for -1+x-x^2+x^4+x^5+x^6 */
  mpz_set_si (p->coeff[0], -1);
  mpz_set_si (p->coeff[1], 1);
  mpz_set_si (p->coeff[2], -1);
  mpz_set_si (p->coeff[3], 0);
  mpz_set_si (p->coeff[4], 1);
  mpz_set_si (p->coeff[5], 1);
  mpz_set_si (p->coeff[6], 1);
  p->deg = 6;
  s = L2_skewness (p, prec);
  ASSERT_ALWAYS (s == 1.0);

  /* test of degree 2 */
  mpz_set_si (p->coeff[0], 42);
  mpz_set_si (p->coeff[1], 0);
  mpz_set_si (p->coeff[2], 17);
  p->deg = 2;
  s = L2_skewness (p, prec);
  ASSERT_ALWAYS (1.571 <= s && s <= 1.572);

  /* test of degree 1 */
  mpz_set_si (p->coeff[0], 42);
  mpz_set_si (p->coeff[1], 17);
  p->deg = 1;
  s = L2_skewness (p, prec);
  ASSERT_ALWAYS (2.470 <= s && s <= 2.471);

  mpz_poly_clear (p);
}

// TODO: add test based on experiments done with rsa896 polynomials
static void
test_size_optimization (void)
{
  mpz_poly f, g, f_opt, g_opt;
  double n;

  /* check size-optimization of some RSA-1024 polynomial */
  mpz_poly_init (f, 6);
  mpz_poly_init (f_opt, 6);
  mpz_poly_init (g, 1);
  mpz_poly_init (g_opt, 1);
  mpz_set_str (g->coeff[1], "479811439908216194249453", 10);
  mpz_set_str (g->coeff[0],
                    "-1817512374557883972669278009547068402044185664937", 10);
  g->deg = 1;
  mpz_set_str (f->coeff[6], "3746994889972677420", 10);
  mpz_set_str (f->coeff[5], "11057269082141058123", 10);
  mpz_set_str (f->coeff[4], "-476255458524556917851536204692614007", 10);
  mpz_set_str (f->coeff[3], "571753863783554395130237455946021515899994741285",
               10);
  mpz_set_str (f->coeff[2],
               "-119215284531776978224224998903910352781670761348", 10);
  mpz_set_str (f->coeff[1],
               "776740771444910094386701822648133452751778906680", 10);
  mpz_set_str (f->coeff[0],
               "-690339709954574842053460672118938320489859967550", 10);
  f->deg = 6;
  n = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
  ASSERT_ALWAYS(106.895 <= n && n <= 106.905);

  size_optimization (f_opt, g_opt, f, g, SOPT_DEFAULT_EFFORT, 0);
  n = L2_skew_lognorm (f_opt, SKEWNESS_DEFAULT_PREC);
  ASSERT_ALWAYS(n <= 87.415);

  size_optimization (f_opt, g_opt, f, g, 3, 0);
  n = L2_skew_lognorm (f_opt, SKEWNESS_DEFAULT_PREC);
  /* note: now size_optimization optimizes the sum of the lognorm and of the
     expected alpha value, thus on this particular example we do no longer get
     a better lognorm with sopt_effort=3 */
  ASSERT_ALWAYS(n <= 87.415);

  mpz_poly_clear (f);
  mpz_poly_clear (g);
  mpz_poly_clear (f_opt);
  mpz_poly_clear (g_opt);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  test_L2_lognorm ();
  test_L2_skewness (0);
  test_L2_skewness (1);
  test_L2_skewness (2);
  test_size_optimization ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
