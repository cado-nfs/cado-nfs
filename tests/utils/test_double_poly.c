#include "cado.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "macros.h"
#include "double_poly.h"
#include "tests_common.h"




/* This computes only roots in [0,s] */
void 
test_double_poly_compute_roots1(const char *poly_str, const char *roots_str, 
                               const double err_margin, const double s,
                               const int verbose)
{
  double_poly poly, roots1, roots2;
  int nr_roots, i;

  double_poly_init(poly, -1);
  double_poly_init(roots1, -1);

  double_poly_set_string(poly, poly_str);
  double_poly_set_string(roots1, roots_str);

  double_poly_init(roots2, poly->deg);

  if (verbose) {
      char * s;
      double_poly_asprint (&s, poly, "x");
      printf("Testing polynomial %s\n", s);
      free(s);
  }

  nr_roots = double_poly_compute_roots(roots2->coeff, poly, s);
  if (nr_roots != (int) roots1->deg + 1) {
    fprintf (stderr, "double_poly_compute_roots() produced wrong number of roots %d, reference has %d\n",
             nr_roots, roots1->deg + 1);
    double_poly_print (stderr, poly, "Polynomial was: ");
    fprintf (stderr, "bound is s=%.16e\n", s);
    abort();
  }
  for (i = 0; i < nr_roots; i++) {
    if (!cmp_double(roots1->coeff[i], roots2->coeff[i], err_margin)) {
      fprintf (stderr, "double_poly_compute_roots() produced wrong roots %f, reference has %f\n",
               roots2->coeff[i], roots1->coeff[i]);
    
      abort();
    }
  }

  double_poly_clear(poly);
  double_poly_clear(roots1);
  double_poly_clear(roots2);
}

static void
test_double_poly_compute_roots(const int verbose)
{
  test_double_poly_compute_roots1("1", "", 1e-9, 3., verbose);

  /* A few roots of 2 */
  test_double_poly_compute_roots1("-2 1", "2", 1e-9, 3., verbose);
  test_double_poly_compute_roots1("-2 0 1", "1.41421356237310", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 1", "1.25992104989487", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 0 1", "1.18920711500272", 1e-6, 3., verbose);
  test_double_poly_compute_roots1("-2 0 0 0 0 1", "1.14869835499704", 1e-6, 3., verbose);

  /* (x-1)*(x-2) */
  test_double_poly_compute_roots1("2 -3 1", "1 2", 1e-6, 3., verbose);
  /* (x-1)*(x-2)*(x-3) */
  test_double_poly_compute_roots1("-6 11 -6 1", "1 2 3", 1e-6, 4., verbose);
  /* (x-1)*(x-2)*(x-3)*(x-4) */
  test_double_poly_compute_roots1("24 -50 35 -10 1", "1 2 3 4", 1e-6, 5., verbose);
  /* (x-1)*(x-2)*(x-3)*(x-4)*(x-5) */
  test_double_poly_compute_roots1("-120 274 -225 85 -15 1", "1 2 3 4 5", 1e-6, 6., verbose);
  
  /* Let f(x+1/x) * x^6 == (x^13-1)/(x-1). Test the positive roots */
  test_double_poly_compute_roots1("-1 3 6 -4 -5 1 1", "0.241073360510646 1.13612949346231 1.77091205130642", 1e-6, 2., verbose);
  /* and the negative ones */
  test_double_poly_compute_roots1("-1 3 6 -4 -5 1 1", "-0.709209774085071 -1.49702149634220 -1.94188363485210", 1e-6, -2., verbose);
  
  /* this is f(x-2). 6 real roots, 0 rational */
  test_double_poly_compute_roots1("1 -21 70 -84 45 -11 1", "0.0581163651478959 0.502978503657798 1.29079022591493 2.24107336051065 3.13612949346231 3.77091205130642", 1e-6, 4., verbose);

  /* dichotomy with too few iterations fails */
  test_double_poly_compute_roots1 ("37718991555021708231785218373859670336563892134301066596876783017325698882362933248 -58290863912589135997939905826055669574526326226812326870330700385351723122688 4974019329969663881845375223408004510305936664806589566321472066027520 0 -2765086156017059372041966183496747450660061680253272064 17192931019341412634837118288585351300728750080 3396573496846254368813196771328", "687391.0 1.1967277e7 4.2033871e7 1.48782523e8", 10.0, 9007199254740992.0, verbose);

  /* false position needs weighting */
  test_double_poly_compute_roots1("1.1771253911282693e+36 1.0293658912886811e+39 1.9888797504712385e+38 -9.7240762665983556e+38 3.5375219923207381e+37 -1.8078625258002374e+40 -1.6114410711830154e+39", "0.469426633865124252781568271387", 0.001, 1, verbose);

  /* false position reaches s==a || s==b very early, need countermeasure */
  test_double_poly_compute_roots1("-4067876448477449548563432969586031412941700056587601428084520911569151852544 299197339446978982940341333993590148850152155031542808049001636857517178880 -88455602043598760130397057784414636527233402886601450561774258429246832640 0 552084484207525010843995104449133201610495729247299431552778240 386336020186016733435543279553138102808608768 38629366113825858026209280", "400274.630694569284177284060753", 0.1, 2199023255552, verbose);

  /* false position needs many iterations */
  test_double_poly_compute_roots1("416305583514625790805142742552103399071483895746667754265057484511150866432 202399505763732049099628933992141454064175135567823362216516742332052144128 0 -2677221347026437285957968988912544408687885411868999680 21555240319368651153052935288520704 5123362746908340224", "8694859813.27109836187832047186 720776737597677797.826174347519", 10, 1e20, verbose);

  /* false position produces b=NaN */
  test_double_poly_compute_roots1 ("-5.1229871591623088e+251 4.8231399628079727e+240 -7.5683722678735590e+228 -1.8935837380523070e+224 -3.4853123818766583e+152", "", 1e60, 1e72, verbose);
}

void
test_double_poly_set (void)
{
  double_poly s, r;

  double_poly_init (s, 2);
  double_poly_init (r, 3);
  s->coeff[0] = -1.0;
  s->coeff[1] = 17.0;
  s->coeff[2] = 42.0;
  double_poly_cleandeg(s, 2);
  double_poly_set (r, s);
  ASSERT_ALWAYS (r->deg == 2);
  ASSERT_ALWAYS (r->coeff[0] == -1.0);
  ASSERT_ALWAYS (r->coeff[1] == 17.0);
  ASSERT_ALWAYS (r->coeff[2] == 42.0);
  double_poly_clear (s);
  double_poly_clear (r);
}

#if GNUC_VERSION_ATLEAST(4,4,0)
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
void
test_double_poly_eval (void)
{
  double_poly s;
  unsigned int deg;
  double v, w;

  double_poly_init (s, 52);
  for (deg = 0, w = 0.0; deg <= 52; deg++)
    {
      s->coeff[deg] = 1.0;
      s->deg = deg;
      v = double_poly_eval (s, 2.0);
      w = 2.0 * w + s->coeff[deg];
      ASSERT_ALWAYS (v == w);
    }
  double_poly_clear (s);
}
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic pop
#endif

void
test_double_poly_derivative (void)
{
  double_poly f, df;

  double_poly_init (f, 2);
  double_poly_init (df, 1);

  f->coeff[0] = 17.0;
  f->coeff[1] = 42.0;
  double_poly_cleandeg(f, 1);
  double_poly_derivative (df, f);
  ASSERT_ALWAYS (df->deg == 0 && df->coeff[0] == 42.0);

  /* Test in-place operation */
  f->coeff[0] = 17.0;
  f->coeff[1] = 42.0;
  f->coeff[2] = 1728.0;
  double_poly_cleandeg(f, 2);
  double_poly_derivative (f, f);
  ASSERT_ALWAYS (f->deg == 1 && f->coeff[0] == 42.0 && f->coeff[1] == 3456);

  f->coeff[0] = 17.0;
  double_poly_cleandeg(f, 0);
  double_poly_derivative (df, f);
  ASSERT_ALWAYS (df->deg == -1);

  double_poly_clear (f);
  double_poly_clear (df);
}

void
test_double_poly_revert (void)
{
  double_poly f;

  double_poly_init (f, 2);

  /* try with degree 2 */
  f->coeff[0] = 1.0;
  f->coeff[1] = 2.0;
  f->coeff[2] = 3.0;
  double_poly_cleandeg(f, 2);
  double_poly_revert (f, f);
  ASSERT_ALWAYS (f->coeff[0] == 3.0 && f->coeff[1] == 2.0 && f->coeff[2] == 1.0);

  /* now with degree 1 */
  double_poly_cleandeg(f, 1);
  double_poly_revert (f, f);
  ASSERT_ALWAYS (f->coeff[0] == 2.0 && f->coeff[1] == 3.0);

  double_poly_clear (f);
}

void
test_double_poly_print ()
{
    char *t;
    int rc;
    double_poly poly;
    double_poly_init(poly, -1);

    double_poly_set_string(poly, "17");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "17");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_set_string(poly, "17 42");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "17+42*x");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_set_string(poly, "17 42 53");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "17+42*x+53*x^2");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    /* We add some examples where non-integer floating point values are
     * printed. Currently this uses the default printing mode of the c++
     * printers, there's no way to change it.
     */
    double_poly_set_string(poly, "17.001e8 0 53");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "1.7001e+09+53*x^2");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_set_string(poly, "17 0 -53.2 99");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "17-53.2*x^2+99*x^3");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_set_string(poly, "1 0 -1 99");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "1-x^2+99*x^3");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_set_string(poly, "0 -1 1");
    double_poly_asprint(&t, poly, "x");
    rc = strcmp(t, "-x+x^2");
    ASSERT_ALWAYS(rc == 0);
    free(t);

    double_poly_clear(poly);
}

void
test_double_poly_set_mpz_poly (void)
{
  double_poly p;
  mpz_poly q;

  mpz_poly_init (q, 2);
  double_poly_init (p, 2);
  mpz_set_ui (q->coeff[2], 17);
  mpz_set_si (q->coeff[1], -42);
  mpz_set_si (q->coeff[0], -3);
  q->deg = 2;
  double_poly_set_mpz_poly (p, q);
  ASSERT_ALWAYS (p->deg == 2 && p->coeff[2] == 17.0 && p->coeff[1] == -42.0 &&
          p->coeff[0] == -3.0);
  double_poly_clear (p);
  mpz_poly_clear (q);
}

unsigned int check_absolute_error(double value, double expected,
        unsigned int error MAYBE_UNUSED)
{
    if (value == 0.0)
        value = 1.0;
    if (expected == 0.0)
        expected = 1.0;
    if ((int)(fabs(log2(fabs(value)) - log2(fabs(expected)))) <= (int) error)
        return 1;
    return 0;
}

MAYBE_UNUSED void test_double_poly_resultant() {
  double_poly f, g;

  double_poly_init(f, -1);
  double_poly_init(g, -1);
  double res = 0.0;
  /*f=x^6+13*x^5+13*x^4+9*x^3+7*x+6*/
  /*g=128*x^2+128*x+128*/
  double_poly_set_string(f, "6 7 0 9 13 13 1");
  double_poly_set_string(g, "128 128 128");
  res = double_poly_resultant(f, g);
  double val = 162727720910848.000000;
  ASSERT_ALWAYS(check_absolute_error(val, res, 0));

  /*f=-3-15*x^1-9*x^2+3*x^3-12*x^4-12*x^5-3*x^6-3*x^7-12*x^8-15*x^9+6*x^10*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7*/
  double_poly_set_string(f, "-3 -15 -9 3 -12 -12 -3 -3 -12 -15 6");
  double_poly_set_string(g, "-6 -13 9 7 -5 -5 11 2");
  res = double_poly_resultant(f, g);
  val = -61519394185549840384.000000;
  ASSERT_ALWAYS(check_absolute_error(val, res, 0));

  /*f=7917871+7917871*x-7916275*x^2-7916275*x^3-7916275*x^4+7917871*x^5+15834944*x^6*/
  /*g=128*x^2+128*x+128*/
  double_poly_set_string(f, "7917871 7917871 -7916275 -7916275 -7916275 7917871 15834944");
  double_poly_set_string(g, "128 128 128");
  res = double_poly_resultant(f, g);
  val = 1102790158070603587092742144.000000;
  ASSERT_ALWAYS(check_absolute_error(val, res, 0));

  /*f=1365*x^6 + 1366*x^5+1368*x^4+1368*x^3+1368*x^2+1366*x+1366*/
  /*g=0*x^3+8320*x^2-50560*x-896*/
  double_poly_set_string(f, "1366 1366 1368 1368 1368 1366 1365");
  double_poly_set_string(g, "-896 -50560 8320 0");
  res = double_poly_resultant(f, g);
  val = 37263864605996575174727132124282880.000000;
  ASSERT_ALWAYS(check_absolute_error(val, res, 0));

  /*f=1365*x^6 + 1366*x^5+1368*x^4+1368*x^3+1368*x^2+1366*x+1366*/
  /*g=15*x^2-43368*x-4753*/

  double_poly_set_string(f, "1366 1366 1368 1368 1368 1366 1365");
  double_poly_set_string(g, "-4753 -43368 15");
  res = double_poly_resultant(f, g);
  val = 11186466747618860118741892404376785.000000;
  ASSERT_ALWAYS(check_absolute_error(val, res, 0));

  double_poly_clear(f);
  double_poly_clear(g);
}

int main()
{
  test_double_poly_compute_roots(0);
  test_double_poly_set ();
  test_double_poly_eval ();
  test_double_poly_derivative ();
  test_double_poly_revert ();
  test_double_poly_print ();
  test_double_poly_set_mpz_poly ();
  test_double_poly_resultant();
  exit(EXIT_SUCCESS);
}
