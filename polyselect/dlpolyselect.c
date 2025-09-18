/* dlpolyselect - discrete logarithm polynomial selection via Joux-Lercier's
   algorithm

Example with Apache 512-bit key (Table 1 of
https://weakdh.org/imperfect-forward-secrecy-ccs15.pdf):

# we first search a degree (3,2) pair
$ ./dlpolyselect -df 3 -dg 2 -N 8372421755538377327377912526045445423027732035562313241965800453667849685158691589507936013805295187219621475007123900107532269487803598942841993804845107 -bound 8
...
c3: 4
c2: 3
c1: -6
c0: -8
Y2: 430506467754036192750431535685941771886700493958222
Y1: 687078639517265884307374391507473735402330150842399
Y0: -533201855944366475689054404823391134825953676216635
skew: 1.23
# f lognorm 1.28, alpha 0.62, score 1.90
# g lognorm 116.59, alpha -2.18, score 114.41
# f+g score 116.31

# then we compare with a degree (4,3) pair
$ ./dlpolyselect -df 4 -dg 3 -N 8372421755538377327377912526045445423027732035562313241965800453667849685158691589507936013805295187219621475007123900107532269487803598942841993804845107 -bound 8
...
c4: 6
c3: 1
c2: -8
c1: -3
c0: 6
Y3: -64071264306884991611859009153886700616
Y2: -115884379190374676852348454783130883186
Y1: 382823080720299801084267253734861739469
Y0: -14529492984288436819691699253895818591
skew: 1.05
# f lognorm 1.15, alpha -0.04, score 1.12
# g lognorm 87.27, alpha -4.07, score 83.20
# f+g score 84.32

To recover the polynomial pair used for the DLP768 record
(https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;a0c66b63.1606 and
http://eprint.iacr.org/2017/067.pdf):

# -modr r -modm m enables one to only consider polynomials f
# with index r mod m in Stage 1 (default is r=0 and m=1).
# This is useful to distribute the work among several machines.
$ ./dlpolyselect -N 1219344858334286932696341909195796109526657386154251328029273656175766870980306505584577389125860826715201547225794072935883258868036433287217994721542199148182841505800433148410869683590659346847659519108393837414567892730579162319 -df 4 -dg 3 -bound 140 -modr 66234274 -modm 100000000
...
c4: 55
c3: 5
c2: -86
c1: 34
c0: -140
Y3: 277260730400349522890422618473498148528706115003337935150
Y2: 217583293626947899787577441128333027617541095004734736415
Y1: -1937981312833038778565617469829395544065255938015920309679
Y0: -370863403886416141150505523919527677231932618184100095924
skew: 1.37
# f lognorm 3.93, alpha -1.96, score 1.97
# g lognorm 130.19, alpha -5.47, score 124.73
# f+g score 126.69

*/

#include "cado.h" // IWYU pragma: keep
#include <float.h> // for DBL_MAX
#include <stdlib.h>
#include <math.h> // pow
#include <time.h>
#include <limits.h>      // for ULONG_MAX
#include <stdio.h>       // for fprintf, printf, stderr, fflush, stdout
#include <string.h>      // for strcmp
#include <gmp.h>         // for mpz_t, mpz_clear, mpz_init, gmp_printf, mpz_...
#include "cado_poly.h"   // for cado_poly_fprintf_MurphyE, cado_poly
#include "macros.h"      // for ASSERT_ALWAYS, ASSERT
#include "omp_proxy.h" // IWYU pragma: keep
#include "auxiliary.h"
#include "gcd.h"        // gcd_uint64
#include "lll.h"        // mat_Z, LLL
#include "mpz_poly.h"
#include "murphyE.h"
#include "rootfinder.h"
#include "polyselect_norms.h"
#include "polyselect_alpha.h"
#include "timing.h"             // for seconds
#include "mpz_mat.h"

/* We assume a difference <= ALPHA_BOUND_GUARD between alpha computed
   with ALPHA_BOUND_SMALL and ALPHA_BOUND. In practice the largest value
   observed is 0.79. */
#define ALPHA_BOUND_GUARD 1.0

/* global variables */
double best_score_f = DBL_MAX, worst_score_f = DBL_MIN, sum_score_f = 0;
unsigned long f_candidate = 0;   /* number of irreducibility tests */
unsigned long f_irreducible = 0; /* number of irreducible polynomials */
double max_guard = DBL_MIN;
int skewed = 0;                 /* boolean (use skewed polynomial) */
unsigned long *count = NULL;    /* use only with -skewed */
int easySM = 0;                 /* see the -easySM option */
int rrf = -1;                   /* see the -rrf option */
int rrg = -1;                   /* see the -rrg option */

/* MPZ_POLY_TIMINGS is defined (or maybe not) in utils/mpz_poly.h
 */
#ifndef CADO_MPZ_POLY_H
#error "please include mpz_poly.h first"
#endif
#ifdef MPZ_POLY_TIMINGS
static double timer[4] = {0.0, };
#endif

double Bf = 0.0, Bg = 0.0, Area = 0.0;
double bestE = 0.0; /* best Murphy-E so far */
int opt_flag = 0; /* 0: optimize "simple" E
                     1: optimize Muprhy E */

#define TIMER_ROOTS 0
#define TIMER_IRRED 1
#define TIMER_LLL   2
#define TIMER_MURPHYE 3

#ifdef MPZ_POLY_TIMINGS
#define START_TIMER double t = seconds_thread ()
#define END_TIMER(x) add_timer (x, seconds_thread () - t)

/* flag=0: computing roots of f mod p
        1: checking irreducibility of f
        2: LLL reduction
        3: computing MurphyE */
static void
add_timer (int flag, double t)
{
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  timer[flag] += t;
}
#else
#define START_TIMER
#define END_TIMER(x)
#endif


static int
check_SM (mpz_poly_srcptr ff, mpz_srcptr ell)
{
    if (ff->deg <= 2) {
        return 1;
    }
    if (ff->deg > 4) {
        fprintf(stderr, "Not implemented\n");
        ASSERT_ALWAYS(0);
    }
    // in degree 3 and 4, the minimum number of SMs is 1. We
    // check that we have at least one root mod ell.
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    int nr = mpz_poly_roots_mpz(NULL, ff, ell, rstate);
    gmp_randclear(rstate);
    return (nr >= 1);
}


/*
  Print two nonlinear poly info. Return non-zero for record polynomials.
*/
static int
print_nonlinear_poly_info (mpz_poly_srcptr ff, double alpha_f, mpz_poly_srcptr gg,
                           int format,  mpz_srcptr n, mpz_srcptr ell)
{
    double skew, logmu[2], alpha_g_approx, alpha_g, score, score_approx;
    static double best_score = DBL_MAX;
    /* the coefficients of g are O(n^(1/df)) */

    /* we use the skewness minimizing the sum of lognorms */
    skew = L2_combined_skewness2 (ff, gg);
    logmu[1] = L2_lognorm (gg, skew);
    logmu[0] = L2_lognorm (ff, skew);
    /* first estimate alpha with a small bound */
    alpha_g_approx = get_alpha (gg, ALPHA_BOUND_SMALL);

    score_approx = logmu[1] + alpha_g_approx + logmu[0] + alpha_f;

    int better;

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    better = (score_approx >= best_score + ALPHA_BOUND_GUARD);

    if (better)
        return 0;

    /* now get a more precise alpha value */
    alpha_g = get_alpha (gg, get_alpha_bound ());

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
        score = logmu[1] + alpha_g + logmu[0] + alpha_f;
        if (score_approx - score > max_guard)
            max_guard = score_approx - score;
    }

    double E = 0.0;
    if (opt_flag == 0)
      {
          int better;

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          better = (score >= best_score);

          if (better)
              return 0; /* only print record scores */
      }
    else /* optimize Murphy-E */
      {
          int better;

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          better = (score >= best_score + 1.0); /* the guard 1.0 seems
                                                   good in practice */
          if (better)
                  return 0;

        /* compute Murphy-E */
        cado_poly p;
	START_TIMER;
        cado_poly_init(p);
        cado_poly_provision_new_poly(p);
        cado_poly_provision_new_poly(p);
        mpz_poly_set(p->pols[ALG_SIDE], ff);
        mpz_poly_set(p->pols[RAT_SIDE], gg);
        p->skew = skew;
        E = MurphyE (p, Bf, Bg, Area, MURPHY_K, get_alpha_bound ());
        cado_poly_clear(p);
	END_TIMER (TIMER_MURPHYE);

        {
            int better;
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
            better = E > bestE;

            if (!better)
                return 0;
        }
      }

    /* Possibly check the number of roots mod ell of f and g, assuming that
     * they have the minimum number of real roots. */
    if (easySM) {
        if (! check_SM(ff, ell))
            return 0;
        if (! check_SM(gg, ell))
            return 0;
    }

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
        bestE = E;

      if (score < best_score)
        best_score = score;

      if (format == 1)
	gmp_printf ("n: %Zd\n", n);
      else
	gmp_printf ("N %Zd\n", n);

      if (format == 1) {
        for (int i = ff->deg ; i >= 0 ; i--)
	  gmp_printf ("c%u: %Zd\n", i, mpz_poly_coeff_const(ff, i));
      }
      else {
        for (int i = ff->deg ; i >= 0 ; i--)
	  gmp_printf ("X%u %Zd\n", i, mpz_poly_coeff_const(ff, i));
      }
      if (format == 1) {
        for (int i = gg->deg ; i >= 0 ; i--)
	  gmp_printf ("Y%u: %Zd\n", i, mpz_poly_coeff_const(gg, i));
      }
      else {
        for (int i = gg->deg ; i >= 0 ; i--)
	  gmp_printf ("Y%u %Zd\n", i, mpz_poly_coeff_const(gg, i));
      }
      printf ("skew: %1.2f\n", skew);
      int nr = mpz_poly_number_of_real_roots (ff);
      printf ("# f lognorm %1.2f, alpha %1.2f, score %1.2f, %d rroot(s)\n",
	      logmu[0], alpha_f, logmu[0] + alpha_f, nr);
      nr = mpz_poly_number_of_real_roots(gg);
      printf ("# g lognorm %1.2f, alpha %1.2f, score %1.2f, %d rroot(s)\n",
	      logmu[1], alpha_g, logmu[1] + alpha_g, nr);
      printf ("# f+g score %1.2f\n", score);
      if (opt_flag)
        cado_poly_fprintf_MurphyE (stdout, NULL, 1, E, Bf, Bg, Area);

      printf ("\n");
      fflush (stdout);
    }
    return 1;
}

/* return the number of polynomials we look for with -skewed:
 * we compute the maximal skewness s = B^(-2/d)
 * the coefficients of degree >= d/2 are bounded by B
 * the coefficients of degree i < d/2 are bounded by B*s^(d/2-i) */
static unsigned long
get_maxtries (unsigned int B, unsigned int d)
{
  unsigned long maxtries = 1, c;
  unsigned int i;
  double skew;

  ASSERT_ALWAYS(d >= 2);

  skew = pow ((double) B, 2.0 / (double) d);
  count = malloc ((d + 1) * sizeof (unsigned long));
  for (i = 0; i <= d; i++)
    {
      if (i == d)
        c = B;         /* coefficient in 1..B */
      else if (i == d - 1)
        c = B + 1;     /* coefficient in 0..B */
      else if (2 * i >= d)
        c = 2 * B + 1; /* coefficient in -B..B */
      else /* i < d/2 */
        {
          c = B * (unsigned long) pow (skew, (double) (d - 2 * i) / 2.0);
          if (i > 0)
            c = 2 * c + 1;
          else
            c = 2 * c;
        }
      count[i] = c;
      if (maxtries >= ULONG_MAX / count[i])
        {
          fprintf (stderr, "Error, too large -bound option\n");
          exit (1);
        }
      maxtries *= count[i];
    }
  return maxtries;
}

static int
generate_f (mpz_poly_ptr ff, unsigned int d, unsigned long idx, unsigned int bound)
{
  unsigned int i, j;
  int ok = 1;

  long * a = malloc ((d + 1) * sizeof (long));

  /* the coefficient of degree j can take count[j] different values */
  for (j = 0; j <= d; j++)
    {
      long k;
      a[j] = idx % count[j];
      idx = idx / count[j];
      if (j == d)
        a[j] ++;      /* coefficient in [1..B] */
      else if (j == 0) /* coefficient in [-k..k] except 0 */
	{
          k = count[j] / 2; /* count[j]=2*k */
          /* 0..k-1 -> -k..-1
             k..2k-1 -> 1..k */
          a[j] = (a[j] < k) ? a[j] - k : a[j] - (k - 1);
	}
      else if (j < d - 1)
	{
          k = count[j] / 2; /* count[j]=2*k+1 */
          a[j] -= k;
	}
    }
  ASSERT_ALWAYS(idx == 0);

  /* Check if polynomial agrees with maximal skewness: for i > d/2 and j < d/2,
     the line going through |a[i]| and |a[j]| should not exceed bound at d/2 */
  for (i = d / 2 + 1; i <= d && ok; i++)
    {
      unsigned long ai = labs (a[i]);
      if (ai == 0)
        continue;
      for (j = 0; j < (d + 1) / 2 && ok; j++)
        {
          unsigned long aj = labs (a[j]);
          double s = pow ((double) aj / (double) ai, 1.0 / (double) (i - j));
          double mid = (double) ai * pow (s, (double) (2 * i - d) / 2.0);
          ok = mid <= (double) bound;
        }
    }

  /* Check if the reverse polynomial has smaller rank.
     This test discards about 7.7% of the polynomials for d=4 and bound=6. */
  if (ok)
    for (j = 0; 2 * j < d; j++)
      if (labs (a[d-j]) != labs(a[j]))
        {
          ok = labs (a[d-j]) < labs(a[j]);
          break;
        }

  /* Since f(x) is equivalent to f(-x), if a[d-1]=0, then the largest a[d-3],
     a[d-5], ... that is non zero should be positive. Discards 13.5% of the
     remaining polynomials. */
  if (ok)
    for (j = d + 1; j >= 2;)
      {
	j -= 2;
	if (a[j] != 0)
	  {
	    ok = a[j] > 0;
	    break;
	  }
      }

  /* Check if +-1 is a root of f. Discards 4.3% of the remaining polynomials. */
  if (ok)
    {
      int value_one = 0, value_minus_one = 0;
      for (j = 0; j <= d; j++)
	{
	  value_one += a[j];
	  value_minus_one += (j & 1) ? -a[j] : a[j];
	}
      ok = value_one != 0 && value_minus_one != 0;
    }

  /* Content test. Discards 2.3% of the remaining polynomials. */
  if (ok)
    {
      unsigned long g = a[d];
      for (j = 0; j < d; j++)
	g = gcd_int64 (g, a[j]);
      ok = g == 1;
    }

  if (ok)
      mpz_poly_setcoeffs_si(ff, a, d);

  free (a);

  return ok;
}

/*
  Generate polynomial f(x) of degree d with rank 'idx',
  with coefficients in [-bound, bound].
  The coefficient of degree d can be taken in [1, bound],
  due to the symmetry f(x) -> -f(x).
  The coefficient of degree d-1 can be taken in [0, bound],
  due to the symmetry f(x) -> f(-x).
  The coefficient of degree 0 should not be 0.
  Thus there are bound*(bound+1)*(2*bound+1)^(d-2)*(2*bound) possible values.
  Return 0 if the corresponding polynomial is not irreducible.
  Return !=0 if the poly is valid.
*/
static int
polygen_JL_f (int d, unsigned int bound, mpz_poly_ptr ff, unsigned long idx)
{
    int ok = 1;
    /* compute polynomial of index idx and check it is irreducible */
    //mpz_poly_init(ff, d);
    if (!skewed) {
        int max_abs_coeffs;
        unsigned long next_counter;
        ok = mpz_poly_setcoeffs_counter(ff, &max_abs_coeffs, &next_counter,
                d, idx, bound);
    } else {
        ok = generate_f (ff, d, idx, bound);
    }

    // to be compatible with the previous version, the count on the number of valid polys was before
    // the content test.
    // Now the content test is included in mpz_poly_setcoeffs_counter()
    // so the number of candidates will be lower than before.
    // to get the previous value, use
    // #define POLY_CONTENT -6
    // if ((ok == 1) || ((ok == 0) && (max_abs_coeffs == POLY_CONTENT))){
    // #undef POLY_CONTENT
    if (ok)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      f_candidate ++;

    /* irreducibility test */
    if (ok)
      {
	START_TIMER;
        ok = mpz_poly_squarefree_p (ff);
        /* check number of real roots */
        if (ok && (easySM || rrf != -1))
          {
            int nr = mpz_poly_number_of_real_roots(ff);
            if (easySM)
              /* check that the number of real roots is minimal */
              ok = nr == (ff->deg & 1);
            else
              ok = nr == rrf;
          }
        if (ok)
          ok = mpz_poly_is_irreducible_z (ff);
	END_TIMER (TIMER_IRRED);
      }

    return ok;
}

/* Generate polynomial g(x) of degree dg, given root 'root' of f mod N.
   It might be better to take into account the skewness of f in the LLL
   lattice, but experimentally this does not give better results (probably
   because LLL is not very sensible to a small change of the skewness).
   kN is the product k*N, where k is the multiplier. */
static void
polygen_JL_g (mpz_t kN, int dg, mpz_mat_ptr g, mpz_t root, double skew_f)
{
    int i, j;
    mpz_t a, b, det, r;
    unsigned long skew;

    skew = skew_f < 0.5 ? 1 : round (skew_f);

    mpz_init (det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    mpz_init_set (r, root);

    mpz_mat_set_ui(g, 0);

    unsigned long skew_powi = 1;
    for (int i = 0; i <= dg ; i++) {
        for (int j = 0; j <= dg ; j++) {
            if (i == 0)
              {
                if (j == 0)
                    mpz_set (mpz_mat_entry(g, j, i), kN);
                else
                  {
                    mpz_neg (mpz_mat_entry(g, j, i), r);
                    mpz_mul (r, r, root);
                  }
              }
            else if (i == j)
	      {
		ASSERT_ALWAYS((double) skew_powi * (double) skew < (double) ULONG_MAX);
		skew_powi *= skew;
		mpz_set_ui (mpz_mat_entry(g, j, i), skew_powi);
	      }
        }
    }

    START_TIMER;
    mpz_mat_LLL(det, g, NULL, a, b);
    END_TIMER (TIMER_LLL);

    /* divide row i back by skew^i */
    skew_powi = 1;
    for (i = 1;  i <= dg; i++)
      {
	skew_powi *= skew;
	for (j = 0; j <= dg; j++)
	  {
	    ASSERT_ALWAYS (mpz_divisible_ui_p (mpz_mat_entry_const(g, j, i), skew_powi));
	    mpz_divexact_ui (mpz_mat_entry(g, j, i), mpz_mat_entry_const(g, j, i), skew_powi);
	  }
      }

    mpz_clear (det);
    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (r);
}

/* JL method to generate degree d and d-1 polynomials.
   Given irreducible polynomial f of degree df, find roots of f mod n,
   and for each root, use Joux-Lercier method to find good polynomials g. */
static void
polygen_JL2 (mpz_t n,
             unsigned int df, unsigned int dg,
	     unsigned long nb_comb, mpz_poly_srcptr f, long bound2, mpz_t ell)
{
    unsigned int nr, format = 1;
    mpz_t *rf, c;
    mpz_poly u;
    long *a;
    double alpha_f;
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    ASSERT_ALWAYS (df >= 3);
    mpz_init (c);
    rf = (mpz_t *) malloc (df * sizeof(mpz_t));
    for (unsigned int i = 0; i < df; i ++)
      mpz_init (rf[i]);

    mpz_mat g;
    mpz_mat_init(g, dg + 1, dg + 1);


    mpz_poly_init (u, dg);
    u->deg = dg;
    a = malloc ((dg + 1) * sizeof (long));

    /* compute roots of the polynomial f */
    START_TIMER;
    nr = mpz_poly_roots_mpz (rf, f, n, rstate);
    END_TIMER (TIMER_ROOTS);
    ASSERT(nr <= df);

    /* update the best and worst score for f (FIXME: even if f has no roots?) */
    double skew_f, lognorm_f, score_f;
    alpha_f = get_alpha (f, get_alpha_bound ());
    skew_f = L2_skewness (f);
    lognorm_f = L2_lognorm (f, skew_f);
    score_f = lognorm_f + alpha_f;

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
        if (score_f < best_score_f)
            best_score_f = score_f;
    }

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
        if (score_f > worst_score_f)
            worst_score_f = score_f;
    }

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
      f_irreducible ++;
      sum_score_f += score_f;
    }

    /* for each root of f mod n, generate the corresponding g */
    for (unsigned int i = 0; i < nr; i ++) {
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i], skew_f);

        /* we skip idx = 0 which should correspond to c[0] = ... = c[dg] = 0 */
        for (unsigned long idx = 1; idx < nb_comb; idx ++)
          {
            unsigned long k = idx;

            /* compute first index */
            a[0] = k % (bound2 + 1);
            k = k / (bound2 + 1);
            for (unsigned int j = 1; j <= dg; j++)
              {
                a[j] = k % (2 * bound2 + 1);
                k = k / (2 * bound2 + 1);
                a[j] = (a[j] <= bound2) ? a[j] : a[j] - (2 * bound2 + 1);
              }
            ASSERT_ALWAYS(k == 0);

            /* compute the linear combination of the rows
             * \sum_{i=0}^{i=dg} a[i] * row[i]
             */
            mpz_poly_set_zero(u);
            for(unsigned int i = 0 ; i <= dg ; i++) {
                for(unsigned int j = 0 ; j <= dg ; j++) {
                    mpz_addmul_si(mpz_poly_coeff(u, j),
                            mpz_mat_entry_const(g, i, j),
                            a[i]);
                }
            }

            mpz_poly_cleandeg(u, dg);

#if 0
            /* If u is not square-free or irreducible, skip it. However, this
               test is very expensive, and non-irreducible polynomials should
               not happen in practice for large input N, thus we disable. */
            if (mpz_cmp_ui (u->coeff[0], 0) == 0 || !mpz_poly_squarefree_p (u)
                || !mpz_poly_is_irreducible_z (u))
              continue;
#endif
            /* check the number of real roots of g */
            if (easySM || rrg != -1)
              {
                int nr = mpz_poly_number_of_real_roots(u);
                int ok;
                if (easySM)
                  ok = nr == (u->deg & 1);
                else
                  ok = nr == rrg;
                if (! ok) 
                    continue; // skip this g
              }

            if (print_nonlinear_poly_info (f, alpha_f, u, format, n, ell))
              {
#if 0 /* print coefficients of record combination */
                for (j = 0; j <= dg; j++)
                  printf ("%ld ", a[j]);
                printf ("\n");
#endif
              }
        }
    }

    /* clear */
    gmp_randclear(rstate);
    free (a);

    mpz_poly_clear (u);

    for (unsigned int i = 0; i < df; i ++)
      mpz_clear (rf[i]);
    free (rf);

    mpz_mat_clear(g);

    mpz_clear (c);
}

/* JL method to generate d and d-1 polynomial.
   Generate polynomial f of degree df with |f[i]| <= bound and index 'idx'. */
static void
polygen_JL1 (mpz_t n,
             unsigned int df, unsigned int dg, unsigned int bound,
             unsigned long idx, unsigned long nb_comb, unsigned int bound2,
             mpz_t ell)
{
    mpz_poly ff;
    int irred;

    mpz_poly_init(ff, df);

    ASSERT_ALWAYS (df >= 3);

    /* generate f of degree d with small coefficients */
    irred = polygen_JL_f (df, bound, ff, idx);
    if (irred)
      polygen_JL2 (n, df, dg, nb_comb, ff, bound2, ell);

    mpz_poly_clear(ff);
}

static void
usage ()
{
    fprintf (stderr, "./dlpolyselect -N xxx -df xxx -dg xxx -bound xxx [-modr xxx] [-modm xxx] [-t xxx] [-easySM <ell>] [-skewed] [-rrf nnn]\n");
    fprintf (stderr, "Mandatory parameters:\n");
    fprintf (stderr, "   -N xxx            input number\n");
    fprintf (stderr, "   -df xxx           degree of polynomial f\n");
    fprintf (stderr, "   -dg xxx           degree of polynomial g\n");
    fprintf (stderr, "   -bound xxx        bound for absolute value of coefficients of f\n");
    fprintf (stderr, "Optional parameters:\n");
    fprintf (stderr, "   -modr r -modm m   processes only polynomials of index r mod m\n");
    fprintf (stderr, "   -t n              uses n threads\n");
    fprintf (stderr, "   -easySM ell       generates polynomials with minimal number of SMs mod ell\n");
    fprintf (stderr, "   -skewed           search for skewed polynomials\n");
    fprintf (stderr, "   -rrf nnn          f should have nnn real roots\n");
    fprintf (stderr, "   -rrg nnn          g should have nnn real roots\n");
    fprintf (stderr, "   -Bf nnn           sieving bound for f\n");
    fprintf (stderr, "   -Bg nnn           sieving bound for g\n");
    fprintf (stderr, "   -area nnn         sieving area\n");
    exit (1);
}

int main(int argc, char const * argv[])
{
    int i;
    mpz_t N;
    mpz_t ell;
    unsigned int df = 0, dg = 0;
    int nthreads = 1;
    unsigned int bound = 4; /* bound on the coefficients of f */
    unsigned long maxtries;
    double t;
    unsigned long modr = 0, modm = 1;

    t = seconds ();
    mpz_init (N);
    mpz_init (ell);

    /* printf command-line */
    printf ("#");
    for (i = 0; i < argc; i++)
        printf (" %s", argv[i]);
    printf ("\n");
    fflush (stdout);

    /* parsing */
    while (argc >= 2 && argv[1][0] == '-')
    {
        if (argc >= 3 && strcmp (argv[1], "-N") == 0) {
            mpz_set_str (N, argv[2], 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-easySM") == 0) {
            mpz_set_str (ell, argv[2], 10);
            easySM = 1;
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-df") == 0) {
            df = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-dg") == 0) {
            dg = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bound") == 0) {
            bound = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-modm") == 0) {
	    modm = strtoul (argv[2], NULL, 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-modr") == 0) {
	    modr = strtoul (argv[2], NULL, 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-t") == 0) {
            nthreads = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        /* if rrf = -1 (default), f might have any number of real roots,
           otherwise it should have exactly 'rrf' real roots */
        else if (argc >= 3 && strcmp (argv[1], "-rrf") == 0) {
            rrf = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        /* if rrg = -1 (default), g might have any number of real roots,
           otherwise it should have exactly 'rrg' real roots */
        else if (argc >= 3 && strcmp (argv[1], "-rrg") == 0) {
            rrg = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 2 && strcmp (argv[1], "-skewed") == 0) {
            skewed = 1;
            argv += 1;
            argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bf") == 0) {
            Bf = atof (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bg") == 0) {
            Bg = atof (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-area") == 0) {
            Area = atof (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else {
            fprintf (stderr, "Invalid option: %s\n", argv[1]);
            usage();
            exit (1);
        }
    }

    if (mpz_cmp_ui (N, 0) <= 0) {
        fprintf (stderr, "Error, missing input number (-N option)\n");
        usage ();
    }

    if (df == 0) {
        fprintf (stderr, "Error, missing degree (-df option)\n");
        usage ();
    }

    if (dg == 0 || dg >= df) {
        fprintf (stderr, "Error, missing or erroneous degree (-dg option): ");
        fprintf (stderr, "one should have dg < df.\n");
        usage ();
    }

    opt_flag = Bf != 0 && Bg != 0 && Area != 0;

    srand (time (NULL));

    ASSERT_ALWAYS (bound >= 1);

    if (skewed)
      maxtries = get_maxtries (bound, df);
    else {
      /* check modm has no common factor with B, B+1, 2B+1 and 2B to avoid
         a bias between classes mod 'modm' */
        if (gcd_uint64 (modm, 2 * bound) != 1 ||
          gcd_uint64 (modm, bound + 1) != 1 ||
          gcd_uint64 (modm, 2 * bound + 1) != 1)
          {
            fprintf (stderr, "Error, modm should be coprime to "
                     "2*bound*(bound+1)*(2*bound+1)\n");
            exit (1);
          }
        double maxtries_double = (double) bound;
        maxtries_double *= (double) (bound + 1);
        maxtries_double *= pow ((double) (2 * bound + 1), (double) (df - 2));
        maxtries_double *= (double) (2 * bound);
        if (maxtries_double >= (double) ULONG_MAX)
            maxtries = ULONG_MAX;
        else
            maxtries = (unsigned long) maxtries_double;
    }

    unsigned int bound2 = 1; /* bound on the coefficients of linear
                                combinations from the LLL short vectors */
    double nb_comb_f;
    unsigned long nb_comb;
    /* we try all combinations u = c[0]*v[0] + ... + c[dg]*v[dg] with
       -bound2 <= c[j] <= bound2, except 0 <= c[0] <= bound2 since u and -u
       are equivalent, and except c[0] = ... = c[dg] = 0, this gives a total
       of (bound2+1)*(2*bound2+1)^dg-1 */
    nb_comb_f = (double) (bound2 + 1) * pow ((double) (2 * bound2 + 1),
                                            (double) dg);
    nb_comb = (nb_comb_f > (double) ULONG_MAX) ? ULONG_MAX
      : (unsigned long) nb_comb_f;

    printf ("# will generate about %lu polynomials\n", maxtries / modm);

#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(dynamic)
#else
    if (nthreads > 1) {
        fprintf(stderr, "Warning: openmp unavailable, -t ignored\n");
    }
#endif
    for (unsigned long c = modr; c < maxtries; c += modm)
      polygen_JL1 (N, df, dg, bound, c, nb_comb, bound2, ell);

    t = seconds () - t;

    printf ("# found %lu irreducible f out of %lu candidates out of %lu\n",
            f_irreducible, f_candidate, maxtries / modm);
    printf ("# best f-score %1.2f, av. %1.2f, worst %1.2f, max alpha-guard %1.2f\n",
            best_score_f, sum_score_f / f_irreducible, worst_score_f,
            max_guard);
    if (max_guard > ALPHA_BOUND_GUARD)
      printf ("# Warning: max_guard > ALPHA_BOUND_GUARD, might "
              "have missed some polynomials\n");
    printf ("# Time %.2fs", t);
#ifdef MPZ_POLY_TIMINGS
    printf (" (roots %.2fs, irred %.2fs, lll %.2fs, MurphyE %.2fs)",
            timer[TIMER_ROOTS], timer[TIMER_IRRED], timer[TIMER_LLL],
	    timer[TIMER_MURPHYE]);
    printf ("\n#");
    print_timings_pow_mod_f_mod_p();
#endif
    printf ("\n");

    mpz_clear (N);
    mpz_clear (ell);
    if (skewed)
        free (count);

    return 0;
}
