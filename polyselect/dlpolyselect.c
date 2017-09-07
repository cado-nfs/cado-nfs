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

#include "cado.h"
#include "auxiliary.h"
#include "utils.h"
#include "mpz_poly.h"
#include "portability.h"
#include "murphyE.h"
#include "ropt_param.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <ctype.h>
#include <stdlib.h>
#include <time.h>

/* We assume a difference <= ALPHA_BOUND_GUARD between alpha computed
   with ALPHA_BOUND_SMALL and ALPHA_BOUND. In practice the largest value
   observed is 0.79. */
#define ALPHA_BOUND_GUARD 1.0

double best_score_f = DBL_MAX, worst_score_f = DBL_MIN;
unsigned long f_candidate = 0;   /* number of irreducibility tests */
unsigned long f_irreducible = 0; /* number of irreducible polynomials */
double max_guard = DBL_MIN;
// #define TIMINGS
#ifdef TIMINGS
double t_roots = 0.0, t_irred = 0.0, t_lll = 0.0;
#endif

/*
  Print two nonlinear poly info. Return non-zero for record polynomials.
*/
static int
print_nonlinear_poly_info (mpz_poly ff, double alpha_f, mpz_poly gg,
                           int format,  mpz_t n)
{
    unsigned int i;
    double skew, logmu[2], alpha_g_approx, alpha_g, score, score_approx;
    int df = ff->deg;
    mpz_t *f = ff->coeff;
    int dg = gg->deg;
    mpz_t *g = gg->coeff;
    static double best_score = DBL_MAX;
    /* the coefficients of g are O(n^(1/df)) */

    /* we use the skewness minimizing the sum of lognorms */
    skew = L2_combined_skewness2 (ff, gg, SKEWNESS_DEFAULT_PREC);
    logmu[1] = L2_lognorm (gg, skew);
    logmu[0] = L2_lognorm (ff, skew);
    /* first estimate alpha with a small bound */
    alpha_g_approx = get_alpha (gg, ALPHA_BOUND_SMALL);

    score_approx = logmu[1] + alpha_g_approx + logmu[0] + alpha_f;

    if (score_approx >= best_score + ALPHA_BOUND_GUARD)
      return 0;

    /* now get a more precise alpha value */
    alpha_g = get_alpha (gg, ALPHA_BOUND);

    score = logmu[1] + alpha_g + logmu[0] + alpha_f;
    if (score_approx - score > max_guard)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      max_guard = score_approx - score;

    if (score >= best_score)
      return 0; /* only print record scores */

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
      best_score = score;

      if (format == 1)
	gmp_printf ("n: %Zd\n", n);
      else
	gmp_printf ("N %Zd\n", n);

      if (format == 1) {
        for (i = df + 1; i -- != 0; )
	  gmp_printf ("c%u: %Zd\n", i, f[i]);
      }
      else {
        for (i = df + 1; i -- != 0; )
	  gmp_printf ("X%u %Zd\n", i, f[i]);
      }
      if (format == 1) {
        for (i = dg + 1; i -- != 0; )
	  gmp_printf ("Y%u: %Zd\n", i, g[i]);
      }
      else {
        for (i = dg + 1; i -- != 0; )
	  gmp_printf ("Y%u %Zd\n", i, g[i]);
      }
      printf ("skew: %1.2f\n", skew);
      printf ("# f lognorm %1.2f, alpha %1.2f, score %1.2f\n",
	      logmu[0], alpha_f, logmu[0] + alpha_f);
      printf ("# g lognorm %1.2f, alpha %1.2f, score %1.2f\n",
	      logmu[1], alpha_g, logmu[1] + alpha_g);
      printf ("# f+g score %1.2f\n", score);
      printf ("\n");
      fflush (stdout);
    }
    return 1;
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
polygen_JL_f (int d, unsigned int bound, mpz_t *f, unsigned long idx)
{
    int ok = 1;
    /* compute polynomial of index idx and check it is irreducible */
    mpz_poly ff;
    ff->deg = d;
    ff->coeff = f;
    //mpz_poly_init(ff, d);
    int max_abs_coeffs;
    unsigned long next_counter;
    ok = mpz_poly_setcoeffs_counter(ff, &max_abs_coeffs, &next_counter, d, idx, bound);

    // to be compatible with the previous version, the count on the number of valid polys was before
    // the content test.
    // Now the content test is included in mpz_poly_setcoeffs_counter()
    // so the number of candidates will be lower than before.
    // to get the previous value, use
    // #define POLY_CONTENT -6
    // if ((ok == 1) || ((ok == 0) && (max_abs_coeffs == POLY_CONTENT))){
    // #undef POLY_CONTENT
    if (ok == 1){
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      f_candidate ++;
    }
    
    /* irreducibility test */
    if(ok){
      ok = mpz_poly_squarefree_p (ff);
#ifdef TIMINGS
      t_irred -= seconds_thread ();
#endif
      if(ok){
	ok = mpz_poly_is_irreducible_z (ff);
#ifdef TIMINGS
	t_irred += seconds_thread ();
#endif
	if(ok){
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
	  f_irreducible ++;
	}
      }
    }
    return ok;
}

/* Generate polynomial g(x) of degree dg, given root 'root' of f.
   It might be better to take into account the skewness of f in the LLL
   lattice, but experimentally this does not give better results (probably
   because LLL is not very sensible to a small change of the skewness). */
static void
polygen_JL_g (mpz_t N, int dg, mat_Z g, mpz_t root)
{
    int i, j;
    mpz_t a, b, det, r;

    mpz_init (det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    mpz_init_set (r, root);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++) {
            mpz_set_ui (g.coeff[i][j], 0);
        }
    }

    for (i = 1;  i <= dg + 1; i++) {
        for (j = 1; j <= dg + 1; j++) {
            if (i == 1)
              {
                if (j == 1)
                    mpz_set (g.coeff[j][i], N);
                else
                  {
                    mpz_neg (g.coeff[j][i], r);
                    mpz_mul (r, r, root);
                  }
              }
            else
              mpz_set_ui (g.coeff[j][i], i==j);
        }
    }

#ifdef TIMINGS
    t_lll -= seconds_thread ();
#endif
    LLL (det, g, NULL, a, b);
#ifdef TIMINGS
    t_lll += seconds_thread ();
#endif

    mpz_clear (det);
    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (r);
}

/* JL method to generate d and d-1 polynomial.
   Given irreducible polynomial f of degree df, find roots of f mod n,
   and for each root, use Joux-Lercier method to find good polynomials g. */
static void
polygen_JL2 (mpz_t n, unsigned int df, unsigned int dg,
	     unsigned long nb_comb, mpz_poly f)
{
    unsigned int i, j, nr, format = 1;
    mpz_t *rf, c;
    mat_Z g;
    mpz_poly *v, u;
    long *a;
    double alpha_f;
    long bound2 = 1;

    ASSERT_ALWAYS (df >= 3);
    mpz_init (c);
    rf = (mpz_t *) malloc (df * sizeof(mpz_t));
    for (i = 0; i < df; i ++)
      mpz_init (rf[i]);
    g.coeff = (mpz_t **) malloc ((dg + 2)*sizeof(mpz_t*));
    g.NumRows = g.NumCols = dg + 1;
    for (i = 0; i <= dg + 1; i ++) {
        g.coeff[i] = (mpz_t *) malloc ((dg + 2)*sizeof(mpz_t));
        for (j = 0; j <= dg + 1; j ++) {
            mpz_init (g.coeff[i][j]);
        }
    }
    v = malloc ((dg + 1) * sizeof (mpz_poly));
    for (j = 0; j <= dg; j++)
      {
        v[j]->deg = dg;
        v[j]->coeff = g.coeff[j + 1] + 1;
      }
    mpz_poly_init (u, dg);
    u->deg = dg;
    a = malloc ((dg + 1) * sizeof (long));

    /* compute roots of the polynomial f */
#ifdef TIMINGS
    t_roots -= seconds_thread ();
#endif
    nr = mpz_poly_roots_mpz (rf, f, n);
#ifdef TIMINGS
    t_roots += seconds_thread ();
#endif
    ASSERT(nr <= df);

    // if (nr > 0)
      {
        double skew_f, lognorm_f, score_f;
        alpha_f = get_alpha (f, ALPHA_BOUND);
        skew_f = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
        lognorm_f = L2_lognorm (f, skew_f);
        score_f = lognorm_f + alpha_f;

        if (score_f < best_score_f)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          best_score_f = score_f;

        if (score_f > worst_score_f)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
          worst_score_f = score_f;
      }

    for (i = 0; i < nr; i ++) {
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i]);

        /* we skip idx = 0 which should correspond to c[0] = ... = c[dg] = 0 */
        for (unsigned long idx = 1; idx < nb_comb; idx ++)
          {
            unsigned long k = idx;

            /* compute first index */
            a[0] = k % (bound2 + 1);
            k = k / (bound2 + 1);
            for (j = 1; j <= dg; j++)
              {
                a[j] = k % (2 * bound2 + 1);
                k = k / (2 * bound2 + 1);
                a[j] = (a[j] <= bound2) ? a[j] : a[j] - (2 * bound2 + 1);
              }
            ASSERT_ALWAYS(k == 0);
            mpz_poly_mul_si (u, v[0], a[0]);
	    for (j = 1; j <= dg; j++)
	      mpz_poly_addmul_si (u, v[j], a[j]);

            /* adjust degree of u */
            for (u->deg = dg; u->deg >= 0 && mpz_cmp_ui (u->coeff[u->deg], 0)
                   == 0; u->deg--);

            if (print_nonlinear_poly_info (f, alpha_f, u, format, n))
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
    free (a);
    mpz_poly_clear (u);
    free (v);
    for (i = 0; i < df; i ++)
      mpz_clear (rf[i]);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++)
            mpz_clear (g.coeff[i][j]);
        free(g.coeff[i]);
    }
    free (g.coeff);
    free (rf);
    mpz_clear (c);
}

/* JL method to generate d and d-1 polynomial.
   Generate polynomial f of degree df with |f[i]| <= bound and index 'idx'. */
static void
polygen_JL1 (mpz_t n, unsigned int df, unsigned int dg, unsigned int bound,
             unsigned long idx, unsigned long nb_comb)
{
    unsigned int i;
    mpz_t *f;
    mpz_poly ff;
    int irred;

    ASSERT_ALWAYS (df >= 3);
    f = (mpz_t *) malloc ((df + 1)*sizeof(mpz_t));
    for (i = 0; i <= df; i ++)
      mpz_init (f[i]);

    ff->deg = df;
    ff->coeff = f;

    /* generate f of degree d with small coefficients */
    irred = polygen_JL_f (df, bound, f, idx);
    if (irred)
      polygen_JL2 (n, df, dg, nb_comb, ff);
    /* clear */
    for (i = 0; i <= df; i ++)
      mpz_clear (f[i]);
    free (f);
}

static void
usage ()
{
    fprintf (stderr, "./dlpolyselect -N xxx -df xxx -dg xxx -bound xxx [-modr xxx] [-modm xxx] [-t xxx]\n");
    exit (1);
}


int
main (int argc, char *argv[])
{
    int i;
    mpz_t N;
    unsigned int df = 0, dg = 0;
    int nthreads = 1;
    unsigned int bound = 4; /* bound on the coefficients of f */
    unsigned long maxtries;
    double t;
    unsigned long modr = 0, modm = 1;

    t = seconds ();
    mpz_init (N);

    /* printf command-line */
    printf ("#");
    for (i = 0; i < argc; i++)
        printf (" %s", argv[i]);
    printf ("\n");
    fflush (stdout);

    /* parsing */
    while (argc >= 3 && argv[1][0] == '-')
    {
        if (argc >= 3 && strcmp (argv[1], "-N") == 0) {
            mpz_set_str (N, argv[2], 10);
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
        fprintf (stderr, "Error, error degree (-df option)\n");
        usage ();
    }

    if (dg == 0 || dg >= df) {
        fprintf (stderr, "Error, missing or error degree (-dg option)\n");
        fprintf (stderr, "       only support dg < df.\n");
        usage ();
    }

    srand (time (NULL));

    ASSERT_ALWAYS (bound >= 1);

    double maxtries_double = (double) bound;
    maxtries_double *= (double) (bound + 1);
    maxtries_double *= pow ((double) (2 * bound + 1), (double) (df - 2));
    maxtries_double *= (double) (2 * bound);
    /* since each coefficient of f of degree 0 to df-1 is chosen randomly
       in [-bound, bound-1], we have (2*bound)^df possible values for f */
    if (maxtries_double >= (double) ULONG_MAX)
      maxtries = ULONG_MAX;
    else
      maxtries = (unsigned long) maxtries_double;

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

#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(dynamic)
#else
    if (nthreads > 1) {
        fprintf(stderr, "Warning: openmp unavailable, -t ignored\n");
    }
#endif
    for (unsigned long c = modr; c < maxtries; c += modm)
      polygen_JL1 (N, df, dg, bound, c, nb_comb);

    t = seconds () - t;

    printf ("found %lu irreducible f out of %lu candidates out of %lu\n",
            f_irreducible, f_candidate, maxtries / modm);
    printf ("best score %1.2f, worst %1.2f, max guard %1.2f\n",
            best_score_f, worst_score_f, max_guard);
    if (max_guard > ALPHA_BOUND_GUARD)
      printf ("Warning: max_guard > ALPHA_BOUND_GUARD, might "
              "have missed some polynomials\n");
    printf ("Time %.2fs", t);
#ifdef TIMINGS
    printf (" (roots %.2fs, irred %.2fs, lll %.2fs)", t_roots, t_irred, t_lll);
#endif
    printf ("\n");

    mpz_clear (N);

    return 0;
}
