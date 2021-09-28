/**
 * @file ropt_io.c
 * Ease input and output effort for ropt_main.c.
 */


#include "cado.h" // IWYU pragma: keep
#include <stdlib.h> // exit
#include <string.h> // strncmp strlen
#include <math.h> // log
#include <stdio.h>
#include <gmp.h>
#include "cado_poly.h"
#include "mpz_poly.h"
#include "ropt_io.h"
#include "ropt.h"
#include "ropt_str.h" // ropt_poly_t ...
#include "size_optimization.h"
#include "usp.h"
#include "murphyE.h"    // MurphyE MURPHY_K
#include "ropt_param.h"    // L1_cachesize MAX_LINE_LENGTH
#include "cachesize_cpuid.h" // cachesize_cpuid
#include "auxiliary.h" // print_poly_fg ALG_SIDE ALPHA_BOUND
#include "ropt_arith.h" // Lemma21
#include "polyselect_norms.h"
#include "polyselect_alpha.h"

/**
 * Get L1 cache size in the beginning.
 */
void
ropt_L1_cachesize ()
{
  /* L1 cache size */
  unsigned int ret = cachesize_cpuid (0);
  if ((2048 <= ret) && (ret <= (1U << 20)))
    L1_cachesize = ret;
#if 0
  else {
    ret = cachesize_guess (0);
    if ((2048 <= ret) && (ret <= (1 << 20))) {
      L1_cachesize = ret;
    }
  }
#endif
  /* since we sieve from -size_tune_sievearray to size_tune_sievearray,
     we set size_tune_sievearray to half the L1 cache size */
  size_tune_sievearray = L1_cachesize / 2;
}


/**
 * Ropt on a single polynomials, called by ropt_on_*() functions.
 */
static void
ropt_common ( ropt_poly_t poly,
              ropt_param_t param,
              int fm)
{
  /* bestpoly, info */
  ropt_info_t info;
  ropt_info_init (info);

  /* set input polynomial as the initial bestpoly,
     assuming that its content is 1. */
  ropt_bestpoly_t bestpoly;
  ropt_poly_setup (poly);
  ropt_bestpoly_init (bestpoly, poly->d);
  ropt_bestpoly_setup (bestpoly, poly->f, poly->g, poly->d);

  /* reducing c5 to raw and skip ropt */
  if (param->gen_raw) {
    mpz_t l, m, res;
    mpz_poly F;
    mpz_init (l);
    mpz_init (m);
    mpz_init (res);
    mpz_poly_init (F, poly->d);
    F->deg = poly->d;
    for (int j = 0; j <= poly->d; j++) {
      mpz_set (F->coeff[j], poly->f[j]);
    }
    /* original polynomial */
    print_poly_fg (F, poly->g, poly->n, 1);
    mpz_set (l, poly->g[1]);
    mpz_neg (m, poly->g[0]);
    Lemma21 (poly->f, poly->n, poly->d, l, m, res);
    mpz_div (poly->f[poly->d], poly->f[poly->d], res);
    Lemma21 (poly->f, poly->n, poly->d, l, m, res);
    ropt_regen_raw (poly->f, poly->d, poly->g);
    mpz_set (l, poly->g[1]);
    mpz_neg (m, poly->g[0]);
    Lemma21 (poly->f, poly->n, poly->d, l, m, res);
    mpz_poly_clear (F);
    mpz_clear (l);
    mpz_clear (m);
    mpz_clear (res);
  }

  mpz_poly F;
  F->coeff = poly->f;
  F->deg = poly->d;
  /* print original or reduced polynomial */
  if (param->gen_raw)
    fprintf (stderr, "\n# Reduced polynomial.\n");
  if (fm==0)
    print_poly_fg (F, poly->g, poly->n, 1);

  /* if sopt and not from -fm */
  if (fm==0) {
    if (param->sopt) {
      mpz_poly F2, G2;
      mpz_poly_init (F2, poly->d);
      mpz_poly_init (G2, 1);
      mpz_poly_setcoeffs (F2, poly->f, poly->d);
      mpz_poly_setcoeffs (G2, poly->g, 1);
      size_optimization (F2, G2, F2, G2, SOPT_DEFAULT_EFFORT, 0);
      fprintf (stderr, "\n# Size-optimized polynomial.\n");
      print_poly_fg (F2, G2->coeff, poly->n, 1);
      for (int j = 0; j <= poly->d; j++)
        mpz_set (poly->f[j], F2->coeff[j]);
      mpz_set (poly->g[0], G2->coeff[0]);
      mpz_set (poly->g[1], G2->coeff[1]);
      mpz_poly_clear (F2);
      mpz_poly_clear (G2);
    }
  }
  
  if (!param->gen_raw) {

    if (param->verbose) {
      fprintf(stderr, "# Info: verbose level: %d\n", param->verbose);
      fprintf(stderr, "# Info: sieving effort: %.0f\n", param->effort);
    }
    
    /* call ropt */
    ropt (poly, bestpoly, param, info);
    fprintf (stderr, "\n# Info: Best E is:\n");
    F->coeff = bestpoly->f;
    print_poly_fg (F, bestpoly->g, poly->n, 1);
  }

  ropt_info_free (info);
  ropt_bestpoly_free (bestpoly, poly->d);
}


/**
 * Read a line.
 */
static inline int
ropt_on_stdin_readpoly_readline ( char *s )
{
  int c;

  while ((c = getchar ()) != '\n' && c != EOF)
    *s++ = c;
  *s = '\0';
  fprintf (stderr, "%s", s);
  return c;
}


/**
 * Read a polynomial from stdin.
 */
static void
ropt_on_stdin_readpoly ( mpz_t N,
                         mpz_t *f,
                         mpz_t *g,
                         mpz_t M,
                         int *degree )
{
  fflush(stdin);
  int i, ret, d;

  /* input buffer */
  char s[MAX_LINE_LENGTH];

  while (feof (stdin) == 0) {
    ret = ropt_on_stdin_readpoly_readline(s);
    if (ret == EOF)
      break;
    if (strlen (s) + 1 >= MAX_LINE_LENGTH) {
      fprintf (stderr, "Error, buffer overflow\n");
      exit (1);
    }
    /* n: input number */
    if (strncmp (s, "n:", 2) == 0) {
      if (mpz_set_str (N, s + 2, 0) != 0) {
        fprintf (stderr, "Error while reading N\n`");
        exit (1);
      }
    }
    /* ci: coeff of degree i */
    else if (sscanf (s, "c%d:", &i) == 1) {
      if (i > MAX_DEGREE) {
        fprintf (stderr, "Error, too large degree %d\n", i);
        exit (1);
      }
      if (mpz_set_str (f[i], s + 3, 0) != 0) {
        fprintf (stderr, "Error while reading f[%d]", i);
        exit (1);
      }
    }
    else if (strncmp (s, "Y1:", 3) == 0) {
      if (mpz_set_str (g[1], s + 3, 0) != 0) {
        fprintf (stderr, "Error while reading Y1");
        exit (1);
      }
    }
    else if (strncmp (s, "Y0:", 3) == 0) {
      if (mpz_set_str (g[0], s + 3, 0) != 0)
      {
        fprintf (stderr, "Error while reading Y0");
        exit (1);
      }
    }
    else if (strncmp (s, "# M", 3) == 0 || strncmp (s, "m:", 2) == 0) {
      if (mpz_set_str (M, s + 2 + (s[0] == '#'), 0) != 0) {
        fprintf (stderr, "Error while reading M or m");
        exit (1);
      }
    }
  }

  for (d = MAX_DEGREE; d > 0 && mpz_cmp_ui (f[d], 0) == 0; d --);
  *degree = d;
  if (mpz_cmp_ui (M, 0) == 0) {
    mpz_t t;
    /* M = -Y0/Y1 mod N */
    mpz_invert (M, g[1], N);
    mpz_neg (M, M);
    mpz_mul (M, M, g[0]);
    mpz_mod (M, M, N);
    /* check M is also a root of the algebraic polynomial mod N */
    mpz_init_set (t, f[d]);
    for (i = d - 1; i >= 0; i --) {
      mpz_mul (t, t, M);
      mpz_add (t, t, f[i]);
      mpz_mod (t, t, N);
    }
    if (mpz_cmp_ui (t, 0) != 0) {
      fprintf (stderr, "Polynomials have no common root mod N\n");
      exit (1);
    }
    mpz_clear (t);
  }
}


/**
 * Ropt on one polynomial from stdin.
 */
void
ropt_on_stdin ( ropt_param_t param )
{

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* read poly from stdin */
  ropt_on_stdin_readpoly (poly->n, poly->f, poly->g, poly->m, &(poly->d));
  param->d = poly->d;

  /* run ropt */
  fprintf (stderr, "\n# Polynomial (# 0).\n");
  ropt_common (poly, param, 0);

  ropt_poly_free (poly);
}


/**
 * Ropt on all CADO-format polynomials in the file.
 */
void
ropt_on_cadopoly ( FILE *file,
                   ropt_param_t param )
{
  unsigned int flag = 0U, count = 0U;
  char str[MAX_LINE_LENGTH];

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* parse each line */
  while (1) {

    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;

    if ( str[0] == 'Y' ) {
      if ( str[1] == '1' ) {
        gmp_sscanf (str, "Y1: %Zd\n", poly->g[1]);
        (flag) ^= (1<<9);
      }
      else if ( str[1] == '0' ) {
        gmp_sscanf (str, "Y0: %Zd\n", poly->g[0]);
        (flag) ^= (1<<8);
      }
      else {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'c') {
      if ( str[1] == '0' ) {
        gmp_sscanf (str, "c0: %Zd\n", poly->f[0]);
        (flag) ^= 1;
      }
      else if ( str[1] == '1' ) {
        gmp_sscanf (str, "c1: %Zd\n", poly->f[1]);
        (flag) ^= (1<<1);
      }
      else if ( str[1] == '2' ) {
        gmp_sscanf (str, "c2: %Zd\n", poly->f[2]);
        (flag) ^= (1<<2);
      }
      else if ( str[1] == '3' ) {
        gmp_sscanf (str, "c3: %Zd\n", poly->f[3]);
        (flag) ^= (1<<3);
      }
      else if ( str[1] == '4' ) {
        gmp_sscanf (str, "c4: %Zd\n", poly->f[4]);
        (flag) ^= (1<<4);
      }
      else if ( str[1] == '5' ) {
        gmp_sscanf (str, "c5: %Zd\n", poly->f[5]);
        (flag) ^= (1<<5);
      }
      else if ( str[1] == '6' ) {
        gmp_sscanf (str, "c6: %Zd\n", poly->f[6]);
        (flag) ^= (1<<6);
      }
      else if ( str[1] == '7' ) {
        gmp_sscanf (str, "c7: %Zd\n", poly->f[6]);
        (flag) ^= (1<<7);
      }
      else {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'n') {
      gmp_sscanf (str, "n: %Zd\n", poly->n);
      (flag) ^= (1<<10);
    }
    else
      continue;

    if ( flag == 2047 || // deg 7
         flag == 1919 || // deg 6
         flag == 1855 || // deg 5
         flag == 1823 || // deg 4
         flag == 1807 ) // deg 3
    {
      if (ropt_poly_setup_check (poly)) {
        fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
        ropt_common (poly, param, 0);
        count ++;
        flag = 0U;
        ropt_poly_refresh (poly);
      }
    }
  } // next line

  ropt_poly_free (poly);
}


/**
 * Regenerate raw polyomial with small c_{d-1}.
 */
void
ropt_regen_raw (mpz_t *f, int d, mpz_t *g)
{
  mpz_t t, d_ad;
  mpz_init (t);
  mpz_init (d_ad);
  
  mpz_mul_si (d_ad, f[d], d);
  mpz_neg (t, f[d-1]);
  mpz_tdiv_q (t, t, d_ad);
  mpz_mul (t, t, g[1]);
  mpz_add (g[0], g[0], t);

  mpz_clear (t);
  mpz_clear (d_ad);
}


/**
 * Ropt on all Msieve-format polynomials in the file.
 */
void
ropt_on_msievepoly ( FILE *file,
                     ropt_param_t param )
{
  unsigned int count = 0U;
  char str[MAX_LINE_LENGTH];

  /* poly struct */
  ropt_poly_t poly;
  ropt_poly_init (poly);

  mpz_t ad, l, m, res;
  mpz_init (ad);
  mpz_init (l);
  mpz_init (m);
  mpz_init (res);

  /* parse each line */
  while (1) {

    /* read ad, l, m */   
    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;
    int ret = gmp_sscanf (str, "%Zd %Zd %Zd", ad, l, m);
    if (ret != 3)
      continue;

    /* set d, n, ad, l, m */
    poly->d = param->d;
    mpz_set (poly->n, param->n);
    mpz_set (poly->f[poly->d], ad);
    mpz_set (poly->g[1], l);
    mpz_neg (poly->g[0], m);

    /* generate polynomial */
    Lemma21 (poly->f, poly->n, poly->d, l, m, res);

    /* scale ad back and generate poly with small c5 */
    if (param->gen_raw) {
      mpz_div (poly->f[poly->d], poly->f[poly->d], res);
      Lemma21 (poly->f, poly->n, poly->d, l, m, res);
      ropt_regen_raw (poly->f, poly->d, poly->g);
      mpz_neg (m, poly->g[0]);
      Lemma21 (poly->f, poly->n, poly->d, l, m, res);
    }

    mpz_poly F, G;
    mpz_poly_init (F, poly->d);
    mpz_poly_init (G, 1);
    mpz_poly_setcoeffs (F, poly->f, poly->d);
    mpz_poly_setcoeffs (G, poly->g, 1);

    /* print before size optimization */
    fprintf (stderr, "\n# Polynomial (# %5d).", count);
    print_poly_fg (F, G->coeff, poly->n, 1);

    /* optimize size */
    if (param->sopt) {
      fprintf (stderr, "\n# Size-optimize only (# %5d).", count);
      size_optimization (F, G, F, G, SOPT_DEFAULT_EFFORT, 0);

      /* print in Msieve format */
      // print_poly_info_short (poly->f, poly->g, poly->d, poly->n);

      /* print in CADO format */
      print_poly_fg (F, G->coeff, poly->n, 1);
    }
    
    /* skip or do ropt */
    if (param->skip_ropt == 0 && param->gen_raw == 0) {
      for (int j = 0; j <= poly->d; j++)
        mpz_set (poly->f[j], F->coeff[j]);
      mpz_set (poly->g[0], G->coeff[0]);
      mpz_set (poly->g[1], G->coeff[1]);
      ropt_common (poly, param, 1);
    }

    mpz_poly_clear (F);
    mpz_poly_clear (G);
    count ++;
  }

  /* free */
  mpz_clear (ad);
  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (res);
  ropt_poly_free (poly);
}



/**
 * Print ad, l, m and polynomial information.
 */
double
print_poly_info_short ( mpz_t *f,
                        mpz_t *g,
                        int d,
                        mpz_t N )
{
  /* print info about the polynomial */
  unsigned int nroots = 0;
  double skew, logmu, alpha, e, alpha_proj;
  int i;
  double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
  mpz_poly F;
  F->coeff = f;
  F->deg = d;

  cado_poly cpoly;
  cado_poly_init(cpoly);
  for (i = 0; i < (d + 1); i++) {
    mpz_set(cpoly->pols[ALG_SIDE]->coeff[i], f[i]);
  }
  for (i = 0; i < 2; i++) {
    mpz_set(cpoly->pols[RAT_SIDE]->coeff[i], g[i]);
  }
  mpz_set (cpoly->n, N);
  cpoly->pols[ALG_SIDE]->deg = d;
  cpoly->pols[RAT_SIDE]->deg = 1;

  /* output original poly */
  gmp_printf ("%Zd ", f[d]);
  mpz_neg (g[0], g[0]);
  for (i = 1; i >= 0; i --) {
    gmp_printf ("%Zd ", g[i]);
  }
  printf ("\n");
  mpz_neg (g[0], g[0]);
  
  /* compute skew, logmu, nroots */
  nroots = numberOfRealRoots ((const mpz_t *) cpoly->pols[ALG_SIDE]->coeff, d, 0, 0, NULL);
  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
  cpoly->skew = skew;
  logmu = L2_lognorm (F, skew);
  alpha = get_alpha (F, get_alpha_bound ());
  alpha_proj = get_alpha_projective (F, get_alpha_bound ());
  e = MurphyE (cpoly, bound_f, bound_g, area, MURPHY_K, get_alpha_bound ());

  printf ("# lognorm: %.2f, alpha: %.2f, (proj: %.2f) E: %.2f, nr: %u, exp_E: %1.2f, MurphyE: %1.2e\n",
          logmu,
          alpha,
          alpha_proj,
          logmu + alpha,
          nroots,
          logmu + expected_alpha(exp_rot[d]*log(skew)),
          e );

  fflush( stdout );
  cado_poly_clear (cpoly);

  return e;
}
