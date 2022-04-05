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
#include "cado_poly.h"
#include "mpz_poly.h"

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
 * Regenerate raw polyomial with small c_{d-1}.
 */
static void
ropt_regen_raw (mpz_poly_srcptr F, mpz_poly_ptr G)
{
  mpz_t t, d_ad;
  mpz_init (t);
  mpz_init (d_ad);
  
  int d = mpz_poly_degree(F);
  mpz_mul_si (d_ad, F->coeff[d], d);
  mpz_neg (t, F->coeff[d-1]);
  mpz_tdiv_q (t, t, d_ad);
  mpz_mul (t, t, G->coeff[1]);
  mpz_add (G->coeff[0], G->coeff[0], t);

  mpz_clear (t);
  mpz_clear (d_ad);
}


/**
 * Ropt on a single polynomials, called by ropt_on_*() functions.
 */
static void
ropt_common ( ropt_poly_ptr poly,
              ropt_param_ptr param,
              int fm)
{
  /* bestpoly, info */
  ropt_info info;
  ropt_info_init (info);

  /* set input polynomial as the initial bestpoly,
     assuming that its content is 1. */
  ropt_bestpoly bestpoly;
  ropt_poly_setup (poly);
  ropt_bestpoly_init (bestpoly);
  ropt_bestpoly_setup (bestpoly, poly->pols[1], poly->pols[0]);

  /* reducing c5 to raw and skip ropt */
  if (param->gen_raw) {
    mpz_t res;
    mpz_poly_ptr F = poly->pols[1];
    mpz_poly_ptr G = poly->pols[0];
    mpz_init (res);

    /* original polynomial */
    print_poly_fg (F, G, poly->n, 1);

    ASSERT_ALWAYS(mpz_poly_degree(G) == 1);

    Lemma21 (F, poly->n, G, res);
    mpz_ptr lc = F->coeff[mpz_poly_degree(F)];
    mpz_div (lc, lc, res);

    Lemma21 (F, poly->n, G, res);
    ropt_regen_raw (F, G);

    Lemma21 (F, poly->n, G, res);

    mpz_clear (res);
  }

  /* print original or reduced polynomial */
  if (param->gen_raw)
    fprintf (stderr, "\n# Reduced polynomial.\n");
  if (fm==0)
    print_poly_fg (poly->pols[1], poly->pols[0], poly->n, 1);

  /* if sopt and not from -fm */
  if (fm==0) {
    if (param->sopt) {
      mpz_poly_ptr F = poly->pols[1];
      mpz_poly_ptr G = poly->pols[0];
      size_optimization (F, G, F, G, SOPT_DEFAULT_EFFORT, 0);
      fprintf (stderr, "\n# Size-optimized polynomial.\n");
      print_poly_fg (F, G, poly->n, 1);
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
    print_poly_fg (bestpoly->pols[1], bestpoly->pols[0], poly->n, 1);
  }
  ropt_info_clear (info);
  ropt_bestpoly_clear (bestpoly);
}



/**
 * Ropt on one polynomial from stdin.
 */
void
ropt_on_stdin ( ropt_param_ptr param )
{
    ropt_on_cadopoly(stdin, param);
}

/**
 * Ropt on all CADO-format polynomials in the file.
 */
void
ropt_on_cadopoly ( FILE *file,
                   ropt_param_ptr param )
{
  cado_poly cpoly;
  cado_poly_init(cpoly);

  /* poly struct */
  ropt_poly poly;
  ropt_poly_init (poly);

  for(int count = 0; cado_poly_read_next_poly_from_stream(cpoly, file); count++) {
      ropt_poly_refresh (poly);

      if (!ropt_poly_setup_check (poly)) 
          continue;

      fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
      ropt_common (poly, param, 0);
  } 

  cado_poly_clear(cpoly);
  ropt_poly_clear (poly);
}



/**
 * Ropt on all Msieve-format polynomials in the file.
 */
void
ropt_on_msievepoly ( FILE *file,
                     ropt_param_ptr param )
{
  unsigned int count = 0U;
  char str[MAX_LINE_LENGTH];

  /* poly struct */
  ropt_poly poly;
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
    mpz_set (poly->n, param->n);

    mpz_poly_ptr F = poly->pols[1];
    mpz_poly_ptr G = poly->pols[0];

    mpz_poly_set_zero(F);
    mpz_poly_setcoeff(F, param->d, ad);

    mpz_poly_set_mpz_ab(G, m, l); /* sets to m-lx */
    mpz_poly_neg(G, G); /* we want lx-m */

    /* generate polynomial */
    Lemma21 (F, poly->n, G, res);

    /* scale ad back and generate poly with small c5 */
    if (param->gen_raw) {
        mpz_ptr lc = F->coeff[mpz_poly_degree(F)];
        mpz_div (lc, lc, res);
        Lemma21 (F, poly->n, G, res);
        ropt_regen_raw (F, G);
        Lemma21 (F, poly->n, G, res);
    }

    /* print before size optimization */
    fprintf (stderr, "\n# Polynomial (# %5d).", count);
    print_poly_fg (F, G, poly->n, 1);

    /* optimize size */
    if (param->sopt) {
      fprintf (stderr, "\n# Size-optimize only (# %5d).", count);
      size_optimization (F, G, F, G, SOPT_DEFAULT_EFFORT, 0);

      /* print in Msieve format */
      // print_poly_info_short (poly->f, poly->g, poly->d, poly->n);

      /* print in CADO format */
      print_poly_fg (F, G, poly->n, 1);
    }
    
    /* skip or do ropt */
    if (param->skip_ropt == 0 && param->gen_raw == 0)
      ropt_common (poly, param, 1);

    count ++;
  }

  /* free */
  mpz_clear (ad);
  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (res);
  ropt_poly_clear (poly);
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
