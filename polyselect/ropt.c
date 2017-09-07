/**
 * @file ropt.c
 * This is the main interface, which is called from
 * either polyselect2*.c or ropt.c. There are three
 * main functions:
 * - ropt_do_both_stages()
 * - ropt_do_stage2()
 * - ropt_polyselect()
 */


#include "cado.h"
#include "ropt.h"
#include "portability.h"
#include "area.h"
#include "size_optimization.h"


/**
 * Find best poly. This is somehow redundant.
 */
void
ropt_get_bestpoly ( ropt_poly_t poly,
                    MurphyE_pq *global_E_pqueue,
                    ropt_bestpoly_t bestpoly )
{
  double ave_MurphyE = 0.0, best_E = 0.0;
  int i, old_i, k;
  mpz_t t, *fuv, *guv;
  mpz_poly Fuv, Guv;

  mpz_init (t);

  /* var for computing E */
  mpz_poly_init (Fuv, poly->d);
  mpz_poly_init (Guv, 1);
  Fuv->deg = poly->d;
  Guv->deg = 1;
  fuv = Fuv->coeff;
  guv = Guv->coeff;
  for (i = 0; i <= poly->d; i++)
    mpz_set (fuv[i], poly->f[i]);
  for (i = 0; i < 2; i++)
    mpz_set (guv[i], poly->g[i]);

  /* output all polys in the global queue */
  for (i = 1; i < global_E_pqueue->used; i ++) {

    old_i = 0;
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i,
                        global_E_pqueue->w[i], 2);

    for (k = 0; k <= poly->d; k++)
      mpz_set (fuv[k], poly->f[k]);

    for (k = 0; k < 2; k++)
      mpz_set (guv[k], poly->g[k]);

    compute_fuv_mp (fuv, poly->f, poly->g, poly->d, global_E_pqueue->u[i],
                    global_E_pqueue->v[i]);

    sopt_local_descent (Fuv, Guv, Fuv, Guv, 1, -1, SOPT_DEFAULT_MAX_STEPS, 0);
    fuv = Fuv->coeff;
    guv = Guv->coeff;

    ave_MurphyE = print_poly_fg (Fuv, guv, poly->n, 0);

    mpz_poly_content (t, Fuv);
    for (k = 0; k <= poly->d; k++) {
      mpz_div (fuv[k], fuv[k], t);
    }
    mpz_poly_content (t, Fuv);

    if ( (ave_MurphyE > best_E) && (mpz_cmp_ui (t, 1) == 0) ) {
      best_E = ave_MurphyE;
      for (k = 0; k <= poly->d; k++)
        mpz_set (bestpoly->f[k], fuv[k]);
      for (k = 0; k < 2; k++)
        mpz_set (bestpoly->g[k], guv[k]);
    }
    rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  }

  mpz_poly_clear (Fuv);
  mpz_poly_clear (Guv);
  mpz_clear (t);
}


/**
 * Root sieve only.
 */
static void
ropt_do_stage2 (ropt_poly_t poly,
                ropt_bestpoly_t bestpoly,
                ropt_param_t param,
                ropt_info_t info )
{

  int i, old_i;
  double alpha_lat;
  mpz_t *fuv, *guv;
  ropt_bound_t bound;
  ropt_s2param_t s2param;
  MurphyE_pq *global_E_pqueue;
  mpz_poly Fuv;

  ropt_bound_init (bound);
  new_MurphyE_pq (&global_E_pqueue, 4);

  mpz_poly_init (Fuv, poly->d);
  Fuv->deg = poly->d;
  fuv = Fuv->coeff;
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in "
             "ropt_do_stage2().\n");
    exit (1);
  }
  for (i = 0; i <= poly->d; i++)
    mpz_set (fuv[i], poly->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], poly->g[i]);

  /* rotate polynomial by f + rot*x^2 */
  old_i = 0;
  old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, param->s2_w, 2);
  //print_poly_fg (poly->f, poly->g, poly->d, poly->n, 1);

  /* reset after rotation */
  ropt_poly_setup (poly);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);

  /* print some basic information */
  compute_fuv_mp (fuv, poly->f, poly->g, poly->d, param->s2_u, param->s2_v);
  alpha_lat = get_alpha (Fuv, ALPHA_BOUND);
  gmp_fprintf ( stderr,
                "\n# Info: Sieve on sublattice, (w, u, v): (%d, %Zd, %Zd) "
                "(mod %Zd)\n"
                "# Info: alpha: %.2f, proj_alpha: %.2f, "
                "exp_alpha: %.2f\n",
                param->s2_w,
                param->s2_u,
                param->s2_v,
                param->s2_mod,
                alpha_lat,
                poly->alpha_proj,
                bound->exp_min_alpha );

  /* root sieve */
  ropt_s2param_init (poly, s2param);
  ropt_s2param_setup_stage2_only (bound, s2param, param,
                                  param->s2_u, param->s2_v, param->s2_mod);
  info->w = param->s2_w;
  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, param->s2_w);

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  old_i = 0;

  /* return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  ropt_bound_free (bound);
  ropt_s2param_free (poly, s2param);
  mpz_poly_clear (Fuv);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  free (guv);

}


/**
 * Ropt linear or quadratic.
 */
static void
ropt_do_both_stages ( ropt_poly_t poly,
                      ropt_bestpoly_t bestpoly,
                      ropt_param_t param,
                      ropt_info_t info )
{
  if (poly->d == 5 || poly->d == 4 || poly->d == 3)
    ropt_linear (poly, bestpoly, param, info);
  else if (poly->d == 6 || poly->d == 7)
    ropt_quadratic (poly, bestpoly, param, info);
  else {
    fprintf (stderr, "Error: only support deg 3, 4, 5, 6 and 7.\n");
    exit(1);
  }
}


/**
 * Start the two-stage root optimization on poly.
 * @param poly contains polynomial.
 * @param param input parameters.
 * @return top polynomial in bestpoly.
 */
void
ropt ( ropt_poly_t poly,
       ropt_bestpoly_t bestpoly,
       ropt_param_t param,
       ropt_info_t info )
{

  /* print cache size */
  if (param->verbose == 2)
    fprintf ( stderr, "# Info: L1_cachesize: %d, "
              "size_tune_sievearray: %d\n",
              L1_cachesize, size_tune_sievearray );

  if (param->stage_flag == 2)
    ropt_do_stage2 (poly, bestpoly, param, info);
  else
    ropt_do_both_stages (poly, bestpoly, param, info);

}


/**
 * Called by polyselect_ropt.
 */
void
ropt_polyselect (cado_poly_ptr output_poly, cado_poly_ptr input_poly,
                 ropt_param_t param, ropt_time_t eacht)
{
  int i;
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* setup poly */
  for (i = 0; i <= input_poly->pols[RAT_SIDE]->deg; i++)
    mpz_set (poly->g[i], input_poly->pols[RAT_SIDE]->coeff[i]);
  for (i = 0; i <= input_poly->pols[ALG_SIDE]->deg; i++)
    mpz_set (poly->f[i], input_poly->pols[ALG_SIDE]->coeff[i]);
  mpz_set (poly->n, input_poly->n);
  ropt_poly_setup (poly);

  ropt_info_t info;
  ropt_info_init (info);


  ropt_bestpoly_t bestpoly;
  ropt_bestpoly_init (bestpoly, poly->d);
  ropt_bestpoly_setup (bestpoly, poly->f, poly->g, poly->d);

  /* call main function */
  ropt_do_both_stages (poly, bestpoly, param, info);
  
  /* bring bestpoly back to polyselect_ropt */
  for (i = 0; i <= input_poly->pols[RAT_SIDE]->deg; i++)
    mpz_set (output_poly->pols[RAT_SIDE]->coeff[i], bestpoly->g[i]);
  mpz_poly_cleandeg (output_poly->pols[RAT_SIDE], input_poly->pols[RAT_SIDE]->deg);
  for (i = 0; i <= input_poly->pols[ALG_SIDE]->deg; i++)
    mpz_set (output_poly->pols[ALG_SIDE]->coeff[i], bestpoly->f[i]);
  mpz_poly_cleandeg (output_poly->pols[ALG_SIDE], input_poly->pols[ALG_SIDE]->deg);
  mpz_set (output_poly->n, input_poly->n);

  /* get time passed from info, use info to keep interface unchanged */
  eacht->ropt_time_stage1 = info->ropt_time_stage1;
  eacht->ropt_time_tuning = info->ropt_time_tuning;
  eacht->ropt_time_stage2 = info->ropt_time_stage2;
  
  /* free */
  ropt_bestpoly_free (bestpoly, poly->d);
  ropt_info_free (info);
  ropt_poly_free (poly);
}
