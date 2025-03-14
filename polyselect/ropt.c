/**
 * @file ropt.c
 * This is the main interface, which is called from
 * either polyselect2*.c or ropt.c. There are three
 * main functions:
 * - ropt_do_both_stages()
 * - ropt_do_stage2()
 * - ropt_polyselect()
 */


#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "cado_poly.h"
#include "macros.h"
#include "mpz_poly.h"
#include "ropt.h"
#include "ropt_arith.h"
#include "ropt_linear.h"
#include "ropt_quadratic.h"
#include "ropt_stage2.h"
#include "ropt_str.h"
#include "ropt_tree.h"
#include "ropt_param.h"
#include "size_optimization.h"
#include "auxiliary.h"


/**
 * Find best poly. This is somehow redundant.
 */
void
ropt_get_bestpoly ( ropt_poly_srcptr poly,
                    MurphyE_pq *global_E_pqueue,
                    ropt_bestpoly_ptr bestpoly
                    )
{
  double ave_MurphyE = 0.0, best_E = 0.0;

  /* var for computing E */
  mpz_poly Fuv, Guv;
  mpz_poly_init (Fuv, -1);
  mpz_poly_init (Guv, -1);

  /* output all polys in the global queue */
  /* WTF starting with 1? Why? XXX */
  for (int i = 1; i < global_E_pqueue->used; i ++) {
      mpz_poly_set(Fuv, poly->cpoly->pols[1]);
      mpz_poly_set(Guv, poly->cpoly->pols[0]);

      /* rotation is by u*x+v+w*x^2 (yes, interesting naming...)
       * Note that u, v, w have different types, and interestingly we
       * have different functions to do this...
       */
      compute_fuv_mp (Fuv, poly->cpoly->pols[1], poly->cpoly->pols[0],
              global_E_pqueue->u[i], global_E_pqueue->v[i]);

      if (global_E_pqueue->w[i])
          rotate_aux (Fuv, Guv, 0, global_E_pqueue->w[i], 2);

      /* This modifies Fuv and Guv */
      sopt_local_descent (Fuv, Guv, Fuv, Guv, 1, -1, SOPT_DEFAULT_MAX_STEPS, 0);

      ave_MurphyE = print_poly_fg (Fuv, Guv, poly->cpoly->n, 0);

      /* Does it even make sense ? After all, if we end up with something
       * with non-trivial content, we would expect that we encountered it
       * earlier, right ?
       */
      mpz_poly_divide_by_content (Fuv);

      if (ave_MurphyE > best_E) {
          best_E = ave_MurphyE;
          mpz_poly_set (bestpoly->f, Fuv);
          mpz_poly_set (bestpoly->g, Guv);
      }
  }

  mpz_poly_clear (Fuv);
  mpz_poly_clear (Guv);
}


/**
 * Root sieve only.
 */
static void
ropt_do_stage2 (ropt_poly_ptr poly,
                ropt_bestpoly_ptr bestpoly,
                ropt_param_ptr param,
                ropt_info_ptr info)
{
  double alpha_lat;

  // int old_i;
  // double alpha_lat;
  // ropt_bound bound;
  // ropt_s2param s2param;
  // MurphyE_pq *global_E_pqueue;
  // mpz_poly Fuv, Guv;

  mpz_poly Fuv, Guv;
  mpz_poly_init (Fuv, -1);
  mpz_poly_init (Guv, -1);

  ropt_bound bound;
  ropt_bound_init (bound);

  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  mpz_poly_set(Fuv, poly->cpoly->pols[1]);
  mpz_poly_set(Guv, poly->cpoly->pols[0]);

  /* rotate polynomial by f + rot*x^2 */
  rotate_aux (Fuv, Guv, 0, param->s2_w, 2);

  /* reset after rotation */

  ropt_poly_setup (poly);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);

  /* print some basic information */
  compute_fuv_mp (Fuv, Fuv, Guv, param->s2_u, param->s2_v);

  alpha_lat = get_alpha (Fuv, get_alpha_bound ());
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
  ropt_s2param s2param;

  ropt_s2param_init (s2param);
  ropt_s2param_setup_stage2_only (bound, s2param, param,
                                  param->s2_u, param->s2_v, param->s2_mod);
  info->w = param->s2_w;
  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, param->s2_w);

  /* return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  ropt_bound_clear (bound);
  ropt_s2param_clear(s2param);
  mpz_poly_clear (Fuv);
  mpz_poly_clear (Guv);
}


/**
 * Ropt linear or quadratic.
 */
static void
ropt_do_both_stages ( ropt_poly_ptr poly,
        ropt_bestpoly_ptr bestpoly,
        ropt_param_ptr param,
        ropt_info_ptr info)
{
     int d = mpz_poly_degree(poly->cpoly->pols[1]);
     if (d == 5 || d == 4 || d == 3)
         ropt_linear (poly, bestpoly, param, info);
     else if (d == 6 || d == 7)
         ropt_quadratic (poly, bestpoly, param, info);
     else {
         fprintf (stderr, "Error: only support deg 3, 4, 5, 6 and 7.\n");
         exit(EXIT_FAILURE);
     }
}


/**
 * Start the two-stage root optimization on poly.
 * @param poly contains polynomial.
 * @param param input parameters.
 * @return top polynomial in bestpoly.
 */
void
ropt ( ropt_poly_ptr poly,
       ropt_bestpoly_ptr bestpoly,
       ropt_param_ptr param,
       ropt_info_ptr info)
{

  /* print cache size */
  if (param->verbose == 2)
    fprintf ( stderr, "# Info: L1_cachesize: %zu, "
              "size_tune_sievearray: %zu\n",
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
ropt_polyselect (cado_poly_ptr output_poly, cado_poly_srcptr input_poly,
                 ropt_param_ptr param, ropt_time_ptr eacht)
{
  ropt_poly poly;
  ropt_poly_init (poly);

  ASSERT_ALWAYS(input_poly->nb_polys == 2);
  /* setup poly */

  cado_poly_set(poly->cpoly, input_poly);

  ropt_poly_setup (poly);

  ropt_info info;
  ropt_info_init (info);


  ropt_bestpoly bestpoly;
  ropt_bestpoly_init (bestpoly);
  ropt_bestpoly_set (bestpoly, poly->cpoly->pols[1], poly->cpoly->pols[0]);

  /* call main function */
  ropt_do_both_stages (poly, bestpoly, param, info);
  
  /* bring bestpoly back to polyselect_ropt */
  for( ; output_poly->nb_polys < 2 ; )
      cado_poly_provision_new_poly(output_poly);
  mpz_poly_set(output_poly->pols[0], bestpoly->g);
  mpz_poly_set(output_poly->pols[1], bestpoly->f);
  mpz_set (output_poly->n, input_poly->n);

  /* get time passed from info, use info to keep interface unchanged */
  eacht->ropt_time_stage1 = info->ropt_time_stage1;
  eacht->ropt_time_tuning = info->ropt_time_tuning;
  eacht->ropt_time_stage2 = info->ropt_time_stage2;
  
  /* free */
  ropt_bestpoly_clear (bestpoly);
  ropt_info_clear (info);
  ropt_poly_clear (poly);
}
