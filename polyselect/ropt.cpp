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

#include <cstdio>
#include <cstdlib>

#include <gmp.h>

#include "cado_poly.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "ropt.hpp"
#include "ropt_arith.hpp"
#include "ropt_linear.hpp"
#include "ropt_quadratic.hpp"
#include "ropt_stage2.hpp"
#include "ropt_str.hpp"
#include "ropt_tree.h"
#include "ropt_param.h"
#include "size_optimization.hpp"
#include "auxiliary.hpp"
#include "polyselect_alpha.h"


/**
 * Find best poly. This is somehow redundant.
 */
void
ropt_get_bestpoly ( ropt_poly const & poly,
                    MurphyE_pq *global_E_pqueue,
                    ropt_bestpoly & bestpoly
                    )
{
  double ave_MurphyE = 0.0, best_E = 0.0;

  /* var for computing E */
  cxx_mpz_poly Fuv, Guv;

  /* output all polys in the global queue */
  /* WTF starting with 1? Why? XXX */
  for (int i = 1; i < global_E_pqueue->used; i ++) {
      Fuv = poly.cpoly[1];
      Guv = poly.cpoly[0];

      /* rotation is by u*x+v+w*x^2 (yes, interesting naming...)
       * Note that u, v, w have different types, and interestingly we
       * have different functions to do this...
       */
      compute_fuv_mp (Fuv, poly.cpoly[1], poly.cpoly[0],
              global_E_pqueue->u[i], global_E_pqueue->v[i]);

      if (global_E_pqueue->w[i])
          rotate_aux (Fuv, Guv, 0, global_E_pqueue->w[i], 2);

      /* This modifies Fuv and Guv */
      sopt_local_descent (Fuv, Guv, Fuv, Guv, 1, -1, SOPT_DEFAULT_MAX_STEPS, 0);

      ave_MurphyE = print_poly_fg (Fuv, Guv, poly.cpoly.n, 0);

      /* Does it even make sense ? After all, if we end up with something
       * with non-trivial content, we would expect that we encountered it
       * earlier, right ?
       */
      mpz_poly_divide_by_content (Fuv);

      if (ave_MurphyE > best_E) {
          best_E = ave_MurphyE;
          bestpoly.f = Fuv;
          bestpoly.g = Guv;
      }
  }
}


/**
 * Root sieve only.
 */
static void
ropt_do_stage2 (ropt_poly & poly,
                ropt_bestpoly & bestpoly,
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

  cxx_mpz_poly Fuv, Guv;

  ropt_bound bound;
  ropt_bound_init (bound);

  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  Fuv = poly.cpoly[1];
  Guv = poly.cpoly[0];

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
                poly.alpha_proj,
                bound->exp_min_alpha );

  /* root sieve */
  ropt_s2param s2param;

  ropt_s2param_setup_stage2_only (bound, s2param, param,
                                  param->s2_u, param->s2_v, param->s2_mod);
  info->w = param->s2_w;
  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, param->s2_w);

  /* return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  ropt_bound_clear (bound);
}


/**
 * Ropt linear or quadratic.
 */
static void
ropt_do_both_stages ( ropt_poly & poly,
        ropt_bestpoly & bestpoly,
        ropt_param_ptr param,
        ropt_info_ptr info)
{
     int d = mpz_poly_degree(poly.cpoly[1]);
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
ropt ( ropt_poly & poly,
       ropt_bestpoly & bestpoly,
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
ropt_polyselect (cxx_cado_poly & output_poly, cxx_cado_poly const & input_poly,
                 ropt_param_ptr param, ropt_time_ptr eacht)
{
  ropt_poly poly;

  ASSERT_ALWAYS(input_poly.nsides() == 2);
  /* setup poly */

  poly.cpoly = input_poly;

  ropt_poly_setup (poly);

  ropt_info info;
  ropt_info_init (info);


  ropt_bestpoly bestpoly;
  ropt_bestpoly_set (bestpoly, poly.cpoly[1], poly.cpoly[0]);

  /* call main function */
  ropt_do_both_stages (poly, bestpoly, param, info);
  
  /* bring bestpoly back to polyselect_ropt */
  for( ; output_poly.nsides() < 2 ; )
      output_poly.provision_new_poly();
  output_poly[0] = bestpoly.g;
  output_poly[1] = bestpoly.f;
  output_poly.n = input_poly.n;
  /* TODO: skewness!! */

  /* get time passed from info, use info to keep interface unchanged */
  eacht->ropt_time_stage1 = info->ropt_time_stage1;
  eacht->ropt_time_tuning = info->ropt_time_tuning;
  eacht->ropt_time_stage2 = info->ropt_time_stage2;
  
  /* free */
  ropt_info_clear (info);
}
