/**
 * @file ropt_quadratic.c
 * Main function for quadratic rotation.
 * Called by ropt.c and then call stage1.c and stage2.c.
 *
 * This is description of the flow (for deg-6 polys):
 * - first ropt_quadratic_tune_stage1() tunes quad rot w and lognorm_incr;
 * -  then ropt_stage1() gets N lats based only on alpha;
 * -  then ropt_quadratic_tune() gets N' << N good lats based on crude test-sieve;
 * -  then ropt_call_sieve() the top few candidates from above.
 *
 * For the above figures:
 *  N   = size_total_sublattices[][1]*ropteffort*TUNE_NUM_SUBLATTICE_STAGE1
 *  N'  = size_total_sublattices[][1]*ropteffort (plus those from stage1_tune)
 *  N'' = size_total_sublattices[][2] (no ropteffort here)
 */


#include "cado.h"
#include "ropt_quadratic.h"
#include "ropt.h"
#include "portability.h"


/**
 * Tune stage 1 for quadratic rotation.
 * - tune qudratic rotation w
 *  -- tune lognorm_incr
 *  -- call ropt_stage1()
 */
static void
ropt_quadratic_tune_stage1 ( ropt_poly_t poly,
                             ropt_bound_t bound,
                             ropt_s1param_t s1param,
                             ropt_param_t param,
                             alpha_pq *alpha_pqueue
#if TUNE_LOGNORM_INCR
                             , alpha_pq *tune_E_pqueue,
                             ropt_info_t info,
                             MurphyE_pq *global_E_pqueue
#endif
  )
{
  const int numw = 3;
  int i, j, k, r, w, old_verbose, old_i, used, *w_good, w_top[numw];
  unsigned int size_alpha_pqueue_all_w =
    s1param->nbest_sl * TUNE_RATIO_STAGE1_FULL_ALPHA * param->effort;
  /* new_alpha_pq requires len >= 2 */
  if (size_alpha_pqueue_all_w < 2)
    size_alpha_pqueue_all_w = 2;
  double score;
#if TUNE_LOGNORM_INCR
  double incr;
#endif
  mpz_t u, v, mod;
  alpha_pq *alpha_pqueue_all_w;
  
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  w_good = (int *) malloc (size_alpha_pqueue_all_w * sizeof(int));

  /* tmp queue for tuning stage 1 */
  new_alpha_pq (&alpha_pqueue_all_w, size_alpha_pqueue_all_w);

  /* parameters for tuning quadratic rotation */
  old_verbose = param->verbose;
  param->verbose = 0;

  /* used in ropt_stage1(), the larger the slower for each
     quadratic rotation w. */
  s1param->nbest_sl = 16;

  /* used in ropt_stage1(), uses this to reduce the length of 
     individual sublattice queue */
  s1param->nbest_sl_tunemode = 1;

  /* Step 1: for each quadratic rotation, find sublattices in a
     quick, but approximated manner */
  k = 0;
  old_i = 0;
  for (i = bound->global_w_boundl; i <= bound->global_w_boundr; i++) {
    if (old_verbose >= 2 && k % 10 == 0)
      fprintf (stderr, "# Info: quadratic rotation range %d*x^2\n", i);
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, i, 2);
    ropt_poly_setup (poly);
    r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue_all_w, i);
    k ++;
  }

  /* get/rotate back */
  param->verbose = old_verbose;
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* fetch all w */
  used = alpha_pqueue_all_w->used - 1;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue_all_w, &w, u, v, mod, &score);
    w_good[i] = w;

#if RANK_SUBLATTICE_BY_E
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: got %4d lat (%d, %Zd, %Zd) "
                   "(mod %Zd), E %.2f\n", i + 1, w, u, v, mod, score);
    }
#else
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: got %4d lat (%d, %Zd, %Zd) "
                   "(mod %Zd), alpha: %.2f\n", i + 1, w, u, v, mod, score);
    }
#endif
  }

  /* Step 2: fetch top three w */
  for (i = 0; i < numw; i ++)
    w_top[i] = bound->global_w_boundr + 1;
  k = 0;
  for (i = used - 1; i > 0; i --) {
    if (k > 2)
      break;    
    r = 0;
    for (j = 0; j < numw; j ++) {
      if (w_top[j] == w_good[i])
        r = 1;
    }
    if (r != 1)
      w_top[k++] = w_good[i];
  }

  
  /* Step 3: for each top w, re-detect (u, v) harder */

  /* consider more sublattices in stage 1 and leave them for tuning */
  old_i = 0;
  for (i = 0; i < numw; i ++) {
    w = w_top[i];
    /* w_top[3] is not full? */
    if (w >= bound->global_w_boundr + 1) continue;
    if (param->verbose >= 2)
      fprintf (stderr, "# Info: find lat on quadratic rotation "
               "by %d*x^2\n", w);

    /* quadratic rotation and setup */
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, w, 2);
    ropt_poly_setup (poly);

    /* tune incr */
#if TUNE_LOGNORM_INCR
    incr = ropt_linear_tune_stage1 (poly, s1param, param, tune_E_pqueue, alpha_pqueue,
                             info, global_E_pqueue, w);
    ropt_bound_setup_incr (poly, bound, param, incr);
    ropt_s1param_setup (poly, s1param, bound, param);
#endif

    /* call ropt_stage1() to find good sublattices */
    s1param->nbest_sl_tunemode = 0; // Note: used in ropt_stage1()
    r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, w);
    remove_rep_alpha (alpha_pqueue);
  }
  if (param->verbose >= 2)
    fprintf (stderr, "\n");

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* clear */
  free_alpha_pq (&alpha_pqueue_all_w);
  free (w_good);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
}


/**
 * Quadratic rotation.
 */
void
ropt_quadratic ( ropt_poly_t poly,
                 ropt_bestpoly_t bestpoly,
                 ropt_param_t param,
                 ropt_info_t info )
{
  double t1, t2, t3;
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  alpha_pq *alpha_pqueue;
#if TUNE_LOGNORM_INCR
  alpha_pq *tune_E_pqueue;
#endif
  MurphyE_pq *global_E_pqueue;

  /* setup bound, s1param, alpha_pqueue, tsieve_E_pqueue */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);
#if TUNE_LOGNORM_INCR
  new_alpha_pq (&tune_E_pqueue, s1param->nbest_sl);
#endif  
  /* here we use some larger value since alpha_pqueue records
     more sublattices to be tuned */
  unsigned long len_full_alpha = s1param->nbest_sl *
    TUNE_RATIO_STAGE1_FULL_ALPHA * param->effort;
  /* new_alpha_pq requires len_full_alpha >= 2 */
  if (len_full_alpha < 2)
    len_full_alpha = 2;
  new_alpha_pq (&alpha_pqueue, len_full_alpha);

  /* [Step 1] find w first and do ropt_stage1() there */
  t1 = seconds_thread ();
#if TUNE_LOGNORM_INCR
  ropt_quadratic_tune_stage1 (poly, bound, s1param, param, alpha_pqueue,
                         tune_E_pqueue, info, global_E_pqueue);
#else
  ropt_quadratic_tune_stage1 (poly, bound, s1param, param, alpha_pqueue);
#endif  
  t1 = seconds_thread () - t1;

  /* [Step 2] rank/tune above found sublattices by short sieving */
  t2 = seconds_thread ();
#if TUNE_LOGNORM_INCR
  ropt_tune_stage2 (poly, bound, s1param, param, info, alpha_pqueue,
                    tune_E_pqueue, global_E_pqueue);
#else
  ropt_tune_stage2 (poly, bound, s1param, param, info, alpha_pqueue,
                    global_E_pqueue);
#endif
  t2 = seconds_thread () - t2;
  
  /* [Step 3] final root sieve */
  t3 = seconds_thread ();
  ropt_call_sieve (poly, bound, s1param, param, info, alpha_pqueue,
                   global_E_pqueue);
  t3 = seconds_thread () - t3;
  
  info->ropt_time_stage1 = t1;
  info->ropt_time_tuning = t2;
  info->ropt_time_stage2 = t3;
  
  if (param->verbose >= 2) {
    fprintf ( stderr, "# Stat: this polynomial (stage 1) took %.2fs\n", t1 );
    fprintf ( stderr, "# Stat: this polynomial (tuning ) took %.2fs\n", t2 );
    fprintf ( stderr, "# Stat: this polynomial (stage 2) took %.2fs\n", t3 );
  }

  /* Step 4, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);
  
  /* free */
  free_MurphyE_pq (&global_E_pqueue);
#if TUNE_LOGNORM_INCR
  free_alpha_pq (&tune_E_pqueue);
#endif
  free_alpha_pq (&alpha_pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
}
