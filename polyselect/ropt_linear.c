/**
 * @file ropt_linear.c
 * Main function for linear rotation.
 * Called by ropt.c and call stage1.c and stage2.c.
 *
 * This is description of the flow (for deg-5 polys):
 * - first ropt_linear_tune_stage1() tunes lognorm_incr and hence bounds;
 * - then  ropt_stage1() gets N lats based only on alpha;
 * - then  ropt_tune_stage2_fast() gets N' << N good lats based on
       crude test-sieve;
 * - then  ropt_tune_stage2_slow() gets N'' ~ N lats based on finer
       and recursive test-sieve;
 * - then  ropt_call_sieve() the top few candidates from above.
 *
 * For the above figures:
 *  N   = size_total_sublattices[][1]*ropteffort*TUNE_RATIO_STAGE1_FULL_ALPHA
 *  N'  = size_total_sublattices[][1]*ropteffort (plus those from stage1_tune)
 *  N'' = size_total_sublattices[][2] (no ropteffort here)
 */ 


#include "cado.h"
#include "ropt_linear.h"
#include "portability.h"
//#define ROPT_LINEAR_TUNE_HARDER


#if TUNE_LOGNORM_INCR
/**
 * Tune the parameters for ropt_stage1().
 *  note even ropt_stage1() only computes sublattice, it depends on
 *  the parameters in `bound'. Here we use test-sieve to tune good
 *  parameters for `bound'. Here, the parameter is
 *   BOUND_LOGNORM_INCR_MAX which affects bound.
 */
double
ropt_linear_tune_stage1 ( ropt_poly_t poly,
                          ropt_s1param_t s1param,
                          ropt_param_t param,
                          alpha_pq *tune_E_pqueue,
                          alpha_pq *alpha_pqueue,
                          ropt_info_t info,
                          MurphyE_pq *global_E_pqueue,
                          unsigned long quad )
{
  info->mode = 0;

  /* setup bound, tmporary pqueue */
  unsigned int s1_size, pqueue_size;
  int i, k, kk, r, w, used, steps = 0;
  ropt_bound_t bound;
  alpha_pq *pqueue;
  ropt_s1param_t s1param_tune;
  mpz_t u, v, mod;

  s1_size = s1param->nbest_sl_tune * param->effort;
  if (s1_size < 2)
    s1_size = 2; /* required by new_alpha_pq() */
  pqueue_size = s1_size;

  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  ropt_bound_init (bound);
  new_alpha_pq (&pqueue, pqueue_size);
  ropt_s1param_init (s1param_tune);

  double score, maxscore,
    start_incr = BOUND_LOGNORM_INCR_MAX,
    best_incr = BOUND_LOGNORM_INCR_MAX,
    upper_incr = start_incr + BOUND_LOGNORM_INCR_MAX_TUNESTEP*param->effort,
    incr_step;
  maxscore = 0.0;
  unsigned long start_modbound, upper_modbound;
  
  /* compute initial bound */
  ropt_bound_setup_incr (poly, bound, param, start_incr);
  ropt_s1param_resetup (poly, s1param_tune, bound, param, s1_size);
  start_modbound = s1param_tune->modbound;
  
  /* compute upper bound modulus for stage1 */
  ropt_bound_setup_incr (poly, bound, param, upper_incr);
  ropt_s1param_resetup (poly, s1param_tune, bound, param, s1_size);
  upper_modbound = s1param_tune->modbound;

  if (param->verbose >= 1) {
    gmp_fprintf (stderr, "# Info: start_modbound : %lu\n",
                 start_modbound);
    gmp_fprintf (stderr, "# Info: end_modbound : %lu\n",
                 upper_modbound);
    gmp_fprintf (stderr, "# Info: tune (incr), modbound: %lu to %lu\n",
                 start_modbound, upper_modbound);
  }

  for (k = 0; k < NUM_DEFAULT_SUBLATTICE; k ++)
    if (default_sublattice_prod[k] == start_modbound)
      break;
  kk = k;
  while (default_sublattice_prod[k] <= upper_modbound) {
    k++;
    steps++;
  }
  incr_step = (upper_incr-start_incr)/steps;

  /* compute incr step */
  while (default_sublattice_prod[kk] <= upper_modbound) {

    /* change bound and test-sieve */
    start_modbound = default_sublattice_prod[kk];
    ropt_bound_setup_incr (poly, bound, param, start_incr);
    ropt_s1param_resetup_modbound (poly, s1param_tune, bound, param,
                                   s1_size, start_modbound); 
    /* num. lat (for each p) bounded by s1_size_each_sublattice_tune[] resp. */
    s1param_tune->nbest_sl_tunemode = 1;
    r = ropt_stage1 (poly, bound, s1param_tune, param, pqueue, quad);
    s1param_tune->nbest_sl_tunemode = 0;
    if (r == -1)
      break;

    /* better not wasting sublattice */
    used = pqueue->used;
    for (i = 1; i < used; i ++) {
      /*
      gmp_fprintf (stderr, "# Info: inserting %.2f on "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   pqueue->alpha[i], pqueue->w[i], pqueue->u[i], pqueue->v[i],
                   pqueue->modulus[i]);
      */
      insert_alpha_pq (alpha_pqueue, pqueue->w[i], pqueue->u[i], pqueue->v[i],
                       pqueue->modulus[i], pqueue->alpha[i]);
    }
    
    /* test root sieve, pqueue changed to have MurphyE scores */
    score = ropt_tune_stage2_fast (poly, s1param_tune, param, info, pqueue,
                                   global_E_pqueue, size_tune_sievearray);
    /* update maxscore */
    if (maxscore < score) {
      maxscore = score;
      best_incr = start_incr;
    }

   /* also better not wasting sieved slots, save to tune_E_pqueue */
    used = pqueue->used - 1;
    for (i = 0; i < used; i ++) {
      extract_alpha_pq (pqueue, &w, u, v, mod, &score);
      insert_alpha_pq (tune_E_pqueue, w, u, v, mod, score);
    }

    /* reset queue */
    reset_alpha_pq (pqueue);
    start_incr += incr_step;
    kk++;
  }
  
  if (param->verbose >= 1) {
    gmp_fprintf (stderr, "# Info: tune (incr), best_lognorm_incr: %f\n",
                 best_incr);
  }

  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  free_alpha_pq (&pqueue);
  ropt_s1param_free (s1param_tune);
  ropt_bound_free (bound);

  info->mode = 0;
  return best_incr;
}
#endif


/**
 * Aux function that copy all E_pqueue to alpha_pqueue
 *  -- use it with caution as it assumes alpha_pqueue now have MurphyE scores
 */
void
ropt_MurphyE_to_alpha ( MurphyE_pq *E_pqueue,
                        alpha_pq *alpha_pqueue )
{
  int i, used = E_pqueue->used;
  for (i = 1; i < used; i ++) {
    /*
      gmp_fprintf (stderr, "# got %.2g"
      "(%d, %Zd, %Zd) (mod %Zd), %d\n",
      E_pqueue->E[i], E_pqueue->w[i], E_pqueue->u[i], E_pqueue->v[i],
      E_pqueue->modulus[i], alpha_pqueue->w[0]);
    */
    insert_alpha_pq (alpha_pqueue, E_pqueue->w[i],
                     E_pqueue->u[i], E_pqueue->v[i], 
                     E_pqueue->modulus[i], -E_pqueue->E[i]);
  }
}


/**
 * This is a fast tuning process which deal with lats from ropt_stage1()
 */
double
ropt_tune_stage2_fast ( ropt_poly_t poly,
                        ropt_s1param_t s1param,
                        ropt_param_t param,
                        ropt_info_t info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune )
{
  /* tune mode, prevent from fprint */
  info->mode = 1;

  int i, w, used, old_i;
  double score;
  double max_score;
  mpz_t u, v, mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue;
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  ropt_s2param_init (poly, s2param);
  new_alpha_pq (&tmp_alpha_pqueue, s1param->nbest_sl);
#if RANK_SUBLATTICE_BY_E
    char stmp[] = "E";  /* RANK_SUBLATTICE_BY_E */
#else
    char stmp[] = "alpha";
#endif

  /* The following is done:
     - test sieve on alpha_pqueue (many) returned from ropt_stage1();
     - reduce the sublattice number to the s1param->nbest_sl;
     - record them to tmp_alpha_pqueue. */
  used = alpha_pqueue->used - 1;
  old_i = 0;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, w, 2);
    ropt_poly_setup (poly);
    //ropt_bound_reset (poly, bound, param); // not really necessary

    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: tune [%4d], %s: %.2f, lat (%d, %Zd, %Zd) "
                   "(mod %Zd)\n", i+1, stmp, score, w, u, v, mod);
    }

    ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                             0, curr_size_tune, NP-1);
    ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);
    insert_alpha_pq (tmp_alpha_pqueue, w, u, v, mod, -info->ave_MurphyE);
    if (param->verbose >= 4)
      gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   info->ave_MurphyE, info->best_MurphyE,
                   w, u, v, mod);

  }
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* save result to alpha_pqueue */
  reset_alpha_pq (alpha_pqueue);
  used =  tmp_alpha_pqueue->used - 1;
  max_score = 0;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);
    if (param->verbose >= 3)
      gmp_fprintf (stderr, "# Info:  got [%4d], E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   i+1, -score, w, u, v, mod);
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
    if (max_score < -score)
      max_score = score;
  }
  max_score = -max_score;

  /* also save best score so far to alpha_pqueue, E_queue 
     is smaller than alpha_queue. Do not seem good in tests */
  //ropt_MurphyE_to_alpha (global_E_pqueue, alpha_pqueue);
  
  /* free s2param */
  free_alpha_pq (&tmp_alpha_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);

  /* normal mode (only control stderr) */
  info->mode = 0;
  return max_score;
}


/* slow but finer tune */
void
ropt_tune_stage2_slow ( ropt_poly_t poly,
                        ropt_bound_t bound,
                        ropt_s1param_t s1param,
                        ropt_param_t param,
                        ropt_info_t info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune,
                        int tune_round,
                        unsigned int curr_nbest )
{
  /* tune mode */
  info->mode = 1;

  int i, j, w, used, old_i;
  double score, old_MurphyE;
  mpz_t u, tmpu, v, mod, old_mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue, *tmp2_alpha_pqueue;
#if TUNE_EARLY_ABORT
  double ave_E, ave_bestE, new_ave_E, new_ave_bestE;
  int acc;
#endif  
  
#ifdef ROPT_LINEAR_TUNE_HARDER
  int k;
#endif
  mpz_init (u);
  mpz_init (tmpu);
  mpz_init (v);
  mpz_init (mod);
  mpz_init (old_mod);
  ropt_s2param_init (poly, s2param);

  if (curr_nbest > (unsigned) alpha_pqueue->used-1)
    curr_nbest = alpha_pqueue->used-1;
  if (curr_nbest < 2) curr_nbest = 2;
  new_alpha_pq (&tmp_alpha_pqueue, curr_nbest);
  new_alpha_pq (&tmp2_alpha_pqueue, curr_nbest);

  /* Step 1: test sieve on alpha_pqueue
     - for each sublattice in s1param->nbest_sl
     - do an early abort based on best_E score */
  used = alpha_pqueue->used - 1;
  old_i = 0;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, w, 2);
    ropt_poly_setup (poly);
    // ropt_bound_reset (poly, bound, param); // not necessary
#if RANK_SUBLATTICE_BY_E
    char stmp[] = "E";
#else
    char stmp[] = "alpha";
#endif
    if (param->verbose >= 1) {
      /* first tune round, score are in E, not murphyE */
      if (tune_round==1) {
        gmp_fprintf (stderr, "# Info: tune [%4d], %s: %.2f, lat (%d, %Zd, %Zd) "
                     "(mod %Zd)\n", i+1,  stmp, score, w, u, v, mod);
      }
      else {
        gmp_fprintf (stderr, "# Info: tune [%4d], E: %.2e, "
                     "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                     -score, i + 1, w, u, v, mod);
      }
    }
    
    /* detect positive good u (may also detect good mod with more time) */
    mpz_set (tmpu, u);
    j = 0;
#if TUNE_EARLY_ABORT
    ave_E = 0.0;
    ave_bestE = 0.0;
    new_ave_E = 0.0;
    new_ave_bestE = 0.0;
    acc = 0;
#endif
    
    while (mpz_cmp_si(tmpu, bound->global_u_boundr) <= 0) {
      
      if (j > TUNE_BOUND_ON_UV_TRIALS)
        break;

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tuning process */
      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, curr_size_tune, NP-1);

        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "lat (%d, %Zd, %Zd) (mod %Zd) "
                       "(%d, %d, %d)\n",
                       info->ave_MurphyE, info->best_MurphyE,
                       w, tmpu, v, mod, i + 1, j, k);
        }

        mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + k]);

        k ++;
      }
      mpz_set (mod, old_mod);

#else

      ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                               0, curr_size_tune, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "lat (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }

#endif

      mpz_add (tmpu, tmpu, mod); /* consider original u itself */

#if TUNE_EARLY_ABORT
      new_ave_E = ave_E + info->ave_MurphyE;
      new_ave_bestE = ave_bestE + info->best_MurphyE;
      if (j>0) {
        if ((ave_E/j>new_ave_E/(j+1)) && (ave_bestE/j > new_ave_bestE/(j+1))) {
          acc ++;
          if (acc >= TUNE_EARLY_ABORT_THR)
            break;
        }
      }
      ave_E = new_ave_E;
      ave_bestE = new_ave_bestE;
#endif      
      j ++;
    }

    /* detect negative good u */
    mpz_set (tmpu, u);
    j = 0;
#if TUNE_EARLY_ABORT
    ave_E = 0.0;
    ave_bestE = 0.0;
    new_ave_E = 0.0;
    new_ave_bestE = 0.0;
    acc = 0;
#endif

    while (mpz_cmp_si (tmpu, bound->global_u_boundl) > 0) {

      if (j > TUNE_BOUND_ON_UV_TRIALS)
        break;
      mpz_sub (tmpu, tmpu, mod);

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tuning process */

      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, curr_size_tune, NP-1);
        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "lat (%d, %Zd, %Zd) (mod %Zd) "
                       "(%d, %d, %d)\n",
                       info->ave_MurphyE, info->best_MurphyE,
                       w, tmpu, v, mod, i + 1, j, k);
        }

        mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + k]);

        k ++;
      }
      mpz_set (mod, old_mod);

#else

      ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                               0, curr_size_tune, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "lat (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }

#endif

#if TUNE_EARLY_ABORT
      new_ave_E = ave_E + info->ave_MurphyE;
      new_ave_bestE = ave_bestE + info->best_MurphyE;
      if (j>0) {
        if ((ave_E/j>new_ave_E/(j+1)) && (ave_bestE/j > new_ave_bestE/(j+1))) {
          acc ++;
          if (acc >= TUNE_EARLY_ABORT_THR)
            break;
        }
      }
      ave_E = new_ave_E;
      ave_bestE = new_ave_bestE;
#endif

      j ++;

    }
  }

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* Step 2: slight larger range sieve on best sublattices */
  used =  tmp_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);

    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, w, 2);

    ropt_poly_setup (poly);

    if (param->verbose >= 1) {
      gmp_fprintf (stderr, "# Info: tune [%4d], E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   -score, i + 1, w, u, v, mod);
    }

    insert_alpha_pq (tmp2_alpha_pqueue, w, u, v, mod, score);

    j = 0;
    old_MurphyE = -score;
    while (j < TUNE_BOUND_ON_MOD_TRIALS) {

      mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + j]);

      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, curr_size_tune, 20);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      /* abort if E is becoming worse */
      if (old_MurphyE > info->ave_MurphyE)
        break;

      old_MurphyE = info->ave_MurphyE;
      insert_alpha_pq (tmp2_alpha_pqueue, w, u, v, mod, -info->ave_MurphyE);
      j ++;
    }
  }

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* End: save to alpha_pqueue */
  reset_alpha_pq (alpha_pqueue);
  used =  tmp2_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tmp2_alpha_pqueue, &w, u, v, mod, &score);
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
  }
  
  /* free s2param */
  free_alpha_pq (&tmp_alpha_pqueue);
  free_alpha_pq (&tmp2_alpha_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (u);
  mpz_clear (tmpu);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (old_mod);

  /* normal mode (only control stderr) */
  info->mode = 0;
}


/**
 * Rank found sublattices in alpha_pqueue by test root sieve.
 * Note: alpha_pqueue contains much more sublattices than 
 * the final sublattices to be root-sieved. In the function,
 * - the 1st pass of tuning reduces #alpha_pqueue to nbest_sl.
 * - the 2nd pass of tuning identifies good u and v.
 * - the 3rd pass of tuning identifies good mod.
 *
 * Note: The 1st is fast so we should do more; 2nd pass is slow;
 *       3rd pass is fast.
 */
void
ropt_tune_stage2 ( ropt_poly_t poly,
                   ropt_bound_t bound,
                   ropt_s1param_t s1param,
                   ropt_param_t param,
                   ropt_info_t info,
                   alpha_pq *alpha_pqueue,
#if TUNE_LOGNORM_INCR
                   alpha_pq *tune_E_pqueue,
#endif
                   MurphyE_pq *global_E_pqueue )
{
#if TUNE_LOGNORM_INCR
  int i, w, used;
  double score;
  mpz_t u, v, mod;
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
#endif
  unsigned int old_nbest = s1param->nbest_sl;
  unsigned int curr_size_tune = size_tune_sievearray;
  int r = 1;

  /* step 2: fast tune for lats from stage 1 */
  if (param->verbose >= 1) {
    printf ("# Info: tune (round %d), tunesize: %u, queuesize: %u to %u\n",
            r, curr_size_tune, alpha_pqueue->used-1, s1param->nbest_sl);
  }
  ropt_tune_stage2_fast (poly, s1param, param, info, alpha_pqueue,
                         global_E_pqueue, curr_size_tune);
  r ++;


  /* aux: in case tune_E_pqueue not empty, combine it */
#if TUNE_LOGNORM_INCR
  used = tune_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tune_E_pqueue, &w, u, v, mod, &score);
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: pre-tune [%4d], E: %.2e, lat (%d, %Zd, %Zd) "
                   "(mod %Zd)\n", i+1, -score, w, u, v, mod);
    }
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
   }
  remove_rep_alpha (alpha_pqueue);
#endif

  /* step 3: recursieve (slower but finer) tune */
  while (1) {
    if (param->verbose >= 1) {
      printf ("# Info: tune (round %d), tunesize: %u, nbest: %u\n",
              r, curr_size_tune, s1param->nbest_sl);
    }
    curr_size_tune = curr_size_tune*4;
    s1param->nbest_sl = s1param->nbest_sl/2;
    ropt_tune_stage2_slow (poly, bound, s1param, param, info, alpha_pqueue,
      global_E_pqueue, curr_size_tune, r, s1param->nbest_sl);
    r ++;
    remove_rep_alpha (alpha_pqueue);
    /* doing 2 rounds of tuning for now */
    if (s1param->nbest_sl < old_nbest/2)
      break;
  }
  s1param->nbest_sl = old_nbest;

#if TUNE_LOGNORM_INCR
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
#endif
  return;
}



/**
 * Call root sieve. It will root-sieve two passes of the following
 * two queues:
 * - alpha_pqueue which records 'average best E' sublattices.
 * - global_E_pqueue which records 'top E' sublattices.
 */
void
ropt_call_sieve ( ropt_poly_t poly,
                  ropt_bound_t bound,
                  ropt_s1param_t s1param,
                  ropt_param_t param,
                  ropt_info_t info,
                  alpha_pq *alpha_pqueue,
                  MurphyE_pq *global_E_pqueue )
{
  int i, w, used;
  double score;
  mpz_t u, v, mod;
  ropt_s2param_t s2param;
  MurphyE_pq *tmp_E_pqueue;
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);

  new_MurphyE_pq (&tmp_E_pqueue, s1param->nbest_sieve);
  ropt_s2param_init (poly, s2param);

  /* remove repetition in global_E_pqueue */
  remove_rep_MurphyE (global_E_pqueue);

  /* remove repetition in alpha_pqueue */
  remove_rep_alpha (alpha_pqueue);

  /* alpha_pqueue contains best ave. E sublattices.
     global_E_pqueue contains best E polynomials. We
     form a final queue using both, with a focus on
     the second. */

  /* Step 1, extract global_E_pqueue --> tmp_E_pqueue */
  used = global_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_MurphyE_pq (global_E_pqueue, &w, u, v, mod, &score);
    insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, score);
    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: best      E: %.2e on [%4d], "
                   "(%d, %Zd, %Zd) (mod %Zd)\n",
                   score, i + 1, w, u, v, mod);
    }
  }

  /* re-insert global_E_pqueue as it may contains some top polynomials found
     in sieving test, which may not be found again due to size property. */
  used = tmp_E_pqueue->used - 1;
  for (i = 1; i < used; i ++) {
    insert_MurphyE_pq (global_E_pqueue, tmp_E_pqueue->w[i],
                       tmp_E_pqueue->u[i], tmp_E_pqueue->v[i], 
                       tmp_E_pqueue->modulus[i], tmp_E_pqueue->E[i]);
  }

  /* Step 2, extract alpha_pqueue */
  used = alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);
    /* score here is negative*/
    insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, -score);
    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: best ave. E: %.2e on [%4d], "
                   "(%d, %Zd, %Zd) (mod %Zd)\n",
                   -score, i + 1, w, u, v, mod);
    }
  }

  /* Step3, final root sieve */
  int old_i = 0;
  used = tmp_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_MurphyE_pq (tmp_E_pqueue, &w, u, v, mod,
                        &score);

    /* rotate */
    old_i = rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, w, 2);
    ropt_poly_setup (poly);

    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: sieve [%4d], E: %.2e, "
                   "(%d, %Zd, %Zd) (mod %Zd)\n",
                   i + 1, score, w, u, v, mod);
    }

    ropt_s2param_setup (bound, s1param, s2param, param, u, v, mod);

    ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

    /* extra sieve for small skewness */
    if (score > info->best_MurphyE) {
      info->mode = 0; // normal mode

      if (param->verbose >= 3) {
        fprintf (stderr, "# Info: resieve [%4d] in "
                 "range %u.\n", i + 1, size_tune_sievearray * 10);
      }
     
      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, size_tune_sievearray * 10, NP - 1);
      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);
    }
  }
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], poly->g[0], old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* free */
  free_MurphyE_pq (&tmp_E_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
}


/**
 * Linear rotation: deg 5.
 */
void
ropt_linear_deg5 ( ropt_poly_t poly,
                   ropt_bestpoly_t bestpoly,
                   ropt_param_t param,
                   ropt_info_t info )
{

  /* setup bound, s1param, alpha_pqueue, global_E_pqueue */
  int r;
  double t1, t2, t3;
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  alpha_pq *alpha_pqueue;
#if TUNE_LOGNORM_INCR
  alpha_pq *tune_E_pqueue;
#endif
  MurphyE_pq *global_E_pqueue;
#if TUNE_LOGNORM_INCR
  double incr;
#endif
  
  /* setup bound, s1param, alpha_pqueue, tsieve_E_pqueue */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl); /* always used */
#if TUNE_LOGNORM_INCR
  new_alpha_pq (&tune_E_pqueue, s1param->nbest_sl);
#endif  
  /* Here we use some larger value since alpha_pqueue records
     more sublattices to be tuned by test root sieve */
  unsigned long len_full_alpha = s1param->nbest_sl *
    TUNE_RATIO_STAGE1_FULL_ALPHA * param->effort;
  if (len_full_alpha < 2)
    len_full_alpha = 2; /* required by new_alpha_pq */
  new_alpha_pq (&alpha_pqueue, len_full_alpha);

  /* [Step 1] tune ropt_bound first */
  t2 = seconds_thread ();
#if TUNE_LOGNORM_INCR
  if (param->s1_num_e_sl==0) {
    incr = ropt_linear_tune_stage1 (poly, s1param, param, tune_E_pqueue,
                                    alpha_pqueue, info, global_E_pqueue, 0);
    ropt_bound_setup_incr (poly, bound, param, incr);
    ropt_s1param_resetup (poly, s1param, bound, param, s1param->nbest_sl);
  }
#endif
  t2 = seconds_thread () - t2;

  
  /* [Step 2] ropt_stage1() find good sublattices */
  t1 = seconds_thread ();
  r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, 0);
  remove_rep_alpha (alpha_pqueue);
  t1 = seconds_thread () - t1;
  if (r == -1) return;
  
  /* [Step 3] rank/tune above found sublattices by short sieving */
  t3 = seconds_thread ();
#if TUNE_LOGNORM_INCR
  ropt_tune_stage2 (poly, bound, s1param, param, info, alpha_pqueue,
                    tune_E_pqueue, global_E_pqueue);
#else
  ropt_tune_stage2 (poly, bound, s1param, param, info, alpha_pqueue,
                    global_E_pqueue);
#endif
  t2 = t2 + (seconds_thread () - t3);

  /* [Step 4] final root sieve */
  t3 = seconds_thread ();
  ropt_call_sieve (poly, bound, s1param, param, info, alpha_pqueue,
                     global_E_pqueue);
  t3 = seconds_thread () - t3;

  /* [Step 5] return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  info->ropt_time_stage1 = t1;
  info->ropt_time_tuning = t2;
  info->ropt_time_stage2 = t3;
  if (param->verbose >= 1) {
    fprintf ( stderr, "# Stat: this polynomial (stage 1) took %.2fs\n", t1 );
    fprintf ( stderr, "# Stat: this polynomial (tuning ) took %.2fs\n", t2 );
    fprintf ( stderr, "# Stat: this polynomial (stage 2) took %.2fs\n", t3 );
  }

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
#if TUNE_LOGNORM_INCR
  free_alpha_pq (&tune_E_pqueue);
#endif
  free_alpha_pq (&alpha_pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
}


/**
 * Linear rotation: deg 4.
 * No tuning at all.
 */
void
ropt_linear_deg34 ( ropt_poly_t poly,
                    ropt_bestpoly_t bestpoly,
                    ropt_param_t param,
                    ropt_info_t info )
{
  unsigned long ub, vb;
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  ropt_s2param_t s2param;
  MurphyE_pq *global_E_pqueue;
  mpz_t u, v, mod;

  mpz_init_set_ui (u, 0);
  mpz_init_set_ui (v, 0);
  mpz_init_set_ui (mod, 1);
  
  /* setup bound, s1param */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);
  ropt_s2param_init (poly, s2param);
  
  /* Step 1, set up lattice (0, 0) (mod 1) and run root sieve */
  ub = (bound->global_u_boundr > 0) ?
    bound->global_u_boundr : (-bound->global_u_boundr);
  vb = mpz_get_ui (bound->global_v_boundr);
  if (vb < 128) vb = 128;
  
  ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                           ub, vb, NP - 1);

  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, 0);

  /* Step 2, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);  
  free_MurphyE_pq (&global_E_pqueue);
  ropt_s1param_free (s1param);
  ropt_s2param_free (poly, s2param);
  ropt_bound_free (bound);
}


/**
 * Linear rotation: for deg 4 and deg 5
 */
void
ropt_linear ( ropt_poly_t poly,
              ropt_bestpoly_t bestpoly,
              ropt_param_t param,
              ropt_info_t info )
{
  if (poly->d == 3)
    ropt_linear_deg34 (poly, bestpoly, param, info);
  else if (poly->d == 5 || poly->d == 4)
    ropt_linear_deg5 (poly, bestpoly, param, info);
  else 
    fprintf (stderr, "Error: ropt_linear() only supports degrees 3-5.");
}
