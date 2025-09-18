/**
 * @file ropt_stage1.c
 * Called by ropt.c to find congruence classes.
 */


#include "cado.h" // IWYU pragma: keep

#include <stdio.h>      // fprintf stderr
#include <stdlib.h> // free malloc exit
#include <limits.h>     // LONG_MAX


#include <gmp.h>

#include "auxiliary.h"  // get_alpha
#include "cado_poly.h"
#include "macros.h"
#include "mpz_poly.h"
#include "polyselect_norms.h"
#include "ropt_arith.h" // reduce_poly_ul ...
#include "ropt_param.h" // TUNE_RATIO_STAGE1_FULL_ALPHA ...
#include "ropt_single_sublattice_priority_queue.h"
#include "ropt_stage1.h"
#include "ropt_str.h" // ropt_poly_t
#include "ropt_sublattice_crt.h"
#include "ropt_sublattice_priority_queue.h"
#include "ropt_tree.h" // node ...
#include "timing.h"             // for milliseconds


/**
 * Find good sublattice, the liftings.
 */
static inline void
find_sublattice_lift ( node *firstchild,
                       single_sublattice_priority_queue_ptr top,
                       unsigned int * f_ui,
                       unsigned int * g_ui,
                       unsigned int * fuv_ui,
                       int d,
                       unsigned int p,
                       char e,
                       char curr_e )
{
  /* recursion end */
  if (firstchild == NULL || curr_e > e)
    return;

  char l;
  unsigned int i, j, k, nroots, pe, pem1, fr, gr;
  node *currnode, *tmpnode = NULL, *tmpnode2 = NULL;

  /* compute p^e */
  pem1 = 1;
  for (l = 0; l < curr_e - 1; l ++)
    pem1 = pem1 * p;
  pe = pem1 * p;

  /* loop until all siblings are checked. */
  currnode = firstchild;
  while (currnode != NULL) {

    /*
      printf("-----\n");
      printf("p: %u, e: %d -> %d\n", p, curr_e - 1, curr_e);
      printf("(u, v) pair: (%u, %u) has roots: \n",
      currnode->u, currnode->v);
      for (i = 0; i < (currnode->nr); i++)
      printf("\tr[%u]: %u\n", i, currnode->r[i]);
      printf("-----\n");
    */

    /* compute f_uv(x) and then evaluate it at r. */
    compute_fuv_ui (fuv_ui, f_ui, g_ui, d, currnode->u, currnode->v, pe);

    /* loop all roots */
    for (nroots = 0; nroots < currnode->nr; nroots++) {

      gr = eval_poly_ui_mod (g_ui, 1, currnode->r[nroots], pe);
      if (gr % p == 0)
        continue;

      /* If the root is multiple */
      if (currnode->roottype[nroots] == 2) {

        fr = eval_poly_ui_mod (fuv_ui, d, currnode->r[nroots], pe);
        fr = fr / pem1;

        /* solve on fr + gr*A = 0 (mod p), where A = i*r + j. */
        fr = (unsigned int) solve_lineq (fr, gr, 0, p);

        /* we want to solve (i, j) in  fr = i*r + j (mod p).
           - if r is not invertible, loop i with fixed j;
           - otherwise, loop j to solve i. */
        if (currnode->r[nroots] % p == 0) {

          for (i = 0; i < p; i ++) {

#if DEBUG_FIND_SUBLATTICE
            fprintf (stderr, "fr: %u, r: %u,  (i, j): (%u, %u) "
                     "-> (%u, %u) (non-invertible, multiple)\n",
                     fr, currnode->r[nroots], i,
                     fr, currnode->u + pem1 * i,
                     currnode->v + pem1 * fr);
#endif
            /* r is a multiple root, add r + k * p^{e-1}. */
            for (k = 0; k < p; k ++) {
              insert_node ( currnode,
                            &tmpnode,
                            currnode->u + pem1 * i,
                            currnode->v + pem1 * fr,
                            currnode->r[nroots] + k * pem1,
                            curr_e, p, pe, 2 );
            }

            /* count the lifted single roots for any lifted
               pairs (u, v). Note the lifted single roots
               will not be computed actually. */
            for (k = 0; k < (currnode->nr); k++) {
              if (currnode->roottype[k] == 1) {
                insert_node ( currnode,
                              &tmpnode,
                              currnode->u + pem1 * i,
                              currnode->v + pem1 * fr,
                              currnode->r[k], curr_e, p, pe, 1 );
              }
            }
          }
        }
        else {
          for (j = 0; j < p; j ++) {

            /* given j, solve i in  fr = i*r + j (mod p). */
            i = solve_lineq (j, currnode->r[nroots], fr, p);

#if DEBUG_FIND_SUBLATTICE
            fprintf (stderr, "fr: %u, r: %u,  (i, j): (%u, %u) "
                     "-> (%u, %u) (invertible, multiple)\n",
                     fr, currnode->r[nroots], i,
                     j, currnode->u + pem1 * i,
                     currnode->v + pem1 * j);
#endif
            /* r is a multiple root, add r + k * p^{e-1}. */
            for (k = 0; k < p; k ++) {
              insert_node ( currnode,
                            &tmpnode,
                            currnode->u + pem1 * i,
                            currnode->v + pem1 * j,
                            currnode->r[nroots] + k * pem1,
                            curr_e, p, pe, 2 );
            }

            /* count the lifted single roots for any lifted
               pairs (u, v). Note the lifted single roots will
               not be computed actually. */
            for (k = 0; k < (currnode->nr); k++) {
              if (currnode->roottype[k] == 1) {
                insert_node ( currnode,
                              &tmpnode,
                              currnode->u + pem1 * i,
                              currnode->v + pem1 * j,
                              currnode->r[k], curr_e, p, pe, 1 );
              }
            }
          }
        }
      }
    }  // next root of current (u, v)

    /* recursieve to next level, curr_e + 1 */
    find_sublattice_lift ( currnode->firstchild,
                           top,
                           f_ui,
                           g_ui,
                           fuv_ui,
                           d,
                           p,
                           e,
                           curr_e + 1 );

    /* If current node is the 2nd bottom leave, add the bottom 
       level leaves with highest valuations to the list and
       delete them. */
    if (curr_e == e) {
      tmpnode = currnode->firstchild;
      while (tmpnode != NULL) {
        struct single_sublattice_info S[1];
        S->u = tmpnode->u;
        S->v = tmpnode->v;
        S->e = e;
        S->val = tmpnode->val;
        single_sublattice_priority_queue_push(top, S);
        tmpnode2 = tmpnode;
        tmpnode = tmpnode->nextsibling;

#if DEBUG_FIND_SUBLATTICE
        fprintf (stderr, "DEBUG_FIND_SUBLATTICE (2nd bottom): "
                 "p: %u, (%u, %u), val: %f, e: %d, max_e: %d\n",
                 p, tmpnode2->u, tmpnode2->v, tmpnode2->val, curr_e, e);
#endif

        free_node (&tmpnode2);
      }
    }

    /* delete current node and move to next sibling. */
    tmpnode = currnode;
    currnode = currnode->nextsibling;
    if (currnode != NULL)
      (currnode->parent)->firstchild = currnode;

#if DEBUG_FIND_SUBLATTICE
    fprintf (stderr, "DEBUG_FIND_SUBLATTICE (bottom): p: %u, "
             "(%u, %u), val: %f, e: %d, max_e: %d\n",
             p, tmpnode->u, tmpnode->v, tmpnode->val, curr_e, e);
#endif

    free_node (&tmpnode);
  } // next sibling of current node
}


/**
 * Find sublattices, the base case. Note, for the (u, v)
 * pairs, we consider all the simple + multi roots.
 */
static inline void
find_sublattice ( single_sublattice_priority_queue_ptr top,
                  ropt_poly_srcptr poly,
                  unsigned int p,
                  char e )
{
  unsigned int pe, r, u, v, fx_ui, gx_ui;
  char i;
  node *currnode, *root;
  mpz_t tmp;
  mpz_init (tmp);

  /* compute p^e */
  pe = 1;
  for (i = 0; i < e; i ++)
    pe = pe * p;

  /* new (u, v, val) tree */
  new_tree (&root);
  root = new_node ();

  /* for each root 0 <= r < p  */
  for (r = 0; r < p; r ++) {

    /* skip these */
    if (mpz_divisible_ui_p(poly->gx[r], p) != 0)
      continue;

    /* use single precision */
    fx_ui = (unsigned int) mpz_fdiv_ui (poly->fx[r], p);
    gx_ui = (unsigned int) mpz_fdiv_ui (poly->gx[r], p);

    for (u = 0; u < p; u ++) {

      /* u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
      mpz_mul (tmp, poly->gx[r], poly->gx[r]);
      mpz_mul_ui (tmp, tmp, u);
      mpz_sub (tmp, tmp, poly->numerator[r]);

      /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
      v =  compute_v_ui (fx_ui, gx_ui, r, u, p);

      /* simple root */
      if (mpz_divisible_ui_p(tmp, p) == 0) {
        insert_node (root, &currnode, u, v, r, 1, p, p, 1);
      }
      /* multiple root */
      else {
        insert_node (root, &currnode, u, v, r, 1, p, p, 2);
      }
    }
  }

  /* If e == 1, stop */
  if (e == 1) {
    node *tmpnode;
    currnode = root->firstchild;
    while (currnode != NULL) {
      struct single_sublattice_info S[1];
      S->u = currnode->u;
      S->v = currnode->v;
      S->val = currnode->val;
      S->e = 1;

      single_sublattice_priority_queue_push(top, S);
      tmpnode = currnode;
      currnode = currnode->nextsibling;
      free_node (&tmpnode);
    }
    tmpnode = NULL;
  }
  /* If e > 1, lift to higher p^e */
  else {
    /* find_sublattice_lift() only consider those pairs with at 
       least one multiple root, hence we need to consider those
       (u,v) which solely have single roots. */
    node *tmpnode = NULL, *lastnode = NULL;
    unsigned int j, c;
    currnode = root->firstchild;
    while (currnode != NULL) {

      c = 1;
      for (j = 0; j < currnode->nr; j ++)
        if (currnode->roottype[j] == 2)
          c = 0;

      /* case when (u, v) only has single roots */
      if (c == 1) {
        struct single_sublattice_info S[1];
        S->u = currnode->u;
        S->v = currnode->v;
        S->val = currnode->nr / (float) (p-1);
        S->e = 1;
        single_sublattice_priority_queue_push(top, S);

        /* delete this node */
        tmpnode = currnode;
        currnode = currnode->nextsibling;
        if (lastnode != NULL)
          lastnode->nextsibling = currnode;
        else
          (currnode->parent)->firstchild = currnode;
        free_node (&tmpnode);
      }
      else {
        lastnode = currnode;
        currnode = currnode->nextsibling;
      }
    }

    /* data struct for the lift */
    unsigned int *f_ui, *g_ui, *fuv_ui;
    int d = mpz_poly_degree(poly->cpoly->pols[1]);
    f_ui = (unsigned int*) malloc ( (d + 1) * sizeof (unsigned int) );
    fuv_ui = (unsigned int*) malloc ( (d + 1) * sizeof (unsigned int) );
    g_ui = (unsigned int*) malloc ( 2 * sizeof (unsigned int) );

    if ( (f_ui == NULL) ||
         (g_ui == NULL) ||
         (fuv_ui == NULL) ) {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit (1);
    }

    /* compute f (mod pe) */
    reduce_poly_uint (f_ui, poly->cpoly->pols[1], pe);
    reduce_poly_uint (g_ui, poly->cpoly->pols[0], pe);

    find_sublattice_lift ( root->firstchild,
                           top,
                           f_ui,
                           g_ui,
                           fuv_ui,
                           d,
                           p,
                           e,
                           2 );
    /* clear */
    free (f_ui);
    free (fuv_ui);
    free (g_ui);
  }

  mpz_clear (tmp);
  free_node(&root);
}


/**
 * Actual length check function.
 */
static inline void
return_combined_sublattice_check_tlen ( ropt_s1param_ptr s1param )
{
  unsigned int i;
  
  /* find the actual s1param->len_e_sl, excluding those 0's */
  s1param->tlen_e_sl = s1param->len_e_sl;
  for (i = 0; i < s1param->len_e_sl; i ++) {
    if (s1param->e_sl[i] == 0)
      s1param->tlen_e_sl --;
  }

  /* check */
  if (s1param->tlen_e_sl < 1) {
    fprintf ( stderr, "Error, number of primes in \"-e\" (len_e_sl) "
              "should be between 1 and 10\n" );
    exit(1);
  }
}

/**
 * Return all sublattices by calling CRT, where for each sublattice,
 * the separate (mod p) valuations are the best ones.
 */
static inline int
return_combined_sublattice ( ropt_poly_srcptr poly,
                             ropt_s1param_ptr s1param,
                             ropt_bound_srcptr bound,
                             sublattice_priority_queue_ptr pqueue,
                             int verbose )
{
  return_combined_sublattice_check_tlen (s1param);

  size_t t_nprimes = s1param->tlen_e_sl;

  unsigned int t_primes[t_nprimes], t_e_sl[t_nprimes];

  /* get t_primes and t_e_sl. These are subsets of ropt_primes and
   * s1param->e_sl (restricting to the positive e_sl values).
   *
   * XXX it is probably feasible to simplify the code path by adopting a
   * "full lattice" approach for e_sl[i]=0 below.
   *
   * XXX I think some of the s1param fields could/should be dropped.
   */
  {
      unsigned int j = 0;
      for (unsigned int i = 0; i < s1param->len_e_sl; i ++) {
          if (s1param->e_sl[i] != 0) {
              t_primes[j] = ropt_primes[i];
              t_e_sl[j++] = s1param->e_sl[i];
          }
      }
      /* j should equal the true length here */
      ASSERT_ALWAYS (j == (unsigned int) t_nprimes);
  }

  /* decide the number of top individual sublattices. Note that
     the s1param->tlen_e_sl must be already set */
  if (s1param->nbest_sl_tunemode == 0)
    ropt_s1param_setup_individual_nbest_sl (s1param);
  else
    ropt_s1param_setup_individual_nbest_sl_tune (s1param);


  single_sublattice_priority_queue * tops;

  tops = (single_sublattice_priority_queue *) malloc(t_nprimes * sizeof(single_sublattice_priority_queue));


  /* for each prime[i], lift the roots */
  for (unsigned i = 0; i < t_nprimes; i ++) {

    single_sublattice_priority_queue_init(tops[i], s1param->individual_nbest_sl[i]);

    /* find individual sublattices */
    find_sublattice (tops[i], poly, t_primes[i], t_e_sl[i]);

    /* if length zero, this set of parameters fails. */
    if (single_sublattice_priority_queue_empty(tops[i])) {
      for (unsigned int j = 0; j <= i ; j++) {
          single_sublattice_priority_queue_clear(tops[j]);
      }
      free(tops);
      return 0;
    }

    if (verbose >= 2)
      fprintf ( stderr,
                "# Info: p: %2u, max_e: %2u, list_size: %6zu\n",
                t_primes[i], t_e_sl[i],
                single_sublattice_priority_queue_size(tops[i]));


  }

  unsigned int count = ropt_sublattice_combine_all_crt(t_nprimes, t_primes, tops, bound, pqueue);

  /* info */
  if (verbose >= 2) {
    fprintf (stderr, "# Info: computed %8u CRTs\n", count);
  }

  for (unsigned int j = 0; j < t_nprimes ; j++) {
      single_sublattice_priority_queue_clear(tops[j]);
  }
  free(tops);

  return 1;
}


/**
 * Call return_all_sublattice() to return good sublattices.
 * Choose the "s1param->nbest_sl" best ones among them if there
 * are more quantities than it. The "best" property is ranked
 * by the size of u.
 */
static inline int
return_best_sublattice ( ropt_poly_srcptr poly,
                         ropt_s1param_ptr s1param,
                         ropt_bound_ptr bound,
                         sublattice_priority_queue_ptr pqueue,
                         int verbose )
{
  unsigned long tmp = bound->global_u_boundr;

  int ret = return_combined_sublattice ( poly,
                                         s1param,
                                         bound,
                                         pqueue,
                                         verbose );

  /* If failed, return. Some individual sublattices has length 0 */
  if (ret == -1) {
    return -1;
  }

  /* If no sublattice is found with u < bound->global_u_bound,
     then we try to enlarge u bound. */
  int count = 1;
  while (sublattice_priority_queue_empty(pqueue)) {
    if (bound->global_u_boundr < (LONG_MAX>>1)) {
      bound->global_u_boundr *= 2;
      if (verbose >= 2) {
        fprintf ( stderr,
                  "# Warn: not enough lat classes. "
                  "Reset \"bound->global_u_bound = %lu\" (#%d)\n",
                  bound->global_u_boundr, count );
      }
    }
    else
      return -1;

    return_combined_sublattice ( poly,
                                 s1param,
                                 bound,
                                 pqueue,
                                 verbose );
    count ++;
  }

  /* info */
  if (verbose >= 2) {
    fprintf ( stderr,
              "# Info: found    %8zu lat (\"size_total_sublattices\") "
              "where |u| < %lu\n",
              sublattice_priority_queue_size(pqueue), bound->global_u_boundr);
  }

  /* recover, for deg 6 poly */
  bound->global_u_boundr = tmp;

  return 1;
}

/* This is used by sublattice_priority_queue_do()
 */
struct transfer_to_alpha_priority_queue_arg {
    alpha_pq * alpha_pqueue;
    cado_poly_srcptr cpoly;
    int current_w;
};

void transfer_to_alpha_priority_queue(mpz_srcptr u, mpz_srcptr v, mpz_srcptr modulus, void * arg)
{
    struct transfer_to_alpha_priority_queue_arg * A = arg;

    /* fuv is f+(u*x+v)*g */
    mpz_poly Fuv;
    mpz_poly_init (Fuv, -1);

    compute_fuv_mp ( Fuv, A->cpoly->pols[1], A->cpoly->pols[0], u, v);

    double alpha_lat;

#if RANK_SUBLATTICE_BY_E
    /* use exp_E as benchmark instead of alpha. */
    alpha_lat = L2_skew_lognorm (Fuv);
    alpha_lat += get_alpha (Fuv, get_alpha_bound ());
#else
    //alpha_lat = get_alpha (fuv, poly->d, primes[s1param->tlen_e_sl-1]);
    alpha_lat = get_alpha (Fuv, get_alpha_bound ());
#endif

#if DEBUG_ROPT_STAGE1
    double logmu = L2_skew_lognorm (Fuv);
    gmp_fprintf ( stderr, "# Info: insert lat #%4d, (w, u, v): "
            "(%d, %Zd, %Zd) (mod %Zd), partial_alpha: %.2f,"
            "lognorm: %.2f\n",
            i,
            current_w,
            pqueue->u[i],
            pqueue->v[i],
            pqueue->modulus[i],
            alpha_lat,
            logmu );
#endif

    /* insert to a global priority queue */
    insert_alpha_pq (A->alpha_pqueue,
            A->current_w,
            u, v, modulus, alpha_lat);
    mpz_poly_clear (Fuv);
}

/**
 * Stage 1: record good sublattices to "alpha_pqueue".
 */
int
ropt_stage1 ( ropt_poly_srcptr poly,
              ropt_bound_ptr bound,
              ropt_s1param_ptr s1param,
              ropt_param_srcptr param,
              alpha_pq *alpha_pqueue,
              int current_w )
{
  int st = 0, re;
  sublattice_priority_queue pqueue;
  size_t len_part_alpha = s1param->nbest_sl *
    TUNE_RATIO_STAGE1_PART_ALPHA * param->effort;
  size_t len_full_alpha = s1param->nbest_sl *
    TUNE_RATIO_STAGE1_FULL_ALPHA * param->effort;

  /* size-cutoff of top sublattices based on partial alpha */
  sublattice_priority_queue_init (pqueue, len_part_alpha);
  
  /* return the nbest sublattices to pqueue ranked by the size of u */
  if (param->verbose >= 2)
    st = milliseconds ();

  re = return_best_sublattice ( poly,
                                s1param,
                                bound,
                                pqueue,
                                param->verbose );

  if (re == -1) {
    sublattice_priority_queue_clear (pqueue);
    return -1;
  }

  if (param->verbose >= 2) {
    gmp_fprintf ( stderr, "# Info: found    %8lu lat with part alpha\n",
            sublattice_priority_queue_size(pqueue));
    gmp_fprintf ( stderr, "# Info: ranked   %8zu lat with full alpha\n",
                  MIN(len_full_alpha, sublattice_priority_queue_size(pqueue)));
    gmp_fprintf ( stderr, "# Info: find best lat took %lums\n",
                  milliseconds () - st );
  }

  if (param->verbose >= 2)
    st = milliseconds ();
  
  struct transfer_to_alpha_priority_queue_arg arg[1];
  arg->alpha_pqueue = alpha_pqueue;
  arg->cpoly = poly->cpoly;
  arg->current_w = current_w;

  sublattice_priority_queue_do(pqueue, transfer_to_alpha_priority_queue, arg);
 
  if (param->verbose >= 2)
    gmp_fprintf ( stderr, "# Info: rank lat took %lums\n",
                  milliseconds () - st );

  /* free priority queue */
  sublattice_priority_queue_clear(pqueue);

  return 0;
}
