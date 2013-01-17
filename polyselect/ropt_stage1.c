/**
 * @file ropt_stage1.c
 * Called by ropt.c to find congruence classes.
 */


#include "ropt_stage1.h"


/**
 * Find good sublattice, the liftings.
 */
static inline void
find_sublattice_lift ( node *firstchild,
                       single_sublattice_pq *top,
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
        insert_single_sublattice_pq ( top,
                                      tmpnode->u,
                                      tmpnode->v,
                                      tmpnode->val,
                                      e );
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

  return;
}


/**
 * Find sublattices, the base case. Note, for the (u, v)
 * pairs, we consider all the simple + multi roots.
 */
static inline void
find_sublattice ( single_sublattice_pq *top,
                  ropt_poly_t poly,
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
      insert_single_sublattice_pq ( top,
                                    currnode->u,
                                    currnode->v,
                                    currnode->val,
                                    1 );
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

        insert_single_sublattice_pq ( top,
                                      currnode->u,
                                      currnode->v,
                                      currnode->nr / ((float)(p - 1)),
                                      1 );

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
    f_ui = (unsigned int*) malloc (
      (poly->d + 1) * sizeof (unsigned int) );
    fuv_ui = (unsigned int*) malloc (
      (poly->d + 1) * sizeof (unsigned int) );
    g_ui = (unsigned int*) malloc (
      2 * sizeof (unsigned int) );

    if ( (f_ui == NULL) ||
         (g_ui == NULL) ||
         (fuv_ui == NULL) ) {
      fprintf (stderr, "Error, cannot allocate memory in "
               "find_sublattice(). \n");
      exit (1);
    }

    /* compute f (mod pe) */
    reduce_poly_ul (f_ui, poly->f, poly->d, pe);
    reduce_poly_ul (g_ui, poly->g, 1, pe);

    find_sublattice_lift ( root->firstchild,
                           top,
                           f_ui,
                           g_ui,
                           fuv_ui,
                           poly->d,
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
 * Compute crt and add (u, v) to queue.
 */
static inline void
return_combined_sublattice_crt ( ropt_s1param_t s1param,
                                 ropt_bound_t bound,
                                 unsigned int *ind,
                                 unsigned int ***individual_sublattices,
                                 float **individual_sublattices_weighted_val,
                                 unsigned int *tprimes,
                                 sublattice_pq *pqueue )
{
  unsigned int i, e, pe[s1param->tlen_e_sl];
  mpz_t sum, inv, re, pe_z, tmpu1, tmpu2;
  float val = 0.0;
  
  mpz_init (sum);
  mpz_init (inv);
  mpz_init (re);
  mpz_init (pe_z);
  mpz_init (tmpu1);
  mpz_init (tmpu2);

  /* compute pe[] */
  for (i = 0; i < s1param->tlen_e_sl; i ++) {
    pe[i] = 1;
    /* note, the true e can be different from the s1param->e_sl[] */
    for (e = 0; e < individual_sublattices[i][ind[i]][2]; e ++)
      pe[i] *= tprimes[i];
  }

  /* compute modulus based on pe[], put it 
     temperorily to s1param->modulus */
  mpz_set_ui (s1param->modulus, 1UL);
  for (i = 0; i < s1param->tlen_e_sl; i ++) {
    mpz_mul_ui (s1param->modulus, s1param->modulus, pe[i]);
  }

  /* compute u */
  mpz_set_ui (sum, 0);
  for (i = 0; i < s1param->tlen_e_sl; i ++) {
    mpz_divexact_ui (re, s1param->modulus, pe[i]);
    mpz_set_ui (pe_z, pe[i]);
    mpz_invert (inv, re, pe_z);
    mpz_mul_ui (inv, inv, individual_sublattices[i][ind[i]][0]);
    mpz_mul (inv, inv, re);
    mpz_add (sum, sum, inv);
    val += log( (double) tprimes[i] ) * 
      individual_sublattices_weighted_val[i][ind[i]];
  }
  mpz_mod (tmpu1, sum, s1param->modulus);
  mpz_sub (tmpu2, s1param->modulus, tmpu1); // tmpu2 > 0
  
  /* if u is good, compute v */
  if ( mpz_cmp_si (tmpu1, bound->global_u_boundr) <= 0 ||
       mpz_cmp_si (tmpu2, -bound->global_u_boundl) <= 0 ) {

    /* compute v */
    mpz_set_ui (sum, 0);
    for (i = 0; i < s1param->tlen_e_sl; i ++) {
      mpz_divexact_ui (re, s1param->modulus, pe[i]);
      mpz_set_ui (pe_z, pe[i]);
      mpz_invert (inv, re, pe_z);
      mpz_mul_ui (inv, inv, individual_sublattices[i][ind[i]][1]);
      mpz_mul (inv, inv, re);
      mpz_add (sum, sum, inv);
    }
    mpz_mod (tmpu2, sum, s1param->modulus);

    /* insert this node */
    insert_sublattice_pq ( pqueue, tmpu1, tmpu2, s1param->modulus, val );
  }

  mpz_clear (sum);
  mpz_clear (inv);
  mpz_clear (re);
  mpz_clear (pe_z);
  mpz_clear (tmpu1);
  mpz_clear (tmpu2);
}


/**
 * Actual length check function.
 */
static inline void
return_combined_sublattice_check_tlen ( ropt_s1param_t s1param )
{
  unsigned int i;
  
  /* find the actual s1param->len_e_sl, excluding those 0's */
  s1param->tlen_e_sl = s1param->len_e_sl;
  for (i = 0; i < s1param->len_e_sl; i ++) {
    if (s1param->e_sl[i] == 0)
      s1param->tlen_e_sl --;
  }

  /* check */
  if (s1param->tlen_e_sl < 1 || s1param->tlen_e_sl > 10) {
    fprintf ( stderr, "Error, number of primes in \"-e\" (len_e_sl) "
              "should be between 1 and 10\n" );
    exit(1);
  }
}


/**
 * Return all sublattices by calling CRT, where for each sublattice,
 * the seperate (mod p) valuations are the best ones.
 */
static inline int
return_combined_sublattice ( ropt_poly_t poly,
                             ropt_s1param_t s1param,
                             ropt_bound_t bound,
                             sublattice_pq *pqueue,
                             int verbose )
{
  /* get s1param->tlen_e_sl */
  return_combined_sublattice_check_tlen (s1param);

  unsigned int i, j, k, count,
    t_primes[s1param->tlen_e_sl],
    t_e_sl[s1param->tlen_e_sl],
    t_size[s1param->tlen_e_sl],
    ind[s1param->tlen_e_sl],
    ***individual_sublattices;
  float **individual_sublattices_weighted_val;
  single_sublattice_pq *top;

  /* sublattice[i][length][] save (u, v) for prime[i] */
  individual_sublattices = (unsigned int ***) malloc (
    s1param->tlen_e_sl * sizeof (unsigned int **) );

  individual_sublattices_weighted_val = (float **) malloc (
    s1param->tlen_e_sl * sizeof (float *) );

  if ( individual_sublattices == NULL ||
       individual_sublattices_weighted_val == NULL ) {
    fprintf ( stderr,
              "Error, cannot allocate memory in "
              "return_combined_sublattice(). \n" );
    exit (1);
  }

  /* get t_primes and t_e_sl */
  for (i = 0, j = 0; i < s1param->len_e_sl; i ++) {
    if (s1param->e_sl[i] != 0) {
      t_primes[j] = primes[i];
      t_e_sl[j++] = s1param->e_sl[i];
    }
  }

  /* decide the number of top individual sublattices. Note that
     the s1param->tlen_e_sl must be already set */
  if (s1param->nbest_sl_tunemode == 0)
    ropt_s1param_setup_individual_nbest_sl (s1param);
  else
    ropt_s1param_setup_individual_nbest_sl_tune (s1param);

  /* for each prime[i], lift the roots */
  for (i = 0; i < s1param->tlen_e_sl; i ++) {

    new_single_sublattice_pq (&top, s1param->individual_nbest_sl[i]);

    /* find individual sublattices */
    find_sublattice (top, poly, t_primes[i], t_e_sl[i]);
    t_size[i] = top->used - 1;

    /* if length zero, this set of parameters fails. */
    if (t_size[i] == 0) {
      for (k = 0; k < s1param->tlen_e_sl; k ++) {
        for (j = 0; j < t_size[i]; j ++) {
          if (individual_sublattices[k][j] != NULL)
            free (individual_sublattices[k][j]);
        }
        if (individual_sublattices[k] != NULL) {
          free (individual_sublattices[k]);
        }
        if (individual_sublattices_weighted_val[k] != NULL) {
          free (individual_sublattices_weighted_val[k]);
        }
      }
      if (individual_sublattices != NULL)
        free (individual_sublattices);
      if (individual_sublattices_weighted_val != NULL)
        free (individual_sublattices_weighted_val);
      return 0;
    }

    /* allocate array for the i-th prime. */
    individual_sublattices[i] = (unsigned int **) malloc (
      t_size[i] * sizeof(unsigned int *) );
    individual_sublattices_weighted_val[i] = (float *) malloc (
      t_size[i] * sizeof(float) );
    if ( (individual_sublattices)[i] == NULL || 
         (individual_sublattices_weighted_val)[i] == NULL ) {
      fprintf ( stderr, "Error, cannot allocate memory in "
                "return_combined_sublattice(). \n" );
      exit (1);
    }

    for (j = 0; j < t_size[i]; j ++) {

      (individual_sublattices)[i][j] = (unsigned int *) 
        malloc ( 3 * sizeof(unsigned int) );
      if ( individual_sublattices[i][j] == NULL ) {
        fprintf ( stderr, "Error, cannot allocate memory in "
                  "return_combined_sublattice(). \n" );
        exit (1);
      }

      individual_sublattices[i][j][0] = top->u[j+1];
      individual_sublattices[i][j][1] = top->v[j+1];
      individual_sublattices[i][j][2] = (unsigned int) top->e[j+1];
      individual_sublattices_weighted_val[i][j] = top->val[j+1];
      /*
        fprintf ( stderr, "SUBLATTICE: #%u, (%u, %u), e: %u, val: %.2f\n",
        j, 
        individual_sublattices[i][j][0],
        individual_sublattices[i][j][1],
        individual_sublattices[i][j][2],
        individual_sublattices_weighted_val[i][j] );
      */
    }

    free_single_sublattice_pq (&top);

    if (verbose >= 2)
      fprintf ( stderr,
                "# Info: p: %2u, max_e: %2u, list_size: %6u\n",
                t_primes[i], t_e_sl[i], t_size[i] );

  }

  /* Loop over combinations of all arrays. This is awkward.
     We could map 0 ... \prod pe[i] to the indices of the
     arrays in the price of using more arithmetic. */

  /* 2 */
  if (s1param->tlen_e_sl == 1) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      return_combined_sublattice_crt ( s1param,
                                       bound,
                                       ind,
                                       individual_sublattices,
                                       individual_sublattices_weighted_val,
                                       t_primes,
                                       pqueue );
  }
  /* 2, 3 */
  else if (s1param->tlen_e_sl == 2) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        return_combined_sublattice_crt ( s1param,
                                         bound,
                                         ind,
                                         individual_sublattices,
                                         individual_sublattices_weighted_val,
                                         t_primes,
                                         pqueue );
  }
  /* 2, 3, 5 */
  else if (s1param->tlen_e_sl == 3) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          return_combined_sublattice_crt ( s1param,
                                           bound,
                                           ind,
                                           individual_sublattices,
                                           individual_sublattices_weighted_val,
                                           t_primes,
                                           pqueue );
  }
  /* 2, 3, 5, 7 */
  else if (s1param->tlen_e_sl == 4) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            return_combined_sublattice_crt ( s1param,
                                             bound,
                                             ind,
                                             individual_sublattices,
                                             individual_sublattices_weighted_val,
                                             t_primes,
                                             pqueue );
  }
  /* 2, 3, 5, 7, 11 */
  else if (s1param->tlen_e_sl == 5) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              return_combined_sublattice_crt ( s1param,
                                               bound,
                                               ind,
                                               individual_sublattices,
                                               individual_sublattices_weighted_val,
                                               t_primes,
                                               pqueue );
  }
  /* 2, 3, 5, 7, 11, 13 */
  else if (s1param->tlen_e_sl == 6) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < t_size[5]; ind[5] ++)
                return_combined_sublattice_crt ( s1param,
                                                 bound,
                                                 ind,
                                                 individual_sublattices,
                                                 individual_sublattices_weighted_val,
                                                 t_primes,
                                                 pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17 */
  else if (s1param->tlen_e_sl == 7) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < t_size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < t_size[6]; ind[6] ++)
                  return_combined_sublattice_crt ( s1param,
                                                   bound,
                                                   ind,
                                                   individual_sublattices,
                                                   individual_sublattices_weighted_val,
                                                   t_primes,
                                                   pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19 */
  else if (s1param->tlen_e_sl == 8) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < t_size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < t_size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < t_size[7]; ind[7] ++)
                    return_combined_sublattice_crt ( s1param,
                                                     bound,
                                                     ind,
                                                     individual_sublattices,
                                                     individual_sublattices_weighted_val,
                                                     t_primes,
                                                     pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23 */
  else if (s1param->tlen_e_sl == 9) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < t_size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < t_size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < t_size[7]; ind[7] ++)
                    for (ind[8] = 0; ind[8] < t_size[8]; ind[8] ++)
                      return_combined_sublattice_crt ( s1param,
                                                       bound,
                                                       ind,
                                                       individual_sublattices,
                                                       individual_sublattices_weighted_val,
                                                       t_primes,
                                                       pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 */
  else if (s1param->tlen_e_sl == 10) {
    for (ind[0] = 0; ind[0] < t_size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < t_size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < t_size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < t_size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < t_size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < t_size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < t_size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < t_size[7]; ind[7] ++)
                    for (ind[8] = 0; ind[8] < t_size[8]; ind[8] ++)
                      for (ind[9] = 0; ind[9] < t_size[9]; ind[9] ++)
                        return_combined_sublattice_crt ( s1param,
                                                         bound,
                                                         ind,
                                                         individual_sublattices,
                                                         individual_sublattices_weighted_val,
                                                         t_primes,
                                                         pqueue );
  }
  /* too aggressive */
  else {
    fprintf ( stderr, "Error, number of primes in \"-e\" (len_e_sl) "
              "should be between 1 and 10\n" );
    exit (1);
  }

  /* info */
  if (verbose >= 2) {
    /* Compute s1param->modulus */
    count = 1;
    for (i = 0; i < s1param->tlen_e_sl; i ++) {
      count *= t_size[i];
    }
    fprintf (stderr, "# Info: computed %u CRTs\n", count);
  }

  /* clearence */
  for (i = 0; i < s1param->tlen_e_sl; i ++) {
    for (j = 0; j < t_size[i]; j ++) {
      free (individual_sublattices[i][j]);
    }
    free (individual_sublattices[i]);
    free (individual_sublattices_weighted_val[i]);
  }
  free (individual_sublattices);
  free (individual_sublattices_weighted_val);

  return 1;
}


/**
 * Call return_all_sublattice() to return good sublattices.
 * Choose the "s1param->nbest_sl" best ones among them if there
 * are more quantities than it. The "best" property is ranked
 * by the size of u.
 */
static inline int
return_best_sublattice ( ropt_poly_t poly,
                         ropt_s1param_t s1param,
                         ropt_bound_t bound,
                         sublattice_pq *pqueue,
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
  while (pqueue->used == 1) {
    if (bound->global_u_boundr < (LONG_MAX>>1)) {
      bound->global_u_boundr *= 2;
      if (verbose >= 2) {
        fprintf ( stderr,
                  "# Warn: not enough sublattice classes. "
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
              "# Info: best %d sublattices (\"size_total_sublattices\") "
              "where |u| < %lu\n",
              pqueue->used - 1, bound->global_u_boundr);
  }

  /* recover, for deg 6 poly */
  bound->global_u_boundr = tmp;

  return 1;
}


/**
 * Stage 1: record good sublattices to "alpha_pqueue".
 */
int
ropt_stage1 ( ropt_poly_t poly,
              ropt_bound_t bound,
              ropt_s1param_t s1param,
              ropt_param_t param,
              alpha_pq *alpha_pqueue,
              int current_w )
{
  int st = 0, i, re;
  double alpha_lat;
  mpz_t *fuv;
  sublattice_pq *pqueue;

  /* size-cutoff of top sublattices */
  new_sublattice_pq (&pqueue, s1param->nbest_sl);

  /* return the nbest sublattices to pqueue ranked by the size of u */
  if (param->verbose >= 2)
    st = milliseconds ();

  re = return_best_sublattice ( poly,
                                s1param,
                                bound,
                                pqueue,
                                param->verbose );

  if (re == -1) {
    free_sublattice_pq (&pqueue);
    return -1;
  }

  if (param->verbose >= 2)
    gmp_fprintf ( stderr, "# Info: find best sublattices took %lums\n",
                  milliseconds () - st );

  /* fuv is f+(u*x+v)*g */
  fuv = (mpz_t*) malloc ((poly->d + 1) * sizeof (mpz_t));
  if (fuv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_stage1().\n");
    exit (1);
  }
  for (i = 0; i <= poly->d; i++)
    mpz_init_set (fuv[i], poly->f[i]);

  if (param->verbose >= 2)
    st = milliseconds ();
  
  /* put pqueue into the global alpha_pqueue ranked by parial alpha */
  for (i = 1; i < pqueue->used; i ++) {

    compute_fuv_mp ( fuv, poly->f, poly->g, poly->d,
                     pqueue->u[i], pqueue->v[i] );


#if RANK_SUBLATTICE_BY_E
    /* use exp_E as benchmark instead of alpha. */
    double skew = L2_skewness (fuv, poly->d, SKEWNESS_DEFAULT_PREC, 
                               DEFAULT_L2_METHOD);
    alpha_lat = L2_lognorm (fuv, poly->d, skew, DEFAULT_L2_METHOD);
    alpha_lat += get_alpha (fuv, poly->d, 2000);
#else
    //alpha_lat = get_alpha (fuv, poly->d, primes[s1param->tlen_e_sl-1]);
    alpha_lat = get_alpha (fuv, poly->d, 2000);
#endif

#if DEBUG_ROPT_STAGE1
     skew = L2_skewness ( fuv, poly->d, SKEWNESS_DEFAULT_PREC,
                                DEFAULT_L2_METHOD );
    double logmu = L2_lognorm (fuv, poly->d, skew, DEFAULT_L2_METHOD);
    gmp_fprintf ( stderr, "# Info: insert sublattice #%4d, (w, u, v): "
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
    insert_alpha_pq ( alpha_pqueue,
                      current_w,
                      pqueue->u[i],
                      pqueue->v[i],
                      pqueue->modulus[i],
                      alpha_lat );
  }

  if (param->verbose >= 2)
    gmp_fprintf ( stderr, "# Info: rank above sublattices took %lums\n",
                  milliseconds () - st );

  /* free priority queue */
  free_sublattice_pq (&pqueue);
  for (i = 0; i <= poly->d; i++) {
    mpz_clear (fuv[i]);
  }
  free (fuv);

  return 0;
}
