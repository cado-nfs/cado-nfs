/**
 * @file ropt_stage2.c
 * Root sieve.
 */


#include "ropt_stage2.h"


/**
 * Print sieve array for debugging.
 */
#if 0
void
print_sievearray ( double **A,
                   long A0,
                   long A1,
                   long B0,
                   long B1,
                   unsigned long K_ST,
                   unsigned long J_ST,
                   unsigned long MOD )
{
  long i, j;

  for (i = 0; i < B1 - B0 + 1; i ++) {
    for (j = 0; j < A1 - A0 + 1; j ++) {
      if (j % 16 == 0) {
        printf ("\n");
        printf ("[%3ld*%lu+%lu, %4ld*%lu+%lu]=[%4ld, %4ld]: ",
                i + B0, MOD, K_ST, j + A0, MOD,
                J_ST, (i + B0)*MOD+K_ST, (j + A0)*MOD+J_ST);
      }
      if (j % 4 == 0)
        printf (" ");
      printf ("%1.1f ", A[i][j]);
    }
    printf ("\n");
  }
  printf ("\n");
}
#endif


/**
 * Init the sieve array.
 */
static inline void
sievearray_init ( sievearray_t sievearray,
                  unsigned int i,
                  unsigned int j )
{
  /* sage: sum ([p/(p^2-1)*log(p) for p in prime_range(200)])
     4.2774597204802731 */
  float tmpf = SUP_ALPHA;
  int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
  unsigned long len = i * j, k;

  /* allocate matrix A. */
  sievearray->len_i = i;
  sievearray->len_j = j;
  sievearray->array = (int16_t *) malloc ( len * sizeof (int16_t) );

  if ( sievearray->array != NULL ) {
    for (k = 0; k < len; k++)
      sievearray->array[k] = tmpu;
  }
  else {
    fprintf ( stderr, "Error, cannot allocate memory in "
              "sievearray_init().\n" );
    exit (1);
  }
}


/**
 * Re-set sieve array.
 */
static inline void
sievearray_reset ( sievearray_t sievearray )
{
  float tmpf = SUP_ALPHA;
  int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
  unsigned long k, len = sievearray->len_i * sievearray->len_j;

  if (sievearray->array != NULL) {
    for (k = 0; k < len; k++)
      sievearray->array[k] = tmpu;
  }
  else {
    fprintf (stderr, "Error, null memory in sievearray_reset().\n");
    exit (1);
  }
}


/**
 * Free sieve array.
 */
static inline void
sievearray_free ( sievearray_t sievearray )
{
  free (sievearray->array);
}


/**
 * Root sieve on ARRAY along the j-line.
 */
static inline long
rootsieve_run_line ( int16_t *sa,
                     long j_bound,
                     long j,
                     unsigned int pe,
                     int16_t sub )
{
  long pel = (long) pe;

  /* bounding j with B0 < tmp + p*j < B1 */
  while ( j <= j_bound ) {
    sa[j] = sa[j] - sub;
    j += pel;
  }
  return j;
}


/**
 * Compute the indices for some special cases
 * Note that pl < pe
 */
static inline long
rootsieve_run_multroot_lift_idx ( unsigned int v,
                                  long Bmin,
                                  mpz_t B,
                                  mpz_t MOD,
                                  unsigned int pl,
                                  unsigned int pe )
{
  if (pl >= pe) {
    fprintf ( stderr, "Non-invertible element in "
              "rootsieve_run_multroot_lift_idx().\n" );
    exit(1);
  }
  /* findthe val_p (MOD), which must be < e */
  unsigned int a, b, c;
  long j_idx;

  /* (B-v)/p^l + (MOD/p^l)*k = 0 (mod p^{e-l})
     Note (B-v)/p^l is exact due to our constructions. */
  a = (unsigned int) mpz_fdiv_ui (B, pe);
  a = (a + pe - v) / pl; // +pe ensure it is positive
  b = (unsigned int) mpz_fdiv_ui (MOD, pe);
  b = b /pl;
  pl = pe / pl;  // it is now p^{e-l}
  c = (unsigned int) solve_lineq (a, b, 0, pl);
  j_idx = (long) ceil (((double) Bmin - (double) c) / (double) pl)
    * (long) pl + (long) c;
  j_idx = ab2ij (Bmin, j_idx);
  return j_idx;
}


#define DEBUG_MULTROOT_LIFT 0
/**
 * rootsieve_run_multroot_lift for the multiple root. Since we are
 * sure that level 1 node only contains one multiple root
 * r, we don't need to (and can't, otherwise, may count
 * repeatedly) consider single root at all.
 */
static inline void
rootsieve_run_multroot_lift ( node *currnode,
                              sievearray_t sa,
                              unsigned int *f_ui,
                              unsigned int *g_ui,
                              unsigned int *fuv_ui,
                              int d,
                              ropt_s2param_t s2param,
                              unsigned int p,
                              unsigned int max_e,
                              unsigned int curr_e,
                              int16_t sub )
{
  /* recursion end */
  if (currnode == NULL || curr_e > max_e || sub == 0)
    return;

  /* variables */
  int16_t subtmp;
  unsigned int i, j, k, l, nroots, pe, pl, pem1, fr, gr, step;
  node *tmpnode = NULL, *tmpnode2 = NULL;
  long i_idx, j_idx;

  /* compute p^e */
  pe = pem1 = 1;
  for (k = 0; k < curr_e - 1; k ++)
    pem1 = pem1 * p;
  pe = pem1 * p;

  /* val of p in MOD */
  l = 0;
  pl = p;
  while (mpz_divisible_ui_p(s2param->MOD, pl) != 0) {
    l ++;
    pl *= p;
  } // pl is p^{l+1}
  pl /= p;

  /* loop until all siblings are checked. */
  while ( currnode != NULL ) {

    /* loop all (lifted multiple) roots */
    for (nroots = 0; nroots < currnode->nr; nroots++) {

      /* compute g(r) */
      gr = eval_poly_ui_mod (g_ui, 1, currnode->r[nroots], pe);
      if (gr % p == 0)
        continue;

      /* compute f_uv(x) and then evaluate it at r. */
      compute_fuv_ui ( fuv_ui, f_ui, g_ui, d, 
                       currnode->u, currnode->v, pe );
      fr = eval_poly_ui_mod (fuv_ui, d, currnode->r[nroots], pe);
      fr = fr / pem1;

      /* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
      fr = (unsigned int) solve_lineq (fr, gr, 0, p); //fr = fr % p;

      /* solve (i, j) in  fr = i*r + j (mod p).
         - if r is not invertible, loop i with fixed j;
         - otherwise, loop j to solve i. */
      if (currnode->r[nroots] % p == 0) {

        for (i = 0; i < p; i ++) {

          /* if p|MOD and current pe|MOD:
             B + MOD*k = v (mod pe), so B = v (mod pe)
             A + MOD*k = u (mod pe), so A = u (mod pe)
             otherwise, no need to record this (u, v) at all. */
          if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {
            if (mpz_divisible_ui_p (s2param->MOD, pe) != 0) {
              if ( mpz_fdiv_ui (s2param->B, pe) !=
                   (currnode->v + fr * pem1) ||
                   mpz_fdiv_ui (s2param->A, pe) !=
                   (currnode->u + i * pem1) )
              {
#if DEBUG_MULTROOT_LIFT
                fprintf ( stderr, "skip found (%u, %u) in level %u.\n",
                          currnode->u + i * pem1,
                          currnode->v + fr * pem1,
                          curr_e );
#endif
                continue;
              }
            }
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {

#if DEBUG_MULTROOT_LIFT
            fprintf ( stderr, "level %u, (%u, %u), r: %u\n",
                      curr_e,
                      currnode->u + i * pem1,
                      currnode->v + fr * pem1,
                      currnode->r[nroots] + k * pem1 );
#endif

            insert_node ( currnode, &tmpnode,
                          currnode->u + i * pem1,
                          currnode->v + fr * pem1,
                          currnode->r[nroots] + k * pem1,
                          curr_e, p, pe, 0 );
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {
            insert_node ( currnode, &tmpnode, currnode->u + pem1 * i,
                          currnode->v + pem1 * fr,
                          currnode->r[nroots] + k * pem1,
                          curr_e, p, pe, 2 );
          }
        }
      }
      /* p -|- (currnode->r[nroots]), loop j to solve i. */
      else {
        for (j = 0; j < p; j ++) {

          /* given j, solve i in  fr = i*r + j (mod p). */
          i = solve_lineq (j, currnode->r[nroots], fr, p);

          /* if p|MOD and current pe|MOD:
             B + MOD*k = v (mod pe), so B = v (mod pe)
             A + MOD*k = u (mod pe), so A = u (mod pe)
             otherwise, no need to record this (u, v) at all. */
          if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {
            if (mpz_divisible_ui_p (s2param->MOD, pe) != 0) {
              if ( mpz_fdiv_ui (s2param->B, pe) !=
                   (currnode->v + j * pem1) ||
                   mpz_fdiv_ui (s2param->A, pe) !=
                   (currnode->u + i * pem1) )
              {
#if DEBUG_MULTROOT_LIFT
                fprintf ( stderr, "skip found (%u, %u) in level %u.\n",
                          currnode->u + i * pem1,
                          currnode->v + fr * pem1,
                          curr_e );
#endif
                continue;
              }
            }
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {
#if DEBUG_MULTROOT_LIFT
            fprintf ( stderr, "level %u, (%u, %u), r: %u\n",
                      curr_e, currnode->u + i * pem1,
                      currnode->v + j * pem1,
                      currnode->r[nroots] + k * pem1 );
#endif
            insert_node ( currnode, &tmpnode, currnode->u + pem1 * i,
                          currnode->v + pem1 * j,
                          currnode->r[nroots] + k * pem1,
                          curr_e, p, pe, 2 );
          }
        }
      }
    }  // next root of current (u, v)

    /* recursieve to next level, curr_e + 1 */
    rootsieve_run_multroot_lift ( currnode->firstchild,
                                  sa,
                                  f_ui,
                                  g_ui,
                                  fuv_ui,
                                  d,
                                  s2param,
                                  p,
                                  max_e,
                                  curr_e + 1,
                                  sub );

    /* we are in the second level from bottom, consider all
       children of this node (bottom nodes) and sieve */
    if (curr_e == max_e) {

      tmpnode = currnode->firstchild;

      while (tmpnode != NULL) {

        /* p|MOD, the inverse doesn't exist */
        if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {

          /* This shouldn't be called since we didn't record them
             during the lift. Only for safty purpose. */
          if (mpz_divisible_ui_p (s2param->MOD, pe) != 0) {
            tmpnode2 = tmpnode;
            tmpnode = tmpnode->nextsibling;
#if DEBUG_MULTROOT_LIFT
            fprintf (stderr, "deleting bottomnode (%u, %u) with %u "
                     "roots in level %u.\n",
                     tmpnode2->u, tmpnode2->v, tmpnode2->nr, curr_e);
#endif
            free_node (&tmpnode2);
            continue;
          }
          /* If current pe-|-MOD, then in B + MOD*x = v (mod pe),
             it is still possible to solve such x. */
          else {
            j_idx = rootsieve_run_multroot_lift_idx
              ( tmpnode->v, s2param->Bmin,
                s2param->B, s2param->MOD, pl, pe );

            i_idx = rootsieve_run_multroot_lift_idx
              ( tmpnode->u, s2param->Amin,
                s2param->A, s2param->MOD, pl, pe );
            step = pe / pl;
          }
        }
        /* p-|-MOD, the inverse exists */
        else {
          j_idx = uv2ij_mod (s2param->B, s2param->Bmin,
                             s2param->MOD, tmpnode->v, pe);
          i_idx = uv2ij_mod (s2param->A, s2param->Amin,
                             s2param->MOD, tmpnode->u, pe);
          step = pe;
        }

        /* be careful about this */
        subtmp = (int16_t) ceil ( (double) sub / (double) pem1
                                  * tmpnode->nr );

        /* For each i block */
        long old_j = j_idx;
        while (i_idx < (long) sa->len_i) {
          j_idx = old_j;
          rootsieve_run_line ( sa->array,
                               i_idx * sa->len_j + sa->len_j - 1,
                               i_idx * sa->len_j + j_idx,
                               step,
                               subtmp );
          i_idx += step;
        }

        tmpnode2 = tmpnode;
        tmpnode = tmpnode->nextsibling;

#if DEBUG_MULTROOT_LIFT
        fprintf ( stderr, "deleting bottomnode (%u, %u) with %u "
                  "roots in level %u AND ",
                  tmpnode2->u, tmpnode2->v, tmpnode2->nr, curr_e );
        fprintf ( stderr, "sieving bottomnode %d in steps %u.\n",
                  subtmp, step ); // pe
#endif

        free_node (&tmpnode2);
      }
    }

    /* delete current node and move to next sibling. */
    tmpnode = currnode;
    currnode = currnode->nextsibling;
    if (currnode != NULL)
      (currnode->parent)->firstchild = currnode;

    /* p|MOD, the inverse doesn't exist */
    if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {

      /* This shouldn't be called since we didn't record them
         during the lift. Only for safty purpose. */
      if (mpz_divisible_ui_p (s2param->MOD, pem1) != 0) {
#if DEBUG_MULTROOT_LIFT
        fprintf (stderr, "deleting (%u, %u) with %u roots in level %u.\n",
                 tmpnode->u, tmpnode->v, tmpnode->nr, curr_e - 1);
#endif
        free_node (&tmpnode);
        continue;
      }
      /* If current pem1-|-MOD, then in B + MOD*x = v (mod pem1),
         it is still possible to solve such x. */
      else {
        j_idx = rootsieve_run_multroot_lift_idx
          ( tmpnode->v, s2param->Bmin,
            s2param->B, s2param->MOD, pl, pem1 );
        i_idx = rootsieve_run_multroot_lift_idx
          ( tmpnode->u, s2param->Amin,
            s2param->A, s2param->MOD, pl, pem1 );
        step = pem1 / pl;
      }
    }
    /* p-|-MOD, the inverse exists */
    else {
      j_idx = uv2ij_mod (s2param->B, s2param->Bmin,
                         s2param->MOD, tmpnode->v, pem1);
      i_idx = uv2ij_mod (s2param->A, s2param->Amin,
                         s2param->MOD, tmpnode->u, pem1);
      step = pem1;
    }

    /* be careful about this */
    subtmp = (int16_t) ceil ( (double) sub / (double) pem1 *
                              (double) p  * tmpnode->nr );

    /* For each i block */
    long old_j = j_idx;
    while (i_idx < (long) sa->len_i) {

      j_idx = old_j;
      rootsieve_run_line ( sa->array,
                           i_idx * sa->len_j + sa->len_j - 1,
                           i_idx * sa->len_j + j_idx,
                           step, subtmp );
      i_idx += step;
    }

#if DEBUG_MULTROOT_LIFT
    fprintf (stderr, "deleting (%u, %u) with %u roots in level %u AND ",
             tmpnode->u, tmpnode->v, tmpnode->nr, curr_e - 1);
    fprintf (stderr, "sieving %d in steps %u.\n", subtmp, step); //pem1
#endif

    free_node (&tmpnode);
  }

  return;
}


/**
 * For multiple root r of the (i, j) (mod p).
 */
static inline void
rootsieve_run_multroot ( sievearray_t sa,
                         ropt_poly_t rs,
                         ropt_s2param_t s2param,
                         unsigned int u,
                         unsigned int i,
                         long j,
                         unsigned int r,
                         unsigned int p,
                         unsigned int e,
                         int16_t sub )
{
  /* sieve f(r) + ((u+i*p)*r + (B + MOD * (j+k*p)))*g(r) for i, k,
     we sieve (mod p) case in the beginning since r would
     be a multiple root for all (j + k*p).
     In the rootsieve_v(), it must be that
     roottype_flag = 1, e = 1; We can sieve starting from j.
     Note if roottype_flag = 2, e must > 1. We are not here. */
  if (e <= 1) {

    /* i is kept the same, next_i denotes the u-dim sieve */
    unsigned int next_i = i;
    long old_j = j;

    /* For each i block */
    while (next_i < sa->len_i) {

      j = old_j;
      rootsieve_run_line ( sa->array,
                           next_i * sa->len_j + sa->len_j - 1,
                           next_i * sa->len_j + j, p, sub );
      next_i += p;
    }
    return;
  }

  /* some variables */
  unsigned int v, pe, *f_ui, *g_ui, *fuv_ui;
  mpz_t tmpz;

  mpz_init_set_ui (tmpz, 0UL);

  /* compute p^e */
  pe = 1;
  for (v = 0; v < e ; v ++)
    pe = pe * p;

  /* use s.p instead of m.p */
  f_ui = (unsigned int*) malloc ((rs->d + 1) * sizeof (unsigned int));
  fuv_ui = (unsigned int*) malloc ((rs->d + 1) * sizeof (unsigned int));
  g_ui = (unsigned int*) malloc ((2) * sizeof (unsigned int));
  if ((f_ui == NULL) || (g_ui == NULL) || (fuv_ui == NULL)) {
    fprintf (stderr, "Error, cannot allocate memory in "
             "rootsieve_run_multroot(). \n");
    exit (1);
  }
  reduce_poly_ul (f_ui, rs->f, rs->d, pe);
  reduce_poly_ul (g_ui, rs->g, 1, pe);

  /* j -> v (mod p) */
  ij2uv (s2param->B, s2param->MOD, s2param->Bmin, j, tmpz);
  v = (unsigned int) mpz_fdiv_ui (tmpz, p);
  // v = (unsigned int) ( (tmp >= 0) ? tmp : tmp + (int) p );

#if DEBUG_MULTROOT_LIFT
  fprintf (stderr, "\n f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
           r, u, r, v, r, p);
#endif

  /* we've already known that r is a multiple root for f_{u, v}. */
  node *tmpnode, *root;
  new_tree (&root);
  root = new_node ();
  insert_node (root, &tmpnode, u, v, r, 1, p, p, 2);

  /* lift to higher p^e */
  rootsieve_run_multroot_lift ( root->firstchild, sa, f_ui, g_ui,
                                fuv_ui, rs->d, s2param, p,
                                e, 2, sub );

  /* free, either root itself or root with a level 1 node.
     Note, if e>1, freeing will be done in the lift function;
     However, if e=1, we need to manaully free root->firstchild. */
  free_node (&root);
  tmpnode = NULL;

  free (f_ui);
  free (fuv_ui);
  free (g_ui);
  mpz_clear (tmpz);
}


#define DEBUG_ROOTSIEVE_UV_ONEBLOCK 0


/**
 * Root sieve for f + u*x*g(x) + v*g(x)
 */
static inline void
rootsieve_one_block ( sievearray_t sa,
                      ropt_poly_t poly,
                      ropt_s2param_t s2param )
{
  char *roottype_flag;
  unsigned int p, r, u, v, k, np, nbb, max_e, totnbb, fr_ui, gr_ui, pe = 1;
  int16_t *subsgl, *submul;
  long i, bblock_size, i_idx, *j_idx, *j_idx_i0;
  float subf;
  mpz_t tmp;

  mpz_init (tmp);
  subsgl = (int16_t *) malloc (s2param->len_p_rs * sizeof(int16_t));
  submul = (int16_t *) malloc (s2param->len_p_rs * sizeof(int16_t));
  j_idx = (long *) malloc (primes[s2param->len_p_rs] * sizeof(long));
  j_idx_i0 = (long *) malloc (primes[s2param->len_p_rs] * sizeof(long));
  roottype_flag = (char *) malloc (primes[s2param->len_p_rs]
                                   * sizeof(char));

  if ( subsgl == NULL || submul == NULL || j_idx == NULL ||
       j_idx_i0 == NULL || roottype_flag == NULL ) {
    fprintf ( stderr, "Error, cannot allocate memory in "
              "rootsieve_one_block().\n" );
    exit (1);
  }

  bblock_size = L1_cachesize;

  totnbb = (unsigned int) ((s2param->Bmax - s2param->Bmin + 1)
                           / bblock_size);

#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
  fprintf ( stderr,
            "# Stat: totnbb: %u, bblock_size: %ld, total_size: %ld, "
            "sa->len_i: %u, sa->len_j: %u\n",
            totnbb, bblock_size, (long) totnbb * bblock_size,
            sa->len_i, sa->len_j );
#endif

  /* Init index array */
  for ( np = 0; np < s2param->len_p_rs; np ++ ) {
    p = primes[np];
    subf = (float) p * log ( (float) p) / ((float) p * (float) p - 1.0);
    subsgl[np] = (int16_t) ceil (subf * 1000.0);
    subf = log ( (float) p) / ( (float) p + 1.0);
    submul[np] = (int16_t) ceil (subf * 1000.0);
    roottype_flag[np] = 0;
  }

  /* Idx holders for each r < B */
  for ( r = 0; r < primes[s2param->len_p_rs]; r ++ ) {
    j_idx[r] = 0;
    j_idx_i0[r] = 0;
  }

  /* For each p < p_bound*/
  for (np = 0; np < s2param->len_p_rs; np ++) {

    p = primes[np];

    /* e depends on p */
    max_e = (unsigned int) (log (200.0) / log ((double) p));
    pe = 1;
    for (k = 0; k < max_e; k ++)
      pe = pe * p;

    /* pe | MOD in A + MOD*i = u (mod pe), don't need to sieve */
    if (mpz_divisible_ui_p (s2param->MOD, pe) != 0)
      continue;

    /* For each u < p */
    for (u = 0; u < p; u++) {

      i = 0;

      /* compute u->i in A + MOD*i = u (mod p) */
      if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {
        if (mpz_fdiv_ui (s2param->A, p) != u)
          continue;
        else {
          i = -1;
        }
      }
      else {
        /* i >= 0 */
        i = uv2ij_mod (s2param->A, s2param->Amin, s2param->MOD, u, p);
      }

      /* i is kept the same, i_idx denotes the u-dim sieve */
      i_idx = i;

      /* For each i block */
      while (i_idx < (long) sa->len_i) {
        
        /* For each j block on a i-line */
        for (nbb = 0; nbb < totnbb; nbb ++) {

          /* For each r < p_bound */
          for (r = 0; r < p; r++) {

            /* skip */
            if (mpz_divisible_ui_p(poly->gx[r], p) != 0)
              continue;

            /* u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
            mpz_mul (tmp, poly->gx[r], poly->gx[r]);
            mpz_mul_ui (tmp, tmp, u);
            mpz_sub (tmp, tmp, poly->numerator[r]);

            /* if nbb == 0:
               -- if i_idx == i, compute indices;
               -- elseif i_idx != 0, retrieve;
               if nbb != 0:
               -- next; */
            if (nbb == 0) {
              if (i_idx == i) {

                /* reset for all r < p for the first block */
                roottype_flag[r] = 0;
                j_idx[r] = 0;
                j_idx_i0[r] = 0;

                /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
                fr_ui = mpz_fdiv_ui (poly->fx[r], p);
                gr_ui = mpz_fdiv_ui (poly->gx[r], p);
                v = compute_v_ui (fr_ui, gr_ui, r, u, p);

                /* compute v->j in B + MOD*j = v (mod p) */
                if (mpz_divisible_ui_p (s2param->MOD, p) != 0) {

                  /* now A + MOD*i = u (mod p)
                     and B + MOD*j = v (mod p) */
                  if (mpz_fdiv_ui (s2param->B, p) == v) {

                    /* case where MOD has p-valuation smaller than e.
                       - If r is multiple, do the lift;
                       - If r is simple, all slots on the array 
                       are equivalent - ignore. */
                    if (mpz_divisible_ui_p(tmp, p) == 0)
                      continue;
                    else
                      roottype_flag[r] = 2;
#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
                    fprintf (stderr, "[mult] f(%u) + ( %u * %u + %u ) "
                             "* g(%u) = 0 (mod %u)\n",
                             r, u, r, v, r, p);
                    fprintf (stderr, " j_idx: %ld\n", i_idx, j_idx[r]);
#endif
                  }
                  else
                    continue;
                }
                else {
                  roottype_flag[r] = 1;

                  /* u -> i in A + MOD*i = u (mod p) has been computed. */
                  /* v -> j in B + MOD*j = v (mod p) */
                  j_idx_i0[r] = uv2ij_mod (s2param->B, s2param->Bmin,
                                           s2param->MOD, v, p);
                  j_idx[r] = j_idx_i0[r];

#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
                  fprintf ( stderr, "[simple] f(%u) + ( %u * %u + %u ) "
                            "* g(%u) = 0 (mod %u)\n",
                            r, u, r, v, r, p );
                  fprintf ( stderr, " i_idx: %ld, j_idx: %ld\n",
                            i_idx, j_idx[r] );
#endif
                }
              }
              /* if nbb == 0, but i_idx != i, retrieve info 
                 from nbb = 0 */
              else {
                /* don't change roottype_flag[np] */
                j_idx[r] = j_idx_i0[r];
              }
            }
            /* nbb != 0, need to find j_idx[r]'s j-position
               with respect to a line */
            else {

              /* prevent overflow when len_i = 1 and j wind back,
                 this happends since the last j-block could too 
                 short and the previous sieve jump out of this block
                 and hence % wind back to the beginning */
              if ((j_idx[r] > (long) sa->len_j) && (nbb == (totnbb-1)))
                  j_idx[r] = sa->len_j;
              else
                j_idx[r] %= sa->len_j;

              /*
              if (r == 1)
                printf ("(u=%u, p=%u, r=%u) nbb: %u, "
                        "idx: %u, j_idx: %u, len_i: %u\n",
                        u, p, r, nbb,
                        i_idx, j_idx[r], sa->len_i);
              */
            }

            /* cases where we don't want to sieve at all */
            if (roottype_flag[r] == 0)
              continue;

            /* r is a simple root for current (u, v, p). Then r
               is simple root for (u, v+i*p, p) for all i. */
            if (mpz_divisible_ui_p(tmp, p) == 0) {
              if (nbb != (totnbb-1)) {
                j_idx[r] = rootsieve_run_line (
                  sa->array,
                  i_idx * sa->len_j + bblock_size * (nbb + 1),
                  i_idx * sa->len_j + j_idx[r],
                  p,
                  subsgl[np] );
              }
              else {
                j_idx[r] = rootsieve_run_line (
                  sa->array,
                  i_idx * sa->len_j + (s2param->Bmax - s2param->Bmin),
                  i_idx * sa->len_j + j_idx[r],
                  p,
                  subsgl[np] );
              }
            }
            /* r is multiple root for current u, p */
            else {
              /* don't sieve in block, assume cm << cs */
              if (nbb == 0 && i_idx == i) {

                /* Two cases:
                   If roottype_flag = 1, then MOD != 0 (mod p),
                   everything is normal;
                   If MOD = 0 (mod p), then at nbb=0, roottype_flag=2;
                   Also must use u (mod pe) instead of u (mod p) */
                rootsieve_run_multroot ( sa, poly, s2param, u, i,
                                         j_idx[r], r, p, max_e,
                                         submul[np] );
              }
            }
          } // next r
        } // next j-block
        if (i_idx == -1)
          break;
        else {
          i_idx += p;
        }
      } // next i-block
    } // next u
  } // next p

  free (subsgl);
  free (submul);
  free (j_idx);
  free (j_idx_i0);
  free (roottype_flag);
  mpz_clear (tmp);
}


/**
 * Rootsieve for f + (u*x +v)*g for each sublattice
 */
static inline void
rootsieve_one_sublattice ( ropt_poly_t poly,
                           ropt_s2param_t s2param,
                           ropt_param_t param,
                           ropt_info_t info,
                           MurphyE_pq *global_E_pqueue )
{
  int i;
  unsigned long j, len_A, len_B, size_array_mem, size_B_block;
  long tmpBmax, tmpBmin;
  double MurphyE, best_MurphyE;
  mpz_t tmpv, tmpu;
  MurphyE_pq *local_E_pqueue;
  sievescore_pq *sievescore;
  sievearray_t sa;

  /* sieving length */
  len_A = (unsigned long) (s2param->Amax - s2param->Amin + 1);
  len_B = (unsigned long) (s2param->Bmax - s2param->Bmin + 1);

  /* size of sieveing array that fits in memory */
  size_array_mem = SIZE_SIEVEARRAY;

  if (len_B < (double) size_array_mem / len_A)
    size_array_mem = len_B * len_A;

  if (size_array_mem <= (unsigned long) len_A) {
    fprintf (stderr, "Error: Amax - Amin + 1 = %lu is too long. "
             "This happend when A is too large while B is "
             "too small. Ropt doesn't support this yet.\n",
             s2param->Amax - s2param->Amin + 1);
    exit(1);
  }

  /* size of B block that fits in memory at one time. Note:
     the siever is inefficient when len_A is large since
     the 2d-sieve is in the short L1-cache blocks. */
  size_B_block = (unsigned long) (
    (double) size_array_mem / (double) len_A );

  /* repair */
  size_array_mem = len_A * size_B_block;

  /* init structs */
  new_sievescore_pq (&sievescore, NUM_TOPALPHA_SIEVEARRAY);
  new_MurphyE_pq (&local_E_pqueue, NUM_TOPE_SUBLATTICE);
  sievearray_init (sa, len_A, size_B_block);
  mpz_init (tmpv);
  mpz_init (tmpu);
  tmpBmax = s2param->Bmax;
  tmpBmin = s2param->Bmin;

  /* for each i -> each u = A + MOD * i */
  int st = milliseconds ();
  do {

    /* fo reach block of size size_B_block */
    s2param->Bmax = s2param->Bmin + size_B_block - 1;
    if (s2param->Bmax > tmpBmax)
      s2param->Bmax = tmpBmax;

#if 0
    fprintf (stderr, "[%ld, %ld] ", s2param->Bmin, s2param->Bmax);
#endif
    
    /* root sieve for v. Note, s2param with updated Bmin 
       and Bmax will be used in rootsieve_v(). */
    rootsieve_one_block (sa, poly, s2param);

    /* compute the alpha of top slots */
    for (j = 0; j < size_array_mem; j++) {

      insert_sievescore_pq ( sievescore,
                             (j - j % (size_B_block)) / size_B_block,
                             j % (size_B_block),
                             sa->array[j] );

#if 0
      ij2uv (s2param->A, s2param->MOD, s2param->Amin,
             (j - j % (size_B_block)) / size_B_block,
             tmpu);
      ij2uv (s2param->B, s2param->MOD, s2param->Bmin,
             j % (size_B_block),
             tmpv);
      gmp_fprintf (stderr, "MAT[]: %d, u: %Zd, v: %Zd, i: %lu, "
                   "j: %lu\n", sa->array[j], tmpu, tmpv,
                   (j - j % (size_B_block)) / size_B_block,
                   j % (size_B_block) );
#endif
    }

    /* put sievescore into the MurphyE priority queue */
    for (i = 1; i < sievescore->used; i ++) {

      ij2uv (s2param->A, s2param->MOD, s2param->Amin, sievescore->i[i],
             tmpu);
      ij2uv (s2param->B, s2param->MOD, s2param->Bmin, sievescore->j[i],
             tmpv);

#if 0
      gmp_fprintf (stderr, "(i: %lu, j: %lu) -> (u: %Zd, v: %Zd), "
                   "score: %d\n", sievescore->i[i], sievescore->j[i],
                   tmpu, tmpv, sievescore->alpha[i]);
#endif

      compute_fuv_mp (s2param->f, poly->f, poly->g, poly->d, tmpu, tmpv);

      /* translation-only optimize */
      mpz_set (s2param->g[0], poly->g[0]);
      mpz_set (s2param->g[1], poly->g[1]);

      optimize_aux (s2param->f, poly->d, s2param->g, 0, 0, CIRCULAR);

      MurphyE = print_poly_fg (s2param->f, s2param->g, poly->d, poly->n, 0);

      insert_MurphyE_pq (local_E_pqueue, info->w, tmpu, tmpv,
                         s2param->MOD, MurphyE);

      /* // Further approximate.
      double skew = L2_skewness (s2param->f, poly->d,
                                 SKEWNESS_DEFAULT_PREC, 
                                 DEFAULT_L2_METHOD);
      MurphyE = L2_lognorm (s2param->f, poly->d, skew, DEFAULT_L2_METHOD);
      double alpha = get_alpha (s2param->f, poly->d, ALPHA_BOUND);
      insert_MurphyE_pq (local_E_pqueue, info->w, tmpu, tmpv, 
                         s2param->MOD,
                         - (MurphyE + alpha));
      */

    }

    /* next j */
    s2param->Bmin = s2param->Bmax + 1;

    /* reset sieve array and priority queue */
    sievearray_reset (sa);
    reset_sievescore_pq (sievescore);

  } while (s2param->Bmin < tmpBmax);

  /* recover */
  s2param->Bmax = tmpBmax;
  s2param->Bmin = tmpBmin;

  /* output polynomials */
  MurphyE = 0.0;
  best_MurphyE = 0.0;
  for (i = 1; i < local_E_pqueue->used; i ++) {
    
    if (param->verbose >= 2 && info->mode == 0) {
      gmp_fprintf ( stderr, "# Found (%3dth) (w=%d, u=%Zd, v=%Zd) "
                    "gives E = %1.2e\n",
                    i, 
                    local_E_pqueue->w[i],
                    local_E_pqueue->u[i],
                    local_E_pqueue->v[i],
                    local_E_pqueue->E[i] );

      if (param->verbose >= 3) {

        compute_fuv_mp (s2param->f, poly->f, poly->g, poly->d,
                        local_E_pqueue->u[i], local_E_pqueue->v[i]);

        print_poly_fg (s2param->f, s2param->g, poly->d, poly->n, 1);

        fprintf (stderr, "\n");
      }
    }

    /* insert E scores to a global queue (before optimization) */
    insert_MurphyE_pq ( global_E_pqueue,
                        local_E_pqueue->w[i],
                        local_E_pqueue->u[i],
                        local_E_pqueue->v[i],
                        local_E_pqueue->modulus[i],
                        local_E_pqueue->E[i] );

    MurphyE += local_E_pqueue->E[i];
    if (best_MurphyE < local_E_pqueue->E[i])
      best_MurphyE = local_E_pqueue->E[i];
  }

  /* compute average E*/
  info->ave_MurphyE = MurphyE / (double) (local_E_pqueue->used - 1);
  info->best_MurphyE = best_MurphyE;
  
  /* output stats */
  if (param->verbose >= 2 && info->mode == 0) {
    gmp_fprintf ( stderr,
                  "# Stat: ave. MurphyE of top %ld polynomials: "
                  "%1.2e (on sublattice %Zd, %Zd)\n",
                  local_E_pqueue->used - 1,
                  info->ave_MurphyE,
                  s2param->A, s2param->B );
    fprintf ( stderr, "# Stat: root sieve took %lums\n",
              milliseconds () - st );
  }

  /* free */
  sievearray_free (sa);
  free_sievescore_pq (&sievescore);
  free_MurphyE_pq (&local_E_pqueue);
  mpz_clear (tmpv);
  mpz_clear (tmpu);
}


/**
 * Stage 2: root sieve on the sublattice points;
 * Call rootsieve_uv().
 */
void
ropt_stage2 ( ropt_poly_t poly,
              ropt_s2param_t s2param,
              ropt_param_t param,
              ropt_info_t info,
              MurphyE_pq *global_E_pqueue,
              int w ) // enforce w explicitly.
{
  int i;

  //  if (param->verbose >= 2 && info->mode == 0)
  if (param->verbose >= 2 && info->mode == 0)
    ropt_s2param_print (s2param);

  for (i = 0; i <= poly->d; i++)
    mpz_set (s2param->f[i], poly->f[i]);

  info->w = w;
  compute_fuv_mp (s2param->f, poly->f, poly->g, poly->d,
                  s2param->A, s2param->B);

  /* sieve fuv */
  rootsieve_one_sublattice (poly, s2param, param, info,
                            global_E_pqueue);
}
