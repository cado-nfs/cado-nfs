/* Code to produce a stage 2 plan for P-1, P+1, and ECM.
   For given B1, B2 determines which value of d and pairs (i,j) to use
   so that the id+-j values cover the primes in ]B1, B2] */

#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stage2.h"
#include "gcd.h"       // gcd_ul
#include "gmp_aux.h"    // nbits
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "macros.h"

static inline int
cmp_uint2 (const void *p, const void *q)
{
  unsigned int *x = (unsigned int *) p;
  unsigned int *y = (unsigned int *) q;
  if (x[1] < y[1])
    return -1;
  else if (x[1] > y[1])
    return 1;
  else
    return (x[0] < y[0]) ? -1 : 1;
}

static inline unsigned long 
eulerphi_ul (unsigned long n)
{
  unsigned long p, r = 1UL;
  
  if (n == 0UL) /* Undefined, we return 0 */
    return 0UL;

  if (n % 2UL == 0UL)
    {
      n /= 2UL;
      while (n % 2UL == 0UL)
        {
          n /= 2UL;
          r *= 2UL;
        }
    }

  for (p = 3UL; p*p <= n; p += 2UL)
    {
      if (n % p == 0UL)
        {
          n /= p;
          r *= p - 1UL;
          while (n % p == 0UL)
            {
              n /= p;
              r *= p;
            }
        }
    }
  /* Now n is either 1 or a prime */
  if (n > 1UL)
    r *= n - 1UL;
  
  return r;
}

static inline unsigned int
compute_v (const unsigned int x, const unsigned int w)
{
  return (x+w/2)/w; /* v = round (x/w) */
}

/* Given x and w return u and v such that x = v*w +/- u and u < w/2
 *    v = floor((x+w/2)/w) = round (x/w)
 *    u = |x - w*v|
 */
static inline void
split_uv (unsigned int *u, unsigned int *v, const unsigned int x,
          const unsigned int w)
{
  *v = compute_v (x, w);
  unsigned int m = (*v) * w;
  *u = (x >= m) ? x - m : m - x;
}

/* We can do better if w <= 30 */
static inline double
stage2_baby_step_cost (const stage2_cost_t *opcost, const unsigned int w)
{
  /* We want to compute uP for all u in U and wP, using a Lucas chains.
   * Here we use the fact that 6|w so, for u in U, u = 1 or 5 mod 6
   */
  double cost = .0;
  /* First, we compute the Lucas chains: 1,2,3,5,6,7,11 */
  cost += 4. * opcost->dadd + 2. * opcost->dbl;
  /* Then, we compute 6*k+1,6*k+5 with 2 dADD for k in [2..floor((w-6)/12)] */
  unsigned int nstep = (w > 6) ? (w-6)/12-1 : 0;
  cost += ((double) nstep) * 2. * opcost->dadd;
  /* Finally, we compute wP */
  if (w == 6)
    cost += 0.; /* nothing to do, we already have 6 */
  else if (w == 12)
    cost += opcost->dbl; /* we have 6, double it to obtain 12 */
  else if (w % 12 == 0) /* last two elements of the chain are w/2-5, w/2-1 */
    cost += 2. * opcost->dadd; /* w/2+1 as (w/2-5)+6, w as (w/2+1)+(w/2-1) */
  else /* w % 12 == 6 */ /* last two elements of the chain are w/2-2, w/2+2 */
    cost += 1. * opcost->dadd + 1. * opcost->dbl;
    /* 4 as dbl(2) and w as (w/2-2)+(w/2+2) */
  return cost;
}

static inline double
stage2_giant_step_cost (const stage2_cost_t *opcost, unsigned int vmin,
                        const unsigned int vmax)
{
  if (vmin == 0 && vmax == 0)
    return 0.;

  if (vmin == 0)
    vmin++;

  if (vmin == vmax) /* do one Montgomery ladder for vmin */
  {
    double cost;
    for (cost = 0; vmin % 2 == 0; vmin >>= 1, cost += opcost->dbl) ;
    int b = nbits (vmin);
    return cost + (b-1) * opcost->dadd + (b-1) * opcost->dbl;
  }
  else
  {
    double cost;
    
    /* 1. Do a Montgomery ladder to compute simultaneously vmin and vmin+1 */
    int b = nbits (vmin);
    if (b > 1 && ((vmin >> (b-2)) & 1) == 0) /* 2nd most signigicant bit is 0 */
      cost = (b-1) * opcost->dadd + (b-1) * opcost->dbl;
    else
      cost = (b-1) * opcost->dadd + b * opcost->dbl;

    /* 2. Compute all v in [vmin+2..vmax]. The cost of computing v is 1 DBL if v
     * is even and v/2 >= vmin, 1 dADD otherwise.
     */
    if (2*vmin <= vmax)
    {
      if (vmin > 2) /* dADD until 2*vmin-1 */
        cost += (vmin-2) * opcost->dadd;
      unsigned int ndadd = (vmax-2*vmin+1)/2;
      unsigned int ndbl = (vmax-2*vmin+2)/2;
      if (vmin == 1)
        ndbl--; /* if vmin=1, 2vmin=vmin+1, no need to do a DBL to compute it */
      cost += ((double) ndadd) * opcost->dadd + ((double) ndbl) * opcost->dbl;
    }
    else /* only dADD */
      cost += (vmax-vmin-1) * opcost->dadd;

    /* Note: Not always optimal. For example, vmin=4 and vmax >= 6, we could
     * replace, in the computation of v=6, one dADD by one DBL, because v=3 was
     * computed during the Montgomery Ladder. For ECM on Montgomery curve, it
     * would save one M, it is not worth writing the specific code for it !
     */

    return cost;
  }
}

/* Return the cost of the stage2 with parameters B1, B2, w and given the cost of
 * operations.
 * If plan is not NULL, it is set to the computed plan.
 * If plan is NULL, only the cost is computed and returned.
 *
 * The unit for the costs of operations is a modular multiplication.
 * The comments of this function are written for an additive group.
 * Assume:
 *    w > 0   (there is an assert that checks this condition)
 *    6 | w   (there is an assert that checks this condition)
 *    w is coprime with all primes in ]B1, B2] (no assert checks this condition)
 *
 * P := { p | p prime and B1 < p <= B2 }
 * U := { u | u in Z, 1 <= u <= w/2 and gcd (u, w) == 1 }
 * V := { round(p/w) | p in P }
 * UV := { (p-w*round(p/w), round(p/w)) | p in P }
 *
 * ECM:
 *  - Q is the input point of stage2 (output of stage1)
 *  - Q is of order q=u+/-v*w => X(uQ)Z(vwQ) - X(vwQ)Z(uQ) = 0 mod q (due to
 *    the formula for differential addition for Montgomery curve)
 *  - The goal of step 2 is to compute the gcd with N of the following product:
 *      m := prod_{(u,v) in UV}{ X(uQ)Z(vwQ) - X(vwQ)Z(uQ) }
 *    We modify the points in order for all the points to have the same
 *    Z-coordinate and compute:
 *      m := Z * prod_{(u,v) in UV}{ X(uQ) - X(vwQ) }
 *         where Z is the common Z-coordinate of all points
 *
 * P-/+1:
 *  - X is the input of stage2 (output of stage1)
 *  - X is of order q=u+/-v*w => V_u (X) - V_vw (X) = 0 mod q (see section 4.7
 *    of Alexander Kruppa's thesis)
 *  - The goal of step 2 is to compute the gcd with N of the following product:
 *      m := prod_{(u,v) in UV}{ V_u(X) - V_vw(X) }
 *
 * Improvment: (u,v) and (-u,v) contribute to the same term in the product 
 * UV := { (u,v) in UxV | v*w+u or v*w-u in P }
 *
 * Improvment: let (u,v) in UV, p and p in P such that v*w+/-u = p and
 *             v*w-/+u is a multiple of q. Then we can remove the pair
 *             corresponding to q in UV (if it is not needed for another prime).
 *
 * We use these two ideas to try and minimize the size of the set UV.
 *
 *
 * In this function, we only considered the largest factor of v*w+/-u in the
 * interval ]B1,B2]. If B1^2 >= B2+w, then v*w+/-u can have at most 1 prime
 * factor in this interval and this function is optimal (in the sense that it
 * minimize #UV while maximizing vmin). If B1^2 < B2+w, this function will still
 * compute a correct stage2 plan but it may not be optimal.
 *
 * Note: the returned cost is not optimal for w <= 30, due to the way the baby
 * step is performed. In practice, we will never use such small value of w.
 */
static double
stage2_one_w (stage2_plan_t *plan, const unsigned int B1, const unsigned int B2,
              const unsigned int w, const stage2_cost_t *opcost,
              const unsigned char *primes, const unsigned int nprimes,
              const unsigned int *largest_factor, const int verbose)
{
  ASSERT_ALWAYS (w > 0);
  ASSERT_ALWAYS (w % 6 == 0);
  /* We assume 6|w in the baby step, and 2|w in different part of the code */

  unsigned int U_len, vmin, vmax, UV_len, nwanted;
  unsigned int *U = NULL, *UV = NULL;
  unsigned char *wanted, *nuv_per_p;
  double cost;
  unsigned int i, j; /* loop variables */

  if (verbose)
    printf ("# %s: w=%u\n", __func__, w);

  /* Set the list of wanted primes (copy the list of primes) */
  wanted = (unsigned char *) malloc ((B2+1) * sizeof (unsigned char));
  ASSERT_ALWAYS (wanted != NULL);
  memcpy (wanted, primes, (B2+1) * sizeof (unsigned char));
  nwanted = nprimes;

  /* We need to unsigned int for each (u,v) pair. We will need at most as many
   * (u,v) pairs as the number of primes in ]B1,B2].
   * We do not need to store the (u,v) pairs if we are not going to build the
   * plan. We need the number of (u,v) pairs to compute the cost.
   * The ith (u,v) pair is stored at indexes 2*i and 2*i+1 (u at even index).
   */
  if (plan)
  {
    UV = (unsigned int *) malloc (2 * nprimes * sizeof (unsigned int));
    ASSERT_ALWAYS (UV != NULL);
  }
  UV_len = 0;

  /* Compute U */
  U_len = (unsigned int) eulerphi_ul ((unsigned long) w) / 2;
  U = malloc (U_len * sizeof (unsigned int));
  ASSERT_ALWAYS (U != NULL);
  for (i = 0, j = 1; j < w / 2; j += 2) /* we know that 2 | w */
    if (gcd_ul ((unsigned long) j, (unsigned long) w) == 1UL)
      U[i++] = j;
  ASSERT_ALWAYS (i == U_len);
  /* There _are_ invertibles mod w ! */
  ASSERT_FOR_STATIC_ANALYZER (U_len > 0);
  if (verbose)
  {
    printf ("# %s: U = { %u", __func__, U[0]);
    for (unsigned int k = 1; k < U_len; k++)
      printf (", %u", U[k]);
    printf (" }\n# %s: #U = %u\n", __func__, U_len);
  }

  /* upper bound on V */
  {
      unsigned int p;
      for (p = B2; primes[p] == 0; p--); /* set p to the largest prime <= B2 */
      vmax = compute_v (p, w);
  }

  if (verbose)
    printf ("# %s: upper bound on V: vmax=%u\n", __func__, vmax);

  /* Count for each prime p in ]B1, B2], the number of (u,v) pairs that could be
   * use for this prime. We saturate at 255 because we only need this to
   * identify primes with only one (u,v) pair.
   */
  nuv_per_p = (unsigned char *) malloc ((B2+1) * sizeof (unsigned char));
  ASSERT_ALWAYS (nuv_per_p != NULL);
  memset (nuv_per_p, 0, (B2+1) * sizeof (unsigned char));

  for (unsigned int v = 0; v <= vmax; v++)
  {
    for (unsigned int k = 0; k < U_len; k++)
    {
      unsigned int u = U[k], n, f;
      
      n = v*w+u;
      f = largest_factor[n];
      if (nuv_per_p[f] != 0xff) /* saturate at 255 */
        nuv_per_p[f]++;
      
      n = (v*w > u) ? v*w-u : 0; /* we chech that v*w-u > 0 */
      f = largest_factor[n];
      if (nuv_per_p[f] != 0xff) /* saturate at 255 */
        nuv_per_p[f]++;
    }
  }

  /* pass 1: for all primes covered by one (u,v) pair, add this pair */
  vmin = vmax;
  for (unsigned int p = B1; p <= B2; p++)
  {
    /* check primes[p] because an (u,v) pair can cover up to 2 primes */
    if (wanted[p] && nuv_per_p[p] == 1)
    {
      unsigned int o, f, u, v;

      wanted[p] = 0; /* remove prime from wanted list */
      nwanted--;

      split_uv (&u, &v, p, w);
      if (verbose > 1)
        printf ("# %s: pass 1: (%u,%u) is the only pair covering %u\n",
                __func__, u, v, p);

      o = (p == v*w+u) ? ((v*w > u) ? v*w-u : 0) : v*w+u;

      f = largest_factor[o];
      if (wanted[f])
      {
        wanted[f] = 0; /* remove prime from wanted list */
        nwanted--;
        if (verbose > 1)
          printf ("# %s: pass 1: (%u,%u) also covers %u [as a factor of %u]\n",
                  __func__, u, v, f, o);
      }

      vmin = (v < vmin) ? v : vmin;

      if (UV) /* add the (u,v) pair in the list (if needed) */
      {
        UV[2*UV_len] = u;
        UV[2*UV_len+1] = v;
      }
      UV_len++;
    }
  }

  if (verbose > 1)
  {
    printf ("# %s: after pass 1: #UV=%u\n", __func__, UV_len);
    printf ("# %s: after pass 1: #covered primes=%u (out of %u)\n", __func__,
            nprimes-nwanted, nprimes);
    printf ("# %s: after pass 1: vmin=%u\n", __func__, vmin);
  }

  /* pass 2: add all (u, v) pairs, with v >= vmin, that covered 2 remaining
   *         primes.
   * pass 3: add all (u, v) pairs, with v >= vmin, that covered 1 remaining
   *         prime.
   * Then we decreased vmin by 1, pass 4 and 5 are identical to pass 2 and 3.
   * We continue to decrease vmin, until we have all the primes that we wanted.
   * In practice pass 2 and 3 are enough.
   */
  for (unsigned int pass = 2; nwanted; pass++)
  {
    unsigned int vstart = (pass == 2 || pass == 3) ? vmax : vmin;
    unsigned int s = (pass % 2) ? 1 : 2;

    for (unsigned int v = vstart; v >= vmin && nwanted; v--)
    {
      for (unsigned int k = 0; k < U_len && nwanted; k++)
      {
        unsigned int nc, f1, f2, u = U[k];
        unsigned int n1 = v*w+u, n2 = (v*w > u) ? v*w-u : 0;
        f1 = largest_factor[n1];
        f2 = largest_factor[n2];
        nc = (wanted[f1] && wanted[f2]) ? 2 : (wanted[f1] || wanted[f2]);
        if (nc == s)
        {
          if (verbose > 1 && s == 2)
            printf ("# %s: pass %u: (%u,%u) covers 2 primes: %u [as a factor "
                    "of %u] and %u [as a factor of %u]\n", __func__, pass, u, v,
                    f1, n1, f2, n2);
          else if (verbose > 1) /* s == 1 */
            printf ("# %s: pass %u: (%u,%u) covers 1 prime: %u [as a factor of "
                    "%u]\n", __func__, pass, u, v, wanted[f1] ? f1 : f2,
                    wanted[f1] ? n1 : n2);

          wanted[f1] = wanted[f2] = 0; /* remove prime(s) from wanted list */
          nwanted -= nc;

          if (UV) /* add the (u,v) pair in the list (if needed) */
          {
            UV[2*UV_len] = u;
            UV[2*UV_len+1] = v;
          }
          UV_len++;
        }
      }
      if (v == 0) /* we need this in the case where vmin = 0 */
        break;
    }

    if (verbose > 1)
    {
      printf ("# %s: after pass %u: #UV=%u\n", __func__, pass, UV_len);
      printf ("# %s: after pass %u: #covered primes=%u (out of %u)\n", __func__,
              pass, nprimes-nwanted, nprimes);
    }
    
    if (nwanted && s == 1)
    {
      ASSERT_ALWAYS (vmin > 0);
      vmin--;
      if (verbose)
        printf ("# %s: decreasing vmin to %u\n", __func__, vmin);
    }
  }

  if (verbose)
  {
    printf ("# %s: lower bound on V: vmin=%u\n", __func__, vmin);
    printf ("# %s: #UV=%u\n", __func__, UV_len);
    printf ("# %s: nprimes/#UV = %.2f\n", __func__, nprimes/(double) UV_len);
  }

#ifndef NDEBUG
  for (unsigned int i = B1; i <= B2; i++)
  {
    if (wanted[i])
      fprintf (stderr, "Error, prime %u is still marked as 'wanted'\n", i);
    ASSERT (wanted[i] == 0);
  }
#endif

  /* set plan if not NULL */
  if (plan)
  {
    plan->B2 = B2;
    plan->w = w;
    plan->U_len = U_len;
    plan->vmin = vmin;
    plan->vmax = vmax;
    plan->U = (unsigned int *) malloc (U_len * sizeof (unsigned int));
    ASSERT_ALWAYS (plan->U != NULL);
    memcpy (plan->U, U, U_len * sizeof (unsigned int));

    /* Alloc one unsigned int for each UV pair, one for each NEXT_D marker and
     * one for the final NEXT_PASS marker.
     */
    plan->pairs = (unsigned int *)
                          malloc ((UV_len+vmax-vmin+1) * sizeof (unsigned int));
    ASSERT_ALWAYS (plan->pairs != NULL);
  
    /* sort by increasing v then increasing u */
    qsort (UV, UV_len, 2*sizeof (unsigned int), cmp_uint2);
    ASSERT_ALWAYS (UV[1] == vmin);
    unsigned int idx = 0, vprev = vmin;
    unsigned int u_idx = 0; /* we store the index of u not the value of u */
    for (unsigned int k = 0; k < UV_len; k++)
    {
      unsigned int u = UV[2*k], v = UV[2*k+1];
      if (v != vprev)
      {
        if (verbose > 1)
          printf ("# %s: pairs[%u] = PAIR_INCR_V\n", __func__, idx);
        plan->pairs[idx++] = PAIR_INCR_V;
        vprev = v;
        u_idx = 0;
      }
      for ( ; U[u_idx] != u ; u_idx++) ; /* u values are in increasing order */
      if (verbose > 1)
        printf ("# %s: pairs[%u] = %u [ corresponds to (u,v)=(%u,%u) ]\n",
                __func__, idx, u_idx, u, v);
      plan->pairs[idx++] = u_idx;
    }
    if (verbose > 1)
      printf ("# %s: pairs[%u] = PAIR_END\n", __func__, idx);
    plan->pairs[idx++] = PAIR_END;
  }

  /* Compute cost */
  unsigned int V_len = vmax - vmin + 1;
  double c;
  cost = stage2_baby_step_cost (opcost, w);
  if (verbose > 1)
    printf ("# %s: cost of baby step: %f\n", __func__, cost);

  c = stage2_giant_step_cost (opcost, vmin, vmax);
  cost += c;
  if (verbose > 1)
    printf ("# %s: cost of giant step: %f\n", __func__, c);
  
  if (opcost->is_ecm)
  {
    c = 4*(U_len + V_len) - 6; /* see ecm.c */
    cost += c;
    if (verbose > 1)
      printf ("# %s: cost of homogenization of coordinates: %f\n", __func__, c);
  }

  c = UV_len;
  cost += c;
  if (verbose > 1)
    printf ("# %s: cost of computing UV: %f\n", __func__, c);

  if (verbose)
    printf ("# %s: total cost: %f\n", __func__, cost);
  
  /* Free memory */
  free (wanted);
  free (nuv_per_p);
  free (UV);
  free (U);

  return cost;
}

/* w_list contains 6 and all w < 1000 such that w = 6 * a * b * c with
 * a in {1, 2, 3}, b in {1, 5, 7, 11}, c in {5, 7, 11, 13, 17, 19} and c > b
 */
#define w_list_len 35
static unsigned int w_list[w_list_len] = { 6, 30, 42, 60, 66, 78, 84, 90, 102,
  114, 126, 132, 156, 198, 204, 210, 228, 234, 306, 330, 342, 390, 420, 462,
  510, 546, 570, 630, 660, 714, 780, 798, 858, 924, 990 };

void 
stage2_make_plan (stage2_plan_t *plan, const unsigned int B1, 
                  const unsigned int B2, const stage2_cost_t *opcost,
                  const int verbose)
{
  if (B2 <= B1)
  {
    memset (plan, 0, sizeof (stage2_plan_t));
    plan->B2 = B2;
    return;
  }

  unsigned char *primes = NULL;
  unsigned int p, nprimes, wmax, wmax_idx;
  unsigned int *largest_factor = NULL;

  if (verbose)
    printf ("# %s: B1=%u B2=%u opcost={.dadd=%f, .dbl=%f .is_ecm=%d}\n",
            __func__, B1, B2, opcost->dadd, opcost->dbl, opcost->is_ecm);

  /* Compute wmax:
   *    wmax <= B1^2-B2 (see comments above stage2_one_w)
   *    wmax <= 6*B1 to ensure that w is coprime with all primes in ]B1, B2]
   */
  for (wmax_idx = w_list_len-1 ; wmax_idx > 0 ; wmax_idx--)
    if (w_list[wmax_idx] <= 6*B1 && w_list[wmax_idx] <= B1*B1-B2)
      break;
  wmax = w_list[wmax_idx];
  if (verbose)
    printf ("# %s: wmax=%u\n", __func__, wmax);

  /* Compute list of primes B1 < p <= B2
   * primes[i] = 1 if and only if i is a prime and B1 < i <= B2.
   *
   * Compute largest factor of all integer n such that B1 < n < B2+wmax
   * largest_factor[n] = f if f is the largest factor of n and B1 < f <= B2
   *                   = 0 otherwise
   */
  primes = (unsigned char *) malloc ((B2+1) * sizeof (unsigned char));
  ASSERT_ALWAYS (primes != NULL);
  memset (primes, 0, (B2+1) * sizeof (unsigned char));

  largest_factor = (unsigned int *) malloc ((B2+wmax) * sizeof (unsigned int));
  ASSERT_ALWAYS (largest_factor != NULL);
  memset (largest_factor, 0, (B2+wmax) * sizeof (unsigned int));

  prime_info pi;
  prime_info_init (pi);
  for (p = 2; p <= B1; p = (unsigned int) getprime_mt (pi));
  /* Now p is the smallest prime > B1 */
  for (nprimes = 0 ; p <= B2; p = getprime_mt (pi), nprimes++)
  {
    primes[p] = 1;
    for (unsigned int kp = p; kp < B2+wmax; kp += p)
      largest_factor[kp] = p;
  }
  prime_info_clear (pi);

  if (verbose > 1)
  {
    printf ("# %s: stage2 primes:", __func__);
    for (unsigned int i = B1; i <= B2; i++)
      if (primes[i])
        printf (" %u", i);
    printf ("\n");
  }
  if (verbose)
    printf ("# %s: nprimes=%u\n", __func__, nprimes);

  /* find the w with the minimun cost */
  double best_cost = 0.;
  unsigned int best_w = 0;
  for (unsigned int w_idx = 0; w_idx < wmax_idx; w_idx++)
  {
    double cost = stage2_one_w (NULL, B1, B2, w_list[w_idx], opcost, primes,
                                nprimes, largest_factor, verbose);
    if (!best_w || cost < best_cost)
    {
      best_cost = cost;
      best_w = w_list[w_idx];
    }
  }
  if (verbose)
    printf ("# %s: best cost %f (with w=%u)\n", __func__, best_cost, best_w);
  ASSERT_ALWAYS (best_w != 0);

  /* write the plan */
  stage2_one_w (plan, B1, B2, best_w, opcost, primes, nprimes, largest_factor,
                verbose);

  /* free memory */
  free (largest_factor);
  free (primes);
}

void stage2_clear_plan (stage2_plan_t *plan)
{
  free (plan->pairs);
  plan->pairs = NULL;
  free (plan->U);
  plan->U = NULL;
  plan->U_len = 0;
  plan->vmin = 0;
  plan->vmax = 0;
  plan->w = 0;
}
