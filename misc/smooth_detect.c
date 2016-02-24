#define _DEFAULT_SOURCE  // for random()
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <gmp.h>
#include "ecm.h"

#include "smooth_detect.h"

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif
#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

// For a bug in ecm ?
double default_B1done;

#define B1MIN 200            // Value of B1 to start each new candidate
                             // If you change it, you must update the
                             // table expected_effort[] below.
#define EA_THRESHOLD 0.8     // Heuristic for early-abort: low = keep many

// Expected effort to extract of prime of p bits.
#define MAX_PBIT 90
static const double expected_effort[MAX_PBIT] = {
  200, 200, 200, 200, 200, 200, 200, 200, 200, 200, // 0 - 9
  200, 200, 200, 200, 200, 200, 200, 211, 219, 243, // 10 - 19
  296, 280, 351, 380, 482, 588, 653, 719, 834, 985, // 20 - 29
  1280, 1489, 2352, 2324, 2714, 3813, 5141, 5891, 7986, 8816, // 30 - 39
  10087, 12264, 14191, 18827, 25633, 25491, 37803, 39392, 44290, 51925, // 40-49
  87203, 79943, 110007, 121644, 147602, // 50 - 54
  174995, 199245, 257190, 279228, 345960, // 55 - 59
  351682, 530703, 640140, 759310, 775311, // 60 - 64
  960249, 1267879, 1174122, 1272107, 1589907, // 65 - 69
  2258437, 2004235, 2903959, 3002629, 3888904, // 70 - 74
  4373729, 4899345, 5218152, 6269843, 7063446, // 75 - 79
  9553542, 9891138, 10623352, 13795248, 17574109, // 80 - 84
  18790448, 23529670, 24757303, 30897420, 31188626, // 85 - 89
};


double get_time() {
  return (double)(clock()) / CLOCKS_PER_SEC;
}


///////////////////////////////////////////////////////////////////////////
// Candidate: structure to store a pair of numbers currently being
// factored.
///////////////////////////////////////////////////////////////////////////

void cand_init(cand_t c) {
  mpz_init(c->u0);
  mpz_init(c->v0);
  mpz_init(c->u);
  mpz_init(c->v);
}

void cand_set(cand_t c, const cand_t c0) {
  mpz_set(c->u0, c0->u0); 
  mpz_set(c->v0, c0->v0); 
  mpz_set(c->u, c0->u); 
  mpz_set(c->v, c0->v); 
  c->lpu = c0->lpu;
  c->lpv = c0->lpv;
  c->effort = c0->effort;
  c->id = c0->id;
}

void cand_init_copy(cand_t c, const cand_t c0) {
  cand_init(c);
  cand_set(c, c0);
}
 
void cand_clear(cand_t c) { 
  mpz_clear(c->u0); 
  mpz_clear(c->v0); 
  mpz_clear(c->u); 
  mpz_clear(c->v); 
} 

void cand_swap(cand_t c, cand_t d) {
  mpz_swap(c->u0, d->u0);
  mpz_swap(c->v0, d->v0);
  mpz_swap(c->u, d->u);
  mpz_swap(c->v, d->v);
  unsigned int t;
  t = c->lpu; c->lpu = d->lpu; d->lpu = t;
  t = c->lpv; c->lpv = d->lpv; d->lpv = t;
  double tt;
  tt = c->effort; c->effort = d->effort; d->effort = tt;
  unsigned long id;
  id = c->id; c->id = d->id; d->id = id;
}

void cand_print(const cand_t c) {
  printf("Candidate id = %lu\n", c->id);
  gmp_printf("u0=%Zd\n", c->u0);
  gmp_printf("v0=%Zd\n", c->v0);
  gmp_printf("u=%Zd (%lu bits) largest prime so far of %lu bits\n",
          c->u, mpz_sizeinbase (c->u, 2), c->lpu);
  gmp_printf("v=%Zd (%lu bits) largest prime so far of %lu bits\n",
          c->v, mpz_sizeinbase (c->v, 2), c->lpv);
  printf("effort=%.0f\n", c->effort);
}

bool cand_is_factored(const cand_t c) {
  return (mpz_cmp_ui(c->u, 1) == 0) && (mpz_cmp_ui(c->v, 1) == 0);
}

// if effort is such that all primes up to b-bits have been removed,
// and if an unfactored part has less than 3*b bits, then there are only
// two factors. If size is more than 2*bound, it can not be smooth.
// Heuristic: we consider that if the effort is twice the average for
// detecting p bits, then all p bits prime have been removed.
bool cand_is_probably_not_smooth(const cand_t c, unsigned int bound) {
  double eff = c->effort;
  unsigned int bits = 0;
  while ((bits < MAX_PBIT) && (2*expected_effort[bits] < eff)) {
    bits++;
  }
  if (bits == MAX_PBIT) {
    return false; // can not conclude; no data for this effort
  }
  // bits -= 3; // take some margin

  for (unsigned int k = 2; k < 4; ++k) {
    unsigned long bu = mpz_sizeinbase(c->u, 2);
    if ((bu < (k+1)*bits) && (bu > k*bound)) {
      printf("Probably not smooth, level %d: %lu bits!\n", k, bu);
      return true;
    }
    unsigned long bv = mpz_sizeinbase(c->v, 2);
    if ((bv < (k+1)*bits) && (bv > k*bound)) {
      printf("Probably not smooth, level %d: %lu bits!\n", k, bv);
      return true;
    }
  }
  return false;
}

void cand_update_check_prime(cand_t c) {
  if (mpz_probab_prime_p(c->u, 1)) {
    c->lpu = MAX(c->lpu, mpz_sizeinbase(c->u, 2));
    mpz_set_ui(c->u, 1);
  }
  if (mpz_probab_prime_p(c->v, 1)) {
    c->lpv = MAX(c->lpv, mpz_sizeinbase(c->v, 2));
    mpz_set_ui(c->v, 1);
  }
}

void cand_set_original_values(cand_t c, const mpz_t u0, const mpz_t v0,
    unsigned long id) {
  mpz_set(c->u0, u0);
  mpz_set(c->v0, v0);
  mpz_set(c->u, u0);
  mpz_set(c->v, v0);
  c->lpu = 0;
  c->lpv = 0;
  c->effort = 0.0;
  c->id = id;
  cand_update_check_prime(c);
}

// Cost is the bitsize of the largest unfactored part
// For factored numbers, this is INT_MAX
int cost(const cand_t c) {
  if (cand_is_factored(c))
    return INT_MAX;
  else
    return MAX(mpz_sizeinbase(c->u, 2), mpz_sizeinbase(c->v, 2));
}

// For use in qsort
int cand_cmp(const void *c1, const void *c2) {
  cand_t *C1 = (cand_t *)c1;
  cand_t *C2 = (cand_t *)c2;
  if (cost(*C1) > cost(*C2))
    return 1;
  if (cost(*C1) < cost(*C2))
    return -1;
  return 0;
}


///////////////////////////////////////////////////////////////////////////
// Pool: contains a list of candidates that are still interesting to 
// try to factor
/////////////////////////////////////////////////////////////////////////////

// The list of candidates must always be sorted by increasing cost.
typedef struct {
  cand_t *list;
  unsigned int n;
} pool_s;
typedef pool_s pool_t[1];

void pool_init(pool_t p) {
  p->list = NULL;
  p->n = 0;
}

void pool_clear(pool_t p) {
  for(unsigned int i = 0; i < p->n; i++)
    cand_clear(p->list[i]);
  free(p->list);
}

// if pool has been modified, sort it again
void pool_sort(pool_t p) {
  qsort(p->list, p->n, sizeof(cand_t), cand_cmp);
}

// Insert candidate in pool, keeping it sorted according to cost
void pool_insert(pool_t p, const cand_t c) {
  p->list = (cand_t *)realloc(p->list, (p->n + 1)*sizeof(cand_t));
  cand_init_copy(p->list[p->n], c);
  for(unsigned int i = p->n;
      i > 0 && cost(p->list[i-1]) > cost(p->list[i]);
      i--) {
    cand_swap(p->list[i-1], p->list[i]);
  }
  p->n += 1;
}

// Remove candidates that have a prime more than lmax, and those that
// are fully factored. If there are still more than max_size elements,
// keep those with smallest cost.
void pool_purge(pool_t p, unsigned int max_size, unsigned int lmax) {
  unsigned int i, j;
  j = 0;
  for (i = 0; i < p->n; ++i) {
    if (!cand_is_factored(p->list[i]) &&
        !cand_is_probably_not_smooth(p->list[i], lmax) &&
        (MAX(p->list[i]->lpu, p->list[i]->lpv) <= lmax)) {
      cand_swap(p->list[j], p->list[i]);
      j++;
    }
  }
  unsigned int newsize = MIN(j, max_size);
  for (i = newsize; i < p->n; ++i)
    cand_clear(p->list[i]);
  p->n = newsize;
  p->list = (cand_t *)realloc(p->list, (p->n)*sizeof(cand_t));
}

void pool_print(const pool_t p) {
  for (unsigned int i = 0; i < p->n; ++i) {
    printf("%u: %lu %lu (%u)\n", i,
        mpz_sizeinbase(p->list[i]->u, 2),
        mpz_sizeinbase(p->list[i]->v, 2),
        MAX(p->list[i]->lpu, p->list[i]->lpv));
  }
}

///////////////////////////////////////////////////////////////////////////
// Stats: keep track of statistics about the exepected number of bits
// obtained after running x curves.
///////////////////////////////////////////////////////////////////////////////

#define N 2048
typedef struct {
  double aver_gain[N];      // average number of bits removed after i curves
  unsigned long nb_test[N]; // size of the sample on which this avearge was done
} stats_s;
typedef stats_s stats_t[1];

void stats_init(stats_t S) {
  for(unsigned int i = 0; i < N; ++i) {
    S->aver_gain[i] = 0.0;
    S->nb_test[i] = 0;
  }
}

void stats_clear(stats_t S) { }

void stats_update(stats_t S, double gain, unsigned int i) {
  double newav = S->aver_gain[i]*S->nb_test[i] + gain;
  S->nb_test[i]++;
  newav /= S->nb_test[i];
  S->aver_gain[i] = newav;
}

///////////////////////////////////////////////////////////////////////////
// General data structure for an instance of the search of smooth numbers
///////////////////////////////////////////////////////////////////////////

typedef struct {
  pool_t pool;
  stats_t stats;
  void *param_next_cand;
  void (*next_cand)(cand_t, void *); // new candidate put in first arg.
  unsigned long target;              // smoothness bound (in bits)
  double current_effort;             // current effort per candidate.
  double max_effort;
  unsigned long max_pool_size;
} context_s;
typedef context_s context_t[1];

double remove_small_factors(mpz_t z) {
  double gain = 0.0;
  while(mpz_divisible_ui_p(z, 2)) {
    mpz_divexact_ui(z, z, 2);
    gain += 1.0;
  }
  while(mpz_divisible_ui_p(z, 3)) {
    mpz_divexact_ui(z, z, 3);
    gain += log2(3.0);
  }
  return gain;
}

// get a B1, so that we can quickly cover the target effort
double get_B1_from_effort(double effort) {
  double B1 = B1MIN;
  double S = B1;
  while (S < effort) {
    B1 += sqrt(B1);
    S += B1;
  }
  return B1;
}

void increase_effort(context_t ctx) {
  ctx->current_effort += sqrt(ctx->current_effort);
  ctx->current_effort = MIN(ctx->current_effort, ctx->max_effort);
}

void my_ecm_factor(mpz_t f, mpz_t z, double B1) {
  ecm_params ecm_par;
  ecm_init(ecm_par);
  long sig = random();
  mpz_set_ui(ecm_par->sigma, sig);
  ecm_par->B1done = default_B1done; /* issue with ECM 6.4.x */
  ecm_factor(f, z, B1, ecm_par);
  ecm_clear(ecm_par);
}

// One step of smoothness detection: get a new candidate, run a bunch of
// ECMs, update the pool, and the stats.
bool smooth_detect_one_step(cand_t winner, context_t ctx) {
  cand_t C;
  cand_init(C);
  // Get a new candidate, not obviously worse than current best.
  double gain_u = 0.0;
  double gain_v = 0.0;
  do {
    ctx->next_cand(C, ctx->param_next_cand);
    gain_u += remove_small_factors(C->u);
    gain_v += remove_small_factors(C->v);
    cand_update_check_prime(C);
  } while (C->lpu >= ctx->target || C->lpv >= ctx->target);

  // Start a loop of ECM
  double B1 = B1MIN;   // initial B1
  int cpt = 0;         // number of curves tried on this number
  mpz_t f;             // output of ecm
  mpz_init(f);
  while (!cand_is_probably_not_smooth(C, ctx->target)
      && !cand_is_factored(C)
      && C->effort < ctx->current_effort) {
    cpt++;
    // u-side
    if (mpz_cmp_ui(C->u, 1) != 0) {
      my_ecm_factor(f, C->u, B1);
      gain_u += log2(mpz_get_d(f));
      stats_update(ctx->stats, gain_u, cpt);
      if (mpz_cmp_ui(f, 1) > 0) {
        mpz_divexact(C->u, C->u, f);
        C->lpu = MAX(C->lpu, mpz_sizeinbase(f, 2));
        cand_update_check_prime(C);
        if (C->lpu >= ctx->target)
          break;
      }
    }
    // v-side
    if (mpz_cmp_ui(C->v, 1) != 0) {
      my_ecm_factor(f, C->v, B1);
      gain_v += log2(mpz_get_d(f));
      stats_update(ctx->stats, gain_v, cpt);
      if (mpz_cmp_ui(f, 1) > 0) {
        mpz_divexact(C->v, C->v, f);
        C->lpv = MAX(C->lpv, mpz_sizeinbase(f, 2));
        cand_update_check_prime(C);
        if (C->lpv >= ctx->target)
          break;
      }
    }
    // if both gain are below average, then abort this candidate
    if (gain_u < EA_THRESHOLD*ctx->stats->aver_gain[cpt] &&
        gain_v < EA_THRESHOLD*ctx->stats->aver_gain[cpt])
      break;
    // remember current effort for this number, and increase B1.
    C->effort += B1;
    B1 += sqrt(B1);
  }
  mpz_clear(f);
   
  // If we had a prime larger than expected, then abort.
  unsigned int l = MAX(C->lpu, C->lpv);
  if (l > ctx->target) {
    cand_clear(C);
    increase_effort(ctx);
    return false;
  }

  // Did we factor the candidate completely?
  // If so, print result, otherwise, insert in pool
  if (cand_is_factored(C)) {
    cand_print(C);
    cand_set(winner, C);
    cand_clear(C);
    return true;
  } else {
    pool_insert(ctx->pool, C);
  }

  B1 = get_B1_from_effort(ctx->current_effort);

  // Run ECM on the numbers in pool, so that they all have received more
  // or less the same effort.
  mpz_init(f);
  for (unsigned int i = 0; i < ctx->pool->n; ++i) {
    cand_s * c = &ctx->pool->list[i][0];
    double effort = ctx->current_effort;
    // more effort for the most promising candidates!
    if (ctx->pool->n > 5) {
        if (i == 0) { effort *= 2; }
        if (i == 1) { effort *= 1.6; }
        if (i == 2) { effort *= 1.3; }
    }
    while (!cand_is_factored(c) && c->effort < effort) {
      if (mpz_cmp_ui(c->u, 1) != 0) {
        my_ecm_factor(f, c->u, B1);
        if (mpz_cmp_ui(f, 1) > 0) {
          mpz_divexact(c->u, c->u, f);
          c->lpu = MAX(c->lpu, mpz_sizeinbase(f, 2));
        }
      }
      if (mpz_cmp_ui(c->v, 1) != 0) {
        my_ecm_factor(f, c->v, B1);
        if (mpz_cmp_ui(f, 1) > 0) {
          mpz_divexact(c->v, c->v, f);
          c->lpv = MAX(c->lpv, mpz_sizeinbase(f, 2));
        }
      }
      cand_update_check_prime(c);
      c->effort += B1;
      if (MAX(c->lpu, c->lpv) > ctx->target)
        break;
    } 
  }
  mpz_clear(f);

  // Go through pool: detect winners, purge bad candidates, keep best.
  for (unsigned int i = 0; i < ctx->pool->n; ++i) {
    cand_s * c = &ctx->pool->list[i][0];
    unsigned int l = MAX(c->lpu, c->lpv);
    if (cand_is_factored(c) && l <= ctx->target) {
      cand_print(c);
      cand_set(winner, c);
    }
  }
  pool_sort(ctx->pool);
  pool_purge(ctx->pool, ctx->max_pool_size, ctx->target);

  increase_effort(ctx);
  return false;
}

void smooth_detect(cand_t C, void (*next_cand)(cand_t, void *),
    void *param_next_cand, unsigned long target,
    const smooth_detect_param_s* param) {
  /* fix issue with ECM 6.4.x */
  {
    ecm_params params;
    ecm_init(params);
    default_B1done = params->B1done;
    ecm_clear(params);
  }

  const smooth_detect_param_s default_param = {2000.0, 1e20, 10};
  if (param == NULL) {
    param = &default_param;
  }

  // Create a context.
  context_t ctx;
  pool_init(ctx->pool);
  stats_init(ctx->stats);
  ctx->param_next_cand = param_next_cand;
  ctx->next_cand = next_cand;
  ctx->target = target;
  ctx->current_effort = param->min_effort;
  ctx->max_effort = param->max_effort;
  ctx->max_pool_size = param->max_pool_size;

  double tm = get_time();
  int cpt = 0;
  bool found = false;
  while (!found) {
    found = smooth_detect_one_step(C, ctx);
    cpt++;
    if (cpt % 20 == 0) {
      printf("***** Pool status after %d candidates in %.1fs\n", cpt,
          get_time()-tm);
      printf("current_effort = %.0f\n", ctx->current_effort);
      printf("current max B1 = %.0f\n",get_B1_from_effort(ctx->current_effort));
      pool_print(ctx->pool);
    }
  }
}

