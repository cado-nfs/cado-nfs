/* Experimental program to initialize the individual logarithm ("descent")
   in DLP.

   Input:
   ./descent-init [-seed s] [-target t] p z
   s is the random seed used (if not given, getpid() is used)
   t is the bit-size of wanted large primes (if not given, assume t=0)
   p is the prime defining the DLP group
   z is the target element

   Output:
   integers e, u0, v0 such that z^e = u0/v0 mod p and u0, v0 are smooth.

   Requires libecm.{a,so} in addition to -lgmp and -lm
   (tested with GMP-ECM 6.4.4).

   Example:
   $ p=53236943330228380237618624445646085674945074907141464418703
   $ z=13236943330228380237618624445646085674945074907141464418703
   $ ./descent-init -seed 1 -target 17 $p $z
   ...
   e=733057497
   u0=65757339031131586688817728662
   v0=364193180709513237520172194279
   u=35983 (16 bits)
   v=125707 (17 bits)
 */

#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <limits.h>
#include <unistd.h>     // getpid
#include <gmp.h>
#include "ecm.h"
#include "macros.h"
#include "misc.h"
#include "portability.h"

#define KEEP 10 /* minimal number of elements in pool */

static double default_B1done;

typedef struct {
  unsigned long e; /* exponent of z */
  mpz_t u0, v0;    /* z^e = u0/v0 mod p */
  mpz_t u, v;      /* unfactored parts of u, v */
  int uprime, vprime;
  double B1;       /* current B1 value */
  unsigned long l; /* potential size */
} cand_t;
typedef cand_t cand[1];

unsigned long
bsize (mpz_t u, int uprime)
{
  if (uprime)
    return mpz_sizeinbase (u, 2);
  else /* if u has n bits, we have 2^(n-1) <= u < 2^n: if u = p*q with p <= q,
          then q >= 2^((n-1)/2) and q has at least floor((n-1)/2)+1 bits */
    return (mpz_sizeinbase (u, 2) + 1) / 2;
}


unsigned long
csize (mpz_t u, int uprime)
{
  if (uprime)
    return mpz_sizeinbase (u, 2);
  else /* if u has n bits, we have 2^(n-1) <= u < 2^n: if u = p*q with p <= q,
          then q >= 2^((n-1)/2) and q has at least floor(((n-1)/2)+1 bits */
    return 0;
}

void
cand_init_set (cand c, unsigned long e, mpz_t u0, mpz_t v0, mpz_t u, mpz_t v,
               int uprime, int vprime, double B1, unsigned long l)
{
  c->e = e;
  mpz_init_set (c->u0, u0);
  mpz_init_set (c->v0, v0);
  mpz_init_set (c->u, u);
  mpz_init_set (c->v, v);
  c->uprime = uprime;
  c->vprime = vprime;
  c->B1 = B1;
  c->l = l;
}

void
cand_clear (cand c)
{
  mpz_clear (c->u0);
  mpz_clear (c->v0);
  mpz_clear (c->u);
  mpz_clear (c->v);
}

void
cand_print_raw (unsigned long e, mpz_t u0, mpz_t v0, mpz_t u, mpz_t v)
{
  printf ("e=%lu\n", e);
  gmp_printf ("u0=%Zd\n", u0);
  gmp_printf ("v0=%Zd\n", v0);
  gmp_printf ("u=%Zd (%lu bits)\n", u, mpz_sizeinbase (u, 2));
  gmp_printf ("v=%Zd (%lu bits)\n", v, mpz_sizeinbase (v, 2));
  fflush (stdout);
}

void
cand_print (cand c)
{
  cand_print_raw (c->e, c->u0, c->v0, c->u, c->v);
}

void
cand_info (cand c)
{
  printf ("u:%c%lu v:%c%lu (B1=%.0f)\n",
          (c->uprime) ? 'p' : 'c', mpz_sizeinbase (c->u, 2),
          (c->vprime) ? 'p' : 'c', mpz_sizeinbase (c->v, 2), c->B1);
}

int
cand_factored (cand c)
{
  return c->uprime && c->vprime;
}

void
params_init (ecm_params params, double B1)
{
  mpz_set_d (params->sigma, B1);
  params->B1done = default_B1done; /* issue with ECM 6.4.x */
}

void
print_params (ecm_params params)
{
  printf ("method=%d\n", params->method);
  gmp_printf ("x=%Zd\n", params->x);
  gmp_printf ("sigma=%Zd\n", params->sigma);
  printf ("sigma_is_A=%d\n", params->sigma_is_A);
  gmp_printf ("go=%Zd\n", params->go);
  printf ("B1done=%.16e\n", params->B1done);
  gmp_printf ("B2min=%Zd\n", params->B2min);
  gmp_printf ("B2=%Zd\n", params->B2);
  gmp_printf ("k=%lu\n", params->k);
  printf ("S=%d\n", params->S);
  printf ("repr=%d\n", params->repr);
  printf ("nobase2step2=%d\n", params->nobase2step2);
  printf ("verbose=%d\n", params->verbose);
  printf ("maxmem=%.16e\n", params->maxmem);
  printf ("use_ntt=%d\n", params->use_ntt);
  printf ("batch=%d\n", params->batch);
}

void
cand_scan (cand c, ecm_params params, double *S2)
{
  mpz_t f;
  unsigned long lu, lv;

  mpz_init (f);
  if (c->uprime == 0)
    {
      params_init (params, c->B1);
      *S2 += c->B1;
      ecm_factor (f, c->u, c->B1, params);
      if (mpz_cmp_ui (f, 1) > 0)
        {
          mpz_divexact (c->u, c->u, f);
          c->uprime = mpz_probab_prime_p (c->u, 1);
          lu = bsize (c->u, c->uprime);
          lv = bsize (c->v, c->vprime);
          c->l = (lu > lv) ? lu : lv;
        }
    }
  if (c->vprime == 0)
    {
      params_init (params, c->B1);
      *S2 += c->B1;
      ecm_factor (f, c->v, c->B1, params);
      if (mpz_cmp_ui (f, 1) > 0)
        {
          mpz_divexact (c->v, c->v, f);
          c->vprime = mpz_probab_prime_p (c->v, 1);
          lu = bsize (c->u, c->uprime);
          lv = bsize (c->v, c->vprime);
          c->l = (lu > lv) ? lu : lv;
        }
    }
  c->B1 = c->B1 + sqrt (c->B1);
  mpz_clear (f);
}

void
cand_set (cand c, cand d)
{
  c->e = d->e;
  mpz_swap (c->u0, d->u0);
  mpz_swap (c->v0, d->v0);
  mpz_swap (c->u, d->u);
  mpz_swap (c->v, d->v);
  c->uprime = d->uprime;
  c->vprime = d->vprime;
  c->B1 = d->B1;
  c->l = d->l;
}

void
cand_swap (cand c, cand d)
{
  unsigned long t;
  int i;
  double s;
  t = c->e; c->e = d->e; d->e = t;
  mpz_swap (c->u0, d->u0);
  mpz_swap (c->v0, d->v0);
  mpz_swap (c->u, d->u);
  mpz_swap (c->v, d->v);
  i = c->uprime; c->uprime = d->uprime; d->uprime = i;
  i = c->vprime; c->vprime = d->vprime; d->vprime = i;
  s = c->B1; c->B1 = d->B1; d->B1 = s;
  t = c->l; c->l = d->l; d->l = t;
}

int
cost_uv (mpz_t u, mpz_t v, int uprime, int vprime)
{
  if (uprime == 0)
    {
      if (vprime == 0)
        {
          return (mpz_cmpabs (u, v) < 0) ? mpz_sizeinbase (u, 2)
            : mpz_sizeinbase (v, 2);
        }
      else
        return mpz_sizeinbase (u, 2);
    }
  else /* u is prime */
    {
      if (vprime)
        return INT_MAX;
      else
        return mpz_sizeinbase (v, 2);
    }
}

/* return the bit-size of the smallest composite between u and v,
   and INT_MAX if both are prime */
int
cost (cand c)
{
  return cost_uv (c->u, c->v, c->uprime, c->vprime);
}

/*****************************************************************************/

typedef struct {
  cand *list;
  unsigned long n;
} pool_t;
typedef pool_t pool[1];

void
pool_init (pool p)
{
  p->list = NULL;
  p->n = 0;
}

void
pool_clear (pool p)
{
  unsigned long i;

  for (i = 0; i < p->n; i++)
    cand_clear (p->list[i]);
  free (p->list);
}

void
pool_add (pool p, unsigned long e, mpz_t u0, mpz_t v0, mpz_t u, mpz_t v,
          int uprime, int vprime, double B1, unsigned long l)
{
  unsigned int i;

  p->list = realloc (p->list, (p->n + 1) * sizeof (cand));
  cand_init_set (p->list[p->n], e, u0, v0, u, v, uprime, vprime, B1, l);
  for (i = p->n; i > 0 && cost (p->list[i-1]) > cost (p->list[i]); i--)
    cand_swap (p->list[i-1], p->list[i]);
  p->n += 1;
}

/* perform an ECM curve for each candidate in the pool */
unsigned long
pool_scan (pool p, ecm_params params, unsigned long L, double S1, double *S2)
{
  unsigned long i, j, lu, lv, l, minl = 0, oldn = p->n;

  /* compute minl */
  for (i = 0; i < p->n; i++)
    {
      lu = mpz_sizeinbase (p->list[i]->u, 2);
      lv = mpz_sizeinbase (p->list[i]->v, 2);
      l = (lu > lv) ? lu : lv;
      if (i == 0 || l < minl)
        minl = l;
    }

  for (i = 0; i < p->n && *S2 < S1; i++)
    {
      /* we focus on the best candidate first */
      do {
        cand_scan (p->list[i], params, S2);
      } while (cand_factored (p->list[i]) == 0 && *S2 < S1);
    }


  /* search for better solutions */
  for (i = 0; i < p->n; i++)
    if (cand_factored (p->list[i]))
      {
        if (p->list[i]->l < L) /* found better solution */
          {
            cand_print (p->list[i]);
            L = p->list[i]->l;
          }
      }

  /* purge candidates */
  for (i = j = 0; i < p->n; i++)
    if (cand_factored (p->list[i]) || p->list[i]->l >= L)
      {
      }
    else if (p->list[i]->uprime && mpz_sizeinbase (p->list[i]->u, 2) >= minl)
      { /* if there is a prime of >= minl bits, it cannot win */
      }
    else if (p->list[i]->vprime && mpz_sizeinbase (p->list[i]->v, 2) >= minl)
      { /* if there is a prime of >= minl bits, it cannot win */
      }
    else
      cand_set (p->list[j++], p->list[i]);

  /* sort remaining candidates */
  p->n = j;
  for (i = 1; i < p->n; i++)
    for (j = i; j > 0 && cost (p->list[j-1]) > cost (p->list[j]); j--)
      cand_swap (p->list[j-1], p->list[j]);

#if 0
  if (j < oldn)
    {
      unsigned long imin, imax;
      p->list = realloc (p->list, p->n * sizeof (cand));
      printf ("Remains %lu candidate(s) in pool", p->n);
      if (p->n == 0)
        printf ("\n");
      else
        {
          printf (", first ");
          cand_info (p->list[0]);
        }
      fflush (stdout);
    }
#endif
  return L;
}

/*****************************************************************************/

int
cputime ()
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

void
HalfGcd (mpz_t a, mpz_t b, mpz_t u)
{
  mpz_t v, w, x, q, r;
  unsigned long n;

  n = (mpz_cmpabs (a, b) > 0) ? mpz_sizeinbase (a, 2) : mpz_sizeinbase (b, 2);
  mpz_set_ui (u, 1);
  mpz_init_set_ui (w, 0);
  mpz_init_set_ui (v, 0);
  mpz_init_set_ui (x, 1);
  mpz_init (q);
  mpz_init (r);
  /* invariant: a = u*a0 + v*b0 */
  while (mpz_cmpabs (a, u) > 0)
    {
      mpz_tdiv_qr (q, r, a, b);
      mpz_swap (a, b);
      mpz_swap (b, r);
      mpz_submul (u, q, w);
      mpz_swap (u, w);
      mpz_submul (v, q, x);
      mpz_swap (v, x);
    }
  mpz_clear (v);
  mpz_clear (w);
  mpz_clear (x);
  mpz_clear (q);
  mpz_clear (r);
}

#define N 2048
double aver_gain[N], gain_u[N], gain_v[N];
unsigned long nb_test[N];

int
main (int argc, char *argv[])
{
  mpz_t p, z, ze, t, u, v, f, u0, v0;
  unsigned long e, L, l, lu, lv, i, seed = 0, target = 0;
  double B1, S, Smax = 1000.0;
  int uprime, vprime;
  ecm_params params;
  pool P;
  double S1 = 0.0; /* total ECM effort spent in first stage */
  double S2 = 0.0; /* total ECM effort spent in second stage */

  while (argc > 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-seed") == 0)
        {
          seed = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (strcmp (argv[1], "-target") == 0)
        {
          target = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  ASSERT (argc == 3);

  for (i = 0; i < N; i++)
    {
      aver_gain[i] = 0.0;
      nb_test[i] = 0;
    }

  /* fix issue with ECM 6.4.x */
  ecm_init (params);
  default_B1done = params->B1done;
  ecm_clear (params);

  pool_init (P);
  mpz_init_set_str (p, argv[1], 10);
  mpz_init_set_str (z, argv[2], 10);
  gmp_printf ("p=%Zd\n", p);
  gmp_printf ("z=%Zd\n", z);
  if (seed == 0)
    seed = getpid ();
  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);
  gmp_randseed_ui(rstate, seed);
  printf ("Using seed %lu\n", seed);
  mpz_init (ze);
  mpz_init (t);
  mpz_init (u);
  mpz_init (v);
  mpz_init (u0);
  mpz_init (v0);
  mpz_init (f);
  L = mpz_sizeinbase (p, 2);
  while (L > target)
    {
      e = u64_random (rstate);
      mpz_powm_ui (u, z, e, p);
      mpz_set (t, p);
      HalfGcd (u, t, v);
      mpz_set (u0, u);
      mpz_set (v0, v);
      /* u = v*z^e mod p thus z^e = u/v mod p */
      
      /* remove factors of 2 and 3 */
      gain_u[0] = gain_v[0] = 0;
      while (mpz_divisible_ui_p (u, 2))
        {
          mpz_divexact_ui (u, u, 2);
          gain_u[0] += 1.0;
        }
      while (mpz_divisible_ui_p (u, 3))
        {
          mpz_divexact_ui (u, u, 3);
          gain_u[0] += log2 (3.0);
        }
      while (mpz_divisible_ui_p (v, 2))
        {
          mpz_divexact_ui (v, v, 2);
          gain_v[0] += 1.0;
        }
      while (mpz_divisible_ui_p (v, 3))
        {
          mpz_divexact_ui (v, v, 3);
          gain_v[0] += log2 (3.0);
        }
      B1 = 100.0;
      S = 0.0;
      unsigned int i = 0;
      uprime = mpz_probab_prime_p (u, 1);
      if (uprime && mpz_sizeinbase (u, 2) >= L)
        continue;
      vprime = mpz_probab_prime_p (v, 1);
      if (vprime && mpz_sizeinbase (v, 2) >= L)
        continue;
      ecm_init (params);
      while (S < Smax && (uprime == 0 || vprime == 0))
        {
          i = i + 1;
	  ASSERT (i < N);
          gain_u[i] = gain_u[i-1];
          gain_v[i] = gain_v[i-1];
          if (uprime == 0)
            {
              params_init (params, B1);
              S1 += B1;
              ecm_factor (f, u, B1, params);
	      gain_u[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_u[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (u, u, f);
                  uprime = mpz_probab_prime_p (u, 1);
                  if (uprime && mpz_sizeinbase (u, 2) >= L)
                    break;
                }
            }

          if (vprime == 0)
            {
              params_init (params, B1);
              S1 += B1;
              ecm_factor (f, v, B1, params);
	      gain_v[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_v[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (v, v, f);
                  vprime = mpz_probab_prime_p (v, 1);
                  if (vprime && mpz_sizeinbase (v, 2) >= L)
                    break;
                }
            }

          if (gain_u[i] < aver_gain[i] && gain_v[i] < aver_gain[i])
            break;
          B1 += sqrt (B1);
          S += B1;
        }
      lu = bsize (u, uprime);
      lv = bsize (v, vprime);
      l = (lu > lv) ? lu : lv;
      if (uprime && vprime)
        {
          if (l < L) /* found better relation */
            {
              cand_print_raw (e, u0, v0, u, v);
              L = l;
            }
        }
      else if (csize (u, uprime) < L && csize (v, vprime) < L &&
               (P->n < KEEP || cost_uv (u, v, uprime, vprime) < 
                cost (P->list[P->n - 1])))
        {
          pool_add (P, e, u0, v0, u, v, uprime, vprime, B1, l);
          L = pool_scan (P, params, L, S1, &S2);
        }
      Smax += sqrt (Smax);
    }
  gmp_randclear(rstate);
  mpz_clear (p);
  mpz_clear (z);
  mpz_clear (ze);
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (u0);
  mpz_clear (v0);
  mpz_clear (f);
  ecm_clear (params);
  pool_clear (P);
}
