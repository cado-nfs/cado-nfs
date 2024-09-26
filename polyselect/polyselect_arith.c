#include "cado.h" // IWYU pragma: keep
#include <stdio.h> // fprintf
#include <stdlib.h>     // exit free malloc
#include <gmp.h>
#include "polyselect_arith.h"
#include "polyselect_special_q.h"    // SPECIAL_Q...
#include "roots_mod.h"
#include "gcd.h"       // for invert_ul
#include "gmp_aux.h"       // mpz_set_uint64


/* first combination of k elements among 0, ..., n-1: 0, 1, 2, 3, \cdots */
void
first_comb (unsigned long k, unsigned long *r)
{
  unsigned long i;
  for (i = 0; i < k; ++i)
    r[i] = i;
}


/* next combination of k elements among 0, 1, ..., n-1,
   return the index of the first increased element (k if finished) */
unsigned long
next_comb ( unsigned long n,
            unsigned long k,
            unsigned long *r )
{
  unsigned long j;

  /* if the last combination */
  if (r[0] == n - k) /* we have n-k, n-k+1, ..., n-1 */
    return k;

  /* if r[k-1] is not equal to n-1, just increase it */
  j = k - 1;
  if (r[j] < n - 1) {
    r[j] ++;
    return j;
  }

  /* find which one we should increase */
  while ( r[j] - r[j-1] == 1)
    j --;

  unsigned long ret = j - 1;
  unsigned long z = ++r[j-1];

  while (j < k) {
    r[j] = ++z;
    j ++;
  }
  return ret;
}


/* debug */
void
print_comb ( unsigned long k,
             unsigned long *r )
{
  unsigned long i;
  for (i = 0; i < k; i ++)
    fprintf (stderr, "%lu ", r[i]);
  fprintf (stderr, "\n");
}


/* return number of n choose k
 *
 * [unused anywhere, it seems]
 */
unsigned long
binomial ( unsigned long n,
        unsigned long k )
{
  if (k > n)
    return 0;
  if (k == 0 || k == n)
    return 1;
  if (2*k > n)
    k = n - k;

  unsigned long tot = n - k + 1, f = tot, i;
  for (i = 2; i <= k; i++) {
    f ++;
    tot *= f;
    tot /= i;
  }
  return tot;
}

/* prepare special-q's roots */
void
comp_sq_roots ( polyselect_poly_header_srcptr header,
                polyselect_qroots_ptr SQ_R,
                gmp_randstate_ptr rstate
                )
{
  unsigned long i, q, nrq;
  uint64_t *rq;

  rq = (uint64_t*) malloc (header->d * sizeof (uint64_t));
  if (rq == NULL)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit(1);
    }

  /* prepare the special-q's */
  for (i = 1; (q = SPECIAL_Q[i]) != 0 ; i++)
  {
    if (polyselect_poly_header_skip (header, q))
      continue;

    if ( mpz_fdiv_ui (header->Ntilde, q) == 0 )
      continue;

    nrq = roots_mod_uint64 (rq, mpz_fdiv_ui (header->Ntilde, q), header->d, q, rstate);
    roots_lift (rq, header->Ntilde, header->d, header->m0, q, nrq);

#ifdef DEBUG_POLYSELECT
    unsigned int j = 0;
    mpz_t r1, r2;
    mpz_init (r1);
    mpz_init (r2);
    mpz_set (r1, header->Ntilde);
    mpz_mod_ui (r1, r1, q*q);

    gmp_fprintf (stderr, "Ntilde: %Zd, Ntilde (mod %u): %Zd\n", 
                header->Ntilde, q*q, r1);

    for (j = 0; j < nrq; j ++) {
      mpz_set (r2, header->m0);

      mpz_add_ui (r2, r2, rq[j]);
      mpz_pow_ui (r2, r2, header->d);
      mpz_mod_ui (r2, r2, q*q);

      if (mpz_cmp (r1, r2) != 0) {
        fprintf (stderr, "Root computation wrong in comp_sq_roots().\n");
        fprintf (stderr, "q: %lu, rq: %lu\n", q, rq[j]);
        exit(1);
      }
    }

    mpz_clear (r1);
    mpz_clear (r2);
#endif

    polyselect_qroots_add (SQ_R, q, nrq, rq);
  }

  /* Reorder R entries by decreasing number of roots (nr).
     It is safe to comment it. */
  polyselect_qroots_rearrange (SQ_R);

  free(rq);
  polyselect_qroots_realloc (SQ_R, SQ_R->size); /* free unused space */
}

/* return the maximal number of special-q's with k elements among lq
 *
 * This is the same as the degree k coefficient of the product
 * \prod_{i=1}^{lq} (1-a_i x)
 * with a_i = SQ_R->nr[i]
 * but it is slightly unsatisfactory that we're apparently unable to
 * compute the result in less time than O(binomial(n,k)*k)... (well, to
 * be honest, it's not a big source of trouble either).
 */
unsigned long
number_comb (polyselect_qroots_srcptr SQ_R, unsigned long k, unsigned long lq)
{
  unsigned long s = 0;
  unsigned long idx[k], j;

  first_comb (k, idx);
  while (1)
    {
      unsigned long p = 1;
      for (j = 0; j < k; j++)
        p *= SQ_R->nr[idx[j]];
      s += p;
      if (next_comb (lq, k, idx) == k)
        break;
    }
  return s;
}

/* This does the reconstruction from a set of residues. An integer r is
 * implicitly given by its residues rq[i] modulo q[i]^2 for a set of
 * primes p in q[0]...q[len-1], and this routines computes r, as well as
 * the product of the q[i]^2's.
 *
 * Note that the q[i] must not be larger than half an unsigned long.
 */
void
crt_sq(mpz_ptr qqz,
       mpz_ptr r, unsigned long *q, unsigned long *rq, unsigned long lq)
{
  mpz_t prod, pprod, mod, inv, sum;
  unsigned long qq[lq];

  mpz_init_set_ui(prod, 1);
  mpz_init(pprod);
  mpz_init(mod);
  mpz_init(inv);
  mpz_init_set_ui(sum, 0);

  for (unsigned long i = 0; i < lq; i++)
    {
      qq[i] = q[i] * q[i];	// q small
      mpz_mul_ui(prod, prod, qq[i]);
    }

  for (unsigned long i = 0; i < lq; i++)
    {
      mpz_divexact_ui(pprod, prod, qq[i]);
      mpz_set_ui(mod, qq[i]);
      mpz_invert(inv, pprod, mod);
      mpz_mul_ui(inv, inv, rq[i]);
      mpz_mul(inv, inv, pprod);
      mpz_add(sum, sum, inv);
    }

  mpz_mod(sum, sum, prod);
  mpz_set(r, sum);
  mpz_set(qqz, prod);

  mpz_clear(prod);
  mpz_clear(pprod);
  mpz_clear(mod);
  mpz_clear(inv);
  mpz_clear(sum);
}

/* given individual q's, return \product q, no rq
 *
 * In plain English: returns the products of the k prime numbers that are
 * indexed by idx_q, and picked from the list of primes that are present
 * in polyselect_qroots_srcptr
 *
 * XXX This function belongs to polyselect_qroots.[ch], IMHO
 */
uint64_t
return_q_norq (polyselect_qroots_srcptr SQ_R, unsigned long *idx_q, unsigned long k)
{
  unsigned long i;
  uint64_t q = 1;

  for (i = 0; i < k; i ++)
    q = q * SQ_R->q[idx_q[i]];
  return q;
}

