#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "polyselect_proots.h"
#include "polyselect_thread.h"
#include "gcd.h"
#include "roots_mod.h"

/* init polyselect_proots_t */
void
polyselect_proots_init (polyselect_proots_ptr R,
              int d,
              unsigned long size )
{
  R->size = size;

  /* length of nr&roots is known now. lengths of roots[i] are TBD. */
  /* +1 for R->nr for end guard in collision_on_each_sq */
  R->nr = (uint8_t *) malloc ((size + 1) * sizeof (*(R->nr)));
  R->roots = (uint64_t **) malloc (size * sizeof (*(R->roots)));

  if (R->nr == NULL || R->roots == NULL) {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit(1);
  }
  for (unsigned int i = 0; i < R->size; i++) {
      R->nr[i] = 0;
      R->roots[i] = (uint64_t *) malloc(d * sizeof(uint64_t));
  }
}


/* add a root to polyselect_proots_t */
void
polyselect_proots_add ( polyselect_proots_ptr R,
             unsigned long nr,
             uint64_t *roots,
             unsigned long index )
{
  unsigned int i;
  R->nr[index] = nr;

  for (i = 0; i < nr; i++)
      R->roots[index][i] = roots[i];
}


/* print roots */
void
polyselect_proots_print (polyselect_proots_srcptr R)
{
  unsigned int i, j;
  for (i = 0; i < R->size; i++) {
    if (R->nr[i] == 0) {
      fprintf (stderr, "NULL\n");
    }
    else {
      for (j = 0; j < R->nr[i]; j ++)
        fprintf (stderr, "%" PRIu64 " ", R->roots[i][j]);
      fprintf (stderr, "\n");
    }
  }
}


/* clear roots */
void
polyselect_proots_clear (polyselect_proots_ptr R)
{
  unsigned int i;

  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}

/* Lift the n roots r[0..n-1] of N = x^d (mod p) to roots of
   N = (m0 + r)^d (mod p^2).
   Return the number of lifted roots (might be less than n if some root is 0).
*/
unsigned long
roots_lift (uint64_t *r, mpz_srcptr N, unsigned long d, mpz_srcptr m0,
            unsigned long p, unsigned long n)
{
  uint64_t pp;
  unsigned long i, j, inv;
  mpz_t tmp, lambda;
  mpz_init (tmp);
  mpz_init (lambda);
  pp = (uint64_t) p;
  pp *= (uint64_t) p;

  if (sizeof (unsigned long) == 8) {
    for (i = j = 0; j < n; j++) {
        if (r[j] == 0)
           continue;
	/* we have for r=r[j]: r^d = N (mod p), lift mod p^2:
	   (r+lambda*p)^d = N (mod p^2) implies
	   r^d + d*lambda*p*r^(d-1) = N (mod p^2)
           lambda = (N - r^d)/(p*d*r^(d-1)) mod p */
	mpz_ui_pow_ui (tmp, r[j], d - 1);
	mpz_mul_ui (lambda, tmp, r[j]);    /* lambda = r^d */
	mpz_sub (lambda, N, lambda);
	mpz_divexact_ui (lambda, lambda, p);
	mpz_mul_ui (tmp, tmp, d);         /* tmp = d*r^(d-1) */
	inv = invert_ul (mpz_fdiv_ui (tmp, p), p);
	mpz_mul_ui (lambda, lambda, inv * p); /* inv * p fits in 64 bits if
						 p < 2^32 */
	mpz_add_ui (lambda, lambda, r[j]); /* now lambda^d = N (mod p^2) */

	/* subtract m0 to get roots of (m0+r)^d = N (mod p^2) */
	mpz_sub (lambda, lambda, m0);
	r[i++] = mpz_fdiv_ui (lambda, pp);
      }
  }
  else {
#if 0   
    printf ("p: %lu, ppl %" PRId64 ": ", p, pp);
#endif
    uint64_t tmp1;
    mpz_t ppz, *rz, tmpz;
    rz = (mpz_t*) malloc (n * sizeof (mpz_t));
    mpz_init (ppz);
    mpz_init (tmpz);
    for (j = 0; j < n; j++) {
      mpz_init_set_ui (rz[j], 0UL);
      mpz_set_uint64 (rz[j], r[j]);
#if 0   
      printf (" %" PRIu64 "", r[j]);
#endif
    }

    for (i = j = 0; j < n; j++) {
	mpz_pow_ui (tmp, rz[j], d - 1);
	mpz_mul (lambda, tmp, rz[j]);    /* lambda = r^d */
	mpz_sub (lambda, N, lambda);
	mpz_divexact_ui (lambda, lambda, p);
	mpz_mul_ui (tmp, tmp, d);         /* tmp = d*r^(d-1) */
	inv = invert_ul (mpz_fdiv_ui (tmp, p), p);
	tmp1 = (uint64_t) inv;
	tmp1 *= (uint64_t) p;
	mpz_set_uint64 (tmpz, tmp1);
	mpz_mul (lambda, lambda, tmpz); 
	mpz_add (lambda, lambda, rz[j]); /* now lambda^d = N (mod p^2) */
	/* subtract m0 to get roots of (m0+r)^d = N (mod p^2) */
	mpz_sub (lambda, lambda, m0);
	mpz_set_uint64 (tmpz, pp);
	mpz_fdiv_r (rz[j], lambda, tmpz);
	r[i++] = mpz_get_uint64 (rz[j]);
      }

    for (j = 0; j < n; j++)
      mpz_clear (rz[j]);
    free (rz);
    mpz_clear (ppz);
    mpz_clear (tmpz);
  }

  mpz_clear (tmp);
  mpz_clear (lambda);
  return i;
}

void polyselect_proots_compute_subtask(polyselect_thread_ptr thread)/*{{{*/
{
    polyselect_primes_table_srcptr pt = thread->team->league->pt;
    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter proots");

    unsigned long tot_roots = 0;
    size_t qt = pt->lenPrimes / nt;
    size_t rt = pt->lenPrimes % nt;
    unsigned long i0 = qt * it + MIN(it, rt);
    unsigned long i1 = i0 + qt + (it < rt);

#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d (%d-th in sync group) enters proots with %d threads [%lu .. %lu)\n", thread->thread_index, it, nt, i0, i1);
#endif

    uint64_t * rp = (uint64_t *) malloc(thread->team->header->d * sizeof(uint64_t));
    FATAL_ERROR_CHECK(!rp, "cannot allocate memory");

    for (unsigned long i = i0; i < i1; i++) {
        unsigned long p = pt->Primes[i];

        /* add fake roots to keep indices */
        if (polyselect_poly_header_skip(thread->team->header, p))
        {
            thread->team->R->nr[i] = 0;	// nr = 0.
            // well, no, don't touch that! It's used later on. And in the
            // non-flat model, it's definitely pre-allocated!
            // thread->team->R->roots[i] = NULL;
            continue;
        }

        unsigned long nrp = roots_mod_uint64(rp,
                mpz_fdiv_ui(thread->team->header->Ntilde, p),
                thread->team->header->d,
                p, thread->rstate);
        tot_roots += nrp;
        nrp = roots_lift(rp,
                thread->team->header->Ntilde,
                thread->team->header->d,
                thread->team->header->m0, p, nrp);
        polyselect_proots_add(thread->team->R, nrp, rp, i);
    }

    polyselect_thread_chronogram_chat(thread, "leave proots");
    free(rp);
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);
    *((unsigned long*)thread->team->task->arg) += tot_roots;
}/*}}}*/

unsigned long polyselect_proots_compute_conductor(polyselect_thread_ptr thread)
{
    /* This is called with the team lock held ! */
    unsigned long tot_roots = 0;

    polyselect_thread_team_post_work(thread->team, thread, polyselect_proots_compute_subtask, &tot_roots);

    return tot_roots;
}/*}}}*/


