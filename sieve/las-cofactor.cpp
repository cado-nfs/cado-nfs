#include "cado.h"
#include <pthread.h>
#include <vector>
#include <utility>

#include "utils.h"
#include "las-cofactor.hpp"
#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "las-types.hpp"

// {{{ las_info::{init,clear,print}_cof_stats
cofactorization_statistics::cofactorization_statistics(param_list_ptr pl)
{
    const char * statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (!statsfilename) {
        file = NULL;
        return;
    }
    file = fopen (statsfilename, "w");
    if (file == NULL) {
        fprintf (stderr, "Error, cannot create file %s\n", statsfilename);
        exit (EXIT_FAILURE);
    }
}

void cofactorization_statistics::call(int bits0, int bits1)
{
    if (!file) return;
    static pthread_rwlock_t lock = PTHREAD_RWLOCK_INITIALIZER;
    pthread_rwlock_rdlock(&lock);
    size_t s0 = cof_call.size();
    pthread_rwlock_unlock(&lock);
    if ((size_t) bits0 >= s0) {
        pthread_rwlock_wrlock(&lock);
        size_t news0 = std::max((size_t) bits0+1, s0 + s0/2);
        cof_call.insert(cof_call.end(), news0-s0, std::vector<uint32_t>());
        cof_success.insert(cof_success.end(), news0-s0, std::vector<uint32_t>());
        s0 = news0;
        pthread_rwlock_unlock(&lock);
    }
    pthread_rwlock_rdlock(&lock);
    size_t s1 = cof_call[bits0].size();
    pthread_rwlock_unlock(&lock);
    if ((size_t) bits1 >= s1) {
        pthread_rwlock_wrlock(&lock);
        size_t news1 = std::max((size_t) bits1+1, s1 + s1/2);
        cof_call[bits0].insert(cof_call[bits0].end(), news1-s1, 0);
        cof_success[bits0].insert(cof_success[bits0].end(), news1-s1, 0);
        s1 = news1;
        pthread_rwlock_unlock(&lock);
    }
    /* no need to use a mutex here: either we use one thread only
       to compute the cofactorization data and if several threads
       the order is irrelevant. The only problem that can happen
       is when two threads increase the value at the same time,
       and it is increased by 1 instead of 2, but this should
       happen rarely. */
    cof_call[s0][s1]++;
}

void cofactorization_statistics::print()
{
    if (!file) return;
    for(size_t bits0 = 0 ; bits0 < cof_call.size() ; ++bits0) {
        for(size_t bits1 = 0 ; bits1 < cof_call[bits0].size() ; ++bits1) {
            fprintf (file, "%zu %zu %" PRIu32 " %" PRIu32 "\n",
                    bits0, bits1,
                    cof_call[bits0][bits1],
                    cof_success[bits0][bits1]);
        }
    }
}

cofactorization_statistics::~cofactorization_statistics()
{
    if (!file) return;
    fclose (file);
}
//}}}

/* {{{ factor_leftover_norm */

#define NMILLER_RABIN 1 /* in the worst case, what can happen is that a
                           composite number is declared as prime, thus
                           a relation might be missed, but this will not
                           affect correctness */
#define IS_PROBAB_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
#define NFACTORS        8 /* maximal number of large primes */

/************************ cofactorization ********************************/

/* {{{ cofactoring area */

/* Return 0 if the leftover norm n cannot yield a relation.

   Possible cases, where qj represents a prime in [B,L], and rj a prime > L
   (assuming L < B^2, which might be false for the DLP descent):
   (0) n >= 2^mfb
   (a) n < L:           1 or q1
   (b) L < n < B^2:     r1 -> cannot yield a relation
   (c) B^2 < n < B*L:   r1 or q1*q2
   (d) B*L < n < L^2:   r1 or q1*q2 or q1*r2
   (e) L^2 < n < B^3:   r1 or q1*r2 or r1*r2 -> cannot yield a relation
   (f) B^3 < n < B^2*L: r1 or q1*r2 or r1*r2 or q1*q2*q3
   (g) B^2*L < n < L^3: r1 or q1*r2 or r1*r2
   (h) L^3 < n < B^4:   r1 or q1*r2, r1*r2 or q1*q2*r3 or q1*r2*r3 or r1*r2*r3
                        -> cannot yield a relation
*/
int
check_leftover_norm (cxx_mpz const & n, siever_config::side_config const & scs)
{
  size_t s = mpz_sizeinbase (n, 2);
  unsigned int lpb = scs.lpb;
  unsigned int mfb = scs.mfb;
  unsigned int klpb;
  double nd, kB, B;

  if (s > mfb)
    return 0; /* n has more than mfb bits, which is the given limit */

  if (scs.lim == 0) {
      /* special case when not sieving */
      return 1;
  }

  if (s <= lpb)
    return 1; /* case (a) */
  /* Note that in the case where L > B^2, if we're below L it's still fine of
     course, but we have no guarantee that our cofactor is prime... */

  nd = mpz_get_d (n);
  B = (double) scs.lim;
  kB = B * B;
  for (klpb = lpb; klpb < s; klpb += lpb, kB *= B)
    {
      /* invariant: klpb = k * lpb, kB = B^(k+1) */
      if (nd < kB) /* L^k < n < B^(k+1) */
	return 0;
    }

  /* Here we have L < n < 2^mfb. If n is composite and we wrongly consider
     it prime, we'll return 0, thus we'll potentially miss a relation, but
     we won't output a relation with a composite ideal, thus a base-2 strong
     prime test is enough. */

  // TODO: maybe we should pass the modulus to the facul machinery
  // instead of reconstructing it.
  int prime=0;
  if (s <= MODREDCUL_MAXBITS) {
      modulusredcul_t m;
      ASSERT(mpz_fits_ulong_p(n));
      modredcul_initmod_ul (m, mpz_get_ui(n));
      prime = modredcul_sprp2(m);
      modredcul_clearmod (m);
  } else if (s <= MODREDC15UL_MAXBITS) {
      modulusredc15ul_t m;
      unsigned long t[2];
      modintredc15ul_t nn;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, n);
      ASSERT_ALWAYS(written <= 2);
      modredc15ul_intset_uls (nn, t, written);
      modredc15ul_initmod_int (m, nn);
      prime = modredc15ul_sprp2(m);
      modredc15ul_clearmod (m);
  } else if (s <= MODREDC2UL2_MAXBITS) {
      modulusredc2ul2_t m;
      unsigned long t[2];
      modintredc2ul2_t nn;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, n);
      ASSERT_ALWAYS(written <= 2);
      modredc2ul2_intset_uls (nn, t, written);
      modredc2ul2_initmod_int (m, nn);
      prime = modredc2ul2_sprp2(m);
      modredc2ul2_clearmod (m);
  } else {
      prime = mpz_probab_prime_p (n, 1);
  }
  if (prime)
    return 0; /* n is a pseudo-prime larger than L */
  return 1;
}

/* This is the header-comment for the old factor_leftover_norm()
 * function, that is now deleted */
/* This function was contributed by Jerome Milan (and bugs were introduced
   by Paul Zimmermann :-).
   Input: n - the number to be factored (leftover norm). Must be composite!
              Assumed to have no factor < B (factor base bound).
          L - large prime bound is L=2^l
   Assumes n > 0.
   Return value:
          -1 if n has a prime factor larger than L
          1 if all prime factors of n are < L
          0 if n could not be completely factored
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/

/* This is the same function as factor_leftover_norm() but it works
   with both norms! It is used when we want to factor these norms
   simultaneously and not one after the other.
   Return values:
   -1  one of the cofactors is not smooth
   0   unable to fully factor one of the cofactors
   1   both cofactors are smooth
*/

int factor_both_leftover_norms(
        std::array<cxx_mpz, 2> & n,
        std::array<std::vector<cxx_mpz>, 2> & factors,
        std::array<unsigned long, 2> const & Bs,
        facul_strategies_t const * strat)
{
    int is_smooth[2] = {FACUL_MAYBE, FACUL_MAYBE};
    /* To remember if a cofactor is already factored.*/

    for (int side = 0; side < 2; side++) {
        factors[side].clear();

        double B = (double) Bs[side];
        /* If n < B^2, then n is prime, since all primes < B have been removed */
        if (mpz_get_d (n[side]) < B * B)
            is_smooth[side] = FACUL_SMOOTH;
    }

    /* call the facul library */
    std::array<int, 2> facul_code = facul_both (factors, n, strat, is_smooth);

    if (is_smooth[0] != FACUL_SMOOTH || is_smooth[1] != FACUL_SMOOTH) {
        if (is_smooth[0] == FACUL_NOT_SMOOTH || is_smooth[1] == FACUL_NOT_SMOOTH)
            return -1;
        else
            return 0;
    }

    /* now we know both cofactors are smooth */
    for (int side = 0; side < 2; side++) {
        /* facul_code[side] is the number of found (smooth) factors */
        for (int i = 0; i < facul_code[side]; i++) {
            /* we know that factors found by facul_both() are primes < L */
            mpz_divexact (n[side], n[side], factors[side][i]);
            /* repeated factors should not be a problem, since they will
               be dealt correctly in the filtering */
        }

        /* since the cofactor is smooth, n[side] is a prime < L here */
        if (mpz_cmp_ui (n[side], 1) > 0) {
            /* 1 is special */
            factors[side].push_back(n[side]);
        }
    }
    return 1; /* both cofactors are smooth */
}


/*}}}*/
/*}}}*/


/*  sieve_info::{init,clear,share}_strategies */
void sieve_info::init_strategies(param_list_ptr pl)
{
    /* TODO: We may also build a strategy book, given that several
     * strategies will be similar. Presently we spend some time creating
     * each of them for the descent case.
     */
    const char *cofactfilename = param_list_lookup_string (pl, "file-cofact");
    /*
     * We create our strategy book from the file given by
     * 'file-cofact'. Otherwise, we use a default strategy given by
     * the function facul_make_default_strategy ().
     */

    FILE* file = NULL;
    if (cofactfilename != NULL) /* a file was given */
        file = fopen (cofactfilename, "r");
    double time_strat = seconds();

    strategies = std::shared_ptr<facul_strategies_t>(
        facul_make_strategies (conf.sides[0].lim,
            conf.sides[0].lpb,
            conf.sides[0].mfb,
            conf.sides[1].lim,
            conf.sides[1].lpb,
            conf.sides[1].mfb,
            (conf.sublat.m == 0), // with sublat, some primes are skipped.
            conf.sides[0].ncurves,
            conf.sides[1].ncurves,
            file, 0), facul_clear_strategies);

    verbose_output_print(0, 1, "# Building/reading strategies took %1.1fs\n",
            seconds() - time_strat);

    if (!strategies) {
        fprintf (stderr, "impossible to read %s\n", cofactfilename);
        abort ();
    }
    if (file)
        fclose (file);
}
