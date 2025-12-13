#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <vector>

#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/modredc_ul.h"
#include "ecm/facul.hpp"
#include "ecm/facul_strategies.hpp"
#include "las-cofactor.hpp"
#include "macros.h"
#include "params.h"

void cofactorization_statistics::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "stats-cofact", "write statistics about the cofactorization step in file xxx");
}

//  las_info::{init,clear,print}_cof_stats
cofactorization_statistics::cofactorization_statistics(param_list_ptr pl)
{
    const char * statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (!statsfilename) {
        file = nullptr;
        return;
    }
    file = fopen (statsfilename, "w");
    if (file == nullptr) {
        fprintf (stderr, "Error, cannot create file %s\n", statsfilename);
        exit (EXIT_FAILURE);
    }
}

void cofactorization_statistics::call(int bits0, int bits1)
{
    if (!file) return;
    std::lock_guard<std::mutex> const dummy(lock);
    size_t s0 = cof_call.size();
    if ((size_t) bits0 >= s0) {
        size_t const news0 = std::max((size_t) bits0+1, s0 + s0/2);
        cof_call.insert(cof_call.end(), news0-s0, std::vector<uint32_t>());
        cof_success.insert(cof_success.end(), news0-s0, std::vector<uint32_t>());
        s0 = news0;
    }
    size_t s1 = cof_call[bits0].size();
    if ((size_t) bits1 >= s1) {
        size_t const news1 = std::max((size_t) bits1+1, s1 + s1/2);
        cof_call[bits0].insert(cof_call[bits0].end(), news1-s1, 0);
        cof_success[bits0].insert(cof_success[bits0].end(), news1-s1, 0);
        s1 = news1;
    }
    ASSERT_ALWAYS((size_t) bits0 < s0);
    ASSERT_ALWAYS((size_t) bits1 < s1);
    /* no need to use a mutex here: either we use one thread only
       to compute the cofactorization data and if several threads
       the order is irrelevant. The only problem that can happen
       is when two threads increase the value at the same time,
       and it is increased by 1 instead of 2, but this should
       happen rarely. */
    cof_call[bits0][bits1]++;
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
//

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
check_leftover_norm (cxx_mpz const & n, siever_side_config const & scs)
{
  size_t const s = mpz_sizeinbase (n, 2);
  unsigned int const lpb = scs.lpb;
  unsigned int const mfb = scs.mfb;
  unsigned int klpb;
  double nd, kB, B;

  ASSERT_ALWAYS(mpz_cmp_ui(n, 0) != 0);

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
   with all norms! It is used when we want to factor these norms
   simultaneously and not one after the other.
   Return values:
   -1  one of the cofactors is not smooth
   0   unable to fully factor one of the cofactors
   1   all cofactors are smooth

  Note: for more than two sides, we may still get a relation even if not all
  cofactors were smooth. Currently, it is not taken into account by this method.
*/

facul_status factor_leftover_norms(
        std::vector<cxx_mpz> const & n,
        std::vector<std::vector<cxx_mpz>> & factors,
        std::vector<unsigned long> const & Bs,
        facul_strategies const & strat)
{
    ASSERT_ALWAYS(Bs.size() == strat.B.size());
    ASSERT_ALWAYS(std::ranges::equal(Bs, strat.B));

    /* call the facul library */
    auto fac = facul_all(n, strat);
    for(auto const & f : fac) {
        if (f.status == FACUL_NOT_SMOOTH)
            return FACUL_NOT_SMOOTH;
    }
    for(size_t side = 0 ; side < fac.size() ; side++) {
        auto & f = fac[side];
        if (f.status == FACUL_MAYBE) {
            /* We couldn't factor this number. So we don't know. It
             * happens also for tiny examples, which is a bit of a pity
             * (see full_p30_JL test for example).
             * Maybe it's due to our incomplete backtracking.
             * I don't have a firm opinion as to whether this needs
             * further investigation or not. At any rate, a "MAYBE" on
             * one side means a global "MAYBE", and for consistency with
             * what the code has been doing for quite some time, let's
             * return that.
             */
            return FACUL_MAYBE;
        }
#ifndef NDEBUG
        cxx_mpz z = 1;
        for(auto const & p : f.primes)
            z *= p;
        ASSERT_ALWAYS(z == n[side]);
#endif
        std::swap(factors[side], f.primes);
    }
    return FACUL_SMOOTH;
}


/*}}}*/
/*}}}*/

