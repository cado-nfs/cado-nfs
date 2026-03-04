/*

We want a function that checks whether a relation is, with high probability, 
a duplicate of an earlier relation. 

"Earlier" requires that we define some ordering on relations. Note that we may
sieve special-q on both sides for some factorisations, in particular SNFS.
It is possible that the same prime occurs on both sides if it divides the 
two polynomials' resultant, which likewise may occur for SNFS.

We define the sieving order as sieving special-q on the algebraic side of 
norm n if n is in the special-q range for the algebraic side, 
then on the rational side of norm n if n is in the special-q range for the 
rational side, then increasing n by 1.

Let Q_0 be the primes in the special-q range on the algebraic side, 
and Q_1 on the rational side. If we sieve only one side i \in {0,1}, 
then the set Q_{1-i} is empty.

When we sieve special-q on only one side i, a relation is "earlier" 
simply if it contains a prime on side i that is in Q_i and smaller than 
the current q. (1)

If we sieve both sides, a relation is earlier if
- (1) holds, or
- we sieve q on the algebraic side and it contains a prime q' on the rational 
  side in Q_1 less than q
- we sieve q on the rational side and it contains a prime q' on the algebraic 
  side in Q_0 less than or equal to q

NOTE: this strategy for the case where we sieve two sides is not implemented.
Currently, las has no way to know about it.

Then we must go through the various steps the siever does to identify 
relations and check the conditions that cause the relation to be found or 
missed.

The relation is a dupe if:
- It may have been found by an "earlier" special-q $q'$
- The relation must be in the I,J sieving area that was used when sieving in 
  the $q'$-lattice
- The estimated norm, minus the log of the FB primes, is within the report 
  threshold. This is tricky, as there are many shortcuts to norm estimation,
  and the contribution of FB primes depends on the choice of log scale and 
  on limits on prime powers that are sieved
- The cofactors of the two norms, without q' but with q on the special-q side, 
  are within the cofactor bounds mfb[ar]
- The cofactors on both sides can be factored by the cofactorization routines 
  that were used when sieving q'.

Thus the function to check for duplicates needs the following information:

- The relation with fully factored norms
- The special-q sieve range on the two sides
- For each q' that was sieved "earlier"
  - The parameters used for generating the special-q' lattice (skewness)
  - The I,J values
  - The log scale
  - The cofactor bounds mfb[ar]
  - The cofactorization parameters

*/

#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <climits>
#include <cstdio>
#include <array>
#include <cstdint>
#include <vector>
#include <gmp.h>

#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "las-duplicate.hpp"
#include "ecm/facul.hpp"
#include "fb.hpp"
#include "gmp_aux.h"
#include "las-choose-sieve-area.hpp"
#include "las-cofactor.hpp"
#include "las-coordinates.hpp"
#include "las-norms.hpp"
#include "las-qlattice.hpp"
#include "las-side-config.hpp"
#include "las-siever-config.hpp"
#include "special-q.hpp"
#include "las-special-q-task-simple.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "relation.hpp"
#include "verbose.h"
#include "utils_cxx.hpp"


/* default verbose level of # DUPECHECK lines */
#define VERBOSE_LEVEL 2

struct sq_with_fac {
  uint64_t q;
  std::vector<uint64_t> facq;
};

// Warning: the entry we create here does not know whether sq is prime
// or composite (it assumes prime). This is fragile!!!
//
// XXX errr. I don't get it. If r=a/b mod all divisors of q, then this
// ought to hold mod q as well, right ?
//
// if we want to deal with the case where b is not coprime to q, the
// problem we have is that of doing a CRT on the projective line. Which
// is doable (well, to start with, the projective line element we're
// looking for is (a:b) anyway). But
// down the line, we stumble on the fact that just a plain integer isn't
// appropriate to represent elements of P1(Z/qZ) when q is composite.
//
static special_q special_q_from_ab(const int64_t a, const uint64_t b, sq_with_fac const & sq, int side)
{
    cxx_mpz p, r;

    mpz_set_ui(p, sq.q);
    mpz_set_uint64(r, b);
    int const ret = mpz_invert(r, r, p);
    ASSERT_ALWAYS(ret);
    mpz_mul_int64(r, r, a);
    mpz_mod(r, r, p);

    return { p, r, side };
}

/* If e == 0, returns 1. Otherwise, if b > lim, returns b.
   Otherwise returns the largest b^k with k <= e and b^k <= lim. */
static unsigned long
bounded_pow(const unsigned long b, const unsigned long e, const unsigned long lim)
{
  if (e == 0) {
    return 1;
  }
  if (b == 0 || b > lim) {
    return b;
  }
  unsigned long r = b;
  /* overflow is the largest ulong that can be multiplied by b without
     overflowing */
  unsigned long const overflow = ULONG_MAX / b;
  for (unsigned long i = 1; i < e && r <= overflow; i++) {
    unsigned long const t = r * b;
    if (t > lim) {
      break;
    }
    r = t;
  }
  return r;
}

/*
 * We subtract the sieve contribution of the factor base primes.
 * I.e., for each p^k || F(a,b) with p <= fbb:
 *   if k > 1 then we reduce k if necessary s.t. p^k <= powlim, but never to k < 1
 *   set l := l - log_b(p^k)
*/
static unsigned char
subtract_fb_log(const unsigned char lognorm,
        double scale,
        unsigned long lim, unsigned long powlim,
        relation const& rel,
        int side,
        special_q const & doing)
{
  unsigned char new_lognorm = lognorm;

  for (unsigned int i = 0; i < rel.sides[side].size(); i++) {
    const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
    int e = rel.sides[side][i].e;
    ASSERT_ALWAYS(e > 0);
    if (p >= lim)
      continue;

    // log(sq) is already subtracted from the lognorm we've received on
    // input.
    // For prime factors of sq, we need to count a smaller valuation
    // because powlim doesn't have to be that large.
    // (note that the list is empty if sq is prime).
    if (side == doing.side)
        for (auto const & f : doing.prime_factors)
            e -= (p == f);

    if (e == 0)
      continue;

    const unsigned long p_pow = bounded_pow(p, e, powlim);
    const unsigned char p_pow_log = fb_log(p_pow, scale, 0);

    if (p_pow_log > new_lognorm) {
      new_lognorm = 0;
    } else {
      new_lognorm -= p_pow_log;
    }
  }
  return new_lognorm;
}

/* Return true if the relation is probably a duplicate of a relation found
 * when sieving [doing].
 * Return false if it is probably not a duplicate */
static bool
sq_finds_relation(las_info const & las,
        special_q const & doing,
        siever_config const & conf,
        special_q_data_class auto const & Q,
        uint32_t J,
        relation const& rel,
        bool must,
        bool talk)
{
  int const logI = conf.logI;
  int const nsides = las.cpoly->nb_polys;

  if (talk) {   // Print some info
      verbose_fmt_print(0, VERBOSE_LEVEL,
              "# DUPECHECK Checking if relation"
              " (a,b) = ({})"
              " is a dupe of sieving special-q {}\n",
              rel.ab(), doing);
      // should stay consistent with "Sieving side" line printed in
      // do_one_special_q()
      verbose_fmt_print(0, VERBOSE_LEVEL,
              "# DUPECHECK Trying {}; I={}; J={};\n",
              Q, 1U << conf.logI, J);
  }

  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */
  int i;
  unsigned int j;
  {
    int ok;
    ok = Q.convert_ab_to_ij(i, j, rel.a, rel.b);
    if (!must && !ok) return false;
    ASSERT_ALWAYS(ok);
  }

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in [doing], then it's not a duplicate */
  if (j >= J || (i < -(1L << (logI-1))) || (i >= (1L << (logI-1))))
  {
    if (talk) verbose_fmt_print(0, VERBOSE_LEVEL,
        "# DUPECHECK (i,j) = ({}, {}) is outside sieve region\n", i, j);
    return false;
  }

  std::vector<uint8_t> remaining_lognorm(nsides);
  bool is_dupe = true;
  for (int side = 0; side < nsides; side++) {
    lognorm_smart const L(conf, las.cpoly, side, Q, conf.logI, J);

    /* This lognorm is the evaluation of fij(i,j)/q on the side with the
     * special-q, so we don't have to subtract it.
     */
    const uint8_t lognorm = L.lognorm(i,j, Q);

    remaining_lognorm[side] = subtract_fb_log(lognorm, L.scale,
            conf.sides[side].lim,
            conf.sides[side].powlim,
            rel,
            side,
            doing);

    if (remaining_lognorm[side] > L.bound) {
      if (talk)
          verbose_fmt_print(0, VERBOSE_LEVEL,
                  "# DUPECHECK On side {},"
                  " remaining lognorm = {} > bound = {}\n",
                  side, remaining_lognorm[side], L.bound);
      is_dupe = false;
    }
  }
  if (talk) {
      verbose_fmt_print(0, VERBOSE_LEVEL,
              "# DUPECHECK relation had i={}, j={},"
              " remaining lognorms {}\n",
              i, j, join(remaining_lognorm, " "));
  }

  if (!is_dupe) {
    return false;
  }

  /* Compute the exact cofactor on each side. This is similar to what we
   * do in subtract_fb_log -- except that now we do the part *above* the
   * factor base. */
  std::vector<cxx_mpz> cof(nsides);
  for (int side = 0; side < nsides; side++) {
    mpz_set_ui(cof[side], 1);
    unsigned long const lim = conf.sides[side].lim;

    for (unsigned int i = 0; i < rel.sides[side].size(); i++) {
      const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
      int e = rel.sides[side][i].e;
      ASSERT_ALWAYS(e > 0);
      if (p <= lim) {
        continue;
      }

      // Remove one occurrence of each factor of sq on side sq_side.
      if (side == doing.side)
        for (auto const & f : doing.prime_factors)
            e -= (p == f);

      for (int i = 0; i < e; ++i) {
        mpz_mul_ui(cof[side], cof[side], p);
      }
    }
  }

  /* Check that the cofactors are within the mfb bound */
  for (int side = 0; side < nsides; ++side) {
    if (!check_leftover_norm (cof[side], conf.sides[side])) {
      verbose_fmt_print(0, VERBOSE_LEVEL,
          "# DUPECHECK cofactor {} is outside bounds\n",
          cof[side]);
      return false;
    }
  }

  std::vector<std::vector<cxx_mpz>> f(nsides);
  auto Bs = siever_side_config::collect_lim(conf.sides);

  facul_status const pass = factor_leftover_norms(cof, f, Bs,
          *las.get_strategies(conf));

  if (pass != FACUL_SMOOTH) {
      if (talk) {
          verbose_fmt_print(0, VERBOSE_LEVEL,
                        "# DUPECHECK norms not both smooth,"
                        " leftover factors: {}\n",
                        join(cof, ", "));
      }
      return false;
  }

  return true;
}

bool
sq_finds_relation(las_info const & las,
        special_q const & doing,
        siever_config const & conf,
        special_q_data_class auto const & Q,
        uint32_t J,
        relation const& rel)
{
    return sq_finds_relation(las, doing, conf, Q, J, rel, false, false);
}

template bool sq_finds_relation(las_info const &,
        special_q const &,
        siever_config const &,
        qlattice_basis const &,
        uint32_t,
        relation const &);
template bool sq_finds_relation(las_info const &,
        special_q const &,
        siever_config const &,
        siqs_special_q_data const &,
        uint32_t,
        relation const &);

static bool
sq_finds_relation(las_info const & las,
        special_q const & doing,
        relation const& rel,
        bool must,
        bool talk)
{
    siever_config conf;
    qlattice_basis Q;
    uint32_t J;
    if (!choose_sieve_area(las, special_q_task_simple(doing), conf, Q, J)) {
        if (talk) verbose_fmt_print(0, VERBOSE_LEVEL, "# DUPECHECK q-lattice discarded\n");
        return false;
    }
    return sq_finds_relation(las, doing, conf, Q, J, rel, must, talk);
}


/* This function decides whether the given (sq,side) was previously
 * sieved (compared to the current special-q stored in [doing]).
 * This takes qmin and qmax into account, on the side of the sq. The side
 * doing.side is irrelevant here.
 */
static int
sq_was_previously_sieved (las_info const & las, const uint64_t sq, int side, special_q const & doing)
{
    /* we use <= and not < since this function is also called with the
     * current special-q */
    if (mpz_cmp_uint64(doing.p, sq) <= 0)
        return 0;

    return sq >= las.dupqmin[side] && sq < las.dupqmax[side];
}

// Warning: this function works with side effects:
//   - it is recursive
//   - it destroys its argument along the way
// Keep only products that fit in 64 bits.
static std::vector<sq_with_fac>
all_multiples(std::vector<uint64_t> & prime_list) {
    if (prime_list.empty())
        return { { .q=1, .facq={} } };

    uint64_t const p = prime_list.back();
    prime_list.pop_back();

    /* this if block is useless, I think */
    if (prime_list.empty())
        return {
            { .q=1, .facq={} },
            { .q=p, .facq={p} },
        };

    std::vector<sq_with_fac> res = all_multiples(prime_list);

    size_t const L = res.size();
    const cxx_mpz pp = p;
    for (size_t i = 0; i < L; ++i) {
        /* create a new multiple, and check if it fits */
        std::vector<uint64_t> facq = res[i].facq;
        facq.push_back(p);
        const cxx_mpz prod = pp * res[i].q;

        if (mpz_fits_uint64_p(prod))
            res.emplace_back(mpz_get_uint64(prod), facq);
    }
    return res;
}

/* Return 1 if the relation, obtained from sieving [doing], is probably a
 * duplicate of a relation found "earlier", and 0 if it is probably not a
 * duplicate
 *
 * We may imagine having different cofactoring strategies depending on q,
 * in which case the code below is incorrect. However, for normal
 * (non-dlp-descent) sieving, there does not seem to be a compelling
 * argument for that.
 */
int
relation_is_duplicate(relation const& rel,
        special_q const & doing,
        las_info const& las)
{
    int const nsides = las.cpoly->nb_polys;

    /* If the special-q does not fit in an unsigned long, we assume it's not a
     * duplicate and just move on.
     * Note: these conditions mean that this function will always return false
     * for SIQS. FIXME
     */
    if (doing.is_prime() && !mpz_fits_uint64_p(doing.p)) {
        return false;
    }

    /* If any large prime doesn't fit in a 64-bit integer, then we assume
     * that we're doing the dlp desecent, in which case we couldn't care
     * less about duplicates check anyway.
     */
    for(auto const & s : rel.sides) {
        for(auto const & pr : s) {
            if (!mpz_fits_uint64_p(pr.p))
                return false;
        }
    }

    for(int side = 0 ; side < nsides ; side++) {
        /* It is allowed to have special-q on both sides, and we wish to
         * check "cross" special-q combinations.
         */

        // Step 1: prepare a list of potential special-q on this side.
        // They involve only prime factors within the allowed bounds
        std::vector<uint64_t> prime_list;

        for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
            uint64_t const p = mpz_get_uint64(rel.sides[side][i].p);

            // can this p be part of valid sq ?
            if (!las.tree->todo->allow_composite_special_q()) {
                /* this check is also done in sq_was_previously_sieved,
                 * but it's easy enough to skip the divisibility test
                 * when we know there's no point.
                 */
                if ((p < las.dupqmin[side]) || (p >= las.dupqmax[side])) {
                    continue;
                }
            } else {
                if (!las.tree->todo->is_in_qfac_range(p))
                    continue;
            }

            /* projective primes are currently not allowed for composite
             * special-q */
            if (mpz_divisible_ui_p(mpz_poly_lc(las.cpoly->pols[side]), p))
                continue;

            // push it in the list of potential factors of sq
            prime_list.push_back(p);
        }

        /* FIXME: For hybrid special-q's, this must be changed! */

        for (auto const & sq : all_multiples(prime_list)) {
            // keep sq only if it was sieved before [doing]
            if (!sq_was_previously_sieved(las, sq.q, side, doing))
                continue;

            // emulate sieving for the valid sq, and check if it finds our
            // relation.
            special_q const other = special_q_from_ab(rel.a, rel.b, sq, side);

            bool const is_dupe = sq_finds_relation(las, other, rel, true, true);
            verbose_fmt_print(0, VERBOSE_LEVEL,
                    "# DUPECHECK relation is probably{} a dupe\n",
                    is_dupe ? "" : " not");
            if (is_dupe) return true;
        }
    }

    return false;
}
