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


#include "cado.h"
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include <string.h>
#include "gmp_aux.h"
#include "las-duplicate.hpp"
#include "las-qlattice.hpp"
#include "las-coordinates.hpp"
#include "las-norms.hpp"
#include "las-cofactor.hpp"

/* default verbose level of # DUPECHECK lines */
#define VERBOSE_LEVEL 2

static void
compute_a_over_b_mod_p(mpz_t r, const int64_t a, const uint64_t b, const mpz_t p)
{
  mpz_set_uint64(r, b);
  int ret = mpz_invert(r, r, p);
  ASSERT_ALWAYS(ret);
  mpz_mul_int64(r, r, a);
  mpz_mod(r, r, p);
}

sieve_info *
fill_in_sieve_info(las_todo_entry const & doing,
                   uint32_t I, uint32_t J,
                   cado_poly_ptr cpoly, siever_config const & conf, int nb_threads)
{
  sieve_info * x = new sieve_info;

  sieve_info & new_si(*x);

  new_si.I = I;
  new_si.J = J;
  new_si.cpoly = cpoly;
  new_si.conf = conf;
  new_si.conf.side = doing.side;
  new_si.doing = doing;

  sieve_range_adjust Adj(doing, cpoly, conf, nb_threads);
  Adj.SkewGauss();
  if (!Adj.sieve_info_adjust_IJ()) {
      delete x;
      return NULL;
  }
  Adj.sieve_info_update_norm_data_Jmax();
  new_si.I = 1UL << Adj.logI;
  new_si.J = Adj.J;

  return x;
}

las_todo_entry special_q_from_ab(const int64_t a, const uint64_t b, const unsigned long p, int side)
{
  mpz_t sq, rho;
  mpz_init_set_ui(sq, p);
  mpz_init(rho);
  compute_a_over_b_mod_p(rho, a, b, sq);
  las_todo_entry doing(sq, rho, side);
  mpz_clear(sq);
  mpz_clear(rho);
  return doing;
}


/* If e == 0, returns 1. Otherwise, if b > lim, returns b.
   Otherwise returns the largest b^k with k <= e and b^k <= lim. */
static unsigned long
bounded_pow(const unsigned long b, const unsigned long e, const unsigned long lim)
{
  if (e == 0) {
    return 1;
  }
  if (b > lim) {
    return b;
  }
  unsigned long r = b;
  /* overflow is the largest ulong that can be multiplied by b without
     overflowing */
  unsigned long overflow = ULONG_MAX / b;
  for (unsigned long i = 1; i < e && r <= overflow; i++) {
    unsigned long t = r * b;
    if (t > lim) {
      break;
    }
    r = t;
  }
  return r;
}

/*
   We subtract the sieve contribution of the factor base primes.
   I.e., for each p^k || F(a,b) with p <= fbb:
     if k > 1 then we reduce k if necessary s.t. p^k <= powlim, but never to k < 1
     set l := l - log_b(p^k)
*/
static unsigned char
subtract_fb_log(const unsigned char lognorm, relation const& rel,
                sieve_info const & si, const int side,
                std::vector<uint64_t> const * facq)
{
  const unsigned long fbb = si.conf.sides[side].lim;
  const unsigned int nb_p = rel.sides[side].size();
  unsigned char new_lognorm = lognorm;

  for (unsigned int i = 0; i < nb_p; i++) {
    lognorm_base & L(*si.sides[side].lognorms);
    const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
    int e = rel.sides[side][i].e;
    ASSERT_ALWAYS(e > 0);
    if (p >= fbb)
      continue;

    if (facq != NULL) {
      // Remove one occurrence of each factor of sq.
      for (auto const & f : *facq) {
        if (p == f) {
          e = e - 1;
        }
      }
    }
    if (e == 0)
      continue;

    const unsigned long p_pow = bounded_pow(p, e, si.conf.sides[side].powlim);
    const unsigned char p_pow_log = fb_log(p_pow, L.scale, 0);
    if (p_pow_log > new_lognorm) {
      new_lognorm = 0;
    } else {
      new_lognorm -= p_pow_log;
    }
  }
  return new_lognorm;
}

struct sq_with_fac {
  uint64_t q;
  std::vector<uint64_t> facq;
};

/* Return true if the relation is probably a duplicate of a relation found
 * when sieving the sq described by si. Return false if it is probably not a
 * duplicate */
bool
sq_finds_relation(sq_with_fac const& sq_fac, const int sq_side,
    relation const& rel, const int nb_threads, sieve_info & old_si,
    int adjust_strategy)
{
  uint64_t sq = sq_fac.q;
  
  // Warning: the entry we create here does not know whether sq is
  // prime or composite (it assumes prime). This is fragile!!!
  las_todo_entry doing = special_q_from_ab(rel.a, rel.b, sq, sq_side);

  sieve_range_adjust Adj(doing, old_si.cpoly, old_si.conf, nb_threads);

  if (!Adj.SkewGauss() || !Adj.Q.fits_31bits() || !Adj.sieve_info_adjust_IJ()) {
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK q-lattice discarded\n");
    return false;
  }

  if (adjust_strategy != 1)
    Adj.sieve_info_update_norm_data_Jmax();

  if (adjust_strategy >= 2)
    Adj.adjust_with_estimated_yield();

  if (adjust_strategy >= 3)
    Adj.sieve_info_update_norm_data_Jmax(true);

  siever_config conf = Adj.config();
  conf.logI_adjusted = Adj.logI;
  conf.side = sq_side;

  /* We don't have a constructor which is well adapted to our needs here.
   * We're going to play dirty tricks, and fill out the stuff by
   * ourselves.
   */
  sieve_info si;
  si.cpoly = old_si.cpoly;
  si.conf = conf;
  si.I = 1UL << conf.logI_adjusted;
  si.recover_per_sq_values(Adj);
  si.strategies = old_si.strategies;


  const uint32_t oldI = si.I, oldJ = si.J;
  si.update_norm_data();
  uint32_t I = si.I, J = si.J;

  
  { // Print some info
    const unsigned long r = mpz_get_ui(doing.r);
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK Checking if relation (a,b) = (%" PRId64 ",%" PRIu64 ") is a dupe of sieving special-q -q0 %" PRIu64 " -rho %lu\n", rel.a, rel.b, sq, r);
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK Using special-q basis a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "\n", si.qbasis.a0, si.qbasis.b0, si.qbasis.a1, si.qbasis.b1);
    if (oldI != I || oldJ != J)
      verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK oldI = %u, I = %u, oldJ = %u, J = %u\n", oldI, I, oldJ, J);
  }
  
  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */
  int i;
  unsigned int j;
  {
    int ok;
    ok = ABToIJ(&i, &j, rel.a, rel.b, si);
    ASSERT_ALWAYS(ok);
  }

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in si, then it's not a duplicate */
  if ((i < 0 && (uint32_t)(-i) > I/2) || (i > 0 && (uint32_t)i > I/2-1) || (j >= J)) {
    verbose_output_print(0, VERBOSE_LEVEL,
        "# DUPECHECK (i,j) = (%d, %u) is outside sieve region\n", i, j);
    return false;
  }

  uint8_t remaining_lognorm[2];
  bool is_dupe = true;
  for (int side = 0; side < 2; side++) {
    lognorm_base & L(*si.sides[side].lognorms);
    /* Here, we assume for now that the log norm estimate is exact,
     * i.e., that it produces the correctly rounded l = log_b(F(a,b)).
     * So we're using a straightforward function designed precisely for
     * that in las-norms. Now, in fact, it would be better to call the
     * _real_ function and pick the entry corresponding to this
     * i,j-coordinate. */
    const uint8_t lognorm = L.lognorm(i,j);
    remaining_lognorm[side] = subtract_fb_log(lognorm, rel, si, side,
        (side == sq_side) ? &(sq_fac.facq) : NULL);
    if (remaining_lognorm[side] > L.bound) {
      verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK On side %d, remaining lognorm = %" PRId8 " > bound = %" PRId8 "\n",
          side, remaining_lognorm[side], L.bound);
      is_dupe = false;
    }
  }
  verbose_output_print(0, VERBOSE_LEVEL,
      "# DUPECHECK relation had i=%d, j=%u, remaining lognorms %" PRId8 ", %" PRId8 "\n",
      i, j, remaining_lognorm[0], remaining_lognorm[1]);
  if (!is_dupe) {
    return false;
  }

  /* Compute the exact cofactor on each side */
  std::array<cxx_mpz, 2> cof;
  for (int side = 0; side < 2; side++) {
    mpz_set_ui(cof[side], 1);
    const unsigned long fbb = old_si.conf.sides[side].lim;
    unsigned int nb_p = rel.sides[side].size();
    for (unsigned int i = 0; i < nb_p; i++) {
      const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
      int e = rel.sides[side][i].e;
      ASSERT_ALWAYS(e > 0);
      if (p <= fbb) {
        continue;
      }
      // Remove one occurrence of each factor of sq on side sq_side.
      if (side == sq_side) {
        for (auto const & facq : sq_fac.facq) {
          if (p == facq) {
            e = e - 1;
          }
        }
      }
      for (int i = 0; i < e; ++i) {
        mpz_mul_ui(cof[side], cof[side], p);
      }
    }
  }

  /* Check that the cofactors are within the mfb bound */
  for (int side = 0; side < 2; ++side) {
    if (!check_leftover_norm (cof[side], si, side)) {
      verbose_output_vfprint(0, VERBOSE_LEVEL, gmp_vfprintf,
          "# DUPECHECK cofactor %Zd is outside bounds\n",
          (__mpz_struct*) cof[side]);
      return false;
    }
  }

  std::array<std::vector<cxx_mpz>, 2> f;
  int pass = factor_both_leftover_norms(cof, f, si);

  if (pass <= 0) {
    verbose_output_vfprint(0, VERBOSE_LEVEL, gmp_vfprintf,
        "# DUPECHECK norms not both smooth, left over factors: %Zd, %Zd\n",
        (mpz_srcptr) cof[0], (mpz_srcptr) cof[1]);
    return false;
  }

  return true;
}


/* This function decides whether the given (sq,side) was previously
 * sieved (compared to the current special-q stored in si.doing).
 * This takes qmin and qmax into account.
 */
static int
sq_was_previously_sieved (const uint64_t sq, int side, sieve_info const & si){
  cxx_mpz Sq;
  mpz_set_uint64(Sq, sq);
  if (mpz_cmp(si.doing.p, Sq) <= 0) /* we use <= and not < since this
					   function is also called with the
					   current special-q */
    return 0;

  return (sq >= si.conf.sides[side].qmin)
      && (sq <  si.conf.sides[side].qmax);
}

// Warning: this function works with side effects:
//   - it is recursive
//   - it destroys its argument along the way
//   - it allocates the result it returns; the caller must free it.
// Keep only products that fit in 64 bits.
std::vector<sq_with_fac>
all_multiples(std::vector<uint64_t> & prime_list) {
  if (prime_list.empty()) {
    std::vector<sq_with_fac> res;
    return res;
  }

  uint64_t p = prime_list.back();
  prime_list.pop_back();

  if (prime_list.empty()) {
    std::vector<sq_with_fac> res;
    res.push_back(sq_with_fac { 1, std::vector<uint64_t> () });
    res.push_back(sq_with_fac { p, std::vector<uint64_t> (1, p) });
    return res;
  }

  std::vector<sq_with_fac> res = all_multiples(prime_list);

  size_t L = res.size();
  cxx_mpz pp;
  mpz_set_uint64(pp, p);
  for (size_t i = 0; i < L; ++i) {
    std::vector<uint64_t> facq;
    facq = res[i].facq;
    facq.push_back(p);

    cxx_mpz prod;
    mpz_mul_uint64(prod, pp, res[i].q);
    if (mpz_fits_uint64_p(prod)) {
      uint64_t q = mpz_get_uint64(prod);
      struct sq_with_fac entry = {q, facq};
      res.push_back(entry);
    }
  }
  return res;
}

/* Return 1 if the relation is probably a duplicate of a relation found
   "earlier", and 0 if it is probably not a duplicate */
int
relation_is_duplicate(relation const& rel, las_info const& las,
                      sieve_info & si, int adjust_strategy)
{
  /* If the special-q does not fit in an unsigned long, we assume it's not a
     duplicate and just move on */
  if (!mpz_fits_uint64_p(si.doing.p)) {
    return false;
  }

  /* If any large prime doesn't fit in an unsigned long, then we assume
   * we don't have a duplicate */
  for(int side = 0 ; side < 2 ; side++) {
    for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
      if (!mpz_fits_uint64_p(rel.sides[side][i].p))
        return false;
    }
  }

  for(int side = 0 ; side < 2 ; side++) {
    // Step 1: prepare a list of potential special-q on this side.
    // They involve only prime factors within the allowed bounds
    std::vector<uint64_t> prime_list;
    for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
      uint64_t p = mpz_get_uint64(rel.sides[side][i].p);

      // can this p be part of valid sq ?
      if (! las.allow_composite_q) {
        if ((p < si.conf.sides[side].qmin) || (p >= si.conf.sides[side].qmax)) {
          continue;
        }
      } else {
        if ((p < las.qfac_min) || (p >= las.qfac_max)) {
          continue;
        }
      }

      cxx_mpz aux;
      mpz_poly_getcoeff(aux, si.cpoly->pols[side]->deg, si.cpoly->pols[side]);
      if (mpz_divisible_ui_p(aux, p))
        continue;

      // push it in the list of potential factors of sq
      prime_list.push_back(p);
    }
    std::vector<sq_with_fac> sq_list = all_multiples(prime_list);

    // Step 2: keep only those that have been sieved before current sq.
    std::vector<sq_with_fac> valid_sq;
    for (auto const & sq : sq_list) {
      if (sq_was_previously_sieved(sq.q, side, si)) {
        valid_sq.push_back(sq);
      }
    }

    // Step 3: emulate sieving for the valid sq, and check if they find
    // our relation.
    for (auto const & sq : valid_sq) {
      bool is_dupe = sq_finds_relation(sq, side, rel, las.nb_threads, si, adjust_strategy);
      verbose_output_print(0, VERBOSE_LEVEL,
          "# DUPECHECK relation is probably%s a dupe\n",
          is_dupe ? "" : " not");
      if (is_dupe) {
        return true;
      }
    }
  }

  return false;
}
