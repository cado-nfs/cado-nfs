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
#include "mpz_array.h"
#include "las-duplicate.hpp"
#include "las-qlattice.hpp"
#include "las-coordinates.hpp"
#include "las-norms.hpp"
#include "las-cofactor.hpp"

/* default verbose level of # DUPECHECK lines */
#define VERBOSE_LEVEL 2

static void
compute_minus_a_over_b_mod_p(mpz_t r, cxx_mpz const & a, cxx_mpz const & b, const mpz_t p)
{
  mpz_neg(r, b);
  mpz_invert(r, r, p);
  mpz_mul(r, r, a);
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

las_todo_entry special_q_from_ab(cxx_mpz const& a, cxx_mpz const& b, const unsigned long p, int side)
{
  cxx_mpz sq, rho;
  mpz_init_set_ui(sq, p);
  compute_minus_a_over_b_mod_p(rho, a, b, sq);
  las_todo_entry doing(sq, rho, side);
  return doing;
}


/* Compute the product of the nr_lp large primes in large_primes,
   skipping sq. If sq occurs multiple times in large_primes, it is
   skipped only once. */
static void
compute_cofactor(mpz_t cof, const unsigned long sq, 
                 const unsigned long *large_primes, const int nr_lp)
{
  mpz_set_ui (cof, 1UL);
  int saw_sq = 0;
  for (int i = 0; i < nr_lp; i++) {
    const unsigned long p = large_primes[i];
    if (sq != 0 && !saw_sq && p == sq) {
      saw_sq = 1;
      continue;
    }
    mpz_mul_ui(cof, cof, p);
  }
  ASSERT_ALWAYS(sq == 0 || saw_sq);
}

/* Here, we assume for now that the log norm estimate is exact, i.e., that
  it produces the correctly rounded l = log_b(F(a,b)).
  To make the estimate more precise, we should call the norm calculation
  function for the correct bucket region and pick the entry corresponding to
  this i,j-coordinate. */
static unsigned char
estimate_lognorm(sieve_info const & si, const int i, const unsigned int j,
                  const int side)
{
  mpz_t norm;
  mpz_init (norm);
  mpz_poly_homogeneous_eval_sisi(norm, si.sides[side].fij, i, j);
  return fb_log(mpz_get_d(norm), si.sides[side].scale * LOG_SCALE, 0.) + GUARD;
  mpz_clear (norm);
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
                sieve_info const & si, const int side)
{
  const unsigned long fbb = si.conf.sides[side].lim;
  const unsigned int nb_p = rel.sides[side].size();
  unsigned char new_lognorm = lognorm;

  for (unsigned int i = 0; i < nb_p; i++) {
    const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
    const int e = rel.sides[side][i].e;
    ASSERT_ALWAYS(e > 0);
    if (p <= fbb) {
      const unsigned long p_pow = bounded_pow(p, e, si.conf.sides[side].powlim);
      const unsigned char p_pow_log = fb_log(p_pow, si.sides[side].scale * LOG_SCALE, 0.);
      if (p_pow_log > new_lognorm) {
        if (0)
          gmp_fprintf(stderr, "Warning: lognorm underflow for relation a,b = %Zx,%Zx\n", (mpz_srcptr) rel[0], (mpz_srcptr) rel[1]);
        new_lognorm = 0;
      } else {
        new_lognorm -= p_pow_log;
      }
    }
  }
  return new_lognorm;
}

/* Return 1 if the relation is probably a duplicate of a relation found
 * when sieving the sq described by si. Return 0 if it is probably not a
 * duplicate */
int
sq_finds_relation(const unsigned long sq, const int sq_side,
                  relation const& rel,
                  const int nb_threads, sieve_info & old_si)
{
  mpz_array_t *f[2] = { NULL, }; /* Factors of the relation's norms */
  uint32_array_t *m[2] = { NULL, }; /* corresponding multiplicities */
  cxx_mpz cof[2];
  int is_dupe = 1; /* Assumed dupe until proven innocent */
  const size_t max_large_primes = 10;
  unsigned long large_primes[2][max_large_primes];
  unsigned int nr_lp[2] = {0, 0};
  int i, pass, ok;
  unsigned int j;


  /* Extract the list of large primes for each side */
  for (int side = 0; side < 2; side++) {
    unsigned int nb_p = rel.sides[side].size();
    for (unsigned int i = 0; i < nb_p; i++) {
      const unsigned long fbb = old_si.conf.sides[side].lim;
      const unsigned long p = mpz_get_ui(rel.sides[side][i].p);
      const int e = rel.sides[side][i].e;
      ASSERT_ALWAYS(e > 0);
      if (p > fbb) {
        for (int i = 0; i < e; i++) {
          if (nr_lp[side] < max_large_primes)
            large_primes[side][nr_lp[side]++] = p;
        }
      }
    }
    /* Compute the cofactor, dividing by sq on the special-q side */
    compute_cofactor(cof[side], side == sq_side ? sq : 0, large_primes[side], nr_lp[side]);
  }

  /* TODO: we might consider doing this even for higher-dimensional
   * sieving. Headaches...
   */
  if (rel.ab().size() != 2)
      return 0;

  las_todo_entry doing = special_q_from_ab(rel[0], rel[1], sq, sq_side);

  sieve_range_adjust Adj(doing, old_si.cpoly, old_si.conf, nb_threads);

  if (!Adj.SkewGauss() || !Adj.Q.fits_31bits() || !Adj.sieve_info_adjust_IJ()) {
      verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK q-lattice discarded\n");
      return 0;
  }

  /* We have no info as to which adjustment option is being used for the
   * current project...
   */
  Adj.sieve_info_update_norm_data_Jmax();
  siever_config conf = Adj.config();
  conf.logI_adjusted = Adj.logI;

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

  const unsigned long r = mpz_get_ui(doing.r);

  const uint32_t oldI = si.I, oldJ = si.J;
  uint32_t I, J; /* Can't declare further down due to goto :( */

  si.update_norm_data();

  verbose_output_vfprint(0, VERBOSE_LEVEL, gmp_vfprintf, "# DUPECHECK Checking if relation (a,b) = (%Zx,%Zx) is a dupe of sieving special-q -q0 %lu -rho %lu\n", (mpz_srcptr) rel[0], (mpz_srcptr) rel[1], sq, r);
  verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK Using special-q basis a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "\n", si.qbasis.a0, si.qbasis.b0, si.qbasis.a1, si.qbasis.b1);

  I = si.I;
  J = si.J;

  if (oldI != I || oldJ != J)
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK oldI = %u, I = %u, oldJ = %u, J = %u\n", oldI, I, oldJ, J);
  
  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */
  ok = ABToIJ(&i, &j, mpz_get_int64(rel[0]), mpz_get_int64(rel[1]), si);

  if (!ok)
    abort();

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in si, then it's not a duplicate */
  if ((i < 0 && (uint32_t)(-i) > I/2) || (i > 0 && (uint32_t)i > I/2-1) || (j >= J)) {
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK (i,j) = (%d, %u) is outside sieve region\n", i, j);
    is_dupe = 0;
    goto clear_and_exit;
  }

  uint8_t remaining_lognorm[2];
  for (int side = 0; side < 2; side++) {
    const uint8_t lognorm = estimate_lognorm(si, i, j, side);
    remaining_lognorm[side] = subtract_fb_log(lognorm, rel, si, side);
    if (remaining_lognorm[side] > si.sides[side].bound) {
      verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK On side %d, remaining lognorm = %" PRId8 " > bound = %" PRId8 "\n",
              side, remaining_lognorm[side], si.sides[side].bound);
      is_dupe = 0;
    }
  }
  verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK relation had i=%d, j=%u, remaining lognorms %" PRId8 ", %" PRId8 "\n",
           i, j, remaining_lognorm[0], remaining_lognorm[1]);
  if (!is_dupe) {
    goto clear_and_exit;
  }

  /* Check that the cofactors are within the mfb bound */
  for (int side = 0; side < 2; ++side) {
    if (!check_leftover_norm (cof[side], si, side)) {
      verbose_output_vfprint(0, VERBOSE_LEVEL, gmp_vfprintf, "# DUPECHECK cofactor %Zd is outside bounds\n", cof);
      is_dupe = 0;
      goto clear_and_exit;
    }
  }

  for(int side = 0 ; side < 2 ; side++) {
      f[side] = alloc_mpz_array (1);
      m[side] = alloc_uint32_array (1);
  }

  // ok, the cast is despicable
  pass = factor_both_leftover_norms((mpz_t*)cof, f, m, si);

  if (pass <= 0)
    verbose_output_vfprint(0, VERBOSE_LEVEL, gmp_vfprintf, "# DUPECHECK norms not both smooth, left over factors: %Zd, %Zd\n", (mpz_srcptr) cof[0], (mpz_srcptr) cof[1]);

  for(int side = 0 ; side < 2 ; side++) {
      clear_uint32_array (m[side]);
      clear_mpz_array (f[side]);
  }

  if (pass <= 0) {
    is_dupe = 0;
    goto clear_and_exit;
  }

clear_and_exit:

  return is_dupe;
}


/* This function decides whether the given (sq,side) was previously
 * sieved (compared to the current sq stored in si.doing).
 * This takes qmax into account.
 */
static int
sq_was_previously_sieved(const unsigned long sq, int side, sieve_info const & si){
  if (side == si.doing.side) {
    int cmp = mpz_cmp_ui(si.doing.p, sq);
    if (cmp <= 0)
      return 0;
    return (sq > si.conf.sides[side].lim);
  } else {
    // Did we sieve other sides?
    if (si.conf.sides[side].qmax == 0)
      return 0;
    // Is q smaller than current, and within the bounds for this side?
    int cmp = mpz_cmp_ui(si.doing.p, sq);
    if (cmp <= 0)
      return 0;
    return (sq > si.conf.sides[side].lim && sq < si.conf.sides[side].qmax);
  }
}

/* For one special-q identified by (sq, side) (the root r is given
   implicitly by the a,b-coordinate in relation), check whether the
   relation is a duplicate of sieving that special-q: check whether 
   (sq, side) was sieved before the special-q specified in si, and
   whether sieving (sq, side) can find the relation with the sieving
   parameters specified in si.conf. */ 

int
check_one_prime(mpz_srcptr zsq, const int side,
                relation const & rel, const int nb_threads,
                sieve_info & si)
{
  if (!mpz_fits_ulong_p(zsq)) {
      /* move on. Probably not a duplicate */
      return 0;
  }
  unsigned long sq = mpz_get_ui(zsq);
  int is_dupe = 0;
  if (sq_was_previously_sieved(sq, side, si)) {
    is_dupe = sq_finds_relation(sq, side, rel, nb_threads, si);
    verbose_output_print(0, VERBOSE_LEVEL, "# DUPECHECK relation is probably%s a dupe\n",
			 is_dupe ? "" : " not");
  }
  return is_dupe;
}

/* Return 1 if the relation is probably a duplicate of a relation found
   "earlier", and 0 if it is probably not a duplicate */
int
relation_is_duplicate(relation const& rel, const int nb_threads,
                      sieve_info & si)
{
  /* If the special-q does not fit in an unsigned long, we assume it's not a
     duplicate and just move on */
  if (!mpz_fits_ulong_p(si.doing.p)) {
    return 0;
  }

  /* If any large prime doesn't fit in an unsigned long, then we assume
   * we don't have a duplicate */
  for(int side = 0 ; side < 2 ; side++) {
      for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
          if (!mpz_fits_ulong_p(rel.sides[side][i].p))
              return false;
      }
  }

  for(int side = 0 ; side < 2 ; side++) {
      for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
          if (check_one_prime(rel.sides[side][i].p, side, rel, nb_threads, si))
              return true;
      }
  }

  return false;
}
