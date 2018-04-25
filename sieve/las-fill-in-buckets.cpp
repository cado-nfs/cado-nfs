#include "cado.h"

#include "fb.hpp"
#include "utils.h"           /* lots of stuff */
#include "bucket.hpp"
#include "modredc_ul.h"
#include "modredc_2ul2.h"
#include "threadpool.hpp"
#include "las-config.h"
#include "las-types.hpp"
#include "las-coordinates.hpp"
#include "las-debug.hpp"
#include "las-arith.hpp"
#include "las-qlattice.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-norms.hpp"
#include "las-smallsieve.hpp"
#include "las-plattice.hpp"

#ifdef USE_CACHEBUFFER
#include "cachebuf.h"
#endif

// FIXME: this function of las.cpp should be somewhere
void * process_bucket_region(timetree_t&, thread_data *th);

/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/**************************************************************************
 * Global DEFINEs for fill_in_buckets, fill_in_k_buckets, fill_in_m_buckets 
 **************************************************************************/

#ifdef HAVE_SSE2							
#define FILL_BUCKET_PREFETCH(PT) do {				\
    _mm_prefetch((char *)(PT), _MM_HINT_T0);			\
  } while (0)
#else
#define FILL_BUCKET_PREFETCH(PT)
#endif


#ifdef USE_CACHEBUFFER
DECLARE_CACHE_BUFFER(bucket_update_shorthint_t, 256)
#endif

#if 0

/* Returns -(-1)^neg num/den (mod p). If gcd(den, p) > 1, returns p. */
static inline unsigned long
compute_1_root_ul(const unsigned long p, const unsigned long num,
                  const unsigned long den, const int k, const int neg)
{
  modulusredcul_t m;
  residueredcul_t rn, rd;

  modredcul_initmod_ul(m, p);
  modredcul_init_noset0(rn, m);
  modredcul_init_noset0(rd, m);

  modredcul_set_ul(rn, num, m);
  modredcul_set_ul(rd, den, m);
  int ok = modredcul_inv(rd, rd, m);
  if (!ok) {
    /* FIXME: p could be a prime power */
    return p;
  }
  modredcul_mul(rn, rn, rd, m);
  if (!neg)
    modredcul_neg(rn, rn, m);

  for (int i = 0; i < k; i++)
    modredcul_div2(rn, rn, m);
  unsigned long r = modredcul_get_ul(rn, m);

  modredcul_clear(rn, m);
  modredcul_clear(rd, m);
  modredcul_clearmod(m);
  return r;
}

static inline unsigned long
compute_1_root_mpz(const unsigned long p, const mpz_t f0, const mpz_t f1)
{
  return compute_1_root_ul(p, mpz_fdiv_ui(f0, p), mpz_fdiv_ui(f1, p), 0, 0);
}

struct contexts_s {
  modredc2ul2_batch_Q_to_Fp_context_t *context_2ul2;
  unsigned long k;
  unsigned long num_ul, den_ul;
  int neg;
};

static inline void
modredc2ul2_set_mpz(modintredc2ul2_t r, const mpz_t a)
{
  unsigned long t[MODREDC2UL2_SIZE];
  size_t n;
  mpz_export (t, &n, -1, sizeof(unsigned long), 0, 0, a);
  ASSERT_ALWAYS(n <= MODREDC2UL2_SIZE);
  modredc2ul2_intset_uls(r, t, n);
}

static void
init_Q_to_Fp_context(struct contexts_s *contexts, mpz_poly_srcptr f)
{
  const int use_Q_to_Fp = 0;
  contexts->context_2ul2 = NULL;
  contexts->num_ul = contexts->den_ul = 0;

  if (!use_Q_to_Fp || f->deg != 1)
    return;

  contexts->k = mpz_scan1(f->coeff[1], 0);
  mpz_t odd_den_mpz;
  mpz_init(odd_den_mpz);
  mpz_tdiv_q_2exp(odd_den_mpz, f->coeff[1], contexts->k);

  /* If exactly one of the two coefficients is negative, we need to flip
     the sign of the root (mod p). The odd_den and num variables contain
     absolute values. */
  contexts->neg = (mpz_sgn(f->coeff[0]) < 0) != (mpz_sgn(f->coeff[1]) < 0);

  modintredc2ul2_t num, odd_den;
  modredc2ul2_intinit(num);
  modredc2ul2_intinit(odd_den);
  /* Set num = |f->coeff[0]|, odd_den = |odd_den_mpz| */
  modredc2ul2_set_mpz(num, f->coeff[0]);
  modredc2ul2_set_mpz(odd_den, odd_den_mpz);

  if (mpz_fits_ulong_p(f->coeff[0]) && mpz_fits_ulong_p(odd_den_mpz))
  {
    contexts->num_ul = mpz_get_ui(f->coeff[0]);
    contexts->den_ul = mpz_get_ui(odd_den_mpz);
  } else if (mpz_sizeinbase(f->coeff[0], 2) >= MODREDC2UL2_MINBITS &&
             mpz_sizeinbase(f->coeff[0], 2) <= MODREDC2UL2_MAXBITS &&
             mpz_sizeinbase(odd_den_mpz, 2) >= MODREDC2UL2_MINBITS &&
             mpz_sizeinbase(odd_den_mpz, 2) <= MODREDC2UL2_MAXBITS )
  {
    contexts->context_2ul2 = modredc2ul2_batch_Q_to_Fp_init (num, odd_den);
  }

  mpz_clear(odd_den_mpz);
  modredc2ul2_intclear(num);
  modredc2ul2_intclear(odd_den);
}

static void
clear_Q_to_Fp_context(struct contexts_s *contexts)
{
  if (contexts->context_2ul2 != NULL)
    modredc2ul2_batch_Q_to_Fp_clear(contexts->context_2ul2);
  contexts->context_2ul2 = NULL;
  contexts->num_ul = contexts->den_ul = 0;
}


/* Compute up to n roots of the transformed polynomial, either by tranforming
   the root stored in the factor base, or, if the polynomial has degree 1 and
   its coefficients fit into unsigned long, by computing the root from the
   transformed polynomial.
   For each root, the modulus (i.e., factor base prime) and
   the root are stored in p[i] and r[i], resp.
   Returns the number of transformed roots produced, which may be less than n
   if we reach the end of the factor base. */
static inline size_t
transform_n_roots(unsigned long *p, unsigned long *r, fb_iterator t,
                  const size_t n, const int side, sieve_info const & si,
                  struct contexts_s *contexts)
{
  size_t i;
  int ok = 0;
  if (contexts->den_ul == 0 && contexts->context_2ul2 == NULL) {
    /* Do the old transform with separate inverse modulo each p */
    for (i = 0; !fb_iterator_over(t) && i < n; i++, fb_iterator_next(t)) {
      fbprime_t current_p = fb_iterator_get_p(t);
      ASSERT_ALWAYS (current_p & 1);
      fbroot_t R = fb_iterator_get_r(t);
      p[i] = current_p;
      r[i] = fb_root_in_qlattice(current_p, R, t->fb->invp, si);
    }
    return i;
  }

  /* Collect primes from fb, store up to n namy of them in p, and the actual
     number of stored primes in i */
  for (i = 0; !fb_iterator_over(t) && i < n; i++, fb_iterator_next(t)) {
    p[i] = fb_iterator_get_p(t);
  }

  if (contexts->den_ul != 0) {
    ok = modredcul_batch_Q_to_Fp(r, contexts->num_ul, contexts->den_ul, contexts->k, p, i);
    if (ok && !contexts->neg) {
      /* We computed r[i] = |g_0| / |g_1| (mod p[i]). If both g0, g1 have the same sign,
         then the root is -g_0/g_1, and we need a sign flip. */
      for (size_t j = 0; j < i; j++)
        if (r[j] != 0)
          r[j] = p[j] - r[j];
    }
  } else if (contexts->context_2ul2 != NULL) {
    ok = modredc2ul2_batch_Q_to_Fp(r, contexts->context_2ul2, contexts->k, !contexts->neg, p, i);
  }

  if(UNLIKELY(!ok)) {
    /* Modular inverse did not exists, i.e., at least one of the factor base
       entries was not coprime to den. Do factor base entries one at a time,
       handling the non-coprime case properly. This should happen at most a
       few times per special-q. */
    mpz_poly_srcptr f = si.sides[side].fij;
    for (size_t j = 0; j < i; j++)
      r[j] = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
  }

#if 0
    /* Very expensive but thorough test */
    for (size_t j = 0; j < i; j++) {
      mpz_poly_srcptr f = si.sides[side].fij;
      unsigned long tr = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
      ASSERT_ALWAYS(tr == r[j]);
    }
#endif

  return i;
}
#endif

/* {{{ Big question: shall we enable bucket-sieving for powers ?
 *
 * There are several difficulties, in fact. One rationale that yields a
 * straight "no" answer is that such primes make very little difference
 * to the smooth part, so we'd better skip them anyway.
 *
 * But it's not the hardest thing.
 *
 * For the small sieve, we create the small_sieve_data from the factor
 * base entries, and we compute the logp accordingly, per entry.
 *
 * For the bucket-sieve, we use the fact that the factor base is sorted
 * in increasing log(p) order, and we create slices with ranges of primes
 * that have the same round(scale*log(p)).
 *
 * Currently, the factor base is sorted by q=p^k. A power that makes the
 * p-valuation go from p^k0 to p^k1 contributes
 * round(k1*log(p))-round(k0*log(p)). Therefore, sorting by q does not
 * mean that log(p)'s are sorted, and we're in trouble because when we
 * take powers aboard in a slice, their log(p) value is not correctly
 * represented. 
 * 
 * Previously, we had the behaviour of setting powlim to bucket_thresh-1,
 * effectively preventing powers from appearing in the bucket-sieve.
 *
 * Now powlim is a factor base parameter, and bucket_thresh comes later,
 * so such a default does not work.
 *
 * The strategy we take here is that *if* we see powers here (and we know
 * that will happen only for the fairly rare general entries), then we do
 * something special:
 *  - either we say that we skip over this entry
 *  - or we defer to apply_buckets the expensive computation of a proper
 *    logp value.
 *
 * Currently we do the former. The latter would be viable since only a
 * small fraction of the apply_one_bucket time is devoted to dealing with
 * general entries, so we could imagine having a branch in there for
 * dealing with them. But that would be quite painful. Furthermore, it
 * would then be mandatory to split the entries with same q, same p, but
 * different k0,k1 pairs (we do encounter these), so that the hint would
 * still open the possibility to infer the value of log(p).
 *
 *
 * Note that it would not be possible to sidestep the issue by sorting
 * the vectors of entries by (k1-k0)*log(p) (which would make a
 * difference only for general entries anyway). This is because even
 * sorting by increasing (k1-k0)*log(p) does not guarantee that
 * round(s*k1*log(p))-round(s*k0*log(p)) increases. (counter-example:
 * s=1, k1*log(p)=0.51, k0*log(p)=0.49 diff=0.02 round-round=1
 *      k1*log(p)=1.49, k0*log(p)=0.51 diff=0.98 round-round=0
 * )
 * }}} */
template<class FB_ENTRY_TYPE>
inline bool discard_power_for_bucket_sieving(FB_ENTRY_TYPE const &) {
    return false;
}
template<>
inline bool discard_power_for_bucket_sieving<fb_entry_general>(fb_entry_general const & e) {
    return e.k > 1;
}

template <class FB_ENTRY_TYPE>
plattices_vector_t
make_lattice_bases(fb_slice<FB_ENTRY_TYPE> const & slice,
        const qlattice_basis &basis,
    const int logI, const sublat_t &sublat)
{
  typename FB_ENTRY_TYPE::transformed_entry_t transformed;
  /* Create a transformed vector and store the index of the fb_slice we
   * currently transform */

  plattices_vector_t result(slice.get_index(), slice.get_weight());
  slice_offset_t i_entry = 0;
  for (auto const & e : slice) {
      increment_counter_on_dtor<slice_offset_t> _dummy(i_entry);
    if (!basis.is_coprime_to(e.p))
      continue;
    if (discard_power_for_bucket_sieving(e))
        continue;
    e.transform_roots(transformed, basis);
    for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
      const fbroot_t r = transformed.get_r(i_root);
      const bool proj = transformed.get_proj(i_root);
      /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
         locations with p | gcd(i,j). */
      if (LIKELY(!proj || r == 0)) {
        plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, logI);
        plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, logI, sublat);
        // Skip (0,0) unless we have sublattices.
        if (!sublat.m)
          ple.next();
        if (LIKELY(pli.a0 != 0)) {
          result.push_back(ple);
        }
      }
    }
  }
  /* This is moved, not copied */
  return result;
}

struct helper_functor_make_lattice_bases {
    int side, level;
    sieve_info const & si;
    precomp_plattice_t & precomp_plattice;
    template<typename T>
        void operator()(T const & s) {
            precomp_plattice.push(side, level, make_lattice_bases(s, si.qbasis, si.conf.logI, si.conf.sublat));
        }
};

void fill_in_buckets_prepare_precomp_plattice(
                        int side,
                        int level,
                        sieve_info const & si MAYBE_UNUSED,
                        precomp_plattice_t & precomp_plattice)
{
    helper_functor_make_lattice_bases F { side, level, si, precomp_plattice };
    si.sides[side].fbs->get_part(level).foreach_slice(F);
}

// At top level, the fill-in of the buckets must interleave
// the root transforms and the FK walks, otherwise we spend a lot of time
// doing nothing while waiting for memory.
// Consequence: we duplicate here the code of make_lattice_bases in fb.cpp
// FIXME: find a way to refactor that.
//
// With Sublat, this function can have two modes:
//   - process the given slice, and store the corresponding FK-basis in
//     precomp_slice for later use.
//   - use the pre-processed precomputed FK_basis.
// If (and only if) we are dealing with (i,j) == (0,1) mod m,
// we are in the second mode.
//
//
// FIXME FIXME FIXME: tons of duplicated code, here!!!
// But putting if() in critical loops can kill performance (I tried...)

template <int LEVEL, class FB_ENTRY_TYPE>
void
fill_in_buckets_toplevel_sublat(bucket_array_t<LEVEL, shorthint_t> &orig_BA,
                sieve_info const & si MAYBE_UNUSED,
                fb_slice<FB_ENTRY_TYPE> const & slice,
                plattices_dense_vector_t * precomp_slice,
                where_am_I & w)
{
  ASSERT_ALWAYS(si.conf.sublat.m != 0);
  bool first_sublat = si.conf.sublat.i0 == 0 && si.conf.sublat.j0 == 1;
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  slice_index_t slice_index = slice.get_index();
  if (!first_sublat) {
    ASSERT(slice_index == precomp_slice->get_index());
  }

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  typename FB_ENTRY_TYPE::transformed_entry_t transformed;

  // FIXME: A LOT OF DUPLICATED CODE, HERE!!!
  if (first_sublat) {
    slice_offset_t i_entry = 0;
    for (auto const & e : slice) {
      increment_counter_on_dtor<slice_offset_t> _dummy(i_entry);
      if (!si.qbasis.is_coprime_to(e.p))
        continue;
      if (discard_power_for_bucket_sieving(e))
          continue;
      e.transform_roots(transformed, si.qbasis);
      for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
        const fbroot_t r = transformed.get_r(i_root);
        const bool proj = transformed.get_proj(i_root);
        /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
           locations with p | gcd(i,j). */
        if (LIKELY(!proj || r == 0)) {
          plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, si.conf.logI);
          // In sublat mode, save it for later use
          plattice_info_dense_t plid(pli, i_entry);
          precomp_slice->push_back(plid);

          plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI, si.conf.sublat);

          if (plattice_enumerate_finished<LEVEL>(ple.get_x()))
            continue;
          if (LIKELY(pli.a0 != 0)) {
            const slice_offset_t hint = ple.get_hint();
            ASSERT(hint == i_entry);
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice.get_prime(hint); 
            WHERE_AM_I_UPDATE(w, p, p);
#else
            const fbprime_t p = 0;
#endif
            // Handle the rare special cases
            const uint32_t I = si.I;
            if ((UNLIKELY(ple.get_inc_c() == 1 && ple.get_bound1() == I - 1))
                ||
                (UNLIKELY(ple.get_inc_c() == I && ple.get_bound1() == I))) {
              continue;
            }

            /* Now, do the real work: the filling of the buckets */
            // Without sublattices, we test (very basic) coprimality,
            while (!plattice_enumerate_finished<LEVEL>(ple.get_x())) {
              BA.push_update(ple.get_x(), p, hint, slice_index, w);
              ple.next();
            }
          } 
        }
      }
    }
  } else { // Use precomputed FK-basis
    for (unsigned int i = 0; i < precomp_slice->size(); ++i) {
      plattice_info_t pli = (*precomp_slice)[i].unpack(si.conf.logI);
      slice_offset_t i_entry = (*precomp_slice)[i].get_hint();

      plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI, si.conf.sublat);

      if (plattice_enumerate_finished<LEVEL>(ple.get_x()))
        continue;
      if (LIKELY(pli.a0 != 0)) {
        const slice_offset_t hint = ple.get_hint();
        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        const fbprime_t p = slice.get_prime(hint); 
        WHERE_AM_I_UPDATE(w, p, p);
#else
        const fbprime_t p = 0;
#endif

        // Handle the rare special cases
        const uint32_t I = si.I;
        if ((UNLIKELY(ple.get_inc_c() == 1 && ple.get_bound1() == I - 1))
            ||
            (UNLIKELY(ple.get_inc_c() == I && ple.get_bound1() == I))) {
          continue;
        }

        /* Now, do the real work: the filling of the buckets */
        // Without sublattices, we test (very basic) coprimality,
        // otherwise not atm. FIXME!
        while (!plattice_enumerate_finished<LEVEL>(ple.get_x())) {
          BA.push_update(ple.get_x(), p, hint, slice_index, w);
          ple.next();
        }
      } 
    }
  }
  orig_BA.move(BA);
}

template <int LEVEL, class FB_ENTRY_TYPE>
void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, shorthint_t> &orig_BA,
                sieve_info const & si MAYBE_UNUSED,
                fb_slice<FB_ENTRY_TYPE> const & slice,
                plattices_dense_vector_t * /* unused */,
                where_am_I & w)
{
  ASSERT_ALWAYS(!si.conf.sublat.m);

  bool first_reg = true;
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  slice_index_t slice_index = slice.get_index();

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  typename FB_ENTRY_TYPE::transformed_entry_t transformed;

  slice_offset_t i_entry = 0;

  for (auto const & e : slice) {
    increment_counter_on_dtor<slice_offset_t> _dummy(i_entry);
    if (!si.qbasis.is_coprime_to(e.p))
      continue;
    if (discard_power_for_bucket_sieving(e))
        continue;
    e.transform_roots(transformed, si.qbasis);
    for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
      const fbroot_t r = transformed.get_r(i_root);
      const bool proj = transformed.get_proj(i_root);
      /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
         locations with p | gcd(i,j). */
      if (LIKELY(!proj || r == 0)) {
        plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, si.conf.logI);
  
        plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI);

        // Skip (i,j)=(0,0)
        ple.next();

        // what does pli.a0 == 0 correspond to ?
        if (LIKELY(pli.a0 != 0)) {
          const slice_offset_t hint = ple.get_hint();
          WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
          const fbprime_t p = slice.get_prime(hint); 
          WHERE_AM_I_UPDATE(w, p, p);
#else
          const fbprime_t p = 0;
#endif

          // Handle the rare special cases
          const uint32_t I = si.I;
          if (UNLIKELY(ple.get_inc_c() == 1 && ple.get_bound1() == I - 1)) {
            // Projective root: only update is at (1,0).
            if (first_reg) {
              uint64_t x = 1 + (I >> 1);
              BA.push_update(x, p, hint, slice_index, w);
            }
            continue;
          }
          if (UNLIKELY(ple.get_inc_c() == I && ple.get_bound1() == I)) {
            // Root=0: only update is at (0,1).
            if (first_reg) {
              uint64_t x = I + (I >> 1);
              BA.push_update(x, p, hint, slice_index, w);
            }
            continue;
          }

          /* Now, do the real work: the filling of the buckets */
          while (!plattice_enumerate_finished<LEVEL>(ple.get_x())) {
            if (LIKELY(ple.probably_coprime()))
              BA.push_update(ple.get_x(), p, hint, slice_index, w);
            ple.next();
          }
        }
      } 
    }
  }
  orig_BA.move(BA);
}



/* {{{ */
template <int LEVEL>
  void
fill_in_buckets_lowlevel(
    bucket_array_t<LEVEL, shorthint_t> &orig_BA,
    sieve_info const & si MAYBE_UNUSED,
    plattices_vector_t & plattices_vector,
    bool first_reg,
    where_am_I & w)
{
    /* The timer stuff is dealt with by the caller */
  const slice_index_t slice_index = plattices_vector.get_index();
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  for (auto & ple : plattices_vector) {
    // Work with a copy, otherwise we don't get all optimizations.
    // Maybe with a wise use of the 'restrict' keyword, we might get
    // what we want, but this is C++11, anyway.
    plattice_enumerate_t pl(ple);

    const slice_offset_t hint = pl.get_hint();
    WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
    /* this is a bit expensive, since we're scanning all parts.
     * Fortunately it's only a debug call anyway. */
    fb_slice_interface const & slice = (*si.sides[w.side].fbs)[slice_index];
    const fbprime_t p = slice.get_prime(hint); 
    WHERE_AM_I_UPDATE(w, p, p);
#else
    const fbprime_t p = 0;
#endif

    // Handle the rare special cases
    const uint32_t I = si.I;
    if (UNLIKELY(pl.get_inc_c() == 1 && pl.get_bound1() == I - 1)) {
        // Projective root: only update is at (1,0).
        if (!si.conf.sublat.m && first_reg) {
            uint64_t x = 1 + (I >> 1);
            BA.push_update(x, p, hint, slice_index, w);
        }
        continue;
    }
    if (UNLIKELY(pl.get_inc_c() == I && pl.get_bound1() == I)) {
        // Root=0: only update is at (0,1).
        if (!si.conf.sublat.m && first_reg) {
            uint64_t x = I + (I >> 1);
            BA.push_update(x, p, hint, slice_index, w);
        }
        continue;
    }

    /* Now, do the real work: the filling of the buckets */
    // Without sublattices, we test (very basic) coprimality,
    // otherwise not atm. FIXME!
    if (!si.conf.sublat.m) {
        while (!plattice_enumerate_finished<LEVEL>(pl.get_x())) {
            if (LIKELY(pl.probably_coprime()))
                BA.push_update(pl.get_x(), p, hint, slice_index, w);
            pl.next();
        }
    } else {
        while (!plattice_enumerate_finished<LEVEL>(pl.get_x())) {
            BA.push_update(pl.get_x(), p, hint, slice_index, w);
            pl.next();
        }
    }

    // save current position, and prepare for next area.
    ple.set_x(pl.get_x());
    ple.advance_to_next_area(LEVEL);
  } 
  orig_BA.move(BA);
}

class fill_in_buckets_parameters: public task_parameters {
public:
  thread_workspaces &ws;
  const int side;
  sieve_info const & si;
  const fb_slice_interface * slice;
  plattices_vector_t * plattices_vector; // content changed during fill-in
  plattices_dense_vector_t * plattices_dense_vector; // for sublat
  const uint32_t first_region0_index;
  fill_in_buckets_parameters(thread_workspaces &_ws, const int _side,
          sieve_info const & _si, const fb_slice_interface *_slice,
          plattices_vector_t *_platt, plattices_dense_vector_t *_dplatt,
          const uint32_t _reg0)
  : ws(_ws), side(_side), si(_si), slice(_slice),
    plattices_vector(_platt), plattices_dense_vector(_dplatt),
    first_region0_index(_reg0)
  {}
};

#if __cplusplus >= 201103L
/* short of a better solution. I know some exist, but it seems way
 * overkill to me.
 *
 * This needs constexpr, though... So maybe I could use a more powerful
 * C++11 trick after all.
 */
#define PREPARE_TEMPLATE_INST_NAMES(F, suffix)				\
    template<int>							\
    struct CADO_CONCATENATE(F, _name) {};				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 0);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 1);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 2);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 3);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 4);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 5);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 6);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 7);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 8);				\
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 9)

#define PREPARE_TEMPLATE_INST_NAME(F, suffix, k)			\
    template<>								\
    struct CADO_CONCATENATE(F, _name)<k> {				\
        static constexpr const char * value = #F "<" #k ">" suffix;	\
    }

/* By tweaking the "" argument, it is possible to have these names
 * embody a suffix like " on side ", so that it's possible tu run
 * parametric timer slots.
 */
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice_internal, "");
// PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice, "");
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_side, "");
PREPARE_TEMPLATE_INST_NAMES(downsort_tree, "");

#define TEMPLATE_INST_NAME(x,y) CADO_CONCATENATE(x, _name)<y>::value
#else
#define TEMPLATE_INST_NAME(x,y) #x " (template)"
#endif

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template<int LEVEL>
task_result *
fill_in_buckets_one_slice_internal(const worker_thread * worker, const task_parameters * _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    // ACTIVATE_TIMER(worker->timer);
    ASSERT_ALWAYS(worker->timer.running());
    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(worker->timer, TEMPLATE_INST_NAME(fill_in_buckets_one_slice_internal, LEVEL));

    where_am_I w;
    WHERE_AM_I_UPDATE(w, psi, & param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, param->first_region0_index);

    try {
        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        bucket_array_t<LEVEL, shorthint_t> &BA =
            param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);

        time_bubble_chaser tt(worker->rank(), time_bubble_chaser::FIB,
                {
                param->side,
                LEVEL,
                param->ws.rank_BA(param->side, BA),
                (int) param->plattices_vector->get_index()
                }
            );

        /* Fill the buckets */
        try {
            fill_in_buckets_lowlevel<LEVEL>(BA, param->si,
                    * param->plattices_vector,
                    (param->first_region0_index == 0), w);
        } catch(buckets_are_full & e) {
            worker->push_chart_item(tt.put());
            param->ws.release_BA(param->side, BA);
            throw e;
        }
        /* Release bucket array again */
        worker->push_chart_item(tt.put());

        param->ws.release_BA(param->side, BA);
    } catch(buckets_are_full & e) {
        e.side = param->side;
        delete param;
        throw e;
    }
    delete param;
    return new task_result;
}

/* we really wish to have a single timing slot for all the instantiations
 * of fill_in_buckets_toplevel_wrapper */
static tdict::slot timer_slot_for_fibt_wrapper("fill_in_buckets_toplevel_wrapper");

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template<int LEVEL, class FB_ENTRY_TYPE>
task_result *
fill_in_buckets_toplevel_wrapper(const worker_thread * worker MAYBE_UNUSED, const task_parameters * _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    // ACTIVATE_TIMER(worker->timer);
    ASSERT_ALWAYS(worker->timer.running());
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     *
     */
    // CHILD_TIMER(worker->timer, TEMPLATE_INST_NAME(fill_in_buckets_one_slice, LEVEL));
    

    timetree_t::accounting_child local_timer_sentry(worker->timer, timer_slot_for_fibt_wrapper);
    // TIMER_CATEGORY(worker->timer, sieving(param->side));

    where_am_I w;
    WHERE_AM_I_UPDATE(w, psi, & param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, shorthint_t> &BA = param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);

        ASSERT(param->slice);
        time_bubble_chaser tt(worker->rank(), time_bubble_chaser::FIB,
                {
                param->side,
                LEVEL,
                param->ws.rank_BA(param->side, BA),
                (int) param->slice->get_index()
                }
                );
        fill_in_buckets_toplevel<LEVEL,FB_ENTRY_TYPE>(BA, param->si,
                *dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice),
                param->plattices_dense_vector, w);

        /* Release bucket array again */
        worker->push_chart_item(tt.put());
        param->ws.release_BA(param->side, BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full & e) {
        e.side = param->side;
        /* release_BA not needed, because failures in reserve() release
         * the lock before throwing */
        delete param;
        throw e;
    }
}
/* same for sublat */
template<int LEVEL, class FB_ENTRY_TYPE>
task_result *
fill_in_buckets_toplevel_sublat_wrapper(const worker_thread * worker MAYBE_UNUSED, const task_parameters * _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    ASSERT_ALWAYS(worker->timer.running());
    timetree_t::accounting_child local_timer_sentry(worker->timer, timer_slot_for_fibt_wrapper);

    where_am_I w;
    WHERE_AM_I_UPDATE(w, psi, & param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        bucket_array_t<LEVEL, shorthint_t> &BA = param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
        time_bubble_chaser tt(worker->rank(), time_bubble_chaser::FIB,
                {
                param->side,
                LEVEL,
                param->ws.rank_BA(param->side, BA),
                (int) param->slice->get_index()
                }
                );
        ASSERT(param->slice);
        fill_in_buckets_toplevel_sublat<LEVEL,FB_ENTRY_TYPE>(BA, param->si,
                *dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice),
                param->plattices_dense_vector, w);
        /* Release bucket array again */
        worker->push_chart_item(tt.put());
        param->ws.release_BA(param->side, BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full & e) {
        e.side = param->side;
        delete param;
        throw e;
    }
}

/* Whether or not we want fill_in_buckets_one_slice to be templatized
 * both for LEVEL and n is not clear. At some point, we're doing code
 * bloat for almost nothing.
 *
 * Now given the code below, it's easy enough to arrange so that we go
 * back to the virtual base fb_slice_interface.
 */
template<int LEVEL>
struct push_slice_to_task_list {
    thread_pool & pool;
    fill_in_buckets_parameters model;
    push_slice_to_task_list(thread_pool&pool, fill_in_buckets_parameters const & m) : pool(pool), model(m) {}
    size_t pushed = 0;
    template<typename T>
    void operator()(T const & s) {
        fill_in_buckets_parameters *param = new fill_in_buckets_parameters(model);
        param->slice = &s;
        typedef typename T::entry_t entry_t;
        task_function_t f = fill_in_buckets_toplevel_wrapper<LEVEL, entry_t>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};
template<int LEVEL>
struct push_slice_to_task_list_saving_precomp {
    thread_pool & pool;
    fill_in_buckets_parameters model;
    /* precomp_plattice_dense_t == std::vector<plattices_dense_vector_t *> */
    precomp_plattice_dense_t & Vpre;
    bool is_first;
    size_t pushed = 0;
    push_slice_to_task_list_saving_precomp(thread_pool&pool, fill_in_buckets_parameters const & m, std::vector<plattices_dense_vector_t *> & Vpre, bool is_first) : pool(pool), model(m), Vpre(Vpre), is_first(is_first) {}
    template<typename T>
    void operator()(T const & s) {
        plattices_dense_vector_t * pre;
        if (is_first) {
            ASSERT_ALWAYS(Vpre.size() == pushed);
            pre = new plattices_dense_vector_t(s.get_index());

            /* TODO: This looks pretty awful. We should reserve the array
             * beforehand, this should save a constant factor on the
             * amortized constant cost of storing an element to a vector.
             */
            Vpre.push_back(pre);
        } else {
            pre = Vpre[pushed];
        }
        fill_in_buckets_parameters *param = new fill_in_buckets_parameters(model);
        param->slice = &s;
        param->plattices_dense_vector = pre;
        typedef typename T::entry_t entry_t;
        task_function_t f = fill_in_buckets_toplevel_sublat_wrapper<LEVEL, entry_t>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};

template <int LEVEL>
static void
fill_in_buckets_one_side(timetree_t& timer, thread_pool &pool, thread_workspaces &ws, sieve_info & si, const int side)
{
    const fb_factorbase::slicing * fbs = si.sides[side].fbs;

    CHILD_TIMER_PARAMETRIC(timer, "fill_in_buckets_one_side on side ", side, "");
    TIMER_CATEGORY(timer, sieving(side));

    slice_index_t slices_pushed = 0;

    fill_in_buckets_parameters model(ws, side, si, NULL, NULL, NULL, 0);

    if (!si.conf.sublat.m) {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_wrapper */
        push_slice_to_task_list<LEVEL> F(pool, model);
        fbs->get_part(LEVEL).foreach_slice(F);
        slices_pushed = F.pushed;
    } else {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_sublat_wrapper */
        auto & Vpre(si.sides[side].precomp_plattice_dense);
        bool is_first = si.conf.sublat.i0 == 0 && si.conf.sublat.j0 == 1;
        push_slice_to_task_list_saving_precomp<LEVEL> F(pool, model, Vpre, is_first);
        fbs->get_part(LEVEL).foreach_slice(F);
        slices_pushed = F.pushed;
    }

#if 0
        /* Process all slices in this factor base part */
        const fb_slice_interface *slice;
        for (slice_index_t slice_index = fb->get_first_slice_index();
             (slice = fb->get_slice(slice_index)) != NULL;
             slice_index++) {
            plattices_dense_vector_t * pre = NULL;
            if (si.conf.sublat.m) {
                pre = new plattices_dense_vector_t(slice_index);
                si.sides[side].precomp_plattice_dense.push_back(pre);
            }
            fill_in_buckets_parameters *param = new fill_in_buckets_parameters(ws, side, si, slice, 
                    NULL, pre, 0);
            pool.add_task(fill_in_buckets_one_slice<LEVEL>, param, 0, 0, (double)slice->get_weight());
            slices_pushed++;
        }
    } else {
        /* We are in sublat mode, and we have already done the first
         * congruence. Re-use precomputed FK-basis */
        for (unsigned int i=0; i < si.sides[side].precomp_plattice_dense.size(); i++){
            plattices_dense_vector_t * pre = si.sides[side].precomp_plattice_dense[i];
            const fb_slice_interface * slice = fb->get_slice(pre->get_index());
            fill_in_buckets_parameters *param = new fill_in_buckets_parameters(ws, side, si, slice,
                    NULL, pre, 0);
            pool.add_task(fill_in_buckets_one_slice<LEVEL>, param, 0, 0, (double)slice->get_weight());
            slices_pushed++;
        }
    }
#endif

    /* we need to check for exceptions due to bucket updates. Because
     * we've pushed all slides at this point, we have no option but to
     * wait for all threads to finish. (or maybe pull back some of the
     * tasks before they're grabbed by worker threads -- that could be an
     * optimization).
     */
    std::vector<buckets_are_full> exceptions_to_throw;

    for (slice_index_t slices_completed = 0; slices_completed < slices_pushed; slices_completed++) {
          task_result *result = pool.get_result();
          delete result;
          /* want to check possible exceptions, too */
          buckets_are_full * e = dynamic_cast<buckets_are_full*>(pool.get_exception());
          if (e) {
              exceptions_to_throw.push_back(*e);
              delete e;
          }
    }
    if (!exceptions_to_throw.empty())
        throw *std::max_element(exceptions_to_throw.begin(), exceptions_to_throw.end());

    /* actually we also want to check that not even the last thread has
     * overrun its bucket pointers.
     *
     * the call below throws an exception if we're overfull.
     */
    ws.buckets_max_full<LEVEL, shorthint_t>();

    pool.accumulate_and_clear_active_time(*timer.current);
    SIBLING_TIMER(timer, "worker thread wait time");
    TIMER_CATEGORY(timer, thread_wait());
    pool.accumulate_and_reset_wait_time(*timer.current);
}


void fill_in_buckets(timetree_t& timer, thread_pool &pool, thread_workspaces &ws, sieve_info & si, int side)
{
    // per se, we're not doing anything here.
    // CHILD_TIMER(timer, __func__);
    plattice_enumerate_t::set_masks(si.conf.logI);
    switch (si.toplevel) {
        case 1:
            plattice_enumerate_area<1>::value = plattice_x_t(si.J) << si.conf.logI;
            fill_in_buckets_one_side<1>(timer, pool, ws, si, side);
            break;
        case 2:
            plattice_enumerate_area<2>::value = plattice_x_t(si.J) << si.conf.logI;
            fill_in_buckets_one_side<2>(timer, pool, ws, si, side);
            break;
        case 3:
            plattice_enumerate_area<3>::value = plattice_x_t(si.J) << si.conf.logI;
            fill_in_buckets_one_side<3>(timer, pool, ws, si, side);
            break;
        default:
            ASSERT_ALWAYS(0);
    }
}

/* }}} */


// first_region0_index is a way to remember where we are in the tree.
// The depth-first is a way to process all the the regions of level 0 in
// increasing order of j-value.
// first_region0_index * nb_lines_per_region0 therefore gives the j-line
// where we are. This is what is called N by WHERE_AM_I and friends.
template <int LEVEL>
void
downsort_tree(
    timetree_t& timer,
    uint32_t bucket_index,
    uint32_t first_region0_index,
    thread_workspaces &ws,
    thread_pool &pool,
    sieve_info & si,
    precomp_plattice_t const & precomp_plattice)
{
  CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort_tree, LEVEL));
  TIMER_CATEGORY(timer, sieving_mixed());
  ASSERT_ALWAYS(LEVEL > 0);

  where_am_I w;
  WHERE_AM_I_UPDATE(w, psi, & si);
  WHERE_AM_I_UPDATE(w, N, first_region0_index);

  double max_full MAYBE_UNUSED = 0.0;

  for (int side = 0; side < 2; ++side) {
    if (!si.sides[side].fb) continue;
    WHERE_AM_I_UPDATE(w, side, side);
    CHILD_TIMER_PARAMETRIC(timer, "side ", side, "");
    TIMER_CATEGORY(timer, sieving(side));
    /* FIRST: Downsort what is coming from the level above, for this
     * bucket index */
    // All these BA are global stuff; see reservation_group.
    // We reserve those where we write, and access the ones for
    // reading without reserving. We require that things at level
    // above is finished before entering here.
    bucket_array_t<LEVEL, longhint_t> & BAout =
      ws.reserve_BA<LEVEL, longhint_t>(side);
    BAout.reset_pointers();
    // This is a fake slice_index. For a longhint_t bucket, each update
    // contains its own slice_index, directly used by apply_one_bucket
    // and purge.
    time_bubble_chaser tt(0, time_bubble_chaser::DS, { side, LEVEL, 0, 0});
    BAout.add_slice_index(0);
    // The data that comes from fill-in bucket at level above:
    {
      const bucket_array_t<LEVEL+1,shorthint_t> * BAin_ptr
        = ws.cbegin_BA<LEVEL+1,shorthint_t>(side);
      while (BAin_ptr != ws.cend_BA<LEVEL+1,shorthint_t>(side)) {
        downsort<LEVEL+1>(BAout, *BAin_ptr, bucket_index, w);
        BAin_ptr++;
      }
    }
    pool.push_chart_item(0, tt.put());

    const int toplevel = si.toplevel;
    if (LEVEL < toplevel - 1) {
      // What comes from already downsorted data above:
      const bucket_array_t<LEVEL+1,longhint_t> * BAin_ptr
        = ws.cbegin_BA<LEVEL+1,longhint_t>(side);
      while (BAin_ptr != ws.cend_BA<LEVEL+1,longhint_t>(side)) { 
        downsort<LEVEL+1>(BAout, *BAin_ptr, bucket_index, w);
        BAin_ptr++;
      }
    }
    ws.release_BA<LEVEL,longhint_t>(side, BAout);

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL, longhint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);

    /* SECOND: fill in buckets at this level, for this region. */
    ws.reset_all_pointers<LEVEL,shorthint_t>(side);
    slice_index_t slices_pushed = 0;
    for (auto const & it : precomp_plattice(side, LEVEL)) {
      fill_in_buckets_parameters *param =
        new fill_in_buckets_parameters(ws, side, si,
            (fb_slice_interface *)NULL, it, NULL, first_region0_index);
      double c = it->get_weight();
      pool.add_task(fill_in_buckets_one_slice_internal<LEVEL>, param, 0, 0, c);
      slices_pushed++;
    }

    std::vector<buckets_are_full> exceptions_to_throw;
    for (slice_index_t slices_completed = 0;
        slices_completed < slices_pushed;
        slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
      /* want to check possible exceptions, too */
      buckets_are_full * e = dynamic_cast<buckets_are_full*>(pool.get_exception());
      if (e) {
          exceptions_to_throw.push_back(*e);
          delete e;
      }
    }
    if (!exceptions_to_throw.empty())
        throw *std::max_element(exceptions_to_throw.begin(), exceptions_to_throw.end());

    pool.accumulate_and_clear_active_time(*timer.current);
    SIBLING_TIMER(timer, "worker thread wait time");
    TIMER_CATEGORY(timer, thread_wait());
    pool.accumulate_and_reset_wait_time(*timer.current);

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL,shorthint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);
  }

#if 0
  /* We've filled buckets at level LEVEL: both downsorted updates from
   * upper level as well as updates that must be processed here. Time to
   * report */
  ws.diagnosis(LEVEL, {si.sides[0].fbs, si.sides[1].fbs});
#endif

  /* RECURSE */
  if (LEVEL > 1) {
    for (unsigned int i = 0; i < si.nb_buckets[LEVEL]; ++i) {
      size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
      uint32_t N = first_region0_index + i*(BRS[LEVEL]/BRS[1]);
      downsort_tree<LEVEL-1>(timer, i, N, ws, pool, si, precomp_plattice);
    }
  } else {
    /* PROCESS THE REGIONS AT LEVEL 0 */
    for (int i = 0; i < ws.thrs[0].plas->nb_threads; ++i) {
      ws.thrs[i].first_region0_index = first_region0_index;
    }
    ws.thread_do_using_pool(pool, &process_bucket_region);
  }
  pool.accumulate_and_clear_active_time(*timer.current);
  SIBLING_TIMER(timer, "worker thread wait time");
  TIMER_CATEGORY(timer, thread_wait());
  pool.accumulate_and_reset_wait_time(*timer.current);
}

/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(timetree_t&,
  uint32_t, uint32_t,
  thread_workspaces &,
  thread_pool &,
  sieve_info &,
  precomp_plattice_t const &)
{
    ASSERT_ALWAYS(0);
}

// other fake instances to please level-2 instanciation.
template <>
void downsort<3>(bucket_array_t<2, longhint_t>&,
        bucket_array_t<3, longhint_t> const&, unsigned int, where_am_I &)
{
    ASSERT_ALWAYS(0);
}

template <>
reservation_array<bucket_array_t<3, longhint_t> > const&
reservation_group::cget<3, longhint_t>() const
{
    ASSERT_ALWAYS(0);
}

// Now the two exported instances

template 
void downsort_tree<1>(timetree_t&, uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, thread_pool &pool, sieve_info & si,
  precomp_plattice_t const & precomp_plattice);

template
void downsort_tree<2>(timetree_t&, uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, thread_pool &pool, sieve_info & si,
  precomp_plattice_t const & precomp_plattice);
