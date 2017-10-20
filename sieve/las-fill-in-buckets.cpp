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
                const fb_slice_interface * const slice,
                plattices_dense_vector_t * precomp_slice,
                where_am_I & w)
{
  ASSERT_ALWAYS(si.conf.sublat.m != 0);
  bool first_sublat = si.conf.sublat.i0 == 0 && si.conf.sublat.j0 == 1;
  ASSERT(slice != NULL);
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  slice_index_t slice_index = slice->get_index();
  if (!first_sublat) {
    ASSERT(slice_index == precomp_slice->get_index());
  }

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  typename FB_ENTRY_TYPE::transformed_entry_t transformed;

  // FIXME: A LOT OF DUPLICATED CODE, HERE!!!
  if (first_sublat) {
    slice_offset_t i_entry = 0;
    const fb_slice<FB_ENTRY_TYPE> * const sl = (const fb_slice<FB_ENTRY_TYPE> * const) slice;
    for (const FB_ENTRY_TYPE *it = sl->begin(); it != sl->end(); it++, i_entry++) {
      if (!si.qbasis.is_coprime_to(it->p))
        continue;
      it->transform_roots(transformed, si.qbasis);
      for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
        const fbroot_t r = transformed.get_r(i_root);
        const bool proj = transformed.get_proj(i_root);
        /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
           locations with p | gcd(i,j). */
        if (LIKELY(!proj || r == 0)) {
          plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, si.conf.logI_adjusted);
          // In sublat mode, save it for later use
          plattice_info_dense_t plid(pli, i_entry);
          precomp_slice->push_back(plid);

          plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI_adjusted, si.conf.sublat);

          if (plattice_enumerate_finished<LEVEL>(ple.get_x()))
            continue;
          if (LIKELY(pli.a0 != 0)) {
            const slice_offset_t hint = ple.get_hint();
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice->get_prime(hint); 
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
      plattice_info_t pli = (*precomp_slice)[i].unpack(si.conf.logI_adjusted);
      slice_offset_t i_entry = (*precomp_slice)[i].hint;

      plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI_adjusted, si.conf.sublat);

      if (plattice_enumerate_finished<LEVEL>(ple.get_x()))
        continue;
      if (LIKELY(pli.a0 != 0)) {
        const slice_offset_t hint = ple.get_hint();
        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        const fbprime_t p = slice->get_prime(hint); 
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
  // printf("%.3f\n", BA.max_full());
  orig_BA.move(BA);
}

template <int LEVEL, class FB_ENTRY_TYPE>
void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, shorthint_t> &orig_BA,
                sieve_info const & si MAYBE_UNUSED,
                const fb_slice_interface * const slice,
                plattices_dense_vector_t * precomp_slice,
                where_am_I & w)
{
  if (si.conf.sublat.m) {
    // TODO: two function calls might yield a penalty, due to 
    // instruction cache pollution (due to templating, we are talking
    // about many, many functions, here).
    // Could worth investigating...
    fill_in_buckets_toplevel_sublat<LEVEL,FB_ENTRY_TYPE>(orig_BA, si, slice, 
        precomp_slice, w);
    return;
  }

  ASSERT(slice != NULL);
  bool first_reg = true;
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  slice_index_t slice_index = slice->get_index();

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  typename FB_ENTRY_TYPE::transformed_entry_t transformed;

  slice_offset_t i_entry = 0;
  const fb_slice<FB_ENTRY_TYPE> * const sl = (const fb_slice<FB_ENTRY_TYPE> * const) slice;
  for (const FB_ENTRY_TYPE *it = sl->begin(); it != sl->end(); it++, i_entry++) {
    if (!si.qbasis.is_coprime_to(it->p))
      continue;
    it->transform_roots(transformed, si.qbasis);
    for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
      const fbroot_t r = transformed.get_r(i_root);
      const bool proj = transformed.get_proj(i_root);
      /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
         locations with p | gcd(i,j). */
      if (LIKELY(!proj || r == 0)) {
        plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, si.conf.logI_adjusted);
  
        plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, si.conf.logI_adjusted);

        // Skip (i,j)=(0,0)
        ple.next();
        if (LIKELY(pli.a0 != 0)) {
          const slice_offset_t hint = ple.get_hint();
          WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
          const fbprime_t p = slice->get_prime(hint); 
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
  // printf("%.3f\n", BA.max_full());
  orig_BA.move(BA);
}



/* {{{ */
template <int LEVEL>
  void
fill_in_buckets_lowlevel(
    bucket_array_t<LEVEL, shorthint_t> &orig_BA,
    sieve_info const & si MAYBE_UNUSED,
    plattices_vector_t *plattices_vector,
    bool first_reg,
    where_am_I & w)
{
    /* The timer stuff is dealt with by the caller */
  const slice_index_t slice_index = plattices_vector->get_index();
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);

  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  for (plattices_vector_t::iterator pl_it = plattices_vector->begin();
      pl_it != plattices_vector->end(); pl_it++) {

    // Work with a copy, otherwise we don't get all optimizations.
    // Maybe with a wise use of the 'restrict' keyword, we might get
    // what we want, but this is C++11, anyway.
    plattice_enumerate_t pl(*pl_it);

    const slice_offset_t hint = pl.get_hint();
    WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
    const fb_slice_interface * slice =
      si.sides[w.side].fb->get_slice(slice_index);
    const fbprime_t p = slice->get_prime(hint); 
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
    pl_it->set_x(pl.get_x());
    pl_it->advance_to_next_area(LEVEL);
  } 
  // printf("%.3f\n", BA.max_full());
  orig_BA.move(BA);
}

class fill_in_buckets_parameters: public task_parameters {
public:
  thread_workspaces &ws;
  const int side;
  sieve_info const & si;
  const fb_slice_interface * const slice;
  plattices_vector_t * const plattices_vector; // content changed during fill-in
  plattices_dense_vector_t * const plattices_dense_vector; // for sublat
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
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice, "");
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

    /* Get an unused bucket array that we can write to */
    bucket_array_t<LEVEL, shorthint_t> *BA_ptr =
        param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
    if (BA_ptr) {
        bucket_array_t<LEVEL, shorthint_t> &BA(*BA_ptr);
        /* Fill the buckets */
        fill_in_buckets_lowlevel<LEVEL>(BA, param->si, param->plattices_vector,
                (param->first_region0_index == 0), w);
        /* Release bucket array again */
        param->ws.release_BA(param->side, BA);
    }
    delete param;
    return new task_result;
}


// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template<int LEVEL>
task_result *
fill_in_buckets_one_slice(const worker_thread * worker MAYBE_UNUSED, const task_parameters * _param)
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
    CHILD_TIMER(worker->timer, TEMPLATE_INST_NAME(fill_in_buckets_one_slice, LEVEL));
    // TIMER_CATEGORY(worker->timer, sieving(param->side));

    where_am_I w;
    WHERE_AM_I_UPDATE(w, psi, & param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    /* Get an unused bucket array that we can write to */
    bucket_array_t<LEVEL, shorthint_t> * BA_ptr = param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
    if (BA_ptr) {
        bucket_array_t<LEVEL, shorthint_t> &BA(*BA_ptr);
        /* Fill the buckets */
        if (param->slice->is_general())
          fill_in_buckets_toplevel<LEVEL,fb_general_entry>(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 0)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<0> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 1)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<1> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 2)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<2> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 3)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<3> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 4)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<4> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 5)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<5> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 6)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<6> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 7)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<7> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 8)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<8> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 9)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<9> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else if (param->slice->get_nr_roots() == 10)
          fill_in_buckets_toplevel<LEVEL,fb_entry_x_roots<10> >(BA, param->si, param->slice, param->plattices_dense_vector, w);
        else
          ASSERT_ALWAYS(0);
        /* Release bucket array again */
        param->ws.release_BA(param->side, BA);
    }
    delete param;
    return new task_result;
}

template <int LEVEL>
static void
fill_in_buckets_one_side(timetree_t& timer, thread_pool &pool, thread_workspaces &ws, const fb_part *fb, sieve_info & si, const int side)
{
    CHILD_TIMER_PARAMETRIC(timer, "fill_in_buckets_one_side on side ", side, "");
    TIMER_CATEGORY(timer, sieving(side));

    slice_index_t slices_pushed = 0;

    if (!si.conf.sublat.m || (si.conf.sublat.i0 == 0 && si.conf.sublat.j0 == 1)) {
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

    for (slice_index_t slices_completed = 0; slices_completed < slices_pushed; slices_completed++) {
          task_result *result = pool.get_result();
          delete result;
          /* want to check possible exceptions, too */
    }
    pool.accumulate_and_clear_active_time(*timer.current);
    SIBLING_TIMER(timer, "worker thread wait time");
    TIMER_CATEGORY(timer, thread_wait());
    pool.accumulate_and_reset_wait_time(*timer.current);
}

void fill_in_buckets(timetree_t& timer, thread_pool &pool, thread_workspaces &ws, sieve_info & si, int side)
{
    // per se, we're not doing anything here.
    // CHILD_TIMER(timer, __func__);
    plattice_enumerate_t::set_masks(si.conf.logI_adjusted);
    switch (si.toplevel) {
        case 1:
            plattice_enumerate_area<1>::value = plattice_x_t(si.J) << si.conf.logI_adjusted;
            fill_in_buckets_one_side<1>(timer, pool, ws,
                    si.sides[side].fb->get_part(si.toplevel), si, side);
            break;
        case 2:
            plattice_enumerate_area<2>::value = plattice_x_t(si.J) << si.conf.logI_adjusted;
            fill_in_buckets_one_side<2>(timer, pool, ws,
                    si.sides[side].fb->get_part(si.toplevel), si, side);
            break;
        case 3:
            plattice_enumerate_area<3>::value = plattice_x_t(si.J) << si.conf.logI_adjusted;
            fill_in_buckets_one_side<3>(timer, pool, ws,
                    si.sides[side].fb->get_part(si.toplevel), si, side);
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
    precomp_plattice_t precomp_plattice)
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
    bucket_array_t<LEVEL, longhint_t> * BAout_ptr =
      ws.reserve_BA<LEVEL, longhint_t>(side);
    if (BAout_ptr) {
        bucket_array_t<LEVEL, longhint_t> & BAout(*BAout_ptr);
        BAout.reset_pointers();
        // This is a fake slice_index. For a longhint_t bucket, each update
        // contains its own slice_index, directly used by apply_one_bucket
        // and purge.
        BAout.add_slice_index(0);
        // The data that comes from fill-in bucket at level above:
        {
          const bucket_array_t<LEVEL+1,shorthint_t> * BAin
            = ws.cbegin_BA<LEVEL+1,shorthint_t>(side);
          while (BAin != ws.cend_BA<LEVEL+1,shorthint_t>(side)) {
            downsort<LEVEL+1>(BAout, *BAin, bucket_index, w);
            BAin++;
          }
        }

        const int toplevel = si.toplevel;
        if (LEVEL < toplevel - 1) {
          // What comes from already downsorted data above:
          const bucket_array_t<LEVEL+1,longhint_t> * BAin
            = ws.cbegin_BA<LEVEL+1,longhint_t>(side);
          while (BAin != ws.cend_BA<LEVEL+1,longhint_t>(side)) { 
            downsort<LEVEL+1>(BAout, *BAin, bucket_index, w);
            BAin++;
          }
        }
        ws.release_BA<LEVEL,longhint_t>(side, BAout);
    } else {
        fprintf(stderr, "encountered full buckets while downsorting\n");
    }

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL, longhint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);

    /* SECOND: fill in buckets at this level, for this region. */
    ws.reset_all_pointers<LEVEL,shorthint_t>(side);
    slice_index_t slices_pushed = 0;
    for (typename std::vector<plattices_vector_t *>::iterator pl_it =
            precomp_plattice[side][LEVEL].begin();
        pl_it != precomp_plattice[side][LEVEL].end();
        pl_it++) {
      fill_in_buckets_parameters *param =
        new fill_in_buckets_parameters(ws, side, si,
            (fb_slice_interface *)NULL, *pl_it, NULL, first_region0_index);
      // TODO: shall we give the weight to help scheduling, here?
      pool.add_task(fill_in_buckets_one_slice_internal<LEVEL>, param, 0);
      slices_pushed++;
    }

    for (slice_index_t slices_completed = 0;
        slices_completed < slices_pushed;
        slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
    }

    pool.accumulate_and_clear_active_time(*timer.current);
    SIBLING_TIMER(timer, "worker thread wait time");
    TIMER_CATEGORY(timer, thread_wait());
    pool.accumulate_and_reset_wait_time(*timer.current);

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL,shorthint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);
  }

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
        uint32_t bucket_index MAYBE_UNUSED,
  uint32_t first_region0_index MAYBE_UNUSED,
  thread_workspaces &ws MAYBE_UNUSED,
  thread_pool &pool MAYBE_UNUSED,
  sieve_info & si MAYBE_UNUSED,
  precomp_plattice_t precomp_plattice MAYBE_UNUSED)
{
    ASSERT_ALWAYS(0);
}

// other fake instances to please level-2 instanciation.
template <>
bucket_array_t<3, longhint_t>::bucket_array_t()
{
    ASSERT_ALWAYS(0);
}

template <>
bucket_array_t<3, longhint_t>::~bucket_array_t()
{
    ASSERT_ALWAYS(0);
}

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
  precomp_plattice_t precomp_plattice);

template
void downsort_tree<2>(timetree_t&, uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, thread_pool &pool, sieve_info & si,
  precomp_plattice_t precomp_plattice);
