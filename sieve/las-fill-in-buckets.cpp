#include "cado.h"

#include "fb.h"
#include "utils.h"           /* lots of stuff */
#include "bucket.h"
#include "modredc_ul.h"
#include "modredc_2ul2.h"
#include "threadpool.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "las-debug.h"
#include "las-arith.h"
#include "las-qlattice.h"
#include "las-fill-in-buckets.h"
#include "las-norms.h"
#include "las-smallsieve.h"
#include "las-plattice.h"

#ifdef USE_CACHEBUFFER
#include "cachebuf.h"
#endif

// FIXME: this function of las.cpp should be somewhere
void * process_bucket_region(thread_data *th);


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
                  const size_t n, const int side, sieve_info_srcptr si,
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
    mpz_poly_srcptr f = si->sides[side]->fij;
    for (size_t j = 0; j < i; j++)
      r[j] = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
  }

#if 0
    /* Very expensive but thorough test */
    for (size_t j = 0; j < i; j++) {
      mpz_poly_srcptr f = si->sides[side]->fij;
      unsigned long tr = compute_1_root_mpz(p[j], f->coeff[0], f->coeff[1]);
      ASSERT_ALWAYS(tr == r[j]);
    }
#endif

  return i;
}
#endif


/* {{{ */
template <int LEVEL>
void
fill_in_buckets(bucket_array_t<LEVEL, shorthint_t> &orig_BA,
                sieve_info_srcptr const si MAYBE_UNUSED,
                plattices_vector_t *plattices_vector,
                bool first_reg,
                where_am_I_ptr w)
{
  const slice_index_t slice_index = plattices_vector->get_index();
  bucket_array_t<LEVEL, shorthint_t> BA;  /* local copy. Gain a register + use stack */
  BA.move(orig_BA);
  
  /* Write new set of pointers for the new slice */
  BA.add_slice_index(slice_index);

  for (plattices_vector_t::iterator pl_it = plattices_vector->begin();
       pl_it != plattices_vector->end(); pl_it++) {

    const slice_offset_t hint = pl_it->get_hint();
    WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
    const fb_slice_interface * slice =
        si->sides[w->side]->fb->get_slice(slice_index);
    const fbprime_t p = slice->get_prime(hint); 
    WHERE_AM_I_UPDATE(w, p, p);
#else
    const fbprime_t p = 0;
#endif

    // Handle the rare special cases
    const uint32_t I = si->I;
    if (UNLIKELY(pl_it->get_inc_c() == 1 && pl_it->get_bound1() == I - 1)) {
        // Projective root: only update is at (1,0).
        if (first_reg) {
            uint64_t x = 1 + (I >> 1);
            BA.push_update(x, p, hint, slice_index, w);
        }
        continue;
    }
    if (UNLIKELY(pl_it->get_inc_c() == I && pl_it->get_bound1() == I)) {
        // Root=0: only update is at (0,1).
        if (first_reg) {
            uint64_t x = I + (I >> 1);
            BA.push_update(x, p, hint, slice_index, w);
        }
        continue;
    }

    /* Now, do the real work: the filling of the buckets */
    while (!plattice_enumerate_finished<LEVEL>(pl_it->get_x())) {
      if (LIKELY(pl_it->probably_coprime()))
        BA.push_update(pl_it->get_x(), p, hint, slice_index, w);
      pl_it->next();
    }

    pl_it->advance_to_next_area(LEVEL);
  } 
  orig_BA.move(BA);
}

class fill_in_buckets_parameters: public task_parameters {
public:
  thread_workspaces &ws;
  const int side;
  sieve_info_srcptr const si;
  const fb_slice_interface * const slice;
  plattices_vector_t * const plattices_vector; // content changed during fill-in
  const uint32_t first_region0_index;
  fill_in_buckets_parameters(thread_workspaces &_ws, const int _side,
          sieve_info_srcptr const _si, const fb_slice_interface *_slice,
          plattices_vector_t *_platt, const uint32_t _reg0)
  : ws(_ws), side(_side), si(_si), slice(_slice),
    plattices_vector(_platt), first_region0_index(_reg0)
  {}
};

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template<int LEVEL>
task_result *
fill_in_buckets_one_slice_internal(const task_parameters *const _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    where_am_I w;
    WHERE_AM_I_UPDATE(w, si, param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->plattices_vector->get_index());

    /* Get an unused bucket array that we can write to */
    bucket_array_t<LEVEL, shorthint_t> &BA =
        param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
    /* Fill the buckets */
    fill_in_buckets<LEVEL>(BA, param->si, param->plattices_vector,
            (param->first_region0_index == 0), w);
    /* Release bucket array again */
    param->ws.release_BA(param->side, BA);
    delete param;
    return new task_result;
}


template<int LEVEL>
task_result *
fill_in_buckets_one_slice(const task_parameters *const _param)
{
    const fill_in_buckets_parameters *param = static_cast<const fill_in_buckets_parameters *>(_param);
    where_am_I w;
    WHERE_AM_I_UPDATE(w, si, param->si);
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    /* Do the root transform and lattice basis reduction for this factor base slice */
    plattices_vector_t *plattices_vector =
        param->slice->make_lattice_bases(param->si->qbasis, param->si->conf->logI);
    /* Get an unused bucket array that we can write to */
    bucket_array_t<LEVEL, shorthint_t> &BA = param->ws.reserve_BA<LEVEL, shorthint_t>(param->side);
    /* Fill the buckets */
    fill_in_buckets<LEVEL>(BA, param->si, plattices_vector, true, w);
    /* Release bucket array again */
    param->ws.release_BA(param->side, BA);
    delete plattices_vector;
    delete param;
    return new task_result;
}

template <int LEVEL>
static void
fill_in_buckets_one_side(thread_pool &pool, thread_workspaces &ws, const fb_part *fb, sieve_info_srcptr const si, const int side)
{
    /* Process all slices in this factor base part */
    const fb_slice_interface *slice;
    slice_index_t slices_pushed = 0;
    for (slice_index_t slice_index = fb->get_first_slice_index();
         (slice = fb->get_slice(slice_index)) != NULL;
         slice_index++) {
        fill_in_buckets_parameters *param = new fill_in_buckets_parameters(ws, side, si, slice, NULL, 0);
        pool.add_task(fill_in_buckets_one_slice<LEVEL>, param, 0);
        slices_pushed++;
    }
    for (slice_index_t slices_completed = 0; slices_completed < slices_pushed; slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
    }
}

void fill_in_buckets_both(thread_pool &pool, thread_workspaces &ws, sieve_info_srcptr si)
{
  plattice_enumerate_t::set_masks(si->conf->logI);
  for (int side = 0; side < si->cpoly->nb_polys; ++side) {
    switch (si->toplevel) {
      case 1:
        plattice_enumerate_area<1>::value = plattice_x_t(si->J) << si->conf->logI;
        fill_in_buckets_one_side<1>(pool, ws,
            si->sides[side]->fb->get_part(si->toplevel), si, side);
        break;
      case 2:
        plattice_enumerate_area<2>::value = plattice_x_t(si->J) << si->conf->logI;
        fill_in_buckets_one_side<2>(pool, ws,
            si->sides[side]->fb->get_part(si->toplevel), si, side);
        break;
      case 3:
        plattice_enumerate_area<3>::value = plattice_x_t(si->J) << si->conf->logI;
        fill_in_buckets_one_side<3>(pool, ws,
            si->sides[side]->fb->get_part(si->toplevel), si, side);
        break;
      default:
        ASSERT_ALWAYS(0);
    }
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
downsort_tree(uint32_t bucket_index,
    uint32_t first_region0_index,
    thread_workspaces &ws,
    thread_pool &pool,
    sieve_info_ptr si,
    precomp_plattice_t precomp_plattice)
{
  ASSERT_ALWAYS(LEVEL > 0);

  where_am_I w;
  WHERE_AM_I_UPDATE(w, si, si);
  WHERE_AM_I_UPDATE(w, N, first_region0_index);

  double max_full MAYBE_UNUSED = 0.0;

  for (int side = 0; side < 2; ++side) { // FIXME: MNFS?
    WHERE_AM_I_UPDATE(w, side, side);
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

    const int toplevel = si->toplevel;
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
            (fb_slice_interface *)NULL, *pl_it, first_region0_index);
      pool.add_task(fill_in_buckets_one_slice_internal<LEVEL>, param, 0);
      slices_pushed++;
    }
    for (slice_index_t slices_completed = 0;
        slices_completed < slices_pushed;
        slices_completed++) {
      task_result *result = pool.get_result();
      delete result;
    }

    max_full = std::max(max_full, ws.buckets_max_full<LEVEL,shorthint_t>());
    ASSERT_ALWAYS(max_full <= 1.0);
  }

  /* RECURSE */
  if (LEVEL > 1) {
    for (unsigned int i = 0; i < si->nb_buckets[LEVEL]; ++i) {
      uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
      uint32_t N = first_region0_index + i*(BRS[LEVEL]/BRS[1]);
      downsort_tree<LEVEL-1>(i, N, ws, pool, si, precomp_plattice);
    }
  } else {
    /* PROCESS THE REGIONS AT LEVEL 0 */
    for (int i = 0; i < ws.thrs[0].las->nb_threads; ++i) {
      ws.thrs[i].first_region0_index = first_region0_index;
    }
    ws.thread_do(&process_bucket_region);
  }
}

/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(uint32_t bucket_index MAYBE_UNUSED,
  uint32_t first_region0_index MAYBE_UNUSED,
  thread_workspaces &ws MAYBE_UNUSED,
  thread_pool &pool MAYBE_UNUSED,
  sieve_info_ptr si MAYBE_UNUSED,
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
        bucket_array_t<3, longhint_t> const&, unsigned int, where_am_I_ptr)
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
void downsort_tree<1>(uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, thread_pool &pool, sieve_info_ptr si,
  precomp_plattice_t precomp_plattice);

template
void downsort_tree<2>(uint32_t bucket_index, uint32_t first_region0_index,
  thread_workspaces &ws, thread_pool &pool, sieve_info_ptr si,
  precomp_plattice_t precomp_plattice);
