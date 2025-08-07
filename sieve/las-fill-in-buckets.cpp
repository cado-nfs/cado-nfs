#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstddef> // for size_t, NULL
#include <cstdint> // for uint32_t, uint64_t

#include <algorithm> // for max_element, min
#include <array>     // for array
#include <limits>    // for numeric_limits
#include <memory>    // for shared_ptr, allocator, __sh...
#include <utility>   // for move
#include <vector>    // for vector, vector<>::iterator

#include "bucket-push-update.hpp"        // for bucket_array_t::push_update
#include "bucket.hpp"                    // for longhint_t (ptr only), buck...
#include "fb-types.hpp"
#include "fb.hpp"                        // for fb_slice, fb_factorbase
#include "las-auxiliary-data.hpp"        // for nfs_aux, nfs_aux::thread_data
#include "las-bkmult.hpp"                // for buckets_are_full
#include "las-config.h"                  // for FB_MAX_PARTS, BUCKET_REGIONS
#include "las-fill-in-buckets.hpp"       // for precomp_plattice_t, downsor...
#include "las-plattice.hpp"              // for plattice_info, plattice_e...
#include "las-process-bucket-region.hpp" // for process_many_bucket_regions
#include "las-qlattice.hpp"              // for qlattice_basis
#include "las-report-stats.hpp"          // for TIMER_CATEGORY
#include "las-siever-config.hpp"         // for siever_config, siever_confi...
#include "las-smallsieve.hpp"            // for small_sieve_activate_many_s...
#include "las-threads-work-data.hpp"     // for nfs_work, nfs_work::side_data
#include "las-where-am-i-proxy.hpp"      // for where_am_I
#include "las-where-am-i.hpp"            // for where_am_I, WHERE_AM_I_UPDATE
#include "macros.h"                      // for CADO_CONCATENATE, ASSERT_AL...
#include "multityped_array.hpp"          // for multityped_array (ptr only)
#include "tdict.hpp"                     // for slot, timetree_t, CHILD_TIMER
#include "threadpool.hpp"                // for thread_pool, worker_thread
#include "utils_cxx.hpp"

/* is this in the std library or not ? */
template <typename T> static inline T const & const_ref(T & x)
{
    return x;
}

/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/**************************************************************************
 * Global DEFINEs for fill_in_buckets, fill_in_k_buckets, fill_in_m_buckets
 **************************************************************************/

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
                  const size_t n, const int side,
                  sieve_info const & si,
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
template <class FB_ENTRY_TYPE>
static inline bool discard_power_for_bucket_sieving(FB_ENTRY_TYPE const &)
{
    /* the entry is not a general entry, therefore k is a const thing
     * equal to 1.
     */
    return false;
}
#ifndef BUCKET_SIEVE_POWERS
template <>
inline bool
discard_power_for_bucket_sieving<fb_entry_general>(fb_entry_general const & e)
{
    return e.k > 1;
}
#endif

/***********************************************************************/
/* multithreaded processing of make_lattice_bases (a.k.a
 * precomp_plattices)
 *
 * This creates one control object per slice, with storage ownership of
 * the control object transfered to the called function. Because we
 * depend on the slice, the type of the object is parameterized by the
 * slice type.
 *
 * We may elect to make the "model" a shared_ptr.
 */

template <int LEVEL>
struct make_lattice_bases_parameters_base : public task_parameters {
    int side;
    nfs_work & ws;
    typename precomp_plattice_t<LEVEL>::vec_type & V;
    make_lattice_bases_parameters_base(
        int side, nfs_work & ws,
        typename precomp_plattice_t<LEVEL>::vec_type & V)
        : side(side)
        , ws(ws)
        , V(V)
    {
    }
};
template <int LEVEL, class FB_ENTRY_TYPE>
struct make_lattice_bases_parameters
    : public make_lattice_bases_parameters_base<LEVEL> {
    typedef make_lattice_bases_parameters_base<LEVEL> super;
    fb_slice<FB_ENTRY_TYPE> const & slice;
    make_lattice_bases_parameters(super const & model,
                                  fb_slice<FB_ENTRY_TYPE> const & slice)
        : super(model)
        , slice(slice)
    {
    }
};

template <int LEVEL, class FB_ENTRY_TYPE>
static task_result * make_lattice_bases(worker_thread * worker MAYBE_UNUSED,
                                        task_parameters * _param, int)
{
    make_lattice_bases_parameters<LEVEL, FB_ENTRY_TYPE> const * param =
        static_cast<
            make_lattice_bases_parameters<LEVEL, FB_ENTRY_TYPE> const *>(
            _param);

    nfs_work const & ws(param->ws);
    qlattice_basis const & Q(ws.Q);
    int const logI = ws.conf.logI;
    sublat_t const & sublat(ws.Q.sublat);
    typename precomp_plattice_t<LEVEL>::vec_type & V(param->V);
    fb_slice<FB_ENTRY_TYPE> const & slice(param->slice);

    slice_index_t const index = slice.get_index();

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;
    /* Create a transformed vector and store the index of the fb_slice we
     * currently transform */

    /* we don't really need the fence at this point, except that we do one
     * shot of "next()", and that happens to need the fence. Only logI is
     * needed, though.
     */
    plattice_enumerator::fence const F(ws.conf.logI, 0);

    plattices_vector_t result(index, slice.get_weight());
    slice_offset_t i_entry = 0;
    for (auto const & e: slice) {
        increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;
        e.transform_roots(transformed, Q);
        for (unsigned char i_root = 0; i_root != transformed.nr_roots;
             i_root++) {
            fbroot_t const r = transformed.get_r(i_root);
            bool const proj = transformed.get_proj(i_root);
            plattice_info const pli =
                plattice_info(transformed.get_q(), r, proj, logI);
            plattice_enumerator ple(pli, i_entry, logI, sublat);
            // Skip (0,0) unless we have sublattices.
            if (!sublat.m)
                ple.next(F);
            if (LIKELY(!pli.is_discarded()))
                result.push_back(ple);
        }
    }
    /* This is moved, not copied. Note that V is a reference. */
    V[result.get_index()] = std::move(result);
    delete param;
    return new task_result;
}

template <int LEVEL> struct push_make_bases_to_task_list {
    thread_pool & pool;
    make_lattice_bases_parameters_base<LEVEL> model;
    size_t pushed = 0;
    push_make_bases_to_task_list(
        thread_pool & pool, make_lattice_bases_parameters_base<LEVEL> const & m)
        : pool(pool)
        , model(m)
    {
    }
    template <typename T> void operator()(T const & s)
    {
        typedef typename T::entry_t E;
        auto param = new make_lattice_bases_parameters<LEVEL, E>(model, s);
        task_function_t f = make_lattice_bases<LEVEL, E>;
        pool.add_task(f, param, 0);
        pushed++;
    }
};

struct helper_functor_prepare_precomp_plattice {
    nfs_work & ws;
    thread_pool & pool;
    int side;
    helper_functor_prepare_precomp_plattice(nfs_work & ws, thread_pool & pool,
                                            int side)
        : ws(ws)
        , pool(pool)
        , side(side)
    {
    }

    template <typename T> /* T is precomp_plattice_t<n> for some level n */
    void operator()(T & precomp_plattice)
    {
        if (T::level >= ws.toplevel)
            return;

        nfs_work::side_data const & wss(ws.sides[side]);
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(T::level);
        typename precomp_plattice_t<T::level>::vec_type & V =
            precomp_plattice.v[side];
        /* We pre-assign the result, so that all threads can write to it
         * comfortably.
         *
         * It would be nice to have a way to notify that all threads here are
         * done with their job.
         */
        V.assign(P.nslices(), plattices_vector_t());
        make_lattice_bases_parameters_base<T::level> const model {side, ws, V};
        push_make_bases_to_task_list<T::level> F {pool, model};
        P.foreach_slice(F);
    }
};

void fill_in_buckets_prepare_plattices(
    nfs_work & ws, thread_pool & pool, int side,
    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> & precomp_plattice)
{
    helper_functor_prepare_precomp_plattice H(ws, pool, side);

    /* this will *not* do anything for level==ws.toplevel, by design */
    multityped_array_foreach(H, precomp_plattice);
}

/***********************************************************************/
/* multithreaded processing of fill_in_buckets_toplevel (both with and
 * without sublattices) is more complicated. First because the important
 * functions are not the ones whose prototype is the one we expect most
 * from a multithreade task, second because we strive to manage
 * exceptions properly. So we go through several quirky paths below.
 */

// At top level, the fill-in of the buckets must interleave
// the root transforms and the FK walks, otherwise we spend a lot of time
// doing nothing while waiting for memory.
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

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void fill_in_buckets_toplevel_sublat(
    bucket_array_t<LEVEL, TARGET_HINT> & orig_BA, nfs_work & ws,
    fb_slice<FB_ENTRY_TYPE> const & slice,
    plattices_dense_vector_t * p_precomp_slice, where_am_I & w)
{
    int const logI = ws.conf.logI;
    qlattice_basis const & Q(ws.Q);

    plattices_dense_vector_t & precomp_slice(*p_precomp_slice);

    ASSERT_ALWAYS(Q.sublat.m);
    bool const first_sublat = Q.sublat.i0 == 0 && Q.sublat.j0 == 1;
    bucket_array_t<LEVEL, TARGET_HINT>
        BA; /* local copy. Gain a register + use stack */
    BA.move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;

    /* top level: the fence we care about is the one defined by J */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J);

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    // FIXME: A LOT OF DUPLICATED CODE, HERE!!!
    if (first_sublat) {
        slice_offset_t i_entry = 0;
        for (auto const & e: slice) {
            increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
            if (!Q.is_coprime_to(e.p))
                continue;
#ifdef BUCKET_SIEVE_POWERS
            /* the combination of bucket-sieving powers + sublattices means
             * that powers of the primes that divide the sublattice determinant
             * may be bucket-sieved. And of course, that leads to problems.
             *
             * technically, Q.sublat.m could be composite, in which case we
             * would have a gcd to compute, here. The only really useful case
             * at the moment is m=3 though.
             */
            if (Q.sublat.m == 3) {
                if (e.p == 3)
                    continue;
            } else if (gcd_ul(e.p, Q.sublat.m) > 1) {
                continue;
            }
#endif
            if (discard_power_for_bucket_sieving(e))
                continue;
            e.transform_roots(transformed, Q);
            for (unsigned char i_root = 0; i_root != transformed.nr_roots;
                 i_root++) {
                fbroot_t const r = transformed.get_r(i_root);
                bool const proj = transformed.get_proj(i_root);
                plattice_info const pli(transformed.get_q(), r, proj, logI);
                // In sublat mode, save it for later use
                precomp_slice.push_back(plattice_info_dense_t(pli, i_entry));

                plattice_enumerator ple(pli, i_entry, logI, Q.sublat);

                if (ple.done(F))
                    continue;
                if (pli.is_discarded())
                    continue;
                slice_offset_t const hint = ple.get_hint();
                ASSERT(hint == i_entry);
                WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
                const fbprime_t p = slice.get_prime(hint);
                WHERE_AM_I_UPDATE(w, p, p);
#else
                const fbprime_t p = 0;
#endif

                typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                    0, p, hint, slice_index);

                // Handle the rare special cases
                // XXX Here, we're not bucket-sieving projective primes at
                // all, and neither do we bucket-sieve primes with root equal
                // to zero.
                if (UNLIKELY(pli.is_vertical_line(logI) ||
                             pli.is_projective_like(logI)))
                    continue;

                /* Now, do the real work: the filling of the buckets */
                // Without sublattices, we test (very basic) coprimality,
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.next(F);
                }
            }
        }
    } else { // Use precomputed FK-basis
        for (unsigned int i = 0; i < precomp_slice.size(); ++i) {
            plattice_info const pli(precomp_slice[i].unpack(logI));
            slice_offset_t const i_entry = precomp_slice[i].get_hint();

            plattice_enumerator ple(pli, i_entry, logI, Q.sublat);

            if (ple.done(F))
                continue;
            if (pli.is_discarded())
                continue;
            slice_offset_t const hint = ple.get_hint();
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice.get_prime(hint);
            WHERE_AM_I_UPDATE(w, p, p);
#else
            const fbprime_t p = 0;
#endif

            typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                0, p, hint, slice_index);

            // Handle (well, do not handle, in fact) the rare special cases
            if (UNLIKELY(pli.is_vertical_line(logI) ||
                         pli.is_projective_like(logI)))
                continue;

            /* Now, do the real work: the filling of the buckets */
            // Without sublattices, we test (very basic) coprimality,
            // otherwise not atm. FIXME!
            while (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                BA.push_update(ple.get_x() >> logB, u, w);
                ple.next(F);
            }
        }
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA.move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, fb_slice<FB_ENTRY_TYPE> const & slice,
                         plattices_dense_vector_t * /* unused */,
                         where_am_I & w)
{
    int const logI = ws.conf.logI;
    qlattice_basis const & Q(ws.Q);

    ASSERT_ALWAYS(!Q.sublat.m);

    bucket_array_t<LEVEL, TARGET_HINT>
        BA; /* local copy. Gain a register + use stack */
    BA.move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    typename FB_ENTRY_TYPE::transformed_entry_t transformed;

    /* top level: the fence we care about is the one defined by J */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J);

    slice_offset_t i_entry = 0;

    /* yes, we want the level-1 regions here */
    int logB1 = LOG_BUCKET_REGIONS[1];
    uint32_t maskB1I = (UINT32_C(1) << std::min(logB1, logI)) - 1;

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    for (auto const & e: slice) {
        increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;
        e.transform_roots(transformed, Q);
        for (unsigned char i_root = 0; i_root != transformed.nr_roots;
             i_root++) {
            fbroot_t const r = transformed.get_r(i_root);
            bool const proj = transformed.get_proj(i_root);
            plattice_info const pli(transformed.get_q(), r, proj, logI);

            plattice_enumerator ple(pli, i_entry, logI);

            // Skip (i,j)=(0,0)
            ple.next(F);

            if (pli.is_discarded())
                continue;

            slice_offset_t const hint = ple.get_hint();
            WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
            const fbprime_t p = slice.get_prime(hint);
            WHERE_AM_I_UPDATE(w, p, p);
#else
            const fbprime_t p = 0;
#endif

            typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(
                0, p, hint, slice_index);

            // Handle the rare special cases
            /* projective-like:
             *
             * ple sets its first position in the (i,j) plane to (1,0),
             * which will typically be the _only_ hit in the normal case.
             *
             * there are more subtle cases that can show up though, because
             * of projective powers, and the combination with
             * adjust-strategy 2 (see bug 30012).
             *
             * the first hit (and only hit on the first line) can be (g,0)
             * for any g. but other lines may hit.
             */
            /* vertical:
             *
             * Root=0: only update is at (0,something).
             * note that "something" might be large !
             */
            if (UNLIKELY(ple.is_projective_like(logI))) {
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    int N = ple.get_x() >> logB;
                    int n =
                        ple.advance_to_end_of_row_or_smallest_region(maskB1I);
                    BA.push_row_update(slice_index, ple.get_inc_step(), N, n, u,
                                       w);
                    ple.next(F);
                }
            } else if (UNLIKELY(pli.is_vertical_line(logI))) {
                if (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.finish();
                }
            } else {
                /* Now, do the real work: the filling of the buckets */
                while (!ple.done(F)) {
                    if (LIKELY(ple.probably_coprime(F))) {
                        u.set_x(ple.get_x() & bmask);
                        BA.push_update(ple.get_x() >> logB, u, w);
                    }
                    ple.next(F);
                }
            }
        }
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA.move(BA);
}

/* {{{ */
/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_lowlevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, plattices_vector_t & plattices_vector,
                         bool first_reg MAYBE_UNUSED, where_am_I & w)
{
    int const logI = ws.conf.logI;
    qlattice_basis const & Q(ws.Q);

    /* The timer stuff is dealt with by the caller */
    slice_index_t const slice_index = plattices_vector.get_index();
    bucket_array_t<LEVEL, TARGET_HINT>
        BA; /* local copy. Gain a register + use stack */
    BA.move(orig_BA);

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    /* top level: the fence we care about is the one defined by the number
     * of bucket regions at this level, most probably (but J might
     * conceivably give a lower bound) */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J,
                                       BUCKET_REGIONS[LEVEL + 1]);

    /* yes, we want the level-1 regions here */
    int logB1 = LOG_BUCKET_REGIONS[1];
    uint32_t maskB1I = (UINT32_C(1) << std::min(logB1, logI)) - 1;

    int logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t bmask =
        (1UL << logB) - 1;

    for (auto & ple_orig: plattices_vector) {
        // Work with a copy, otherwise we don't get all optimizations.
        // Maybe with a wise use of the 'restrict' keyword, we might get
        // what we want, but this is C++11, anyway.
        //
        // FIXME we're c++11 now. Look into this.
        plattice_enumerator ple(ple_orig);

        slice_offset_t const hint = ple.get_hint();
        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        /* this is a bit expensive, since we're scanning all parts.
         * Fortunately it's only a debug call anyway. */
        fb_slice_interface const & slice =
            (*w->sides[w->side].fbs)[slice_index];
        fbprime_t const p = slice.get_prime(hint);
        WHERE_AM_I_UPDATE(w, p, p);
#else
        const fbprime_t p = 0;
#endif

        typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(0, p, hint,
                                                                slice_index);

        // Handle the rare special cases
        /* see fill_in_bucket_toplevel. */
        if (UNLIKELY(ple.is_projective_like(logI))) {
            if (Q.sublat.m)
                continue; /* XXX headaches ! */

            while (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                int N = ple.get_x() >> logB;
                int n = ple.advance_to_end_of_row_or_smallest_region(maskB1I);
                BA.push_row_update(slice_index, ple.get_inc_step(), N, n, u, w);
                ple.next(F);
            }
            /* we now do the end of loop normally: store x into ple_orig, and
             * then advance to the next area. This is because more rows can
             * be interesting as we go towards increasing j's
             */
        } else if (UNLIKELY(ple.is_vertical_line(logI))) {
            if (Q.sublat.m)
                continue; /* XXX headaches ! */

            if (!ple.done(F)) {
                u.set_x(ple.get_x() & bmask);
                BA.push_update(ple.get_x() >> logB, u, w);
                // ple.next(F);
                ple.finish();
            }
        } else {
            /* Now, do the real work: the filling of the buckets */
            // Without sublattices, we test (very basic) coprimality,
            // otherwise not atm. FIXME!
            if (!Q.sublat.m) {
                while (!ple.done(F)) {
                    if (LIKELY(ple.probably_coprime(F))) {
                        u.set_x(ple.get_x() & bmask);
                        BA.push_update(ple.get_x() >> logB, u, w);
                    }
                    ple.next(F);
                }
            } else {
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.next(F);
                }
            }
        }
        // save current position, and prepare for next area.
        ple_orig.set_x(ple.get_x());
        ple_orig.advance_to_next_area(F);
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA.move(BA);
}

template <int LEVEL> class fill_in_buckets_parameters : public task_parameters
{
  public:
    nfs_work & ws;
    nfs_aux & aux;
    int const side;
    fb_slice_interface const * slice;
    plattices_vector_t * plattices_vector; // content changed during fill-in
    plattices_dense_vector_t * plattices_dense_vector; // for sublat
    uint32_t const first_region0_index;
    where_am_I w;

    fill_in_buckets_parameters(nfs_work & _ws, nfs_aux & aux, int const _side,
                               fb_slice_interface const * _slice,
                               plattices_vector_t * _platt,
                               plattices_dense_vector_t * _dplatt,
                               uint32_t const _reg0, where_am_I const & w)
        : ws(_ws)
        , aux(aux)
        , side(_side)
        , slice(_slice)
        , plattices_vector(_platt)
        , plattices_dense_vector(_dplatt)
        , first_region0_index(_reg0)
        , w(w)
    {
    }
};

/* short of a better solution. I know some exist, but it seems way
 * overkill to me.
 *
 * This needs constexpr, though... So maybe I could use a more powerful
 * C++11 trick after all.
 */
#define PREPARE_TEMPLATE_INST_NAMES(F, suffix)                                 \
    template <int> struct CADO_CONCATENATE(F, _name) {                         \
    };                                                                         \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 0);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 1);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 2);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 3);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 4);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 5);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 6);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 7);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 8);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 9)

#define PREPARE_TEMPLATE_INST_NAME(F, suffix, k)                               \
    template <> struct CADO_CONCATENATE(F, _name)<k> {                         \
        static constexpr const char * value = #F "<" #k ">" suffix;            \
    }

/* By tweaking the "" argument, it is possible to have these names
 * embody a suffix like " on side ", so that it's possible tu run
 * parametric timer slots.
 */
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice_internal, "");
PREPARE_TEMPLATE_INST_NAMES(downsort, "");
PREPARE_TEMPLATE_INST_NAMES(downsort_tree, " (dispatcher only)");

#define TEMPLATE_INST_NAME(x, y) CADO_CONCATENATE(x, _name)<y>::value

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template <int LEVEL, typename TARGET_HINT>
static task_result *
fill_in_buckets_one_slice_internal(worker_thread * worker,
                                   task_parameters * _param, int)
{
    fill_in_buckets_parameters<LEVEL> * param =
        static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    nfs_work & ws(param->ws);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    MARK_TIMER_FOR_SIDE(timer, side);

    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(timer,
                TEMPLATE_INST_NAME(fill_in_buckets_one_slice_internal, LEVEL));

    w = param->w;
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, param->first_region0_index);

    try {

        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);

        /* Fill the buckets */
        try {
            fill_in_buckets_lowlevel<LEVEL>(BA, ws, *param->plattices_vector,
                                            (param->first_region0_index == 0),
                                            w);
        } catch (buckets_are_full & e) {
            wss.release_BA(BA);
            throw e;
        }
        /* Release bucket array again */
        wss.release_BA(BA);
    } catch (buckets_are_full & e) {
        delete param;
        throw e;
    }
    delete param;
    return new task_result;
}

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_wrapper(worker_thread * worker MAYBE_UNUSED,
                                 task_parameters * _param, int)
{
    fill_in_buckets_parameters<LEVEL> * param =
        static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        fb_slice<FB_ENTRY_TYPE> const * sl =
            dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice);
        ASSERT_ALWAYS(sl != NULL);
        fill_in_buckets_toplevel<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
            BA, ws, *sl, param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
        delete param;
        throw e;
    }
}
/* same for sublat */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_sublat_wrapper(worker_thread * worker MAYBE_UNUSED,
                                        task_parameters * _param, int)
{
    fill_in_buckets_parameters<LEVEL> * param =
        static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        fill_in_buckets_toplevel_sublat<LEVEL, FB_ENTRY_TYPE>(
            BA, ws,
            *dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice),
            param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
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
template <int LEVEL, typename TARGET_HINT> struct push_slice_to_task_list {
    thread_pool & pool;
    fill_in_buckets_parameters<LEVEL> model;
    push_slice_to_task_list(thread_pool & pool,
                            fill_in_buckets_parameters<LEVEL> const & m)
        : pool(pool)
        , model(m)
    {
    }
    size_t pushed = 0;
    template <typename T> void operator()(T const & s)
    {
        fill_in_buckets_parameters<LEVEL> * param =
            new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        typedef typename T::entry_t entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_wrapper<LEVEL, entry_t, TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};
template <int LEVEL, typename TARGET_HINT>
struct push_slice_to_task_list_saving_precomp {
    thread_pool & pool;
    fb_factorbase::slicing::part const & P;
    fill_in_buckets_parameters<LEVEL> model;
    /* precomp_plattice_dense_t == std::vector<plattices_dense_vector_t> */
    typename precomp_plattice_dense_t<LEVEL>::type & Vpre;
    size_t pushed = 0;
    push_slice_to_task_list_saving_precomp(
        thread_pool & pool, fb_factorbase::slicing::part const & P,
        fill_in_buckets_parameters<LEVEL> const & m,
        typename precomp_plattice_dense_t<LEVEL>::type & Vpre)
        : pool(pool)
        , P(P)
        , model(m)
        , Vpre(Vpre)
    {
    }
    template <typename T> void operator()(T const & s)
    {
        /* we're pushing the global index, relative to all fb parts */
        slice_index_t const idx = s.get_index();
        ASSERT_ALWAYS((size_t)idx == pushed);

        plattices_dense_vector_t & pre(Vpre[idx]);

        fill_in_buckets_parameters<LEVEL> * param =
            new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        param->plattices_dense_vector = &pre;

        typedef typename T::entry_t entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_sublat_wrapper<LEVEL, entry_t,
                                                    TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};

template <int LEVEL, typename TARGET_HINT>
static void fill_in_buckets_one_side(nfs_work & ws, nfs_aux & aux,
                                     thread_pool & pool, int const side,
                                     where_am_I & w)
{
    timetree_t & timer(aux.rt.timer);
    nfs_work::side_data & wss(ws.sides[side]);

    /* We're just pushing tasks, here. */
    BOOKKEEPING_TIMER(timer);

    fill_in_buckets_parameters<LEVEL> const model(ws, aux, side, NULL, NULL,
                                                  NULL, 0, w);

    if (!ws.Q.sublat.m) {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_wrapper */
        push_slice_to_task_list<LEVEL, TARGET_HINT> F(pool, model);
        wss.fbs->get_part(LEVEL).foreach_slice(F);
    } else {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_sublat_wrapper */
        auto & Vpre(wss.precomp_plattice_dense.get<LEVEL>());
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(LEVEL);

        /* This way we can spare the need to expose the copy contructor
         * of the container's value_type */
        if (ws.Q.sublat.i0 == 0 && ws.Q.sublat.j0 == 1) {
            /* first sublat */
            Vpre = typename precomp_plattice_dense_t<LEVEL>::type(P.nslices());
        }

        ASSERT_ALWAYS(Vpre.size() == P.nslices());
        push_slice_to_task_list_saving_precomp<LEVEL, TARGET_HINT> F(
            pool, P, model, Vpre);
        P.foreach_slice(F);
    }
}

void fill_in_buckets_toplevel(nfs_work & ws, nfs_aux & aux, thread_pool & pool,
                              int side, where_am_I & w)
{
    // per se, we're not doing anything here.
    // CHILD_TIMER(timer, __func__);

    if (ws.conf.needs_resieving()) {
        switch (ws.toplevel) {
        case 1:
            fill_in_buckets_one_side<1, shorthint_t>(ws, aux, pool, side, w);
            break;
        case 2:
            fill_in_buckets_one_side<2, shorthint_t>(ws, aux, pool, side, w);
            break;
        case 3:
            fill_in_buckets_one_side<3, shorthint_t>(ws, aux, pool, side, w);
            break;
        default:
            ASSERT_ALWAYS(0);
        }
    } else {
        switch (ws.toplevel) {
        case 1:
            fill_in_buckets_one_side<1, emptyhint_t>(ws, aux, pool, side, w);
            break;
        case 2:
            fill_in_buckets_one_side<2, emptyhint_t>(ws, aux, pool, side, w);
            break;
        case 3:
            fill_in_buckets_one_side<3, emptyhint_t>(ws, aux, pool, side, w);
            break;
        default:
            ASSERT_ALWAYS(0);
        }
    }
}

/* }}} */

/* multithreaded implementation of the downsort procedure. It becomes a
 * bottleneck sooner than one might think.
 *
 */

/* This is auxiliary only. We downsort stuff that wa already downsorted.
 * So it applies only if LEVEL+1 is itself not the toplevel.
 * For this reason, we must have a specific instantiation that reduces
 * this to a no-op if LEVEL+1>=3, because there's no longhint_t for level
 * 3 presently.
 */
template <int LEVEL>
static void downsort_aux(fb_factorbase::slicing const & fbs, nfs_work & ws,
                         nfs_aux & aux, thread_pool & pool, int side,
                         uint32_t bucket_index, where_am_I & w)
{
    nfs_work::side_data & wss(ws.sides[side]);
    // What comes from already downsorted data above:
    for (auto const & BAin: wss.bucket_arrays<LEVEL + 1, longhint_t>()) {
        pool.add_task_lambda(
            [&, side, w](worker_thread * worker, int bucket_index) {
                nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                timetree_t & timer(aux.get_timer(worker));
                ENTER_THREAD_TIMER(timer);
                MARK_TIMER_FOR_SIDE(timer, side);
                taux.w = w;
                CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                auto & BAout(
                    wss.reserve_BA<LEVEL, longhint_t>(wss.rank_BA(BAin)));
                downsort<LEVEL + 1>(fbs, BAout, BAin, bucket_index, taux.w);
                wss.template release_BA<LEVEL, longhint_t>(BAout);
            },
            bucket_index, 0);
    }
}
template <>
void downsort_aux<2>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}

// first_region0_index is a way to remember where we are in the tree.
// The depth-first is a way to process all the the regions of level 0 in
// increasing order of j-value.
// first_region0_index * nb_lines_per_region0 therefore gives the j-line
// where we are. This is what is called N by WHERE_AM_I and friends.

template <int LEVEL, bool WITH_HINTS>
static void downsort_tree_inner(
    nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p, thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> & precomp_plattice,
    where_am_I & w)
{
    int nsides = ws.las.cpoly->nb_polys;
    nfs_aux & aux(*aux_p);
    timetree_t & timer(aux.rt.timer);

    typedef typename hints_proxy<WITH_HINTS>::l my_longhint_t;
    typedef typename hints_proxy<WITH_HINTS>::s my_shorthint_t;

    CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort_tree, LEVEL));
    TIMER_CATEGORY(timer, sieving_mixed());
    ASSERT_ALWAYS(LEVEL > 0);

    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    for (int side = 0; side < nsides; ++side) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;

        WHERE_AM_I_UPDATE(w, side, side);
        TIMER_CATEGORY(timer, sieving(side));
        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        {
            /* This is the "dictionary" that maps slice indices to actual fb
             * entries. We rarely need it, except when downsorting short entries
             * in the case where we've eliminated the hint
             */
            fb_factorbase::slicing const & fbs(*wss.fbs);

            size_t const n2s =
                wss.bucket_arrays<LEVEL + 1, my_shorthint_t>().size();
            size_t const n1l = wss.bucket_arrays<LEVEL, my_longhint_t>().size();
            /* otherwise the code here can't work */
            ASSERT_ALWAYS(n2s == n1l);

            /* We create one output array for each input array, and we
             * process them in parallel. There would be various ways to
             * achieve that.
             */
            for (auto const & BAin:
                 wss.bucket_arrays<LEVEL + 1, my_shorthint_t>()) {
                pool.add_task_lambda(
                    [&, side, w](worker_thread * worker, int bucket_index) {
                        nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                        timetree_t & timer(aux.get_timer(worker));
                        taux.w = w;
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                        auto & BAout(wss.reserve_BA<LEVEL, my_longhint_t>(
                            wss.rank_BA(BAin)));
                        // This is a fake slice_index. For a longhint_t bucket,
                        // each update contains its own slice_index, directly
                        // used by apply_one_bucket
                        // and purge.
                        BAout.reset_pointers();
                        BAout.add_slice_index(
                            std::numeric_limits<slice_index_t>::max());
                        downsort<LEVEL + 1>(fbs, BAout, BAin, bucket_index,
                                            taux.w);
                        wss.template release_BA<LEVEL, my_longhint_t>(BAout);
                    },
                    bucket_index, 0);
            }
            // What comes from already downsorted data above. We put this in
            // an external function because we need the code to be elided or
            // LEVEL >= 2.
            if (LEVEL < ws.toplevel - 1)
                downsort_aux<LEVEL>(fbs, ws, aux, pool, side, bucket_index, w);
        }

        /* SECOND: fill in buckets at this level, for this region. */
        wss.reset_all_pointers<LEVEL, my_shorthint_t>();

        for (auto & it: precomp_plattice.get<LEVEL>()(side)) {
            pool.add_task(
                fill_in_buckets_one_slice_internal<LEVEL, my_shorthint_t>,
                new fill_in_buckets_parameters<LEVEL> {
                    ws, aux, side, (fb_slice_interface *)NULL, &it, NULL,
                    first_region0_index, const_ref(w)},
                0, 0, it.get_weight());
        }
    }

    /* RECURSE */
    if (LEVEL > 1) {
        for (int i = 0; i < ws.nb_buckets[LEVEL]; ++i) {
            size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
            uint32_t const N = first_region0_index + i * (BRS[LEVEL] / BRS[1]);
            downsort_tree<LEVEL - 1>(ws, wc_p, aux_p, pool, i, N,
                                     precomp_plattice, w);
        }
    } else {
        /* Prepare for PBR: we need to precompute the small sieve positions
         * for all the small sieved primes.
         *
         * For ws.toplevel==1, we don't reach here, of course, and the
         * corresponding initialization is done with identical code in
         * las.cpp
         */
        ASSERT(ws.toplevel > 1);
        for (int side = 0; side < nsides; side++) {
            nfs_work::side_data const & wss(ws.sides[side]);
            if (wss.no_fb())
                continue;
            pool.add_task_lambda(
                [=, &ws, &aux](worker_thread * worker, int) {
                    timetree_t & timer(aux.get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);
                    SIBLING_TIMER(timer, "prepare small sieve");
                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;
                    SIBLING_TIMER(timer, "small sieve start positions");
                    /* When we're doing 2-level sieving, there is probably
                     * no real point in doing ssdpos initialization in
                     * several passes.
                     */
                    small_sieve_prepare_many_start_positions(
                        wss.ssd, first_region0_index,
                        std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE,
                                 ws.nb_buckets[1]),
                        ws.conf.logI, ws.Q.sublat);
                    small_sieve_activate_many_start_positions(wss.ssd);
                },
                0);
        }

        pool.drain_queue(0);
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();
        auto exc = pool.get_exceptions<buckets_are_full>(0);
        if (!exc.empty()) {
            throw *std::max_element(exc.begin(), exc.end());
        }

        /* PROCESS THE REGIONS AT LEVEL 0 */
        process_many_bucket_regions(ws, wc_p, aux_p, pool, first_region0_index,
                                    w);

        /* We need that, because the next downsort_tree call in the loop
         * above (for LEVEL>1) will reset the pointers while filling the 1l
         * buckets -- and we read the 1l buckets from PBR.
         */
        if (ws.toplevel > 1)
            pool.drain_queue(0);
    }
}

template <int LEVEL>
void downsort_tree(
    nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p, thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> & precomp_plattice,
    where_am_I & w)
{
    if (ws.conf.needs_resieving()) {
        downsort_tree_inner<LEVEL, true>(ws, wc_p, aux_p, pool, bucket_index,
                                         first_region0_index, precomp_plattice,
                                         w);
    } else {
        downsort_tree_inner<LEVEL, false>(ws, wc_p, aux_p, pool, bucket_index,
                                          first_region0_index, precomp_plattice,
                                          w);
    }
}
/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(nfs_work &, std::shared_ptr<nfs_work_cofac>,
                      std::shared_ptr<nfs_aux>, thread_pool &, uint32_t,
                      uint32_t,
                      multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> &,
                      where_am_I &)
{
    ASSERT_ALWAYS(0);
}

// other fake instances to please level-2 instantiation.
template <>
void downsort<3>(fb_factorbase::slicing const &,
                 bucket_array_t<2, longhint_t> &,
                 bucket_array_t<3, longhint_t> const &, unsigned int,
                 where_am_I &)
{
    ASSERT_ALWAYS(0);
    /* Note that the const getter for 3l buckets (which don't exist) is
     * instantiated in bucket.cpp as erroring out, too */
}

// Now the two exported instances

template void downsort_tree<1>(
    nfs_work &, std::shared_ptr<nfs_work_cofac>, std::shared_ptr<nfs_aux> aux_p,
    thread_pool &, uint32_t, uint32_t,
    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> &, where_am_I &);
template void
downsort_tree<2>(nfs_work &, std::shared_ptr<nfs_work_cofac>,
                 std::shared_ptr<nfs_aux>, thread_pool &, uint32_t, uint32_t,
                 multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> &,
                 where_am_I &);
