#ifndef CADO_LAS_FILL_IN_BUCKETS_INL
#define CADO_LAS_FILL_IN_BUCKETS_INL

#include <cstdint>

#include <utility>

#include "bucket.hpp"
#include "bucket-push-update.hpp"       // IWYU pragma: keep
#include "fb-types.hpp"
#include "fb.hpp"
#include "las-plattice.hpp"
#include "las-qlattice.hpp"
#include "threadpool.hpp"
#include "las-threads-work-data.hpp"
#include "las-fill-in-buckets.hpp"
#include "macros.h"
#include "chronograms.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "las-sublat.hpp"

template <int LEVEL, class FB_ENTRY_TYPE, uint32_t M>
void make_lattice_bases(
        worker_thread * worker,
        int side,
        nfs_work & ws,
        qlattice_basis const & Q,
        /* We're passing it by value because it lives as a temporary at
         * the caller site */
        sublat_t<M> sublat,
        precomp_plattice_t<LEVEL> & V,
        fb_slice<FB_ENTRY_TYPE> const & slice)
{
    int const logI = ws.conf.logI;

    auto const index0 = ws.sides[side].fbs->get_part(LEVEL).first_slice_index;
    auto const index = slice.get_index();
    auto const relative_index = index - index0;
    ASSERT_ALWAYS(relative_index < V.size());

    auto tt = worker->trace(chronograms::PCLAT { side, LEVEL, index });

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
#ifdef BUCKET_SIEVE_POWERS
        /* Same guard as in fill_in_buckets_toplevel_impl: a prime power
         * whose prime divides the sublattice modulus cannot be handled by
         * the sublattice machinery, because the congruence it defines
         * cannot be divided through by m. Bucket-sieving powers means
         * that, say, 2^15 shows up here as a bucket-sieved entry with
         * p = 2, and without this test it produced updates that do not
         * divide -- which only showed up once the toplevel was 2 or more,
         * since at toplevel 1 the guard in the toplevel fill covers it. */
        if (sublat_t<M>::not_coprime(e.p))
            continue;
#endif
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
            if constexpr (M == 1)
                ple.next(F);
            if (LIKELY(!pli.is_discarded()))
                result.push_back(ple);
        }
    }
    /* This is moved, not copied. Note that V is a reference. */
    V[relative_index] = std::move(result);
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
// With sublattices, the Franke-Kleinjung bases do not depend on which
// residue class is being sieved, so they are computed once -- while
// processing the class (i,j) == (0,1) mod m -- stored in precomp_slice,
// and replayed for the other classes.
//
// The three ways the top level can be run. Selecting between them with a
// template parameter rather than a runtime flag is what lets the three
// bodies below be a single one: every `if constexpr` in
// fill_in_buckets_toplevel_impl disappears at compile time, so the inner
// enumeration loop is exactly as tight as it was when this file carried
// three hand-specialised copies of it.

template <int LEVEL,
         class FB_ENTRY_TYPE,
         hint_type TARGET_HINT,
         uint32_t M,
         bool is_first_sublat>
static void fill_in_buckets_toplevel_impl(
    bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
    nfs_work & ws,
    fb_slice<FB_ENTRY_TYPE> const & slice,
    qlattice_basis const & Q,
    sublat_t<M> const & sublat,
    plattices_dense_vector_t * p_precomp_slice,
    where_am_I & w)
{
    constexpr bool with_sublat = M > 1;

    int const logI = ws.conf.logI;

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    /* top level: the fence we care about is the one defined by J */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J);

    /* yes, we want the level-1 regions here */
    int const logB1 MAYBE_UNUSED = LOG_BUCKET_REGIONS[1];
    uint32_t const maskB1I MAYBE_UNUSED =
        (UINT32_C(1) << std::min(logB1, logI)) - 1;

    int const logB = LOG_BUCKET_REGIONS[LEVEL];
    typename bucket_array_t<LEVEL, TARGET_HINT>::update_t::br_index_t const
        bmask = (1UL << logB) - 1;

    /* Everything that happens once a p-lattice is known, which is the
     * only part that is performance-critical. */
    auto handle_one_lattice = [&](plattice_info const & pli,
                                  slice_offset_t const i_entry) {
        plattice_enumerator ple = [&]() {
            return plattice_enumerator(pli, i_entry, logI, sublat);
        }();

        if constexpr (with_sublat) {
            if (ple.done(F))
                return;
        } else {
            /* Skip (i,j)=(0,0) */
            ple.next(F);
        }

        if (pli.is_discarded())
            return;

        slice_offset_t const hint = ple.get_hint();
        ASSERT(hint == i_entry);
        WHERE_AM_I_UPDATE(w, h, hint);
#ifdef TRACE_K
        const fbprime_t p = slice.get_prime(hint);
        WHERE_AM_I_UPDATE(w, p, p);
#else
        const fbprime_t p = 0;
#endif

        typename bucket_array_t<LEVEL, TARGET_HINT>::update_t u(0, p, hint,
                                                                slice_index);

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
            if constexpr (with_sublat) {
                /* plattice_enumerator_base::starting_point() bails out
                 * with UMAX whenever i1 or j1 is zero, so under
                 * sublattices these lattices contribute nothing at all
                 * and ple is already done. That is a loss of yield, not
                 * a correctness problem, but it is a real loss and it
                 * should be fixed rather than papered over. */
                return;
            } else {
                while (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    int const N = ple.get_x() >> logB;
                    int const n =
                        ple.advance_to_end_of_row_or_smallest_region(maskB1I);
                    BA.push_row_update(slice_index, ple.get_inc_step(), N, n, u,
                                       w);
                    ple.next(F);
                }
            }
        } else if (UNLIKELY(pli.is_vertical_line(logI))) {
            if constexpr (with_sublat) {
                return; /* same story as just above */
            } else {
                if (!ple.done(F)) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                    ple.finish();
                }
            }
        } else {
            /* Now, do the real work: the filling of the buckets */
            while (!ple.done(F)) {
                /* Without sublattices, we test (very basic) coprimality.
                 * With sublattices we do not: for an even modulus the
                 * test is vacuous anyway, since the parities of i and j
                 * are then fixed over the whole class and the class
                 * where both are even is never sieved. For an odd
                 * modulus it would still be worth something. */
                bool keep = true;
                if constexpr (!with_sublat)
                    keep = LIKELY(ple.probably_coprime(F));
                if (keep) {
                    u.set_x(ple.get_x() & bmask);
                    BA.push_update(ple.get_x() >> logB, u, w);
                }
                ple.next(F);
            }
        }
    };

    if (!is_first_sublat) {
        for (auto const & psl: *p_precomp_slice) {
            slice_offset_t const hint = psl.get_hint();
            plattice_info const pli(psl.unpack(slice.begin()[hint].get_q()));
            handle_one_lattice(pli, hint);
        }
    } else {
        if constexpr (M > 1) {
            /* This is an upper bound (exact for fixed-root entries), and
             * it avoids the slack of a geometric growth, which the
             * memory estimate does not count. */
            p_precomp_slice->reserve(slice.size() *
                    std::max(1, int(FB_ENTRY_TYPE::fixed_nr_roots)));
        }
        typename FB_ENTRY_TYPE::transformed_entry_t transformed;
        slice_offset_t i_entry = 0;
        for (auto const & e: slice) {
            increment_counter_on_dtor<slice_offset_t> const _dummy(i_entry);
            if (!Q.is_coprime_to(e.p))
                continue;
#ifdef BUCKET_SIEVE_POWERS
            if constexpr (with_sublat) {
                /* the combination of bucket-sieving powers + sublattices
                 * means that powers of the primes that divide the
                 * sublattice determinant may be bucket-sieved. And of
                 * course, that leads to problems. */
                if (sublat_t<M>::not_coprime(e.p))
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

                if constexpr (M > 1) {
                    // In sublat mode, save it for later use
                    p_precomp_slice->push_back(
                        plattice_info_dense_t(pli, i_entry,
                                              transformed.get_q()));
                }

                handle_one_lattice(pli, i_entry);
            }
        }
    }
    // printf("%.3f\n", BA.max_full());
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL,
    class FB_ENTRY_TYPE,
    hint_type TARGET_HINT>
static void fill_in_buckets_toplevel(
    bucket_array_t<LEVEL, TARGET_HINT> & orig_BA, nfs_work & ws,
    fb_slice<FB_ENTRY_TYPE> const & slice,
    qlattice_basis const & Q,
    plattices_dense_vector_t * p_precomp_slice,
    where_am_I & w)
{
    static_assert(!TARGET_HINT::is_long_v);
    /* The sublattice we visit first computes the Franke-Kleinjung bases
     * and stores them; the ones after that replay them. This is the only
     * runtime test, and it is outside everything that matters. */
    dispatch_sublat(Q.sublat, [&](auto sublat) {
            constexpr uint32_t M = decltype(sublat)::modulus;
            if (Q.sublat.is_first())
            fill_in_buckets_toplevel_impl<LEVEL, FB_ENTRY_TYPE, TARGET_HINT, M, true>(
                    orig_BA, ws, slice, Q, sublat, p_precomp_slice, w);
            else
            fill_in_buckets_toplevel_impl<LEVEL, FB_ENTRY_TYPE, TARGET_HINT, M, false>(
                    orig_BA, ws, slice, Q, sublat, p_precomp_slice, w);
            });
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, hint_type TARGET_HINT>
static void
fill_in_buckets_lowlevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws,
                         qlattice_basis const & Q,
                         plattices_vector_t & plattices_vector,
                         uint32_t const /*first_region0_index*/,
                         where_am_I & w)
{
    int const logI = ws.conf.logI;

    /* The timer stuff is dealt with by the caller */
    slice_index_t const slice_index = plattices_vector.get_index();

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    /* we used to look up BUCKET_REGIONS[LEVEL + 1] here, which doesn't
     * really seem to make sense. I expect that ws.J actually yields a
     * stricter bound in all cases.
     */
    plattice_enumerator::fence const F(ws.conf.logI, ws.J,
            (LEVEL + 1 < FB_MAX_PARTS ? BUCKET_REGIONS[LEVEL + 1] : SIZE_MAX));
    /* just checking... */
    if (LEVEL == FB_MAX_PARTS - 1)
        ASSERT_ALWAYS((ws.J << ws.conf.logI) < (BUCKET_REGIONS[LEVEL] << 8));

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
            if (Q.sublat.m > 1)
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
            if (Q.sublat.m > 1)
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
            if (Q.sublat.m == 1) {
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
    orig_BA = std::move(BA);
}

#endif /* CADO_LAS_FILL_IN_BUCKETS_INL */
