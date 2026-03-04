#ifndef CADO_LAS_FILL_IN_BUCKETS_INL
#define CADO_LAS_FILL_IN_BUCKETS_INL


template <int LEVEL, class FB_ENTRY_TYPE>
static task_result * make_lattice_bases(worker_thread * worker MAYBE_UNUSED,
                                        task_parameters * _param, int)
{
    auto const * param =
        static_cast<
            make_lattice_bases_parameters<LEVEL, FB_ENTRY_TYPE> const *>(
            _param);

    nfs_work const & ws(param->ws);
    auto const & Q(param->Q);
    int const logI = ws.conf.logI;
    sublat_t const & sublat(Q.sublat);
    auto & V(param->V);
    auto const & slice(param->slice);

    auto const index0 = ws.sides[param->side].fbs->get_part(LEVEL).first_slice_index;
    auto const index = slice.get_index();
    auto const relative_index = index - index0;
    ASSERT_ALWAYS(relative_index < V.size());

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
    V[relative_index] = std::move(result);
    delete param;
    return new task_result;
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
    qlattice_basis const & Q,
    fb_slice<FB_ENTRY_TYPE> const & slice,
    plattices_dense_vector_t * p_precomp_slice, where_am_I & w)
{
    int const logI = ws.conf.logI;

    plattices_dense_vector_t & precomp_slice(*p_precomp_slice);

    ASSERT_ALWAYS(Q.sublat.m);
    bool const first_sublat = Q.sublat.i0 == 0 && Q.sublat.j0 == 1;
    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

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
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, fb_slice<FB_ENTRY_TYPE> const & slice,
                         qlattice_basis const & Q,
                         plattices_dense_vector_t * /* unused */,
                         where_am_I & w)
{
    int const logI = ws.conf.logI;

    ASSERT_ALWAYS(!Q.sublat.m);

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

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
                    BA.push_row_update(slice_index, ple.get_inc_step(),
                            N, n, u, w);
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
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, typename TARGET_HINT>
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
    orig_BA = std::move(BA);
}

#endif /* CADO_LAS_FILL_IN_BUCKETS_INL */
