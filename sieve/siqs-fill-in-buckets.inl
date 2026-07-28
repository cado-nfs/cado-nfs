#ifndef CADO_SIQS_FILL_IN_BUCKETS_INL
#define CADO_SIQS_FILL_IN_BUCKETS_INL


template <int LEVEL, class FB_ENTRY_TYPE>
void make_lattice_bases(worker_thread * worker MAYBE_UNUSED,
        int side,
        nfs_work & ws,
        nfs_aux &,
        siqs_special_q_data const & Q,
        precomp_plattice_t<LEVEL> & V,
        fb_slice<FB_ENTRY_TYPE> const & slice)
{
    int const logI = ws.conf.logI;
    ASSERT_ALWAYS(!Q.sublat.m);
    ASSERT_ALWAYS(side == 0);

    auto const index0 = ws.sides[side].fbs->get_part(LEVEL).first_slice_index;
    auto const index = slice.get_index();
    auto const relative_index = index - index0;
    ASSERT_ALWAYS(relative_index < V.size());

    static_assert(LEVEL != MAX_TOPLEVEL);
    int logB = LOG_BUCKET_REGIONS[LEVEL+1];
    ASSERT_ALWAYS(logB > logI);
    unsigned int logB_minus_logI = logB - logI;

    plattices_vector_t result(index, slice.get_weight());
    for (slice_offset_t i_entry = -1; auto const & e: slice) {
        ++i_entry;
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;

        ASSERT_ALWAYS((e.get_q() >> logI) != 0); /* p must be > I */

        result.emplace_back(Q, e, logB_minus_logI, i_entry);
    }
    /* This is moved, not copied. Note that V is a reference. */
    V[relative_index] = std::move(result);
}

/* T1 and T2 are already computed for all the possible n lowest bit of j.
 * Given fixed values for the highest bits of j, this function computes the hits
 * for all the roots.
 * It iterates over all roots and calls the method prepare_for_root which
 * computes T2' (using precomputed T2).
 * It then looks for pair (t1, d1), (t2, d2) in T1 x T2' such that t1+t2 % p is
 * in [0, I[ (note that t1+t2 is in [0, 2p[). It corresponds to the update
 *  N = (((d1 xor d2) << logI) | (t1 + t2 % p)) >> logB;
 *  x = (((d1 xor d2) << logI) | (t1 + t2 % p)) % 2^logB;
 *
 * Assumes pp > I. It should be checked in make_lattices_base or in
 * fill_in_buckets_toplevel.
 */
template <typename BA_t>
void fill_in_buckets_siqs_compute_hits(
        siqs_largesieve & ple,
        BA_t & BA,
        int const logI,
        int const logB,
        uint32_t const j_high,
        slice_index_t const slice_index,
        MAYBE_UNUSED where_am_I & w)
{
    auto const pp = ple.pp;
    WHERE_AM_I_UPDATE(w, h, ple.hint);
    WHERE_AM_I_UPDATE(w, p, pp);
    typename BA_t::update_t::br_index_t const bmask = (1UL << logB) - 1U;
    typename BA_t::update_t u(0, pp, ple.hint, slice_index);

    uint32_t mask = ple.prepare_for_root_mask(j_high);
    uint32_t s = ple.prepare_for_root_precomp(logI, j_high);

    for (unsigned int i = 0; i < ple.nroots; ++i) {
        auto T2s = ple.prepare_for_root(ple.roots[i], s, mask);

        /* First type of hit: t1+t2 in [0, I[. As T1 and T2 are sorted in
         * increasing order, we loop over T1, then T2 and break the inner
         * loop once t1+t2 >= I. So number of iteration is #hits+#T1.
         */
        for (auto const [t2, d2]: T2s) {
            for (auto const [t1, d1]: ple.T1) {
                uint64_t i = t1+t2;
                if (i >> logI) {
                    break; /* stop when i=t1+t2 becomes larger than I */
                }
                uint64_t j = d1 xor d2;
                uint64_t x = (j << logI) | i;
                u.set_x(x & bmask);
                BA.push_update(x >> logB, u, w);
            }
        }

        /* Second type of hit: t1+t2 in [p, p+I[. As T1 and T2 are sorted in
         * increasing order, we loop over T1, then loop over T2 in reverse
         * order. We look for t2i such that t1+t2i < p+I (and record the
         * position) and then loop over T2 (again in reverse order) starting
         * from this value and break the inner loop once t1+t2 < p.
         * Starting from the second value of t1, we look for t2i using the
         * previous recorded value.
         * So number of iteration is #hits+#T1+(a term bounded by #T2).
         */
        /* T1 does not change during iteration */
        auto it1 = ple.T1.rbegin();
        auto T1rend = ple.T1.rend();
        for (auto const [t2, d2]: T2s) {
            for (auto it = it1; it != T1rend; ++it) {
                auto const [t1, d1] = *it;
                auto t12 = t1+t2;
                if (t12 < pp) {
                    break; /* stop when t1+t2 becomes smaller than p */
                }
                uint64_t i = t12-pp;
                if (i >> logI) {
                    ++it1; /* can be skipped for all next t1 values */
                    continue; /* skip if i=t1+t2-p is too large */
                }
                uint64_t j = d2 xor d1;
                uint64_t x = (j << logI) | i;
                u.set_x(x & bmask);
                BA.push_update(x >> logB, u, w);
            }
        }
    }
}

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void fill_in_buckets_toplevel_sublat(
    bucket_array_t<LEVEL, TARGET_HINT> &,
    nfs_work &,
    siqs_special_q_data const &,
    fb_slice<FB_ENTRY_TYPE> const &,
    plattices_dense_vector_t *,
    where_am_I &)
{
    throw std::runtime_error("sublat is not supported in SIQS");
}

template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
void
fill_in_buckets_toplevel(bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
                         nfs_work & ws, fb_slice<FB_ENTRY_TYPE> const & slice,
                         siqs_special_q_data const & Q,
                         plattices_dense_vector_t * /* unused */,
                         where_am_I & w)
{
    if (LEVEL == 3) {
        throw std::runtime_error("Level 3 bucket sieving is not supported in SIQS");
    }
    int const logI = ws.conf.logI;
    size_t logJ = nbits(ws.J) - 1u; /* 2^m has m+1 bits */
    ASSERT_ALWAYS(ws.J == 1u << logJ); /* J must be a power of 2 */

    ASSERT_ALWAYS(!Q.sublat.m);

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    slice_index_t const slice_index = slice.get_index();

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    int logB = LOG_BUCKET_REGIONS[LEVEL];

    std::size_t memsize = siqs_largesieve::memory_required(Q, logJ);
    std::byte data[memsize];
    std::span<std::byte> mem(data, memsize);

    for (slice_offset_t hint = -1; auto const & e: slice) {
        ++hint;
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;

        ASSERT_ALWAYS((e.get_q() >> logI) != 0); /* p must be > I */
        siqs_largesieve ple(Q, e, logJ, hint, mem);

        fill_in_buckets_siqs_compute_hits(ple, BA, logI, logB, 0U, slice_index,
                                          w);
    }
    orig_BA = std::move(BA);
}

/* TARGET_HINT is shorthint_t or void */
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_lowlevel(
        bucket_array_t<LEVEL, TARGET_HINT> & orig_BA,
        nfs_work & ws,
        siqs_special_q_data const & Q,
        plattices_vector_t & plattices_vector,
        uint32_t const first_region0_index,
        where_am_I & w)
{
    int const logI = ws.conf.logI;
    size_t logJ = nbits(ws.J) - 1u; /* 2^m has m+1 bits */
    ASSERT_ALWAYS(ws.J == 1u << logJ); /* J must be a power of 2 */

    ASSERT_ALWAYS(!Q.sublat.m);

    /* The timer stuff is dealt with by the caller */
    slice_index_t const slice_index = plattices_vector.get_index();

    /* local copy. Gain a register + use stack */
    bucket_array_t<LEVEL, TARGET_HINT> BA = std::move(orig_BA);

    /* Write new set of pointers for the new slice */
    BA.add_slice_index(slice_index);
    WHERE_AM_I_UPDATE(w, i, slice_index);

    int logB = LOG_BUCKET_REGIONS[LEVEL];

    ASSERT_ALWAYS(logB >= logI);
    uint32_t const j_high = first_region0_index << (logB - logI);

    for (auto & ple: plattices_vector) {
        fill_in_buckets_siqs_compute_hits(ple, BA, logI, logB, j_high,
                                          slice_index, w);
    }
    orig_BA = std::move(BA);
}

#endif /* CADO_SIQS_FILL_IN_BUCKETS_INL */
