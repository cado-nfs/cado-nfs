#ifndef CADO_SIQS_FILL_IN_BUCKETS_INL
#define CADO_SIQS_FILL_IN_BUCKETS_INL


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
    ASSERT_ALWAYS(!Q.sublat.m);
    auto & V(param->V);
    auto const & slice(param->slice);

    auto const index0 = ws.sides[param->side].fbs->get_part(LEVEL).first_slice_index;
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

        siqs_largesieve ple(Q, e, logB_minus_logI, i_entry);
        ASSERT_ALWAYS((ple.get_pp() >> logI) != 0); /* p must be > I */

        result.push_back(ple);
    }
    /* This is moved, not copied. Note that V is a reference. */
    V[relative_index] = std::move(result);
    delete param;
    return new task_result;
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
        std::vector<siqs_largesieve::T_elt> & T2scratch,
        MAYBE_UNUSED where_am_I & w)
{
    auto const pp = ple.pp;
    WHERE_AM_I_UPDATE(w, h, ple.hint);
    WHERE_AM_I_UPDATE(w, p, pp);
    typename BA_t::update_t::br_index_t const bmask = (1UL << logB) - 1U;
    typename BA_t::update_t u(0, pp, ple.hint, slice_index);

    for (uint32_t root: ple.roots) {
        ple.prepare_for_root(T2scratch, root, logI, j_high);

        /* First type of hit: t1+t2 in [0, I[. As T1 and T2 are sorted in
         * increasing order, we loop over T1, then T2 and break the inner
         * loop once t1+t2 >= I. So number of iteration is #hits+#T1.
         */
        for (auto const & [t1, d1]: ple.T1) {
            for (auto const & [t2, d2]: T2scratch) {
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
         * So number of iteration is #hits+#T2.
         */
        /* T2scratch does not change during iteration */
        auto it2 = T2scratch.rbegin();
        auto T2rend = T2scratch.rend();
        for (auto const & [t1, d1]: ple.T1) {
           for (auto it = it2; it != T2rend; ++it) {
               auto const [t2, d2] = *it;
               if (t1+t2 < pp) {
                   break; /* stop when t1+t2 becomes smaller than I */
               }
               uint64_t i = t1+t2-pp;
               if (i >> logI) {
                   ++it2; /* can be skipped for all next t1 values */
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

    std::vector<siqs_largesieve::T_elt> T2s;

    for (slice_offset_t hint = -1; auto const & e: slice) {
        ++hint;
        if (!Q.is_coprime_to(e.p))
            continue;
        if (discard_power_for_bucket_sieving(e))
            continue;

        ASSERT_ALWAYS((e.get_q() >> logI) != 0); /* p must be > I */
        siqs_largesieve ple(Q, e, logJ, hint);

        fill_in_buckets_siqs_compute_hits(ple, BA, logI, logB, 0U, slice_index,
                                          T2s, w);
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

    std::vector<siqs_largesieve::T_elt> T2s;

    ASSERT_ALWAYS(logB >= logI);
    uint32_t const j_high = first_region0_index << (logB - logI);

    for (auto & ple: plattices_vector) {
        fill_in_buckets_siqs_compute_hits(ple, BA, logI, logB, j_high,
                                          slice_index, T2s, w);
    }
    orig_BA = std::move(BA);
}

#endif /* CADO_SIQS_FILL_IN_BUCKETS_INL */
