#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <new>
#ifdef TRACE_K
#include <type_traits>
#endif

/* bucket.hpp and bucket-push-update.hpp sort of go hand in hand, but as
 * they're written we can't make two standalone files that include
 * eachother */
#include "bucket-push-update.hpp"

#include "bucket.hpp"
#include "memory.h"
#include "verbose.h"

#include "fb-types.hpp"
#include "fb.hpp"
#include "iqsort.h"
#include "las-config.hpp"
#include "las-memory.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "macros.h"
#ifdef TRACE_K
#include "las-output.hpp"
#endif
#ifdef HAVE_SSE2
#include "smallset.hpp"
#endif
#include "portability.h"

template <int LEVEL, typename HINT> struct bucket_update_t;

static size_t bucket_misalignment(size_t const sz, size_t const sr MAYBE_UNUSED)
{
    size_t size = sz;
#ifdef HAVE_SSE2
    /* Round up to a multiple of CACHELINESIZE to make SSE2 aligned accesses
       work */
    size = iceildiv(size, 16U) * 16U;
#endif
    return size;
}

/* Set the read and write pointers of the buckets back to the respective bucket
   start, and set nr_slices back to 0. */
template <int LEVEL, typename HINT>
void bucket_array_t<LEVEL, HINT>::reset_pointers()
{
    std::copy_n(bucket_start.get(), pointer_pack, bucket_write.get());
    std::copy_n(bucket_start.get(), pointer_pack, bucket_read.get());
    slice_start.clear();
    slice_index.clear();
    for (auto & r: row_updates)
        r.clear();
}

/* Allocate enough memory to be able to store _n_bucket buckets, each of at
   least min_bucket_size entries. If enough (or more) memory was already
   allocated, does not shrink the allocation. */
template <int LEVEL, typename HINT>
void bucket_array_t<LEVEL, HINT>::allocate_memory(
    las_memory_accessor & memory, uint32_t const new_n_bucket,
    double const fill_ratio, int logI, slice_index_t const prealloc_slices)
{
    static_assert(LEVEL < FB_MAX_PARTS);
    /* Don't try to allocate anything, nor print a message, for sieving levels
     * where the corresponding factor base part is empty. We do want to
     * reset the pointers, though!!
     */
    if (fill_ratio == 0.) {
        reset_pointers();
        return;
    }

    /* We'll allocate bucket regions of size 2^LOG_BUCKET_REGIONS[LEVEL],
     * and those will be used to cover lines of size 2^logI.
     *
     * If LOG_BUCKET_REGIONS[LEVEL] > logI, then we'll have both
     * even and odd lines in there, so 75% of the locations will receive
     * updates.
     *
     * If LOG_BUCKET_REGIONS[LEVEL] <= logI, it's different. The
     * line ordinate for bucket region number N will be N, maybe shifted
     * right by a few bits (logI - LOG_BUCKET_REGIONS[LEVEL]).
     * Bucket regions for which this line ordinate is even will receive
     * updates for 50% of the locations only, in contrast to 100% when the
     * ordinate is even.
     */

    size_t const Q = 0.25 * fill_ratio * BUCKET_REGIONS[LEVEL];

    size_t bs_even, bs_odd;
    size_t new_big_size;
    /* The number of updates has a variance which is very often close to
     * its mean. By aiming at a number of updates that is at most
     * N+ndev*sqrt(N), the probability to overrun is at most the probability
     * that this random variable exceeds its mean by ndev standard deviations
     * or more. Now it's a sum of random variables (number of hits per line
     * for p) that are independent but not identically distributed.
     * Probably there are some tail bounds which can tell me how rare this
     * is, depending on ndev.
     */
    int const ndev = NB_DEVIATIONS_BUCKET_REGIONS;

    uint32_t bitmask_line_ordinate = 0;
    if (LOG_BUCKET_REGIONS[LEVEL] <= logI) {
        bitmask_line_ordinate = UINT32_C(1)
                                << (logI - LOG_BUCKET_REGIONS[LEVEL]);
        ASSERT_ALWAYS(new_n_bucket % 2 == 0);
        bs_even = 2 * Q + ndev * sqrt(2 * Q);
        bs_odd = 4 * Q + ndev * sqrt(4 * Q);
        bs_even = bucket_misalignment(bs_even, sizeof(update_t));
        bs_odd = bucket_misalignment(bs_odd, sizeof(update_t));
        new_big_size =
            (bs_even + bs_odd) * (new_n_bucket / 2);
    } else {
        bs_even = bs_odd = 3 * Q + ndev * sqrt(3 * Q);
        bs_even = bucket_misalignment(bs_even, sizeof(update_t));
        bs_odd = bucket_misalignment(bs_odd, sizeof(update_t));
        new_big_size = bs_odd * new_n_bucket;
    }

    /* add 1 megabyte to each bucket array, so as to deal with overflowing
     * buckets */
    new_big_size += (1 << 20) / sizeof(update_t);

    size_t const new_size_b_align =
        ((sizeof(void *) * new_n_bucket + 0x3F) & ~((size_t)0x3F));
    size_t const new_pointer_pack = new_size_b_align / sizeof(update_t *);

    if (new_big_size > big_data.get_deleter().size) {
        if (bitmask_line_ordinate) {
            verbose_fmt_print(
                0, 3,
                "# [{}{}] Allocating {} bytes for {} buckets"
                " of {} to {} update entries of {} bytes each\n",
                LEVEL, HINT::rtti[0], new_big_size * sizeof(update_t),
                new_n_bucket, bs_even, bs_odd, sizeof(update_t));
        } else {
            verbose_fmt_print(
                0, 3,
                "# [{}{}] Allocating {} bytes for {} buckets"
                " of {} update entries of {} bytes each\n",
                LEVEL, HINT::rtti[0], new_big_size * sizeof(update_t),
                new_n_bucket, bs_even, sizeof(update_t));
        }
        big_data = memory.make_unique_physical_array<update_t>(new_big_size, true);
    }

    if (!big_data)
        throw std::bad_alloc();

    // bucket_size = new_bucket_size;
    n_bucket = new_n_bucket;

    if (new_pointer_pack > pointer_pack) {
        if (pointer_pack) {
            verbose_fmt_print(0, 1,
                    "# [{}{}] Changing bucket allocation"
                    " from {} to {} pointers\n",
                    LEVEL, HINT::rtti[0], pointer_pack,
                    new_pointer_pack);
        }
        pointer_pack = new_pointer_pack;

        /* bucket_start is allocated as an array of n_bucket+1 pointers */
        size_t const alloc_bstart = std::max(pointer_pack, n_bucket + 1);
        bucket_start = make_unique_aligned_array<update_t *>(alloc_bstart, 0x40);
        std::fill_n(bucket_start.get(), alloc_bstart, nullptr);

        bucket_read = make_unique_aligned_array<update_t *>(pointer_pack, 0x40);
        std::fill_n(bucket_read.get(), pointer_pack, nullptr);

        /* there is probably _some_ sense in making bucket_write
         * page-aligned, but the exact rationale was never spelled out in
         * clear, I believe. Probably because the set of pointers at
         * bucket_write[] is the one that is constantly accessed during
         * fill-in-buckets, and it would be embarrassing if it had to
         * span several pages.
         */
        bucket_write = make_unique_aligned_array<update_t *>(pointer_pack, pagesize());
        std::fill_n(bucket_write.get(), pointer_pack, nullptr);
    }

    if (prealloc_slices) {
        slice_start.reserve(prealloc_slices);
        slice_index.reserve(prealloc_slices);
    }

    /* Spread bucket_start pointers equidistantly over the big_data array */
    update_t * cur = big_data.get();
    bucket_start[0] = cur;
    for (uint32_t i = 0; bucket_start[i] = cur, i < n_bucket; i++) {
        cur += ((i & bitmask_line_ordinate) != 0) ? bs_odd : bs_even;
    }

    row_updates.assign(n_bucket, typename decltype(row_updates)::value_type());

    reset_pointers();
#ifdef SAFE_BUCKET_ARRAYS
    verbose_fmt_print(0, 0, "# WARNING: SAFE_BUCKET_ARRAYS is on !\n");
#endif
}

/* Returns how full the fullest bucket is, as a fraction of its size */
template <int LEVEL, typename HINT>
double bucket_array_t<LEVEL, HINT>::max_full(unsigned int * fullest_index) const
{
    double max = 0;
    for (unsigned int i = 0; i < n_bucket; ++i) {
        double const j =
            (double)nb_of_updates(i) / room_allocated_for_updates(i);
        if (max < j) {
            max = j;
            if (fullest_index)
                *fullest_index = i;
        }
    }
    return max;
}
template <int LEVEL, typename HINT>
double bucket_array_t<LEVEL, HINT>::average_full() const
{
    size_t a = 0;
    for (unsigned int i = 0; i < n_bucket; ++i)
        a += nb_of_updates(i);
    return (double)a / (bucket_start[n_bucket] - bucket_start[0]);
}

template <int LEVEL, typename HINT>
void bucket_array_t<LEVEL, HINT>::log_this_update(update_t const update
                                                  MAYBE_UNUSED,
                                                  uint64_t const bucket_number
                                                  MAYBE_UNUSED,
                                                  where_am_I & w
                                                  MAYBE_UNUSED) const
{
#if defined(TRACE_K)
    size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
    unsigned int const saveN = w->N;
    /* flatten the (N,x) coordinate as if relative to a unique array of
     * level-1 bucket regions */
    unsigned int const x = update.x % BRS[1];
    unsigned int const N =
        w->N + bucket_number * BRS[LEVEL] / BRS[1] + (update.x / BRS[1]);

    WHERE_AM_I_UPDATE(w, x, x);
    WHERE_AM_I_UPDATE(w, N, N);

    if (trace_on_spot_Nx(w->N, w->x)) {
        verbose_fmt_print(
            TRACE_CHANNEL, 0,
            "# Pushed hit at location (x={}, side {}), from factor base entry "
            "(slice_index={}, slice_offset={}, p={}), "
            "to BA<{}>[{}]\n",
            (unsigned int)w->x, w->side, (unsigned int)w->i, (unsigned int)w->h,
            w->p, LEVEL, (unsigned int)w->N);
        if (std::is_same<HINT, longhint_t>::value) {
            verbose_fmt_print(TRACE_CHANNEL, 0,
                              "# Warning: did not check divisibility"
                              " during downsorting p={}\n",
                              w->p);
        } else {
            ASSERT_ALWAYS(test_divisible(w));
        }
    }
    WHERE_AM_I_UPDATE(w, N, saveN);
#endif
}

#if 0
/* no longer used */
/* A compare function suitable for sorting updates in order of ascending x
   with qsort() */
static inline int
bucket_cmp_x (const bucket_update_t<1, longhint_t>  *a, const bucket_update_t<1, longhint_t>  *b)
{
  if (a->x < b->x)
    return -1;
  if (a->x == b->x)
    return 0;
  return 1;
}
#endif

template <int LEVEL, typename HINT> void bucket_single<LEVEL, HINT>::sort()
{
//  qsort (start, write - start, sizeof (bucket_update_t<1, longhint_t> ),
//	 (int(*)(const void *, const void *)) &bucket_cmp_x);
#define islt(a, b) ((a)->x < (b)->x)
    QSORT(update_t, start.get(), size(), islt);
#undef islt
}

/*
void bucket_primes_t::purge(bucket_array_t<1, shorthint_t> const & BA,
                            int const i, fb_factorbase::slicing const & fb,
                            unsigned char const * S)
{
    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        slice_index_t const slice_index = BA.get_slice_index(i_slice);
        for (auto const & it: BA.slice_range(i, i_slice)) {
            if (UNLIKELY(S[it.x] != 255)) {
                fbprime_t const p = fb[slice_index].get_prime(it.hint);
                push_update(bucket_update_t<1, primehint_t>(it.x, p, 0, 0));
            }
        }
    }
    // read the projective updates, and put the surviving ones to the
    // main array.
    for (auto const & ru: BA.row_updates[i]) {
        fbprime_t p = fb[ru.slice_index].get_prime(ru.hint);
        // longhint_t h(0, ru.hint, ru.slice_index);
        // have to walk the updates to see the ones that match
        bucket_update_t<1, shorthint_t> u = ru;
        for (int nx = ru.n + 1; nx--;) {
            if (UNLIKELY(S[u.x] != 255)) {
                push_update(bucket_update_t<1, primehint_t>(u.x, p, 0, 0));
            }
            u.x += ru.inc;
        }
    }
}
*/

/*
template <typename HINT>
static inline bucket_update_t<1, longhint_t>
to_longhint(const bucket_update_t<1, HINT> &update, slice_index_t slice_index);
*/

static inline bucket_update_t<1, longhint_t>
to_longhint(bucket_update_t<1, shorthint_t> const & update,
            slice_index_t const slice_index)
{
    return { update.x, 0, update.hint, slice_index };
}

static inline bucket_update_t<1, longhint_t>
to_longhint(bucket_update_t<1, longhint_t> const & update,
            slice_index_t const slice_index MAYBE_UNUSED)
{
    return update;
}

template <typename HINT>
void bucket_array_complete::purge(bucket_array_t<1, HINT> const & BA,
                                  int const i, unsigned char const * S)
{
    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        slice_index_t const slice_index = BA.get_slice_index(i_slice);
        for (auto const & it: BA.slice_range(i, i_slice)) {
            if (UNLIKELY(S[it.x] != 255)) {
                push_update(to_longhint(it, slice_index));
            }
        }
    }
    if ((uint32_t)i < BA.n_bucket) {
        for (auto const & ru: BA.row_updates[i]) {
            // longhint_t h(0, ru.hint, ru.slice_index);
            /* have to walk the updates to see the ones that match */
            bucket_update_t<1, HINT> u = ru;
            for (int nx = ru.n + 1; nx--;) {
                if (UNLIKELY(S[u.x] != 255)) {
                    push_update(to_longhint(u, ru.slice_index));
                }
                u.x += ru.inc;
            }
        }
    }
}

template void bucket_array_complete::purge<shorthint_t>(
    bucket_array_t<1, shorthint_t> const & BA, int const i,
    unsigned char const * S);

template void bucket_array_complete::purge<longhint_t>(
    bucket_array_t<1, longhint_t> const & BA, int const i,
    unsigned char const * S);

#if defined(HAVE_SSE2)

template <typename HINT, int SIZE>
void bucket_array_complete::purge_1(
    const bucket_array_t<1, HINT> & BA, const int i,
    const std::vector<typename bucket_update_t<1, HINT>::br_index_t> &
        survivors)
{
    smallset<SIZE, typename bucket_update_t<1, HINT>::br_index_t> surv_set(
        survivors);

    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        slice_index_t const slice_index = BA.get_slice_index(i_slice);
        for (auto const & it: BA.slice_range(i, i_slice)) {
            if (UNLIKELY(surv_set.contains(it.x))) {
                push_update(to_longhint(it, slice_index));
            }
        }
    }
}

template <typename HINT>
void bucket_array_complete::purge(
    bucket_array_t<1, HINT> const & BA, int const i, unsigned char const * S,
    std::vector<typename bucket_update_t<1, HINT>::br_index_t> const &
        survivors)
{
    size_t const items_per_size =
        smallset<1, typename bucket_update_t<1, HINT>::br_index_t>::nr_items;
    if (survivors.size() == 0)
        return;
    switch ((survivors.size() + items_per_size - 1) / items_per_size) {
    case 1:
        purge_1<HINT, 1>(BA, i, survivors);
        break;
    case 2:
        purge_1<HINT, 2>(BA, i, survivors);
        break;
    case 3:
        purge_1<HINT, 3>(BA, i, survivors);
        break;
    case 4:
        purge_1<HINT, 4>(BA, i, survivors);
        break;
    case 5:
        purge_1<HINT, 5>(BA, i, survivors);
        break;
    case 6:
        purge_1<HINT, 6>(BA, i, survivors);
        break;
    case 7:
        purge_1<HINT, 7>(BA, i, survivors);
        break;
    case 8:
        purge_1<HINT, 8>(BA, i, survivors);
        break;
    default:
        purge(BA, i, S);
    }
}

template void bucket_array_complete::purge<shorthint_t>(
    bucket_array_t<1, shorthint_t> const &, int, unsigned char const *,
    std::vector<typename bucket_update_t<1, shorthint_t>::br_index_t> const &);
#endif

template <int INPUT_LEVEL>
void downsort(fb_factorbase::slicing const & fbs MAYBE_UNUSED,
              bucket_array_t<INPUT_LEVEL - 1, longhint_t> & BA_out,
              const bucket_array_t<INPUT_LEVEL, shorthint_t> & BA_in,
              uint32_t bucket_number, where_am_I & w)
{
    BA_out.add_slice_index(0);

    /* Time recording for this function is done by the caller
     * (downsort_wrapper)
     */
    /* Rather similar to purging, except it doesn't purge */
    int logB = LOG_BUCKET_REGIONS[INPUT_LEVEL - 1];
    using BA_out_t = bucket_array_t<INPUT_LEVEL - 1, longhint_t>;
    using lower_update_t = BA_out_t::update_t;
    using lower_row_update_t = BA_out_t::row_update_t;
    decltype(lower_update_t::x) maskB = (1 << logB) - 1;

    for (slice_index_t i_slice = 0; i_slice < BA_in.get_nr_slices();
         i_slice++) {
        slice_index_t const slice_index = BA_in.get_slice_index(i_slice);
        WHERE_AM_I_UPDATE(w, i, slice_index);
        for (auto const & it: BA_in.slice_range(bucket_number, i_slice)) {
            WHERE_AM_I_UPDATE(w, p, fbs[slice_index].get_prime(it.hint));
            WHERE_AM_I_UPDATE(w, h, it.hint);
            longhint_t h(0, it.hint, slice_index);
            lower_update_t u_low(it.x & maskB, h);
            BA_out.push_update(it.x >> logB, u_low, w);
        }
    }

    for (auto const & ru: BA_in.row_updates[bucket_number]) {
        WHERE_AM_I_UPDATE(w, i, ru.slice_index);
        WHERE_AM_I_UPDATE(w, p, fbs[ru.slice_index].get_prime(ru.hint));
        WHERE_AM_I_UPDATE(w, h, ru.hint);
        /* XXX
         * we have a problem here when
         * logB[INPUT_LEVEL] >= logI > logB[INPUT_LEVEL-1]
         * ru.n corresponds to a full line, but that doesn't fit in one
         * region at level INPUT_LEVEL-1.
         */

        /* This assumes that we push separate row updates for each
         * lowest-level bucket.
         */
        ASSERT(((ru.x + ru.n * ru.inc) >> logB) == (ru.x >> logB));
        longhint_t h(0, ru.hint, ru.slice_index);
        lower_update_t u_low(ru.x & maskB, h);
        lower_row_update_t ru_low(u_low, ru.slice_index, ru.inc, ru.n);
        BA_out.row_updates[ru.x >> logB].push_back(ru_low);

        /* I think that the approach below would work as well, but it's
         * not tested. It would also be more costly.
         */
#if 0
        auto x0 = ru.x;
        for(int nx = ru.n + 1 ; nx ; ) {
            auto x = x0;
            int n;
            for(n = 0 ; (x & maskB) == (x0 & maskB) ; x += ru.inc, nx--, n++);
            /* we can fit n updates in this sub-bucket */
            lower_update_t u_low(x0 & maskB, h);
            lower_row_update_t ru_low(u_low, ru.slice_index, ru.inc, ru.n);
            BA_out.row_updates[x0 >> logB].push_back(ru_low);
            x0 = x;
        }
#endif
    }
}

template <int INPUT_LEVEL>
void downsort(fb_factorbase::slicing const & fbs MAYBE_UNUSED /* unused */,
              bucket_array_t<INPUT_LEVEL - 1, longhint_t> & BA_out,
              bucket_array_t<INPUT_LEVEL, longhint_t> const & BA_in,
              uint32_t bucket_number, where_am_I & w)
{
    /* longhint updates don't write slice end pointers, so there must be
       exactly 1 slice per bucket */
    ASSERT_ALWAYS(BA_in.get_nr_slices() == 1);
    ASSERT(BA_in.get_slice_index(0) == 0);

    int logB = LOG_BUCKET_REGIONS[INPUT_LEVEL - 1];
    using BA_out_t = bucket_array_t<INPUT_LEVEL - 1, longhint_t>;
    using lower_update_t = BA_out_t::update_t;
    using lower_row_update_t = BA_out_t::row_update_t;
    decltype(lower_update_t::x) maskB = (1 << logB) - 1;

    for (auto const & it: BA_in.slice_range(bucket_number, 0)) {
        BA_out.push_update(it.x >> logB, lower_update_t(it.x & maskB, it), w);
    }

    for (auto const & ru: BA_in.row_updates[bucket_number]) {
        WHERE_AM_I_UPDATE(w, i, ru.slice_index);
        WHERE_AM_I_UPDATE(w, p, fbs[ru.slice_index].get_prime(ru.hint));
        WHERE_AM_I_UPDATE(w, h, ru.hint);
        lower_update_t u_low(ru.x & maskB, ru);
        lower_row_update_t ru_low(u_low, ru.slice_index, ru.inc, ru.n);
        BA_out.row_updates[ru.x >> logB].push_back(ru_low);
    }
}

template <int INPUT_LEVEL>
void downsort(fb_factorbase::slicing const & fbs,
              bucket_array_t<INPUT_LEVEL - 1, logphint_t> & BA_out,
              bucket_array_t<INPUT_LEVEL, emptyhint_t> const & BA_in,
              uint32_t bucket_number, where_am_I & w)
{
    BA_out.add_slice_index(0);

    /* Time recording for this function is done by the caller
     * (downsort_wrapper)
     */

    int logB = LOG_BUCKET_REGIONS[INPUT_LEVEL - 1];
    using BA_out_t = bucket_array_t<INPUT_LEVEL - 1, logphint_t>;
    using lower_update_t = BA_out_t::update_t;
    using lower_row_update_t = BA_out_t::row_update_t;
    decltype(lower_update_t::x) maskB = (1 << logB) - 1;

    /* Rather similar to purging, except it doesn't purge */
    for (slice_index_t i_slice = 0; i_slice < BA_in.get_nr_slices();
         i_slice++) {
        slice_index_t const slice_index = BA_in.get_slice_index(i_slice);
        // WHERE_AM_I_UPDATE(w, i, slice_index);
        for (auto const & it: BA_in.slice_range(bucket_number, i_slice)) {
            logphint_t h = fbs[slice_index].get_logp();
            lower_update_t u_low(it.x & maskB, h);
            BA_out.push_update(it.x >> logB, u_low, w);
        }
    }
    for (auto const & ru: BA_in.row_updates[bucket_number]) {
        WHERE_AM_I_UPDATE(w, i, ru.slice_index);
        // we have no slice offset here, so no p.
        // WHERE_AM_I_UPDATE(w, p, fbs[ru.slice_index].get_prime(ru.hint));
        // WHERE_AM_I_UPDATE(w, h, ru);
        logphint_t h = fbs[ru.slice_index].get_logp();
        lower_update_t u_low(ru.x & maskB, h);
        lower_row_update_t ru_low(u_low, ru.slice_index, ru.inc, ru.n);
        BA_out.row_updates[ru.x >> logB].push_back(ru_low);
    }
}

template <int INPUT_LEVEL>
void downsort(fb_factorbase::slicing const & fbs MAYBE_UNUSED /* unused */,
              bucket_array_t<INPUT_LEVEL - 1, logphint_t> & BA_out,
              bucket_array_t<INPUT_LEVEL, logphint_t> const & BA_in,
              uint32_t bucket_number, where_am_I & w)
{
    int logB = LOG_BUCKET_REGIONS[INPUT_LEVEL - 1];
    using BA_out_t = bucket_array_t<INPUT_LEVEL - 1, logphint_t>;
    using lower_update_t = BA_out_t::update_t;
    using lower_row_update_t = BA_out_t::row_update_t;
    decltype(lower_update_t::x) maskB = (1 << logB) - 1;

    /* logphint updates don't write slice end pointers, so there must be
       exactly 1 slice per bucket */
    ASSERT_ALWAYS(BA_in.get_nr_slices() == 1);
    ASSERT(BA_in.get_slice_index(0) == 0);

    for (auto const & it: BA_in.slice_range(bucket_number, 0)) {
        BA_out.push_update(it.x >> logB, lower_update_t(it.x & maskB, it), w);
    }
    for (auto const & ru: BA_in.row_updates[bucket_number]) {
        WHERE_AM_I_UPDATE(w, i, ru.slice_index);
        // no slice_offset in logphint either
        // WHERE_AM_I_UPDATE(w, p, fbs[ru.slice_index].get_prime(ru.hint));
        // WHERE_AM_I_UPDATE(w, h, ru);
        lower_update_t u_low(ru.x & maskB, ru);
        lower_row_update_t ru_low(u_low, ru.slice_index, ru.inc, ru.n);
        BA_out.row_updates[ru.x >> logB].push_back(ru_low);
    }
}

/* Explicitly instantiate the versions of downsort() that we'll need:
   downsorting shorthint from level 3 and level 2, and downsorting
   longhint from level 2.
   (for the moment we have no way to expose the level-3 case, so unless
   we find one, let's just drop this code)
 */

#if MAX_TOPLEVEL >= 2
template void downsort<2>(fb_factorbase::slicing const &,
                          bucket_array_t<1, longhint_t> & BA_out,
                          bucket_array_t<2, shorthint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & w);

template void downsort<2>(fb_factorbase::slicing const &,
                          bucket_array_t<1, logphint_t> & BA_out,
                          bucket_array_t<2, emptyhint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & w);
#endif

#if MAX_TOPLEVEL >= 3
template void downsort<3>(fb_factorbase::slicing const &,
                          bucket_array_t<2, longhint_t> & BA_out,
                          bucket_array_t<3, shorthint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & wr);

template void downsort<3>(fb_factorbase::slicing const &,
                          bucket_array_t<2, logphint_t> & BA_out,
                          bucket_array_t<3, emptyhint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & wr);
template void downsort<2>(fb_factorbase::slicing const &,
                          bucket_array_t<1, logphint_t> & BA_out,
                          bucket_array_t<2, logphint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & w);

template void downsort<2>(fb_factorbase::slicing const &,
                          bucket_array_t<1, longhint_t> & BA_out,
                          bucket_array_t<2, longhint_t> const & BA_in,
                          uint32_t bucket_number, where_am_I & w);

#endif
static_assert(MAX_TOPLEVEL == 3);

/* Instantiate concrete classes that we need or some methods do not get
   compiled and cause "undefined reference" errors during linking. */
template class bucket_single<1, primehint_t>;
template class bucket_single<1, longhint_t>;

template class bucket_array_t<1, shorthint_t>;
template class bucket_array_t<1, emptyhint_t>;

#if MAX_TOPLEVEL >= 2
template class bucket_array_t<2, shorthint_t>;
template class bucket_array_t<2, emptyhint_t>;
template class bucket_array_t<1, longhint_t>;
template class bucket_array_t<1, logphint_t>;
#endif

#if MAX_TOPLEVEL >= 3
template class bucket_array_t<3, shorthint_t>;
template class bucket_array_t<3, emptyhint_t>;
template class bucket_array_t<2, longhint_t>;
template class bucket_array_t<2, logphint_t>;
#endif
static_assert(MAX_TOPLEVEL == 3);
