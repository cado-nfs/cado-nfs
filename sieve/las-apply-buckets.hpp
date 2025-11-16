#ifndef CADO_LAS_APPLY_BUCKETS_HPP
#define CADO_LAS_APPLY_BUCKETS_HPP

#include "cado_config.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include "fb-types.hpp"
#include "bucket.hpp"
#include "cado-endian.h"
#include "las-where-am-i.hpp" // for where_am_I, WHERE_AM_I_UPDATE
#include "macros.h"

/* {{{ apply_buckets */
template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
    void
    apply_one_update(unsigned char * const S,
                     bucket_update_t<1, HINT> const & u,
                     const unsigned char logp, where_am_I & w)
{
    WHERE_AM_I_UPDATE(w, h, u.hint_for_where_am_i());
    WHERE_AM_I_UPDATE(w, x, u.x);
    sieve_increase(S + u.x, logp, w);
}

/* 
 * downsorted buckets might come from several parts (part 2 and
 * above only, though). Depending on the type of hint we're considering,
 * the logp lookup might call for a lookup on part 1 only, or possibly on
 * several parts.
 */
template<typename HINT> struct fb_part_lookup {};
template<typename HINT>
requires HINT::allowed_at_toplevel
struct fb_part_lookup<HINT> {
    static size_t min(fb_factorbase::slicing const &) { return 1; }
    static size_t max(fb_factorbase::slicing const &) { return 2; }
};
template<typename HINT>
requires (!HINT::allowed_at_toplevel)
struct fb_part_lookup<HINT> {
    static size_t min(fb_factorbase::slicing const &) { return 2; }
    static size_t max(fb_factorbase::slicing const & fbs) { return fbs.nparts(); }
};

/* TODO: for HINT=shorthint, fb_part_lookup will really be a proxy to
 * fbs.get_part(1), so having it inside the loop is wasteful
 */
template <int LEVEL, typename HINT>
inline void apply_row_updates_for_one_bucket(
    unsigned char * S, bucket_array_t<LEVEL, HINT> const & BA, int const i,
    fb_factorbase::slicing const & fbs, where_am_I & w)
{
    if ((uint32_t)i >= BA.n_bucket)
        return;

    /* p is logged at push_update, not at apply_update */
    WHERE_AM_I_UPDATE(w, p, 0);
    for (auto & ru: BA.row_updates[i]) {
        bucket_update_t<LEVEL, HINT> u = ru;
        auto const * fb_slice = fbs.get(ru.slice_index,
                fb_part_lookup<HINT>::min(fbs),
                fb_part_lookup<HINT>::max(fbs));
        ASSERT_ALWAYS(fb_slice != nullptr);
        unsigned char const logp = fb_slice->get_logp();
        // WHERE_AM_I_UPDATE(w, p, fbp[ru.slice_index].get_prime(u.hint));
        for (size_t n = ru.n + 1; n--;) {
            apply_one_update<HINT>(S, u, logp, w);
            u.x += ru.inc;
        }
    }
}

template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
    void
    apply_one_bucket(unsigned char * S, bucket_array_t<1, HINT> const & BA,
                     const int i, fb_factorbase::slicing const & fbs,
                     where_am_I & w)
{
    WHERE_AM_I_UPDATE(w, p, 0);

    /* We only want 1s buckets here. */
    static_assert(HINT::allowed_at_toplevel);

    /* This code is only used for 1s buckets. We do not have any
     * downsorted buckets here, so all primes must come from part 1
     * of the factor base.
     */
    auto const & fbp = fbs.get_part(1);

    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        auto sl = BA.slice_range(i, i_slice);
        auto it = sl.begin();
        auto it_end = sl.end();

        slice_index_t const slice_index = BA.get_slice_index(i_slice);

        unsigned char const logp = fbp[slice_index].get_logp();

        /* TODO: the code below is quite possibly correct, except perhaps for
         * the treatment of where_am_I stuff. I get inconsistent reports, esp
         * when I vary the number of threads.
         */
#if 1
        const bucket_update_t<1, HINT> * next_align;
        if (sizeof(bucket_update_t<1, HINT>) == 4) {
            next_align = (bucket_update_t<1, HINT> *)(((size_t)it + 0x3F) &
                                                      ~((size_t)0x3F));
            if (UNLIKELY(next_align > it_end))
                next_align = it_end;
        } else {
            next_align = it_end;
        }

        while (it != next_align)
            apply_one_update<HINT>(S, *it++, logp, w);

        while (it + 16 <= it_end) {
            uint64_t x0, x1, x2, x3, x4, x5, x6, x7;
            uint16_t x;
#ifdef HAVE_SSE2
#if defined(__ICC) || defined(__INTEL_COMPILER)
            /* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=45414 */
            _mm_prefetch(((char const *)it) + 256, _MM_HINT_NTA);
#else
            _mm_prefetch(((unsigned char *)it) + 256, _MM_HINT_NTA);
#endif
#endif
            x0 = ((uint64_t *)it)[0];
            x1 = ((uint64_t *)it)[1];
            x2 = ((uint64_t *)it)[2];
            x3 = ((uint64_t *)it)[3];
            x4 = ((uint64_t *)it)[4];
            x5 = ((uint64_t *)it)[5];
            x6 = ((uint64_t *)it)[6];
            x7 = ((uint64_t *)it)[7];
            it += 16;
            __asm__ __volatile__(
                ""); /* To be sure all x? are read together in one operation */
#if defined(CADO_LITTLE_ENDIAN)
#define INSERT_2_VALUES(X)                                                     \
    do {                                                                       \
        (X) >>= 16;                                                            \
        x = (uint16_t)(X);                                                     \
        WHERE_AM_I_UPDATE(w, x, x);                                            \
        sieve_increase(S + x, logp, w);                                        \
        (X) >>= 32;                                                            \
        WHERE_AM_I_UPDATE(w, x, X);                                            \
        sieve_increase(S + (X), logp, w);                                      \
    } while (0);
#elif defined(CADO_BIG_ENDIAN)
#define INSERT_2_VALUES(X)                                                     \
    do {                                                                       \
        x = (uint16_t)(X);                                                     \
        WHERE_AM_I_UPDATE(w, x, x);                                            \
        sieve_increase(S + x, logp, w);                                        \
        (X) >>= 32;                                                            \
        x = (uint16_t)(X);                                                     \
        WHERE_AM_I_UPDATE(w, x, x);                                            \
        sieve_increase(S + x, logp, w);                                        \
    } while (0);
#else
#error "You must #include \"cado-endian.h\""
#endif
            INSERT_2_VALUES(x0);
            INSERT_2_VALUES(x1);
            INSERT_2_VALUES(x2);
            INSERT_2_VALUES(x3);
            INSERT_2_VALUES(x4);
            INSERT_2_VALUES(x5);
            INSERT_2_VALUES(x6);
            INSERT_2_VALUES(x7);
        }
#endif
        while (it != it_end)
            apply_one_update<HINT>(S, *it++, logp, w);
    }
    apply_row_updates_for_one_bucket<1, HINT>(S, BA, i, fbs, w);
}

// Create the four instances, longhint_t and logphint_t are specialized.
template void apply_one_bucket<shorthint_t>(
    unsigned char * S, bucket_array_t<1, shorthint_t> const & BA, int const i,
    fb_factorbase::slicing const & fbs, where_am_I & w);

template void apply_one_bucket<emptyhint_t>(
    unsigned char * S, bucket_array_t<1, emptyhint_t> const & BA, int const i,
    fb_factorbase::slicing const & fbs, where_am_I & w);

template <>
void apply_one_bucket<longhint_t>(unsigned char * S,
                                  bucket_array_t<1, longhint_t> const & BA,
                                  int const i,
                                  fb_factorbase::slicing const & fbs,
                                  where_am_I & w)
{
    /* We need to access the full slicing, not the part, because
     * downsorted buckets might come from several parts (part 2 and
     * above only, though).
     */
    WHERE_AM_I_UPDATE(w, p, 0);
    // There is only one fb_slice. Slice indices are embedded in the hints.
    ASSERT(BA.get_nr_slices() == 1);
    ASSERT(BA.get_slice_index(0) == 0); // std::numeric_limits<slice_index_t>::max());
    for (auto const & it: BA.slice_range(i, 0)) {
        slice_index_t slice_index = it.slice_index;
        auto const * fb_slice = fbs.get(slice_index, 2);
        ASSERT_ALWAYS(fb_slice != nullptr);
        unsigned char const logp = fb_slice->get_logp();
        apply_one_update<longhint_t>(S, it, logp, w);
    }
    apply_row_updates_for_one_bucket<1, longhint_t>(S, BA, i, fbs, w);
}

template <>
void apply_one_bucket<logphint_t>(unsigned char * S,
                                  bucket_array_t<1, logphint_t> const & BA,
                                  int const i,
                                  fb_factorbase::slicing const & fbs,
                                  where_am_I & w)
{
    WHERE_AM_I_UPDATE(w, p, 0);
    // There is only one fb_slice. logp's are embedded in the hints.
    ASSERT(BA.get_nr_slices() == 1);
    ASSERT(BA.get_slice_index(0) == 0); // std::numeric_limits<slice_index_t>::max());
    for (auto const & it: BA.slice_range(i, 0)) {
        apply_one_update<logphint_t>(S, it, it.logp, w);
    }
    apply_row_updates_for_one_bucket<1, logphint_t>(S, BA, i, fbs, w);
}
/* }}} */

#endif /* CADO_LAS_APPLY_BUCKETS_HPP */
