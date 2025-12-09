#ifndef CADO_LAS_CONFIG_HPP
#define CADO_LAS_CONFIG_HPP

#include "cado_config.h" // HAVE_SSE2

#include <cstddef>

/* factor base is split in parts 0, 1, ..., FB_MAX_PARTS-1.
 *
 * The toplevel can thus be at most FB_MAX_PARTS-1.
 *
 * level-0 is for small-sieved primes.
 * level-1 is for primes that are bucket sieved in one go, without
 * downsorting.
 * level-2 is for primes that are bucket sieved in a first phase, then
 * downsorted.
 * for level-3, downsorting is a bit more work, done in downsort_aux
 */
#define FB_MAX_PARTS 4

#define MAX_TOPLEVEL ((FB_MAX_PARTS) - 1)


#ifdef __cplusplus
extern "C" {
#endif

extern void set_LOG_BUCKET_REGION();

extern int LOG_BUCKET_REGION;
extern int LOG_BUCKET_REGION_step;

extern int LOG_BUCKET_REGIONS[FB_MAX_PARTS];

extern size_t BUCKET_REGION;
extern size_t BUCKET_REGIONS[FB_MAX_PARTS];

extern int NB_DEVIATIONS_BUCKET_REGIONS;

extern int las_production_mode;

#ifdef __cplusplus
}
#endif

#define DESCENT_DEFAULT_GRACE_TIME_RATIO 0.2 /* default value */

/* (Re-)define this to support larger q. This is almost mandatory for the
 * descent.
 *
 * Modify the flag here **ONLY** if you intend to produce a custom las
 * siever that supports large q. The "normal" use case for this flag is
 * the las_descent binary, which is anyway _compiled_ with -DDLP_DESCENT
 * -DSUPPORT_LARGE_Q
 *
 * Tinkering with this definition is also the way to go if you want a
 * las_tracek binary that works with large q's similar to the ones we
 * have in the descent.
 */

#ifndef SUPPORT_LARGE_Q
#define xxxSUPPORT_LARGE_Q
#endif

#ifndef BUCKET_SIEVE_POWERS
#define BUCKET_SIEVE_POWERS
#endif

/* Guard for the logarithms of norms, so that the value does not wrap around
   zero due to roundoff errors. */
#define LOGNORM_GUARD_BITS 1

/* See PROFILE flag above */
/* Some functions should not be inlined when we profile or it's hard or
   impossible to tell them apart from the rest in the profiler output */
#ifdef PROFILE
#define NOPROFILE_INLINE
#define NOPROFILE_STATIC
#else
#define NOPROFILE_INLINE static inline
#define NOPROFILE_STATIC static
#endif

/* A memset with less than MEMSET_MIN bytes is slower than a fixed
 * memset (which is inlined with special code). So, if possible, it is
 * worthwhile to write:
 *   if (LIKELY(ts <= MEMSET_MIN)) memset (S, i, MEMSET_MIN); else memset (S, i,
 * ts); Some write-ahead allocation has to be provided for at malloc() time.
 */
#define MEMSET_MIN 64

#ifdef __cplusplus
extern "C" {
#endif
void las_display_config_flags();
#ifdef __cplusplus
}
#endif

#endif /* CADO_LAS_CONFIG_HPP */
