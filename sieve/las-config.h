#ifndef LAS_CONFIG_H_
#define LAS_CONFIG_H_

#include "cado_config.h"  // HAVE_SSE2
#include <stddef.h>

/* un-sieving of locations where gcd(i,j)>1 instead of testing gcd for
 * each survivor. Appears slower than default. This code has always been
 * #ifdef'd out, but maybe can be improved enough to make it worthwhile
 */
#define xxxUNSIEVE_NOT_COPRIME  /* see las-unsieve.c */

#define FB_MAX_PARTS 4

#ifdef __cplusplus
extern "C" {
#endif

#define LOG_NB_BUCKETS_2 8
#define LOG_NB_BUCKETS_3 8

extern void set_LOG_BUCKET_REGION();

extern int LOG_BUCKET_REGION;
extern int LOG_BUCKET_REGIONS[FB_MAX_PARTS];

extern size_t BUCKET_REGION;
extern size_t BUCKET_REGIONS[FB_MAX_PARTS];

extern int NB_DEVIATIONS_BUCKET_REGIONS;

extern int las_production_mode;

#ifdef __cplusplus
}
#endif

#define DESCENT_DEFAULT_GRACE_TIME_RATIO 0.2    /* default value */

/* (Re-)define this to support larger q. This is almost mandatory for the
 * descent.
 *
 * Modify the flag here **ONLY** if you intend to produce a custom las
 * siever that supports large q. The "normal" use case for this flag is
 * the las_descent binary, which is anyway _compiled_ with -DDLP_DESCENT
 * -DSUPPORT_LARGE_Q
 */

#ifndef SUPPORT_LARGE_Q
#define xxxSUPPORT_LARGE_Q
#endif

/* See bug #30012. The fix actually has a very adverse effect in some
 * cases, and we do not want to make this regression stick in the master
 * branch until we merge
 * https://gitlab.inria.fr/cado-nfs/cado-nfs/-/merge_requests/29
 * (which, among other things, has the "right" fix for this bug).
 */

#ifndef FIX_30012
#define xxxFIX_30012
#endif

/* Define SKIP_GCD3 to skip updates where 3 divides gcd(i,j) in the
   bucket sieving phase. Slightly slower than not skipping them
   in single-thread mode, but might be useful for multi-threading,
   or when memory is tight */
// #define SKIP_GCD3

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
 *   if (LIKELY(ts <= MEMSET_MIN)) memset (S, i, MEMSET_MIN); else memset (S, i, ts);
 * Some write-ahead allocation has to be provided for at malloc() time.
 */
#define MEMSET_MIN 64

/* Should we use a cache-line buffer when converting kilo-bucket updates to
   regular bucket updates? Requires SSE2 if enabled. */
#ifdef HAVE_SSE2
// #define USE_CACHEBUFFER 1
#endif 

#ifdef __cplusplus
extern "C" {
#endif
void las_display_config_flags();
#ifdef __cplusplus
}
#endif

#endif	/* LAS_CONFIG_H_ */
