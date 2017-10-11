#ifndef LAS_CONFIG_H_
#define LAS_CONFIG_H_

#ifdef HAVE_SSE2
#define SSE_NORM_INIT
#endif

/* Number of bits used to estimate the norms with the old reference code.
 * Unused otherwise.
 * This should be large enough: it must be such that all norms are
 * smaller than 2^(2^NORM_BITS)
 * This imposes NORM_BITS >= 8, or even >= 9 for large factorizations. */
#define NORM_BITS 10

/* define PROFILE to keep certain functions from being inlined, in order to
   make them show up on profiler output */
//#define PROFILE

/* (for debugging only) define TRACE_K, to something non-zero,
 * in order to get tracing information on a
 * particular relation.  In particular this traces the sieve array entry
 * corresponding to the relation. Upon startup, the three values below
 * are reconciled.
 *
 * This activates new command lines arguments: -traceab, -traceij, -traceNx.
 * (see las-coordinates.cpp for the description of these)
 */
#ifndef TRACE_K
#define xxxTRACE_K
#endif

/* Define CHECK_UNDERFLOW to check for underflow when subtracting
   the rounded log(p) from sieve array locations */
//#define CHECK_UNDERFLOW

/* Define TRACK_CODE_PATH in order to have the where_am_I structures
 * propagate info on the current situation of the data being handled.
 * This more or less makes the variables global, in that every function
 * can then access the totality of the variables. But it's for debug and
 * inspection purposes only.
 *
 * Note that WANT_ASSERT_EXPENSIVE, a flag which exists in broader
 * context, augments the scope of the tracking here by performing a
 * divisibility test on each sieve update. This is obviously very
 * expensive, but provides nice checking.
 *
 * Another useful tool for debugging is the sieve-area checksums that get
 * printed with verbose output (-v) enabled.
 */
#define xxxTRACK_CODE_PATH
#define xxxWANT_ASSERT_EXPENSIVE

/* TRACE_K *requires* TRACK_CODE_PATH -- or it displays rubbish */
#if defined(TRACE_K) && !defined(TRACK_CODE_PATH)
#define TRACK_CODE_PATH
#endif

/* idem for CHECK_UNDERFLOW */
#if defined(CHECK_UNDERFLOW) && !defined(TRACK_CODE_PATH)
#define TRACK_CODE_PATH
#endif

/* un-sieving of locations where gcd(i,j)>1 instead of testing gcd for
 * each survivor. Appears slower than default. This code has always been
 * #ifdef'd out, but maybe can be improved enough to make it worthwhile
 */
#define xxxUNSIEVE_NOT_COPRIME  /* see las-unsieve.c */

/* Optimal bucket region: 2^16 = 64K == close to L1 size.
 * It is possible to put a higher value, in order to set I > 16.
 * However, this will have a bad impact on the memory usage, and on
 * efficiency, due to worse memory access. See bucket.h .
*/
#ifndef LOG_BUCKET_REGION
#define LOG_BUCKET_REGION 16
#endif

#define NB_BUCKETS_2 256
#define NB_BUCKETS_3 256

#define BUCKET_REGION (UINT64_C(1) << LOG_BUCKET_REGION)

#define BUCKET_REGION_1 BUCKET_REGION
#define BUCKET_REGION_2 NB_BUCKETS_2*BUCKET_REGION_1
#define BUCKET_REGION_3 NB_BUCKETS_3*BUCKET_REGION_2

#define BUCKET_REGIONS { 0, BUCKET_REGION_1, BUCKET_REGION_2, BUCKET_REGION_3 }

#define DESCENT_DEFAULT_GRACE_TIME_RATIO 0.2    /* default value */

/* (Re-)define this to support larger q. This is almost mandatory for the
 * descent. */
#ifndef SUPPORT_LARGE_Q
#define xxxSUPPORT_LARGE_Q
#endif

/* Define SKIP_GCD3 to skip updates where 3 divides gcd(i,j) in the
   bucket sieving phase. Slightly slower than not skipping them
   in single-thread mode, but might be useful for multi-threading,
   or when memory is tight */
// #define SKIP_GCD3

/* Guard for the logarithms of norms, so that the value does not wrap around
   zero due to roundoff errors. */
#define GUARD 1

/* GUARD+LOG_MAX should be as near as possible from 256, to get more accuracy
   in the norm computations, but not too much, otherwise a norm might be
   rounded to zero. */
#define LOG_MAX (255.9 - (double) GUARD)

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

/* A memset with less than MEMSET_MIN bytes is slower than an fixed
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

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
void las_display_config_flags();
#ifdef __cplusplus
}
#endif

#endif	/* LAS_CONFIG_H_ */
