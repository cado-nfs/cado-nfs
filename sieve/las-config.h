#ifndef LAS_CONFIG_H_
#define LAS_CONFIG_H_

#include "cado.h"

#ifdef HAVE_SSE2
#define SSE_NORM_INIT
#endif

/* Thresholds of the kilo-buckets & mega-buckets.
 * When the number of buckets I*J/BUCKET_REGION = I*(I/2)/65536
 * is less to THRESHOLD_K_BUCKETS, an one pass sort algorithm is
 * used to computed the buckets : function fill_in_buckets.
 * When the number is >= THRESHOLD_K_BUCKETS and < THRESHOLD_M_BUCKETS,
 * a two passes sort algorithm is used : function fill_in_k_buckets.
 * When the number is >= THRESHOLD_M_BUCKETS, a three pass sort
 * algorithm is used : function fill_in_m_buckets.
 * These threshold are very critical for the performance.
 * The "good" values for THRESHOLD_K_BUCKETS are between 4 096 & 32 768;
 * the "good" values for THRESHOLD_M_BUCKETS are between 32 768 &
 * 1 048 576 (it's possible the 3 passes sort is always slower than
 * the 2 passes).
 * Mandatory: THRESHOLD_K_BUCKETS >= 16;
 * THRESHOLD_M_BUCKETS >= (THRESHOLD_K_BUCKETS) * 4.
*/
#ifndef THRESHOLD_K_BUCKETS
#define THRESHOLD_K_BUCKETS 32768   /* 512 */
#endif
#ifndef THRESHOLD_M_BUCKETS
#define THRESHOLD_M_BUCKETS 1048576 /* 131072 */
#endif

#if THRESHOLD_K_BUCKETS < 16
#error THRESHOLD_K_BUCKETS must be >= 16
#endif
#if THRESHOLD_M_BUCKETS < (THRESHOLD_K_BUCKETS * 4)
#error THRESHOLD_M_BUCKETS must be >= (THRESHOLD_K_BUCKETS * 4)
#endif

/* Number of bits used to estimate the norms
 * This should be large enough: it must be such that all norms are
 * smaller than 2^NORM_BITS.
 * This imposes NORM_BITS >= 8, or even >= 9 for large factorizations. */
#define NORM_BITS 10

/* Smart norm computation. Compute only one norm per n (max HORI_NORM_STRIDE)
 * on a line, and propagate up to VERT_NORM_STRIDE rows above it. 
 * These approximations speed-up the norm computation, but maybe put a little 
 * more pressure on the cofactorisation step, since the useless 
 * cofactorisations might be more frequent. 
 * See the result of test_init_norms_bucket_region for precision.
 * If SMART_NORM is not defined, an exact algorithm with a fast log2 is used;
 * the result is always exact -/+ 1.
 * HORI_NORM_STRIDE & VERT_NORM_STRIDE are useless if SMART_NORM is undefined.
 */ 
#define SMART_NORM 1
#define HORI_NORM_STRIDE 8
#define VERT_NORM_STRIDE 4 

/* define PROFILE to keep certain functions from being inlined, in order to
   make them show up on profiler output */
//#define PROFILE

/* (for debugging only) define TRACE_K, and exactly one of the TRACE_*
 * values to something non-zero, in order to get tracing information on a
 * particular relation.  In particular this traces the sieve array entry
 * corresponding to the relation. Upon startup, the three values below
 * are reconciled.
 *
 * (see also las-coordinates.c)
 */
#define xxxTRACE_K
// #define TRACE_AB { 2039914353344275UL,6656604L }
// #define TRACE_IJ
// #define TRACE_Nx { 0,1655 }

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

/* This triggers code which fills bucket in several passes, one for each
 * congruence class mod 2 (three such, the trivial one leading to
 * spurious reports). It's currently only part of the story, and at the
 * moment it is almost neutral in terms of efficiency (small slowdown
 * because of access pattern).  But it's the way to go if one wants to
 * support I=16. There are many other places where changes must be made.
 * This particular flag affects only the treatment of the ``bucket
 * sieved'' primes, not the pattern-sieved, or small-sieved.
 */
#define MOD2_CLASSES_BS 0       /* define to 0 or 1 */

/* un-sieving of locations where gcd(i,j)>1 instead of testing gcd for
 * each survivor. Appears slower than default. This code has always been
 * #ifdef'd out, but maybe can be improved enough to make it worthwhile
 */
#define xxxUNSIEVE_NOT_COPRIME  /* see las-unsieve.c */

/* Optimal bucket region: 2^16 = 64K == close to L1 size.
   MANDATORY. Don't change this value.
*/
#ifndef LOG_BUCKET_REGION
#define LOG_BUCKET_REGION 16
#endif

#if LOG_BUCKET_REGION != 16
#error "LOG_BUCKET_REGION must (mandatory!) be equal to 16."
#endif

/* This flag is necessary to support I=16. Otherwise it's a useless
 * burden
 */
#define SUPPORT_I16

/* This is currently used to enable some code paths specific to the
 * descent. The mid-term plan is to remove this compile-time flag.
 */
#define xxxDLP_DESCENT
#define DESCENT_GRACE_TIME_RATIO 0.4

/* Define this to support larger q. This is almost mandatory for the
 * descent. */
#define xxxSUPPORT_LARGE_Q

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

/* A memset with less MEMSET_MIN bytes is slower than an fixed memset
   (which is inlined with special code). So, if it's possible, the optimal
   memset is
   if (LIKELY(ts <= MEMSET_MIN)) memset (S, i, MEMSET_MIN); else memset (S, i, ts);
   S += ts;
   So, all S' malloc must be increased of MEMSET_MIN. */
#define MEMSET_MIN 64

/* Should we use a cache-line buffer when converting kilo-bucket updates to
   regular bucket updates? Requires SSE2 if enabled. */
#ifdef HAVE_SSE2
#define USE_CACHEBUFFER 1
#endif 

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void las_display_config_flags(FILE * stream);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_CONFIG_H_ */
