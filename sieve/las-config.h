#ifndef LAS_CONFIG_H_
#define LAS_CONFIG_H_

#include "cado.h"

#ifdef HAVE_SSE2
#define SSE_NORM_INIT
#endif

/* Number of bits used to estimate the norms
 * This should be large enough: it must be such that all norms are
 * smaller than 2^NORM_BITS.
 * This imposes NORM_BITS >= 8, or even >= 9 for large factorizations. */
#define NORM_BITS 10

/* Smart norm computation.
   These are two initializations of the algebraics/rationnals.
   These 2 initializations compute F(i, const j)=f(i) for each line.
   f(i)=log2(abs(sum[k=0...d] Ak i^k)+1.)*scale+GUARD = [GUARD...254],
   where Ak are the coefficients of the polynomial.

   The classical initialization is slow; the only optimization is done
   by the fast computation of the log2. The associated error guarantees
   the precision of the results +/- 1.
   This initialization is done if SMART_NORM is undefined.

   The smart initialization is very fast, and done if SMART_NORM is defined.
   This initialization computes first the roots of F, F', F" which each defined
   a line which intercepts (0, 0.).
   For each j, the roots of f, f' and f" are computed : there are F, F', F" roots * j.
   Because the absolute value, the roots of f have an unstable neighbourhood :
   f "bounces" on the horizontal axis.
   The roots of f" have also an unstable neighbourhood (inflexion points of f).
   The roots of f have a stable neighbourhood.
   So, the neighbourhoods of f(root(f)) and f(root(f")) are computed until
   on each side of the root there are SMART_NORM_STABILITY identical values,
   so until f has a local horizontal stability on the left and on the right of
   the root.
   A root has a maximal influence of SMART_NORM_INFLUENCE values on its
   neighbourhood (on each sides).
   
   These roots defined some segments of contiguous values of f(i). Some of them are
   reduced to a point (f(roots(f')); the lenght of the others is between 
   SMART_NORM_STABILITY * 2 + 1 and SMART_NORM_MAX_INFLUENCE * 2 + 1.

   3 artificials roots are inserted: -I/2 and (I/2)-1 as two roots of f' (two
   one-point segment), and 0.0 as a root of f, because near the neighbourhood of 0.
   f is very unstable.

   To compute all missing values of f(i) between ]-I/2,(I/2)-1[, a polygonal
   approximation is done between all these segments.
   The polygonal approximation between (i0,f(i0)-(i1,f(i1)) has 2 parameters :
   - The minimal lenght between i0 and i1, SMART_NORM_LENGTH. Below this
   value, the values f(i0...i1) are approximated by a line.
   - The maximal acceptable distance between f((i0+i1)/2) and (f(i0)+f(i1))/2,
   SMART_NORM_DISTANCE. Above this value, (i0,f(i0)-(i1,f(i1)) is approximated
   by the polygonal approximations of (i0,f(i0)-((i0+i1)/2,(f(i0)+f(i1))/2) and
   ((i0+i1)/2,(f(i0)+f(i1))/2)-(i1,f(i1)).

   The maximal error of the smart initialization depends on its four parameters,
   but should be less or egal to the double of SMART_NORM_DISTANCE.
   Modify the parameters of the smart norm algorithm is not a good idea. If you try
   it, use test_init_norm_bucket_region in the tests part in order to have an idea
   of the corresponding errors.
*/ 
#define SMART_NORM

#define SMART_NORM_STABILITY 3      /* Min:2; No max. Optimal: 3-4 */
#define SMART_NORM_INFLUENCE 10     /* Min:4; > SMART_NORM_STABILITY. No max. Optimal: 8-12 */
#define SMART_NORM_LENGTH    8      /* Min:3; no max. Optimal: 8-16 */
#define SMART_NORM_DISTANCE  1.     /* Min:1.; no max. Optimal: 1.0-1.5; NB: float. */

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
 * (see las-coordinates.c for the description of these)
 */
#define xxxTRACE_K

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
 * It is possible to put a higher value, in order to set I > 16.
 * However, this will have a bad impact on the memory usage, and on
 * efficiency, due to worse memory access. See bucket.h .
*/
#ifndef LOG_BUCKET_REGION
#define LOG_BUCKET_REGION 16
#endif

#define NB_BUCKETS_2 256
#define NB_BUCKETS_3 256

// Portability of ULL suffix?
#define BUCKET_REGION (UINT64_C(1) << LOG_BUCKET_REGION)

#define BUCKET_REGION_1 BUCKET_REGION
#define BUCKET_REGION_2 NB_BUCKETS_2*BUCKET_REGION_1
#define BUCKET_REGION_3 NB_BUCKETS_3*BUCKET_REGION_2

#define BUCKET_REGIONS { 0, BUCKET_REGION_1, BUCKET_REGION_2, BUCKET_REGION_3 }

/* This is currently used to enable some code paths specific to the
 * descent. The mid-term plan is to remove this compile-time flag.
 */
#define xxxDLP_DESCENT

/* used for a temporary attempt */
#define DLP_MNFSQ

#define DESCENT_DEFAULT_GRACE_TIME_RATIO 0.2    /* default value */

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
// #define USE_CACHEBUFFER 1
#endif 

/* A special ultrafast memset for las. Independant of MEMSET_MIN.
   Only for x86 64. */
// #define LAS_MEMSET

#include <stdio.h>

void las_display_config_flags();

#endif	/* LAS_CONFIG_H_ */
