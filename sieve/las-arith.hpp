#ifndef CADO_LAS_ARITH_HPP
#define CADO_LAS_ARITH_HPP

// #define LAS_ARITH_INVMOD_HISTOGRAM 1

#include "cado_config.h"
#include "las-config.hpp"

#include <cstdint>
#include <cinttypes>

#include "macros.h"
#include "fb-types.hpp"
#include "misc.h"
#include "arith/mod_ul.h"
#include "verbose.h"

#ifdef LAS_ARITH_INVMOD_HISTOGRAM
#include "utils_cxx.hpp"
#endif

/* This header file is also #include'd by tests/sieve/torture-redc.cpp,
 * which (as its name says) checks that redc_32 and redc_u32 hold to
 * their promises.
 */

/* define STAT to get statistics on the frequency of carries in redc_u32 and
   redc_32 */
// #define STAT

/** Returns true if p is small enough so that the first carry in our redc_*()
 *  functions cannot occur, and false otherwise.
 *
 *   If redc_no_carry(p) is true, then redc_*<false>() functions may be used
 *   which are slightly faster, otherwise redc_*<true>() functions must be used.
 *
 *   \param [in] p, the modulus for REDC to check
*/
static inline bool redc_no_carry(const uint32_t p) {
    return p < (UINT32_C(1) << 31);
}

/** Unsigned Redc_32 based on 64-bit arithmetic
 *   \param [in] x unsigned 64-bit integer, we require x < 2^32 * p
 *   \param [in] p is an odd number < 2^32. If CARRYCHECK is false, we require
 *               that p satisfies redc_no_carry(p).
 *   \param [in] invp is -1/p mod 2^32.
 *   \return x/2^32 mod p as an integer in [0, p[
*/

template <bool CARRYCHECK>
static inline uint32_t // NO_INLINE
redc_u32(const uint64_t x, const uint32_t p, const redc_invp_t invp)
{
  ASSERT_EXPENSIVE(CARRYCHECK || redc_no_carry(p));
  uint32_t t = (uint32_t) x * invp;   /* t = x * invp mod 2^32 */
  /* x + t*p is bounded by 2^32*p-1+(2^32-1)*p < 2*2^32*p
   * therefore, we must pay attention to the carry flag while doing the
   * addition.
   */
  uint64_t tp = (uint64_t)t * (uint64_t)p;
  uint64_t xtp = x;

  /* do xtp += tp. If the addition produces a carry, then set t = p,
   * and t = 0 otherwise */

  t = 0;
  if (CARRYCHECK) {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__("addq %[tp],%[xtp]\n" 
          "cmovc %[p], %[t]\n"
          : [xtp]"+r"(xtp), [t]"+r"(t) : [tp]"r"(tp), [p] "r" (p)
          : "cc");
#else
  /* This conditional branch is never taken unless p > 2^31.
     In sieving tests, the branch was a little faster than the asm CMOV.
     With factor base primes > 2^31, this may change. */
  xtp += tp;
  if (xtp < tp) t = p;
#endif
#ifdef STAT
  static int count = 0, carry = 0;
  count ++;
  carry += (xtp < tp) ? 1 : 0;
  if ((count % 1000000) == 0)
    printf ("redc_u32: count=%d carry=%d\n", count, carry);
#endif
  } else {
    xtp += tp;
  }
  uint32_t u = xtp >> 32;
  // by construction, xtp is divisible by 2^32. u is such that
  // 0 <= u < 2*p ; however the representative that we have is capped to
  // 2^32, and may wrap around.
  /* if cf is true, then with u' := 2^32 + u, we have 2^32 <= u' < 2*p,
   * thus u' - p < p, which ensures that (1) a borrow will occur in the
   * subtraction u - p, which will compensate for the carry in x + t*p,
   * and (2) the final result will be < p as wanted */

  /* Alternative: t2=u-t; u-=p; if carry: u=t2. 4 insn instead of 3, dependent
   * chain also has length 3 but with one MOV. Turned out to be slightly slower
   * on i3-6100. */

  /* With -O3, gcc turns this into a conditional branch. :( */
  if (u >= p) t = p;
  return u - t;
}

/** Variable-width unsigned REDC.
 *   \param [in] x unsigned 32-bit integer. We require x < p.
 *   \param [in] p is an odd number < 2^32.
 *   \param [in] invp is -1/p mod 2^32.
 *   \param [in] s is an integer 0 <= s < 32
 *   \return x/2^s mod p as an integer in [0, p[
*/
static inline uint32_t // NO_INLINE
varredc_u32(const uint32_t x, const uint32_t p, const redc_invp_t invp, const uint8_t s)
{
  ASSERT(s < 32); /* shift by word width is undefined. caller should use
                     redc_u32() instead. */
  ASSERT(x < p);
  const uint32_t mask = (UINT32_C(1) << s) - 1;
  uint32_t t = ((uint32_t) x * invp) & mask;   /* t = x * invp mod 2^s, 0 <= t <= 2^s-1 */
  uint64_t tp = (uint64_t)t * (uint64_t)p;     /* 0 <= tp <= (2^s-1)*p */
  uint64_t xtp = x + tp;                       /* 0 <= xtp <= 2^s*p - 1 <= 2^31*p - 1 < 2^64. No carry */

  ASSERT((xtp & mask) == 0);
  uint32_t u = xtp >> s;                       /* 0 <= u < p */
  return u;
}


/** Unsigned modular multiplication with REDC
 *   \param [in] a unsigned 32-bit integer
 *   \param [in] b unsigned 32-bit integer. We require a*b < 2^32*p
 *   \param [in] p is an odd number < 2^32.
 *   \param [in] invp is -1/p mod 2^32.
 *   \return (a*b)/2^32 mod p as an integer in [0, p[
*/

template <bool CARRYCHECK>
static inline uint32_t
mulmodredc_u32(const uint32_t a, const uint32_t b, const uint32_t p, const redc_invp_t invp)
{
  uint64_t x = (uint64_t) a * (uint64_t) b;
  return redc_u32<CARRYCHECK>(x, p, invp);
}

/** Signed redc_32 based on 64-bit arithmetic
 *   \param [in] x is some signed integer in ]-2^32*p, 2^32*p[ (fitting in int64_t)
 *   \param [in] p is an odd number < 2^32. If CARRYCHECK is false, we require
 *               that p satisfies redc_no_carry(p).
 *   \param [in] invp is -1/p mod 2^32.
 *   \return x/2^32 mod p as an integer in [0, p[
*/
template <bool CARRYCHECK>
static inline uint32_t // NO_INLINE
redc_32(const int64_t x, const uint32_t p, const redc_invp_t invp)
{
  ASSERT_EXPENSIVE(x > 0 && (x >> 32) < p || x < 0 && ((-x) >> 32) < p);
  ASSERT_EXPENSIVE(CARRYCHECK || redc_no_carry(p));
  uint32_t t = (uint32_t)x * invp;
  uint64_t tp = (uint64_t)t * (uint64_t)p;
#if 1
  /* gcc 6.4, 7.5, 8.3, 9.3, 10.2 all turn this into test/lea(or add)/cmovns,
   * which seems pretty good. On Skylake i3-6100 it makes root_transform()
   * 1-2% percent faster than the MUL below */
  uint64_t xtp = (x >= 0) ? x : x + ((uint64_t) p << 32);
#else
  uint64_t xtp = x + (uint64_t) p * (uint64_t) ((x & 0x8000000000000000) >> 31);
#endif
  /* now xtp >= 0, and we are in the same case as redc_u32,
     thus the same analysis as redc_u32 applies here */
  t = 0;
  if (CARRYCHECK) {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__("addq %[tp],%[xtp]\n" 
          "cmovc %[p], %[t]\n"
          : [xtp]"+r"(xtp), [t]"+r"(t) : [tp]"r"(tp), [p] "r" (p)
          : "cc");
#else
  /* Same thing as in redc_u32 */
  xtp += tp;
  if (xtp < tp) t = p;
#endif
  } else {
    xtp += tp;
  }
  /* Timings with revision 8de5fc9
     on Intel i5-4590 at 3.3Ghz with turbo-boost disabled:
     torture-redc 10000000:
     assembly: redc_32: 8388608 tests in 0.0668s
     C code  : redc_32: 8388608 tests in 0.0679s
     assembly+(if 0 instead of if 1): redc_32: 8388608 tests in 0.0680s */
  uint32_t u = xtp >> 32;
  if (u >= p) t = p;
  return u - t;
}

/* NOTE: we used to have redc_64, invmod_redc_64, and invmod_64 ; actually we
 * never really needed these functions, and they're gone since b9a4cbb40 ;
 * maybe temporarily.
 */

MAYBE_UNUSED
static inline fbprime_t
invmod_po2 (fbprime_t n)
{
  fbprime_t r;
  
  ASSERT (n & 1);
  r = (3 * n) ^ 2;
  r *= 2 - r * n;
  r *= 2 - r * n;
  r *= 2 - r * n;
  if (sizeof (fbprime_t) > 4)
    r *= 2 - r * n;
  return r;
}

/* put in pa[0] the value of 1/pa[0] mod b and return non-zero if the inverse
   exists, otherwise return 0 */
NOPROFILE_INLINE int
invmod_32 (uint32_t *pa, uint32_t b)
{
  ASSERT (sizeof(unsigned long) >= 4);
  modulusul_t m;
  residueul_t r;
  int rc;
  modul_initmod_ul (m, b);
  modul_init (r, m);
  modul_set_ul (r, *pa, m); /* With mod reduction */
  if ((rc = modul_inv(r, r, m)))
    *pa = modul_get_ul (r, m);
  modul_clear (r, m);
  modul_clearmod (m);
  return rc;
}

/* Requires a < m and b <= m, then r == a+b (mod m) and r < m */
static inline uint32_t
addmod_u32 (const uint32_t a, const uint32_t b, const uint32_t m)
{
#if (defined(__i386__) && defined(__GNUC__)) || defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  {
    uint32_t t = a + b, tr = a - m;

    __asm__ __VOLATILE (
      "add %2, %0\n\t"   /* tr += b */
      "cmovnc %1, %0\n\t"  /* if (!cy) tr = t */
      : "+&r" (tr)
      : "rm" (t), "g" (b)
      : "cc"
    );
    return tr;
  }
#else
  return (b >= m - a) ? (b - (m - a)) : (a + b);
#endif
}

/* TODO: this is a close cousin of modredcul_inv, but the latter does
 * 64-bit redc */

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// return result on succes (new a value), UINT32_MAX on failure

static inline uint32_t // NO_INLINE
invmod_redc_32(uint32_t a, const uint32_t orig_b, const uint32_t invb)
{
  uint32_t b = orig_b;
  ASSERT (a < b);
  if (UNLIKELY(!a)) return a; /* or we get infinite loop */
  if (UNLIKELY(!(b & 1))) {
    uint32_t pa = a;
    invmod_32(&pa, b);
    return pa;
  }
  const uint32_t p = b;
  uint32_t u = 1, v = 0, lsh = cado_ctz(a);
  uint8_t t = lsh;
  // make a odd
  a >>= lsh;
  
  // Here a and b are odd, and a < b
  while (true) {
    uint32_t diff1 = b - a;
    if (diff1 == 0)
      goto done;
    lsh = cado_ctz(diff1);
    diff1 >>= lsh;
    const uint32_t v2 = v << lsh;
    const uint32_t diff2 = (a - b) >> lsh;
    v += u;
    u <<= lsh;
    t += lsh;
#if 1 && (defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) || defined(__i386__) && defined(__GNUC__))
    __asm__(
      "cmp %[b], %[a]\n"          // Compare b, a
      "cmovb %[diff1], %[b]\n"    // if a<b:  b = diff1 (= (b - a) >> lsh)
      "cmovnb %[diff2], %[a]\n"   // if a>=b: a = diff2 (= (a - b) >> lsh)
      "cmovnb %[v], %[u]\n"       // if a>=b: u = uv (= u + v)
      "cmovnb %[v2], %[v]\n"      // if a>=b: v = v2 (= v << lsh)
      : [a] "+r" (a), [b] "+r" (b), [u] "+r" (u), [v] "+r" (v)
      : [diff1] "rm" (diff1), [diff2] "rm" (diff2), [v2] "rm" (v2)
      /* On 32-bit i386 there are not enough registers to use only "r" for
       * all the inputs, hence "rm". */
      : "cc"
    );
#else
    /* Unfortunately, gcc 7.5.0 and gcc 10.2.0 both produce a conditional
     * branch rather than a bunch of CMOVs. The branch is slower. */
    if (a < b) {
      b = diff1;
    } else {
      a = diff2;
      u = v;
      v = v2;
    }
#endif
  };
done:

  if (UNLIKELY(a != 1)) return 0;
  
#ifdef LAS_ARITH_INVMOD_HISTOGRAM
/* Gather and print histogram of t values */
  static StaticHistogram t_hist(256, "# Histogram of t values", true);
  static StaticHistogram calls(1, "# invmod_redc_32 calls");
  t_hist.inc(t);
  calls.inc();
#endif
  // Here, the inverse of a is u/2^t mod b. We want the result to be
  // u/2^32 mod b and divide or multiply by a power of 2 accordingly.
  if (t >= 32) {
    u = varredc_u32(u, orig_b, invb, t - 32);
  } else {
    while (t++ < 32) {
        u = addmod_u32(u, u, p);
    }
  }
#undef T3
#undef T4
  return u;
}


/** Compute r[i] = 2^32*a[i]^(-1) mod p, for 0 <= i < n.
 *
 * The a[i] are non-negative integers and r[i] are integers with 0 <= r[i] < p.
 * r and a must be non-overlapping.
 * \param r [out] The n inverses, if all exist, otherwise undefined
 * \param a [in] The n integers to invert. Need not be reduced modulo p
 * \param n [in] Number of integers to invert
 * \param p [in] The modulus for the inversion. Must be odd
 * \param invp [in] Must satisfy p*invp == -1 (mod 2^32)
 * \return true if all inverses exist, false otherwise returns
 */

template <bool CARRYCHECK>
static inline bool
batchinvredc_u32 (uint32_t *r, const uint32_t *a, const size_t n,
              const uint32_t p, const redc_invp_t invp)
{
  uint32_t R;

  /* For the sake of readability of the comments, let
     $P_i := \prod_{0 \leq j \leq i} a_j$, for $0 \leq i < n$.
     This is just a notational aid for the comments, the variables
     P_i do not occur in the code. */

  /* A weak test that a and r are non-overlapping. Sadly, C++11 leaves
   * ordering comparisons of pointers undefined unless they are both part
   * of the same array, which we can't guarantee here. */
  ASSERT_EXPENSIVE(a != r);

  if (n == 0)
    return 1;

  /* Reduce a[0] % m, and store in r[0]. */
  r[0] = (a[0] >= p) ? a[0] % p : a[0];

  for (size_t i = 1; i < n; i++) {
    r[i] = mulmodredc_u32<CARRYCHECK>(r[i - 1], a[i], p, invp);
  }

  /* Here, $r_i = (\beta^{-i} P_i ) \pmod{p}$, for $0 \leq i < n$,
     with $\beta = 2^{32}$.
     I.e., $r_0 = a_0, r_1 = a_0 * a_1 * \beta^{-1}$, etc., */

  /* R := \beta * r_{n-1}^{-1} */
  R = invmod_redc_32(r[n - 1], p, invp);
  if (R == 0 || R == UINT32_MAX)
    return false;

  /* R = P_{n-1}^{-1} * \beta^{n} mod p */

  for (size_t i = n - 1; i > 0; i--) {
    /* Invariant here: $R = P_i^{-1} * \beta^{i+1} \pmod{p}$ */
    r[i] = mulmodredc_u32<CARRYCHECK>(R, r[i-1], p, invp);
    /* $r_i := R * r_{i-1} * beta^{-1}
             = (P_i^{-1} * \beta^{i+1}) * (\beta^{-i+1} P_{i-1}) * \beta^{-1}
             = \beta * a_i^{-1}$. */
    R = mulmodredc_u32<CARRYCHECK>(R, a[i], p, invp);
    /* Here, $R = R * a_i * \beta^{-1} \pmod{p} = 
       P_{i-1}^{-1} * \beta^{i} \pmod{p}$. */
  }
  r[0] = R;
//#define LAS_ARITH_COMPARE_BATCHINV 1
#ifdef LAS_ARITH_COMPARE_BATCHINV
  for (size_t i = 0; i < n; i++) {
    uint32_t product = mulmodredc_u32<CARRYCHECK>(r[i], a[i], p, invp);
    if (product != 1) {
      verbose_fmt_print(1, 0,
              "batchinv_u32: 1/{} (mod {}) wrong: {}\n",
              a[i], p, r[i]);
      ASSERT(0);
    }
  }
#endif
  return true;
}


fbprime_t is_prime_power(fbprime_t q);

#endif	/* CADO_LAS_ARITH_HPP */
