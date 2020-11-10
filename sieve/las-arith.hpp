#ifndef LAS_ARITH_HPP_
#define LAS_ARITH_HPP_

#include "las-config.h"   // for NOPROFILE_INLINE
#include "cado_config.h"  // for HAVE_GCC_STYLE_AMD64_INLINE_ASM

#include <cstdint>        // for uint32_t, uint64_t, uint8_t, int64_t
#include <cinttypes>      // for PRI* macros

#include "macros.h"       // for ASSERT, UNLIKELY, GNUC_VERSION_ATLEAST, MAY...
#include "fb-types.h"     // for fbprime_t
#include "misc.h"          // cado_ctz
#include "mod_ul.h"        // for modul_clear, modul_clearmod, modul_get_ul
#include "verbose.h"

/* This header file is also #include'd by tests/sieve/torture-redc.cpp,
 * which (as its name says) checks that redc_32 and redc_u32 hold to
 * their promises.
 */

/* If defined to any value, use asm() CMOV in redc_u?32() */
// #define LAS_ARITH_REDC_USE_ASM 1

/* If defined, uses the old invmod code for speed comparison */
// #define LAS_ARITH_HPP_OLD_INVMOD_REDC_32 1

/* define STAT to get statistics on the frequency of carries in redc_u32 and
   redc_32 */
// #define STAT

/** Unsigned Redc_32 based on 64-bit arithmetic
 *   \param [in] x unsigned 64-bit integer, we require x < 2^32 * p
 *   \param [in] p is an odd number < 2^32.
 *   \param [in] invp is -1/p mod 2^32.
 *   \return x/2^32 mod p as an integer in [0, p[
*/

static inline uint32_t // NO_INLINE
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
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
#if defined(LAS_ARITH_REDC_USE_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__("addq %[tp],%[xtp]\n" 
          "cmovc %[p], %[t]\n"
          : [xtp]"+r"(xtp), [t]"+r"(t) : [tp]"r"(tp), [p] "r" (p)
          : "cc");
#else
  /* This conditional branch is rarely taken unless p is close to 2^32.
     In sieving tests, the branch was a little faster than the asm CMOV.
     With special-q close to 2^32, this may change. */
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

  if (u >= p) t = p;
  return u - t;
}

/** Unsigned modular multiplication with REDC
 *   \param [in] a unsigned 32-bit integer
 *   \param [in] b unsigned 32-bit integer. We require a*b < 2^32*p
 *   \param [in] p is an odd number < 2^32.
 *   \param [in] invp is -1/p mod 2^32.
 *   \return (a*b)/2^32 mod p as an integer in [0, p[
*/

static inline uint32_t
mulmodredc_u32(const uint32_t a, const uint32_t b, const uint32_t p, const uint32_t invp)
{
  uint64_t x = (uint64_t) a * (uint64_t) b;
  return redc_u32(x, p, invp);
}

/** Signed redc_32 based on 64-bit arithmetic
 *   \param [in] x is some signed integer in ]-2^32*p, 2^32*p[ (fitting in int64_t)
 *   \param [in] p is an odd number, 0 < p < 2^32.
 *   \param [in] invp is -1/p mod 2^32.
 *   \return x/2^32 mod p as an integer in [0, p[
*/
static inline uint32_t // NO_INLINE
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
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
#if defined(LAS_ARITH_REDC_USE_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__("addq %[tp],%[xtp]\n" 
          "cmovc %[p], %[t]\n"
          : [xtp]"+r"(xtp), [t]"+r"(t) : [tp]"r"(tp), [p] "r" (p)
          : "cc");
#else
  /* Same thing as in redc_u32 */
  xtp += tp;
  if (xtp < tp) t = p;
#endif
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

/* TODO: this is a close cousin of modredcul_inv, but the latter does
 * 64-bit redc */

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// return result on succes (new a value), UINT32_MAX on failure

#ifndef LAS_ARITH_HPP_OLD_INVMOD_REDC_32

static inline uint32_t // NO_INLINE
invmod_redc_32(uint32_t a, uint32_t b) {
  ASSERT (a < b);
  if (UNLIKELY(!a)) return a; /* or we get infinite loop */
  if (UNLIKELY(!(b & 1))) {
    uint32_t pa = a;
    invmod_32(&pa, (uint32_t) b);
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
    uint32_t v2 = v << lsh;
    uint32_t diff2 = (a - b) >> lsh;
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
  const uint32_t fix = (p+1)>>1;
  
  // TODO: can we use variable-width REDC for the division by 2^t here?
  // Here, the inverse of a is u/2^t mod b.
#define T3 do { uint8_t sig = (uint8_t) u; u >>= 1; if (sig & 1) u += fix; } while (0)
#define T4 do { u <<= 1; if (u >= p) u -= p; } while (0)
#if 1 /* Original code */
  if (t > 32)
    do T3; while (--t > 32);
  else
    while (t++ < 32) T4;
#else
  /* Duff's device (cf. Wikipedia) */
  t -= 32;
  if (LIKELY(t)) {
    if (LIKELY((int8_t) t > 0)) {
      uint8_t n = (t + 7) >> 3;
      switch (t & 7) {
          case 0: do { T3; no_break();
                      case 7: T3; no_break();
                      case 6: T3; no_break();
                      case 5: T3; no_break();
                      case 4: T3; no_break();
                      case 3: T3; no_break();
                      case 2: T3; no_break();
                      case 1: T3;
                  } while (--n > 0);
      }
    } else {
      uint8_t n = ((t = -t) + 7) >> 3;
      switch (t & 7) {
            case 0: do { T4; no_break();
                        case 7: T4; no_break();
                        case 6: T4; no_break();
                        case 5: T4; no_break();
                        case 4: T4; no_break();
                        case 3: T4; no_break();
                        case 2: T4; no_break();
                        case 1: T4;
                    } while (--n > 0);
      }
    }
  }
#endif
#undef T3
#undef T4
  return u;
}

#else /* LAS_ARITH_HPP_OLD_INVMOD_REDC_32 */

NOPROFILE_INLINE uint32_t
invmod_redc_32(uint32_t a, uint32_t b) {

  ASSERT (a < b);
  if (UNLIKELY(!a)) return a; /* or we get infinite loop */
  if (UNLIKELY(!(b & 1))) {
    uint32_t pa = a;
    invmod_32(&pa, (uint32_t) b);
    return pa;
  }
  const uint32_t p = b;
  uint32_t u = 1, v = 0, lsh = cado_ctz(a);
  uint8_t t = lsh;
  // make a odd
  a >>= lsh;
  
  // Here a and b are odd, and a < b
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
#define T1 "sub %0,%1\n " /* tzcnt */ "rep; bsf %1,%5\n add %2,%3\n shr %%cl,%1\n add %%cl,%4\n shl %%cl,%2\n cmp %0,%1\n "
#define T2 "sub %1,%0\n " /* tzcnt */ "rep; bsf %0,%5\n add %3,%2\n shr %%cl,%0\n add %%cl,%4\n shl %%cl,%3\n cmp %1,%0\n "
  __asm__ ( ".balign 8\n 0:\n"						\
	    T1 " je 9f\n jb  1f\n"					\
	    T1 " je 9f\n jb  1f\n"					\
	    T1 " je 9f\n jb  1f\n"					\
	    T1 " je 9f\n jb  1f\n"					\
	    T1 " je 9f\n jae 0b\n"					\
	    ".balign 8\n 1:\n"						\
	    T2 " je 9f\n jb  0b\n"					\
	    T2 " je 9f\n jb  0b\n"					\
	    T2 " je 9f\n jb  0b\n"					\
	    T2 " je 9f\n jb  0b\n"					\
	    T2 " jb 0b\n ja  1b\n"					\
	    ".balign 8\n 2:\n"						\
	    "9: \n"							\
	    : "+r" (a), "+r" (b), "+r" (u), "+r" (v), "+r" (t), "+c" (lsh));
#else
#define T1 do { b-=a; lsh=cado_ctz(b); v+=u; b>>=lsh; t+=lsh; u<<=lsh; if (a==b) goto ok; } while (0)
#define T2 do { a-=b; lsh=cado_ctz(a); u+=v; a>>=lsh; t+=lsh; v<<=lsh; if (b==a) goto ok; } while (0)
  {
      for (;;) {
          do {
              T1; if (a > b) break;
              T1; if (a > b) break;
              T1; if (a > b) break;
              T1; if (a > b) break;
              T1;
          } while (a < b);
          do {
              T2; if (b > a) break;
              T2; if (b > a) break;
              T2; if (b > a) break;
              T2; if (b > a) break;
              T2;
          } while (b < a);
      }
    ok: ; /* Need something after the label */
  }
#endif
#undef T1
#undef T2

  if (UNLIKELY(a != 1)) return 0;
  const uint32_t fix = (p+1)>>1;
  
  // Here, the inverse of a is u/2^t mod b.  
#define T3 do { uint8_t sig = (uint8_t) u; u >>= 1; if (sig & 1) u += fix; } while (0)
#define T4 do { u <<= 1; if (u >= p) u -= p; } while (0)
#if 1 /* Original code */
  if (t > 32)
    do T3; while (--t > 32);
  else
    while (t++ < 32) T4;
#else
  /* Duff's device (cf. Wikipedia) */
  t -= 32;
  if (LIKELY(t)) {
    if (LIKELY((int8_t) t > 0)) {
      uint8_t n = (t + 7) >> 3;
      switch (t & 7) {
          case 0: do { T3; no_break();
                      case 7: T3; no_break();
                      case 6: T3; no_break();
                      case 5: T3; no_break();
                      case 4: T3; no_break();
                      case 3: T3; no_break();
                      case 2: T3; no_break();
                      case 1: T3;
                  } while (--n > 0);
      }
    } else {
      uint8_t n = ((t = -t) + 7) >> 3;
      switch (t & 7) {
            case 0: do { T4; no_break();
                        case 7: T4; no_break();
                        case 6: T4; no_break();
                        case 5: T4; no_break();
                        case 4: T4; no_break();
                        case 3: T4; no_break();
                        case 2: T4; no_break();
                        case 1: T4;
                    } while (--n > 0);
      }
    }
  }
#endif
#undef T3
#undef T4
  return u;
}


#endif /* LAS_ARITH_HPP_OLD_INVMOD_REDC_32 */
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

static inline bool
batchinvredc_u32 (uint32_t *r, const uint32_t *a, const size_t n,
              const uint32_t p, const uint32_t invp)
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
    r[i] = mulmodredc_u32(r[i - 1], a[i], p, invp);
  }

  /* Here, $r_i = (\beta^{-i} P_i ) \pmod{p}$, for $0 \leq i < n$,
     with $\beta = 2^{32}$.
     I.e., $r_0 = a_0, r_1 = a_0 * a_1 * \beta^{-1}$, etc., */

  /* R := \beta * r_{n-1}^{-1} */
  R = invmod_redc_32(r[n - 1], p);
  if (R == 0 || R == UINT32_MAX)
    return false;

  /* R = P_{n-1}^{-1} * \beta^{n} mod p */

  for (size_t i = n - 1; i > 0; i--) {
    /* Invariant here: $R = P_i^{-1} * \beta^{i+1} \pmod{p}$ */
    r[i] = mulmodredc_u32(R, r[i-1], p, invp);
    /* $r_i := R * r_{i-1} * beta^{-1}
             = (P_i^{-1} * \beta^{i+1}) * (\beta^{-i+1} P_{i-1}) * \beta^{-1}
             = \beta * a_i^{-1}$. */
    R = mulmodredc_u32(R, a[i], p, invp);
    /* Here, $R = R * a_i * \beta^{-1} \pmod{p} = 
       P_{i-1}^{-1} * \beta^{i} \pmod{p}$. */
  }
  r[0] = R;
//#define LAS_ARITH_COMPARE_BATCHINV 1
#ifdef LAS_ARITH_COMPARE_BATCHINV
  for (size_t i = 0; i < n; i++) {
    uint32_t product = mulmodredc_u32(r[i], a[i], p, invp);
    if (product != 1) {
      verbose_output_print(1, 0, "batchinv_u32: 1/%" PRIu32 " (mod %" PRIu32
        ") wrong: %" PRIu32 "\n", a[i], p, r[i]);
      ASSERT(0);
    }
  }
#endif
  return true;
}


fbprime_t is_prime_power(fbprime_t q);

#endif	/* LAS_ARITH_HPP_ */
