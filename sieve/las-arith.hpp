#ifndef LAS_ARITH_HPP_
#define LAS_ARITH_HPP_

#include <cstdint>
#include <cinttypes>

#include "fb-types.h"
#include "utils.h"
#include "las-config.h"
#include "utils/misc.h" /* cado_ctzl */

/* This header file is also #include'd by tests/sieve/torture-redc.cpp,
 * which (as its name says) checks that redc_32 and redc_u32 hold to
 * their promises.
 */

/* define STAT to get statistics on the frequency of carries in redc_u32 and
   redc_32 */
// #define STAT

// Redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd number < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some integer in [0, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t) x * invp;   /* t = x * invp mod 2^32 */
  /* x + t*p is bounded by 2^32*p-1+(2^32-1)*p < 2*2^32*p
   * therefore, we must pay attention to the carry flag while doing the
   * addition.
   */
  uint64_t tp = (uint64_t)t * (uint64_t)p;
  uint64_t xtp = x;
  int cf;
  /* do xtp += tp, get carry out in cf */
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && GNUC_VERSION_ATLEAST(6,0,0)
  asm("addq %[tp],%[xtp]\n" : [xtp]"+r"(xtp), "=@ccc"(cf) : [tp]"r"(tp));
#else
  /* With GCC 9.2.1 the following code is as fast as the above assembly code.
     Example on Intel i5-4590 at 3.3Ghz with turbo-boost disabled:
     torture-redc 10000000:
     assembly: redc_u32: 8388608 tests in 0.0760s
     C code  : redc_u32: 8388608 tests in 0.0743s
  */
  xtp += tp;
  cf = xtp < tp;
#endif
#ifdef STAT
  static int count = 0, carry = 0;
  count ++;
  carry += cf != 0;
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

  return (cf || u >= p) ? u - p : u;
}

// Signed redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd number < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some signed integer in ]-2^32*p, 2^32*p[ (fitting in int64_t)
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t)x * invp;
  uint64_t tp = (uint64_t)t * (uint64_t)p;
  uint64_t xtp = (x >= 0) ? x : x + ((uint64_t) p << 32);
  /* the following does the same without any branch, but seems slightly
     worse with GCC 9.2.1 */
  // uint64_t xtp = x + (uint64_t) p * (uint64_t) ((x & 0x8000000000000000) >> 31);
  /* now xtp >= 0, and we are in the same case as redc_u32,
     thus the same analysis as redc_u32 applies here */
  int cf;
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && GNUC_VERSION_ATLEAST(6,0,0)
  asm("addq %[tp],%[xtp]\n" : [xtp]"+r"(xtp), "=@ccc"(cf) : [tp]"r"(tp));
#else
  xtp += tp;
  cf = xtp < tp;
#endif
  uint32_t u = xtp >> 32;
  return (cf || u >= p) ? u - p : u;
}

/* NOTE: we used to have redc_64 and invmod_redc_64 ; actually we never
 * really needed these functions, and they're gone since b9a4cbb40 ;
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

// Compute in place 1/x mod p, return non-zero if modular inverse exists,
// for uint64_t.
// Fallback function for 32-bit archis.
static inline int
fallback_invmod_64(uint64_t *x, uint64_t p)
{
    mpz_t xx, pp;
    mpz_init(xx);
    mpz_init(pp);
    mpz_set_uint64(xx, *x);
    mpz_set_uint64(pp, p);
    int rc = mpz_invert(xx, xx, pp);
    *x = mpz_get_uint64(xx);
    mpz_clear(xx);
    mpz_clear(pp);
    return rc;
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

/* put in pa[0] the value of 1/pa[0] mod b and return non-zero if the inverse
   exists, otherwise return 0 */
NOPROFILE_INLINE int
invmod_64 (uint64_t *pa, uint64_t b)
{
#if LONG_BIT == 64
  ASSERT (sizeof(unsigned long) >= 8);
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
#else
  return fallback_invmod_64(pa, b);
#endif
}

/* TODO: this is a close cousin of modredcul_inv, but the latter does
 * 64-bit redc */

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// return result on succes (new a value), UINT32_MAX on failure
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

fbprime_t is_prime_power(fbprime_t q);

#endif	/* LAS_ARITH_HPP_ */
