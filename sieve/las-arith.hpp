#ifndef LAS_ARITH_HPP_
#define LAS_ARITH_HPP_

#include <stdint.h>
#include <inttypes.h>

#include "fb-types.h"
#include "utils.h"
#include "las-config.h"
#include "utils/misc.h" /* cado_ctzl */

// Redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32. FIXME: p < 2^31? (see below)
//   * invp is -1/p mod 2^32.
//   * x is some integer in [0, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t) x * invp;                            /* t = x * invp mod 2^32 */
  uint32_t u = (x + (uint64_t)t * (uint64_t)p) >> 32;
  /* x + t*p is bounded by 2^32*p-1+(2^32-1)*p < 2*2^32*p:
     we might want p < 2^31 so that there is no overflow */
  t = u - p;
  if ((int32_t) t >= 0) u = t;
  return u;
}

// Signed redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some signed integer in ]-2^32*p, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t)x * invp;
  uint32_t u = (x + (uint64_t)t * (uint64_t)p) >> 32;
  // might be too large by p, or too small by p.
  t = u;
  u += p;
  if ((int32_t) t >= 0) u = t;
  t -= p;
  if ((int32_t) t >= 0) u = t;
#ifndef NDEBUG
  if (UNLIKELY(u >= p)) {
    fprintf(stderr, "BUG in redc_32. x = %" PRId64
	    " p = %u, invp = %u, u = %d\n", x, p, invp, u);
    if (x < 0) {
      fprintf(stderr, "x/2^32 = -%" PRId64, (-x)>>32);
      if (((-x)>>32) < p) {
	fprintf(stderr, ", within bounds\n");
      } else {
	fprintf(stderr, ", OUT OF BOUNDS\n");
      }
    } else {
      fprintf(stderr, "x/2^32 = -%" PRId64, x>>32);
      if ((x>>32) < p) {
	fprintf(stderr, ", within bounds\n");
      } else {
	fprintf(stderr, ", OUT OF BOUNDS\n");
      }
    }
    ASSERT(0);
    /* TODO: Fall back to safer code ? */
  }
#endif
  return u;
}

#define HAVE_redc_64
/* This does a full mul, and should grok slightly larger bounds than the
 * version above. Presumably, as long as x fits within 64 bits, (63 bits,
 * with sign), things should be ok. (TODO: check).
 */
static inline uint64_t
redc_64(const int64_t x, const uint32_t p, const uint64_t invp)
{
#if LONG_BIT == 64
  int64_t t = ((uint64_t)x)*invp;
  uint64_t u;

  // First, compute:
  //    u = (p*t + x) / 2^64
  // This requires 128-bit word addition.
  /* Need the high part of a 64x64 mul */
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
			"    mulq    %[p]\n"
			"    addq    %[x], %%rax\n"
			"    adcq    $0, %%rdx\n"
			: "=&d" (u)
			: "a" (t), [x] "rm" (x), [p] "r" ((uint64_t) p));
#else
  ASSERT (sizeof(unsigned long) == 8);
  uint64_t rax, rdx;
  ularith_mul_ul_ul_2ul(&rax, &rdx, p, t);
  ularith_add_ul_2ul(&rax, &rdx, x);
  u = rdx;
#endif
  /* As per the early clobber on rdx, x can't be put in there.
   * Furthermore, since t goes to rax, x doesn't go there either. Thus
   * it is reasonable to assume that it is still unmodified after the
   * asm block */
  u-=x<0;     /* FIXME. Can I get around this ? */
  t = u;
  u += p;
  if ((int64_t) t >= 0) u = t;
  t -= p;
  if ((int64_t) t >= 0) u = t;
  return u;
#else // 32-bit support, via gmp
    uint64_t u;
    mpz_t xx, pp, invpp;
    mpz_t tt, uu;
    mpz_init(xx); mpz_init(pp); mpz_init(invpp);
    mpz_init(tt); mpz_init(uu);
    mpz_set_int64(xx, x);
    mpz_set_ui(pp, p);
    mpz_set_uint64(invpp, invp);

    mpz_mul(tt, xx, invpp);
    mpz_fdiv_r_2exp(tt, tt, 64); // mod 2^64, take non-negative remainder.
    mpz_mul(uu, pp, tt);
    mpz_add(uu, uu, xx);
    ASSERT(mpz_divisible_2exp_p(uu, 64));
    mpz_tdiv_q_2exp(uu, uu, 64);
    u = mpz_get_uint64(uu);
    while (u >= p) {
        u -= p;
    }
    mpz_clear(xx); mpz_clear(pp); mpz_clear(invpp);
    mpz_clear(tt); mpz_clear(uu);
    return u;
#endif
}


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

// Compute in place 1/x mod p, return if modular inverse exists,
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
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
#define T3 "shr %2\n    lea (%2,%3,1),%0\n      cmovc   %0,%2\n "
#define T4 "add %2,%2\n mov %2,%0\n sub %4,%2\n cmovl   %0,%2\n "
  __asm__ ( "    cmp $0x20, %1\n je 30f\n       jb 19f\n"		\
	    ""								\
	    "    sub $0x26, %1\n                je 16f\n   js  17f\n"	\
	    ".balign 8\n 10:\n" T3 T3 T3 T3 T3 T3			\
	    "    sub $0x6,  %1\n ja 10b\n       je 16f\n"		\
	    "17: cmp $0xfc, %1\n je 12f\n       jb 11f\n"		\
	    "    cmp $0xfe, %1\n je 14f\n       jb 13f\n   jmp 15f\n"	\
	    "16: " T3 "15: " T3 "14: " T3 "13: " T3 "12: " T3 "11: " T3 \
	    "    jmp 30f\n"						\
	    ""								\
	    "19: neg %1\n        add $0x1a,%1\n je 26f\n   js  27f\n"	\
	    ".balign 8\n 20:\n" T4 T4 T4 T4 T4 T4			\
	    "    sub $0x6,  %1\n ja 20b\n       je 26f\n"		\
	    "27: cmp $0xfc, %1\n je 22f\n       jb 21f\n"		\
	    "    cmp $0xfe, %1\n je 24f\n       jb 23f\n   jmp 25f\n"	\
	    "26: " T4 "25: " T4 "24: " T4 "23: " T4 "22: " T4 "21: " T4 \
	    ""								\
	    "30:\n"							\
	    : "=&r" (v), "+r" (t), "+r" (u) : "r" (fix), "r" (p));
#else
#define T3 do { uint8_t sig = (uint8_t) u; u >>= 1; if (sig & 1) u += fix; } while (0)
#define T4 do { u <<= 1; if (u >= p) u -= p; } while (0)
#if 0 /* Original code */
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
#endif
#undef T3
#undef T4
  return u;
}

/* Only used for together with redc_64 */
NOPROFILE_INLINE int
invmod_redc_64(uint64_t a, uint64_t b)
{
  ASSERT (a < b);
  if (UNLIKELY(!a)) return a; /* or we get infinite loop */
  if (UNLIKELY(!(b & 1))) {
    invmod_64(&a, b);
    return a;
  }
  const uint64_t p = b;
  uint64_t u = 1, v = 0, lsh = cado_ctz(a);
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
  const uint64_t fix = (p+1)>>1;
  
  // Here, the inverse of a is u/2^t mod b.  
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
#define T3 "shr %2\n    lea (%2,%3,1),%0\n      cmovc   %0,%2\n "
#define T4 "add %2,%2\n mov %2,%0\n sub %4,%2\n cmovl   %0,%2\n "
  __asm__ ( "    cmp $0x40, %1\n je 30f\n       jb 19f\n"		\
	    ""								\
	    "    sub $0x46, %1\n                je 16f\n   js  17f\n"	\
	    ".balign 8\n 10:\n" T3 T3 T3 T3 T3 T3			\
	    "    sub $0x6,  %1\n ja 10b\n       je 16f\n"		\
	    "17: cmp $0xfc, %1\n je 12f\n       jb 11f\n"		\
	    "    cmp $0xfe, %1\n je 14f\n       jb 13f\n   jmp 15f\n"	\
	    "16: " T3 "15: " T3 "14: " T3 "13: " T3 "12: " T3 "11: " T3 \
	    "    jmp 30f\n"						\
	    ""								\
	    "19: neg %1\n        add $0x3a,%1\n je 26f\n   js  27f\n"	\
	    ".balign 8\n 20:\n" T4 T4 T4 T4 T4 T4			\
	    "    sub $0x6,  %1\n ja 20b\n       je 26f\n"		\
	    "27: cmp $0xfc, %1\n je 22f\n       jb 21f\n"		\
	    "    cmp $0xfe, %1\n je 24f\n       jb 23f\n   jmp 25f\n"	\
	    "26: " T4 "25: " T4 "24: " T4 "23: " T4 "22: " T4 "21: " T4 \
	    ""								\
	    "30:\n"							\
	    : "=&r" (v), "+r" (t), "+r" (u) : "r" (fix), "r" (p));
#else
#define T3 do { uint8_t sig = (uint8_t) u; u >>= 1; if (sig & 1) u += fix; } while (0)
#define T4 do { u <<= 1; if (u >= p) u -= p; } while (0)
#if 0 /* Original code */
  if (t > 64)
    do T3; while (--t > 64);
  else
    while (t++ < 64) T4;
#else
  /* Duff's device (cf. Wikipedia) */
  t -= 64;
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
#endif
#undef T3
#undef T4
  return u;
}

fbprime_t is_prime_power(fbprime_t q);

#endif	/* LAS_ARITH_HPP_ */
