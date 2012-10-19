#ifndef LAS_ARITH_H_
#define LAS_ARITH_H_

#include <stdint.h>

#include "fb.h"
#include "utils.h"
#include "las-config.h"

#ifdef __cplusplus
extern "C" {
#endif

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
    uint32_t t = x*invp;                            /* t = x * invp mod 2^32 */
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32; /* x + t*p is bounded
            by 2^32*p-1+(2^32-1)*p < 2*2^32*p: we might want p < 2^31 so that
            there is no overflow */
    if (y >= p)
        y -= p;
    ASSERT(y<p);
    return y;
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
    uint32_t t = ((uint64_t)x)*invp;
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32;
    // might be too large by p, or too small by p.
    if ((int32_t)y < 0)
        y += p;
    else if (y >= p)
        y -= p;
#ifndef NDEBUG
    if (UNLIKELY(y >= p)) {
        fprintf(stderr, "BUG in redc_32. x = %" PRId64
                " p = %u, invp = %u, y = %d\n", x, p, invp, y);
        if (x < 0) {
            fprintf(stderr, "x/2^32 = -%"PRId64, (-x)>>32);
            if (((-x)>>32) < p) {
                fprintf(stderr, ", within bounds\n");
            } else {
                fprintf(stderr, ", OUT OF BOUNDS\n");
            }
        } else {
            fprintf(stderr, "x/2^32 = -%"PRId64, x>>32);
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
    return y;
}

MAYBE_UNUSED
static inline fbprime_t
invmod_po2 (fbprime_t n)
{
  fbprime_t r;

  ASSERT (n % 2 != 0);

  r = (3 * n) ^ 2;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  if (sizeof (fbprime_t) > 4)
    r = 2 * r - r * r * n;

  return r;
}


NOPROFILE_INLINE int
invmod (uint64_t *pa, uint64_t b)
{
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


/* Define to 0 or 1. Table lookup seems to be slightly faster than
   ctzl() on Opteron and slightly slower on Core2, but doesn't really
   make much difference in either case. */
#define LOOKUP_TRAILING_ZEROS 1

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// a is modified in place
// return 1 on succes, 0 on failure
NOPROFILE_INLINE int
invmod_REDC(uint64_t *pa, uint64_t b) {
#if LOOKUP_TRAILING_ZEROS
  static const unsigned char trailing_zeros[256] =
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
#endif
  uint64_t a, u, v, fix, p = b;
  int t, lsh;
  a = *pa;

  if (UNLIKELY(*pa == 0))
    return 0; /* or we get infinite loop */

  if (UNLIKELY(b % 2UL == 0))
    return invmod(pa, b);

  fix = (b+1)>>1;

  ASSERT (a < b);

  u = 1; v = 0; t = 0;

  // make a odd
  lsh = ctzl(a);
  a >>= lsh;
  t += lsh;
  /* v <<= lsh; ??? v is 0 here */

  // Here a and b are odd, and a < b
  do {
    do {
      b -= a; v += u;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b];
	b >>= lsh;
	t += lsh;
	u <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(b);
      b >>= lsh;
      t += lsh;
      u <<= lsh;
#endif
    } while (a<b);
    if (UNLIKELY(a == b))
      break;
    do {
      a -= b; u += v;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a];
	a >>= lsh;
	t += lsh;
	v <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
#endif
    } while (b < a);
  } while (a != b);
  if (a != 1)
    return 0;

  // Here, the inverse of a is u/2^t mod b.
  while (t>32) {
    unsigned long sig = u & 1UL;
    u >>= 1;
    if (sig)
      u += fix;
    --t;
  }
  while (t < 32)
    {
      u <<= 1;
      if (u >= p)
        u -= p;
      t ++;
    }
  *pa = u;
  return 1;
}

fbprime_t is_prime_power(fbprime_t q);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_ARITH_H_ */
