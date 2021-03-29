#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#include "gcd.h"
#include "mod_ul.h"
#include "macros.h"   // for ASSERT
#include "misc.h" /* for cado_ctzl */

uint64_t
gcd_int64 (const int64_t a, const int64_t b)
{
  return gcd_uint64(safe_abs64(a), safe_abs64(b));
}

uint64_t
gcd_uint64 (uint64_t a, uint64_t b)
{
  uint64_t t;

  if (b == 0)
    return a;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a; /* 0 <= t < a */
      b = a;
      a = t; /* 0 <= a < b */
    }

  return b;
}

unsigned long
gcd_ul (unsigned long a, unsigned long b)
{
  unsigned long t;

  if (b == 0)
    return a;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}

/* return g=gcd(a,b), together with xa such that xa * a = g mod b
 *
 * *xa is in [0..b/g-1]
 *
 * *xa is a multiplicative inverse of a/g mod b/g
 *
 */
unsigned long
xgcd_ul (unsigned long * xa, unsigned long a, unsigned long b)
{
    if (!xa) return gcd_ul(a, b);

    /* b == 0 isn't well defined anyway */
    if (b == 0) {
        *xa = 1;
        return a;
    }
    if (b == a) {
        *xa = 0;
        return a;
    }

    unsigned long r0 = b;
    unsigned long r1 = a;
    unsigned long u0 = 0; /* -u0 * a = r0 mod b */
    unsigned long u1 = 1; /*  u1 * a = r1 mod b */
    // unsigned long v0 = 1; /* -u0 * a + v0 * b = r0 */
    // unsigned long v1 = 0; /*  u1 * a - v1 * b = r1 */
    int i = 0;
    for( ; r1 > 0 ; i++) {
        // this invariant works, but of course it can't be checked when
        // u0*a overflows !
        // ASSERT((i & 1) ? ((u0 * a - r0) % b == 0) : ((u0 * a + r0) % b == 0));
        // (bogus?) ASSERT(r0 * u0 <= b);
        // (bogus?) ASSERT(r1 * u1 <= b);
        unsigned long q = r0 / r1;
        unsigned long r2 = r0 - q * r1;
        unsigned long u2 = u0 + q * u1;
        /* Note that we have the following matrix identities:
         *
         * ( u1 v1 )   ( 0 1 )   ( u0 v0 )
         * ( u2 v2 ) = ( 1 q ) * ( u1 v1 )
         *
         * ( r1 )   ( 0  1 )   ( r0 )
         * ( r2 ) = ( 1 -q ) * ( r1 )
         */
        r0 = r1; r1 = r2;
        u0 = u1; u1 = u2;
    } 
    /* Let M=(v,u,V,U) be the product of the matrices (0,1,1,q) that
     * appear in the first matrix identity above (with the q's in reverse
     * order, since we do a product on the left in the iteration).
     *
     * We claim that
     *          0 <= u <= b/2g
     *          0 <= v <= a/2g
     *          U = a/g
     *          V = b/g
     *
     * Assume for a moment that a and b are coprime (g=1)
     *
     * if a>b, M is [0,1,1,0] times the matrix for (b,a), so we're left
     * with the case a<b
     *
     *
     * if i==1, a and b being coprime, we have a=1 and b = a*q.  Note
     * that we have b>=2.  M is [0,1,1=a,q=b], and the claim holds.
     *
     * if i is larger, let (a',b')=(r,a)=(b-aq,a). The sequence for
     * (a',b') has length i-1 and the recurrence applies. M' is
     * [x',y',b-aq,a] with 0 <= x' <= b'/2 and 0<=y'<=a'/2
     * M is therefore [y',x'+qy',a,b], and it is easy to check that the
     * claim holds
     *
     * Now if a and b are *not* coprime, the same holds by observing the
     * bounds given by algorithm with inputs a/g and b/g
     *
     * Note: the assumption that a and b are distinct, and non zero, is
     * essential.
     */
    // ASSERT(r0 * v1 == b);
    if ((i & 1) == 0 && u0) {
        u0 = b/r0 - u0;
    }
    ASSERT((unsigned long) u0 < b/r0);
    *xa = u0;
    return r0;
}

/* return 1/a mod p, assume 0 <= a < p */
unsigned long
invert_ul (unsigned long a, unsigned long p)
{
  modulusul_t q;
  residueul_t b;

  modul_initmod_ul (q, p);
  modul_init (b, q);
  ASSERT (a < p);
  modul_set_ul_reduced (b, a, q);
  modul_inv (b, b, q);
  a = modul_get_ul (b, q);
  modul_clear (b, q);
  modul_clearmod (q);
  return a;
}

/* Binary gcd; any input allowed. */
uint64_t
bin_gcd_int64_safe (const int64_t a, const int64_t b)
{
  uint64_t ua = safe_abs64(a), ub = safe_abs64(b);
  int s, t;

  if (ua == 0)
    return ub;

  if (ub == 0)
    return ua;

  s = cado_ctz64 (ua);
  t = cado_ctz64 (ub);
  ua >>= s;
  ub >>= t;
  if (t < s)
    s = t;
  /* Here ua, ub > 0 and both odd */

  while (1)
    {
      while (ua >= ub) {
        /* Here, ua >= ub > 0, and ua, ub are both odd */
        ua -= ub;
        if (ua == 0)
          return ub << s;
        ua >>= cado_ctz64 (ua);
      }
      while (ub >= ua) {
       /* Here, ub >= ua > 0, and ua, ub are both odd */
       ub -= ua;
       if (ub == 0)
         return ua << s;
       ub >>= cado_ctz64 (ub);
      }
    }
}

