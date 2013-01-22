/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 3 unsigned longs 
   in size.
   Moduli must be odd and have the upper word non-zero. Residues are stored 
   in Montgomery form, reduction after multiplication is done with REDC. 
   Due to inlining, this file must be included in the caller's source code with 
   #include */

/* Naming convention: all function start with modredc3ul, for 
   MODulus REDC 3 Unsigned Longs minus X bits, followed by underscore, 
   functionality of function (add, mul, etc), and possibly underscore and 
   specification of what argument types the function takes (_ul, etc). */

#ifndef MODREDC_3UL_H
#define MODREDC_3UL_H

/**********************************************************************/
#include <assert.h>
#if defined(MODTRACE)
#include <stdio.h>
#endif
#include <limits.h>
#include "macros.h"
#include "ularith.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

/*********************************************************************/
/* Helper macros */

/* A macro for function renaming. All functions here start with 
   modredc3ul_ */
#define MODREDC3UL_RENAME(x) modredc3ul_##x

#define MODREDC3UL_SIZE 3
#define MODREDC3UL_MINBITS (2 * LONG_BIT + 1)
#define MODREDC3UL_MAXBITS (3 * LONG_BIT - 2) 

typedef unsigned long residueredc3ul_t[MODREDC3UL_SIZE];
typedef unsigned long modintredc3ul_t[MODREDC3UL_SIZE];
typedef struct { 
  modintredc3ul_t m;
  residueredc3ul_t one;
  unsigned long invm;
} __modulusredc3ul_t;
typedef __modulusredc3ul_t modulusredc3ul_t[1];


/* ==================== Functions used internally ==================== */

static inline int
modredc3ul_intlt (const modintredc3ul_t a, const modintredc3ul_t b);

static inline void
modredc3ul_add (residueredc3ul_t r, const residueredc3ul_t a, 
		 const residueredc3ul_t b, const modulusredc3ul_t m);
static inline void
modredc3ul_get_uls (modintredc3ul_t r, const residueredc3ul_t s, 
		     const modulusredc3ul_t m MAYBE_UNUSED);

#ifndef umul_ppmm
#define umul_ppmm(xh, xl, a, b) \
    ularith_mul_ul_ul_2ul(&(xl), &(xh), a, b);
#endif

#ifndef ctzl
#ifdef __GNUC__
#define ctzl(x)         __builtin_ctzl(x)
#else
static inline int clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if (GMP_LIMB_BITS == 64)
        if (x >> 32) { a += 32; x >>= 32; }
#endif  
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = GMP_LIMB_BITS - 2 - a + t[x];
        return res;
}
static inline int ctzl(unsigned long x)
{
	return GMP_LIMB_BITS - clzl(x & - x);
}
#endif
#endif

static unsigned long
add_2(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  unsigned long r, s, t, cy, cy1, cy2;
  cy = 0;

  r = x[0];
  s = r + y[0];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[0] = t;
  r = x[1];
  s = r + y[1];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[1] = t;
  return cy;
}

static unsigned long
add_3(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  unsigned long r, s, t, cy, cy1, cy2;
  cy = 0;

  r = x[0];
  s = r + y[0];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[0] = t;
  r = x[1];
  s = r + y[1];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[1] = t;
  r = x[2];
  s = r + y[2];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[2] = t;
  return cy;
}

static unsigned long
sub_3(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  unsigned long r, s, t, cy, cy1, cy2;
  cy = 0;

  r = x[0];
  s = r - y[0];
  cy1 = s > r;
  t = s - cy;
  cy2 = t > s;
  cy = cy1 | cy2;
  z[0] = t;
  r = x[1];
  s = r - y[1];
  cy1 = s > r;
  t = s - cy;
  cy2 = t > s;
  cy = cy1 | cy2;
  z[1] = t;
  r = x[2];
  s = r - y[2];
  cy1 = s > r;
  t = s - cy;
  cy2 = t > s;
  cy = cy1 | cy2;
  z[2] = t;
  return cy;
}

static inline void lshift_3(unsigned long *a, int cnt) {
  int i;
  int dnt = GMP_NUMB_BITS - cnt;
  if (cnt != 0) {
    for (i = 3-1; i>0; --i) {
      a[i] <<= cnt;
      a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
  }
}

static inline void long_lshift_3(unsigned long *a, int off, int cnt) {
  int i;
  int dnt = GMP_NUMB_BITS - cnt;
  assert (off > 0);
  if (cnt != 0) {
    for (i = 3-1; i>off; --i) {
      a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
    }
    a[off] = a[0]<<cnt;
    for (i = off-1; i>=0; --i) {
      a[i] = 0UL;
    }
  } else {
    for (i = 3-1; i >= off; --i)
      a[i] = a[i-off];
    for (i = off-1; i >= 0; --i)
      a[i] = 0;
  }
}


static inline void rshift_3(unsigned long *a, int cnt) {
  int i;
  int dnt = GMP_NUMB_BITS - cnt;
  if (cnt != 0) {
    for (i = 0; i < 3-1; ++i) {
      a[i] >>= cnt;
      a[i] |= (a[i+1] << dnt);
    }
    a[3-1] >>= cnt;
  }
}


static inline void long_rshift_3(unsigned long *a, int off, int cnt) {
  int i;
  int dnt = GMP_NUMB_BITS - cnt;
  assert (off > 0);
  if (cnt != 0) {
    for (i = 0; i < 3 - off - 1; ++i) {
      a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
    }
    a[3-off-1] = a[3-1]>>cnt;
    for (i = 3-off; i < 3; ++i) {
      a[i] = 0UL;
    }
  } else {
    for (i = 0; i < 3-off; ++i)
      a[i] = a[i+off];
    for (i = 3-off; i < 3; ++i)
      a[i] = 0;
  }
}

static unsigned long
addmul1_3(unsigned long *z, const unsigned long *x, const unsigned long c)
{
  unsigned long hi,lo,carry,buf;
  carry = 0;

  umul_ppmm(hi,lo,c,x[0]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[0];
  lo += buf;
  carry += (lo<buf);
  z[0] = lo;

  umul_ppmm(hi,lo,c,x[1]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[1];
  lo += buf;
  carry += (lo<buf);
  z[1] = lo;

  umul_ppmm(hi,lo,c,x[2]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[2];
  lo += buf;
  carry += (lo<buf);
  z[2] = lo;
  z[3] += carry;
  return (z[3]<carry);
}

static void
MAYBE_UNUSED
add_nc_1(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  unsigned long r, s, t, cy, cy1, cy2;
  cy = 0;

  r = x[0];
  s = r + y[0];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[0] = t;
}

static void
MAYBE_UNUSED
add_nc_2(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  unsigned long r, s, t, cy, cy1, cy2;
  cy = 0;

  r = x[0];
  s = r + y[0];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[0] = t;
  r = x[1];
  s = r + y[1];
  cy1 = s < r;
  t = s + cy;
  cy2 = t < s;
  cy = cy1 | cy2;
  z[1] = t;
}

static void
addmul1_nc_1(unsigned long *z, const unsigned long *x, const unsigned long c)
{
  unsigned long hi,lo,carry,buf;
  carry = 0;

  umul_ppmm(hi,lo,c,x[0]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[0];
  lo += buf;
  carry += (lo<buf);
  z[0] = lo;
  z[1] += carry;
}

static void
addmul1_nc_2(unsigned long *z, const unsigned long *x, const unsigned long c)
{
  unsigned long hi,lo,carry,buf;
  carry = 0;

  umul_ppmm(hi,lo,c,x[0]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[0];
  lo += buf;
  carry += (lo<buf);
  z[0] = lo;

  umul_ppmm(hi,lo,c,x[1]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[1];
  lo += buf;
  carry += (lo<buf);
  z[1] = lo;
  z[2] += carry;
}

static void
addmul1_nc_3(unsigned long *z, const unsigned long *x, const unsigned long c)
{
  unsigned long hi,lo,carry,buf;
  carry = 0;

  umul_ppmm(hi,lo,c,x[0]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[0];
  lo += buf;
  carry += (lo<buf);
  z[0] = lo;

  umul_ppmm(hi,lo,c,x[1]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[1];
  lo += buf;
  carry += (lo<buf);
  z[1] = lo;

  umul_ppmm(hi,lo,c,x[2]);
  lo += carry;
  carry = (lo<carry) + hi;
  buf = z[2];
  lo += buf;
  carry += (lo<buf);
  z[2] = lo;
  z[3] += carry;
}
    
static void
mul_3(unsigned long *z, const unsigned long *x, const unsigned long *y)
{
  int i;
  for (i = 0; i < 2*3; ++i) 
    z[i] = 0;
  addmul1_nc_3 (z+0, x, y[0]);
  addmul1_nc_3 (z+1, x, y[1]);
  addmul1_nc_3 (z+2, x, y[2]);
}

static void
sqr_3(unsigned long *z, const unsigned long *x)
{
  unsigned long buf[2*3];
  int i;

  for (i = 0; i < 2*3; ++i)
    buf[i] = 0;
  addmul1_nc_1(buf+1, x, x[1]);
  addmul1_nc_2(buf+2, x, x[2]);

  umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
  umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
  umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
  mpn_lshift(buf, buf, 2*3, 1);
  mpn_add_n(z, z, buf, 2*3);
}

static int
cmp_3(const unsigned long *x, const unsigned long *y)
{
  int i;
  for (i = 3-1; i >= 0; --i) {
    if (x[i] > y[i])
      return 1;
    if (x[i] < y[i])
      return -1;
  }
  return 0;
}

static void
redc_3(unsigned long *z, unsigned long *x, const unsigned long *mip, const unsigned long *p) {

  int i;
  unsigned long cy;
  for (i = 0; i < 3; ++i) {
    unsigned long t = x[i]*mip[0];
    cy = addmul1_3(x+i, p, t);
    assert (x[i] == 0);
    x[i] = cy;
  }
  cy = add_2(x+3+1, x+3+1, x);
  cy += x[2];
  if (cy || cmp_3(x+3, p) >= 0)
    sub_3(z, x+3, p);
  else
    for (i = 0; i < 3; ++i)
      z[i] = x[i+3];
}

static void
mod_3(unsigned long *z, const unsigned long *x, const unsigned long *p)
{
  int i;
  unsigned long q[3+1], r[3];
  ASSERT (p[3-1] != 0);
  mpn_tdiv_qr(q, r, 0, x, 2*3, p, 3);
  for (i = 0; i < 3; ++i)
    z[i] = r[i];
}

static void
MAYBE_UNUSED
mgy_encode_3(unsigned long *z, const unsigned long *x, const unsigned long *p) 
{
  unsigned long t[6];
  int i;
  for (i = 0; i < 3; ++i) {
    t[i] = 0;
    t[i+3] = x[i];
  }
  mod_3(z, t, p);
}

static void
MAYBE_UNUSED
mgy_decode_3(unsigned long *z, const unsigned long *x, const unsigned long *invR, const unsigned long *p) 
{
  unsigned long t[6];
  mul_3(t, x, invR);
  mod_3(z, t, p);
}

static int
invmod_3(mp_limb_t *res, const mp_limb_t *x, const mp_limb_t *p) {
  mp_limb_t u[3], v[3], a[3], b[3], fix[3];
  int i, t, lsh;

  u[0] = 1UL; v[0] = 0UL;
  a[0] = x[0]; b[0] = p[0];
  for (i=1; i < 3; ++i) {
    u[i] = 0UL; v[i] = 0UL;
    a[i] = x[i]; b[i] = p[i];
  }
  
  if (cmp_3(a, v) == 0) {
    res[0]=0; res[1]=0; res[2]=0;
    return 0;
  }

  add_3(fix, b, u);
  rshift_3(fix, 1);

  assert (cmp_3(a,b) < 0);

  t = 0;
  
  if (a[0] != 0) {
    lsh = ctzl(a[0]);
    rshift_3(a, lsh);
    t += lsh;
    lshift_3(v, lsh);
  } else { // rare...
//    fprintf(stderr, "XOURIG\n");
    i = 1;
    while (a[i] == 0)
      ++i;
    assert (i <= 3);
    lsh = ctzl(a[i]);
    long_rshift_3(a, i, lsh);
    t += lsh + i*GMP_NUMB_BITS;
    long_lshift_3(v, i, lsh);
  }

  do {
    do {
      sub_3(b, b, a);
      add_3(v, v, u);
      if (b[0] != 0) {
        lsh = ctzl(b[0]);
        rshift_3(b, lsh);
        t += lsh;
        lshift_3(u, lsh);
      } else {  // Should almost never occur.
 //       fprintf(stderr, "XOURIG\n");
        i = 1;
        while (b[i] == 0)
          ++i;
        assert (i <= 3);
        lsh = ctzl(b[i]);
        long_rshift_3(b, i, lsh);
        t += lsh + i*GMP_NUMB_BITS;
        long_lshift_3(u, i, lsh);
      }
    } while (cmp_3(a,b) < 0);
    if (cmp_3(a, b) == 0)
      break;
    do {
      sub_3(a, a, b);
      add_3(u, u, v);
      if (a[0] != 0) {
        lsh = ctzl(a[0]);
        rshift_3(a, lsh);
        t += lsh;
        lshift_3(v, lsh);
      } else { // rare...
//        fprintf(stderr, "XOURIG\n");
        i = 1;
        while (a[i] == 0)
          ++i;
        assert (i <= 3);
        lsh = ctzl(a[i]);
        long_rshift_3(a, i, lsh);
        t += lsh + i*GMP_NUMB_BITS;
        long_lshift_3(v, i, lsh);
      }
    } while (cmp_3(b,a)<0);
  } while (cmp_3(a,b) != 0);
  {
    int ok = 1;
    if (a[0] != 1)
      ok = 0;
    else {
      for (i = 1; i < 3; ++i) 
        if (a[1] != 0) 
	  ok = 0;
    }
    if (!ok) {
      for (i = 0; i < 3; ++i)
        res[i] = a[i];
      return 0;
    }
  }
  while (t>0) {
    mp_limb_t sig = u[0] & 1UL;
    rshift_3(u, 1);
    if (sig)
      add_3(u, u, fix);
    --t;
  }
  for (i = 0; i < 3; ++i) {
    res[i] = u[i];
  }
  return 1;
}






MAYBE_UNUSED
static inline void
modredc3ul_tomontgomery (residueredc3ul_t r, const residueredc3ul_t s,
			  const modulusredc3ul_t m)
{
  int i;

  r[0] = s[0];
  r[1] = s[1];
  r[2] = s[2];
  /* TODO FIXME: ridiculously slow */
  for (i = 0; i < 3 * LONG_BIT; i++)
    modredc3ul_add (r, r, r, m);
}

/* Converts s out of Montgomery form by dividing by 2^(3*LONG_BIT).
   Requires s < m. */
MAYBE_UNUSED
static inline void
modredc3ul_frommontgomery (residueredc3ul_t r, const residueredc3ul_t s,
			    const modulusredc3ul_t m)
{
  unsigned long t[6];
  t[0] = s[0]; t[1] = s[1]; t[2] = s[2];
  t[3] = 0; t[4] = 0; t[5] = 0;
  redc_3(r, t, &m[0].invm, m[0].m);
}


/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredc3ul_intset (modintredc3ul_t r, const modintredc3ul_t s)
{
  r[0] = s[0];
  r[1] = s[1];
  r[2] = s[2];
}

MAYBE_UNUSED
static inline void
modredc3ul_intset_ul (modintredc3ul_t r, const unsigned long s)
{
  r[0] = s;
  r[1] = 0UL;
  r[2] = 0UL;
}

MAYBE_UNUSED
static inline int
modredc3ul_intequal (const modintredc3ul_t a, const modintredc3ul_t b)
{
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}

MAYBE_UNUSED
static inline int
modredc3ul_intequal_ul (const modintredc3ul_t a, const unsigned long b)
{
  return (a[0] == b && a[1] == 0UL && a[2] == 0UL);
}

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc3ul_intlt (const modintredc3ul_t a, const modintredc3ul_t b)
{
    modintredc3ul_t t;

    modredc3ul_intset (t, a);
    return ularith_sub_3ul_3ul_cy (&(t[0]), &(t[1]), &(t[2]), b[0], b[1], b[2]);
}

MAYBE_UNUSED
static inline int
modredc3ul_intcmp (const modintredc3ul_t a, const modintredc3ul_t b)
{
  if (a[2] < b[2])
    return -1;
  if (a[2] > b[2])
    return 1;
  // here a[2] == b[2]
  if (a[1] < b[1])
    return -1;
  if (a[1] > b[1])
    return 1;
  // here a[1:2] == b[1:2]
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc3ul_intcmp_ul (const modintredc3ul_t a, const unsigned long b)
{
  if (a[2] > 0UL)
    return 1;
  if (a[1] > 0UL)
    return 1;
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc3ul_intfits_ul (const modintredc3ul_t a)
{
  return (a[2] == 0UL && a[1] == 0UL);
}

MAYBE_UNUSED
static inline void
modredc3ul_intadd (modintredc3ul_t r, const modintredc3ul_t a,
		    const modintredc3ul_t b)
{
  modintredc3ul_t t;
  modredc3ul_intset (t, a);
  ularith_add_3ul_3ul (&(t[0]), &(t[1]), &(t[2]), b[0], b[1], b[2]);
  modredc3ul_intset (r, t);
}

MAYBE_UNUSED
static inline void
modredc3ul_intsub (modintredc3ul_t r, const modintredc3ul_t a,
		    const modintredc3ul_t b)
{
  modintredc3ul_t t;
  modredc3ul_intset (t, a);
  ularith_sub_3ul_3ul (&(t[0]), &(t[1]), &(t[2]), b[0], b[1], b[2]);
  modredc3ul_intset (r, t);
}

/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline int
modredc3ul_intbits (const modintredc3ul_t a)
{
  if (a[2] > 0UL)
    return 3*LONG_BIT - ularith_clz (a[2]);

  if (a[1] > 0UL)
    return 2*LONG_BIT - ularith_clz (a[1]);

  if (a[0] > 0UL)
    return LONG_BIT - ularith_clz (a[0]);
  
  return 0;
}


MAYBE_UNUSED
static inline void
modredc3ul_intshr (modintredc3ul_t r, const modintredc3ul_t s, const int i)
{
  r[0] = s[0];
  ularith_shrd (&(r[0]), s[1], i);
  r[1] = s[1];
  ularith_shrd (&(r[1]), s[2], i);
  r[2] = s[2] >> i;
}


MAYBE_UNUSED
static inline void
modredc3ul_intshl (modintredc3ul_t r, const modintredc3ul_t s, const int i)
{
  r[2] = s[2];
  ularith_shld (&(r[2]), s[1], i);
  r[1] = s[1];
  ularith_shld (&(r[1]), s[0], i);
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc3ul_intdivexact (modintredc3ul_t r, const modintredc3ul_t n,
                         const modintredc3ul_t d)
{
  unsigned long q[3], rr[3];
  q[0] = 0; q[1] = 0; q[2] = 0;
  int nd = 3;
  while (d[nd-1] == 0)
      nd--;
  mpn_tdiv_qr(q, rr, 0, n, 3, d, nd);
  r[0] = q[0];
  r[1] = q[1];
  r[2] = q[2];
}


/* r = n%d */
MAYBE_UNUSED
static inline void
modredc3ul_intmod (modintredc3ul_t r, const modintredc3ul_t n,
                    const modintredc3ul_t d)
{
  unsigned long q[3], rr[3];
  rr[0] = 0; rr[1] = 0; rr[2] = 0;
  int nd = 3;
  while (d[nd-1] == 0)
      nd--;
  mpn_tdiv_qr(q, rr, 0, n, 3, d, nd);
  r[0] = rr[0];
  r[1] = rr[1];
  r[2] = rr[2];
}


/* Functions for the modulus */

/* Init the modulus from a multi-word integer. s must point to an array of
   at least two unsigned longs, where s[0] is the low word of the modulus, 
   and s[1] is the high word. */
MAYBE_UNUSED
static inline void
modredc3ul_initmod_uls (modulusredc3ul_t m, const modintredc3ul_t s)
{
  ASSERT (s[2] > 0UL);
//  ASSERT (s[1] < (1UL << (LONG_BIT - 2)));
  modredc3ul_intset (m[0].m, s);
  m[0].invm = -ularith_invmod (s[0]);

  modredc3ul_intset_ul (m[0].one, 0UL);
  modredc3ul_intsub (m[0].one, m[0].one, m[0].m); /* 2^192 - m */
  modredc3ul_intmod (m[0].one, m[0].one, m[0].m);
  
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc3ul_t t;
    modredc3ul_get_uls (t, m[0].one, m);
    ASSERT_EXPENSIVE (modredc3ul_intequal_ul (t, 1UL));
  }
#endif
}


/* Returns the modulus to an array of unsigned longs. */
MAYBE_UNUSED
static inline void
modredc3ul_getmod_uls (modintredc3ul_t r, const modulusredc3ul_t m)
{
  modredc3ul_intset (r, m[0].m);
}


MAYBE_UNUSED
static inline void
modredc3ul_clearmod (modulusredc3ul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residueredc3ul_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc3ul_init (residueredc3ul_t r, const modulusredc3ul_t m MAYBE_UNUSED)
{
  modredc3ul_intset_ul (r, 0UL);
}


/* Initialises a residueredc3ul_t, but does not set it to zero. For fixed 
   length residueredc3ul_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc3ul_init_noset0 (residueredc3ul_t r MAYBE_UNUSED, 
			 const modulusredc3ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc3ul_clear (residueredc3ul_t r MAYBE_UNUSED, 
		   const modulusredc3ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc3ul_set (residueredc3ul_t r, const residueredc3ul_t s, 
		 const modulusredc3ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (s, m[0].m));
  modredc3ul_intset (r, s);
}


MAYBE_UNUSED
static inline void
modredc3ul_set_ul (residueredc3ul_t r, const unsigned long s, 
		    const modulusredc3ul_t m)
{
  modredc3ul_intset_ul (r, s);
  modredc3ul_tomontgomery (r, r, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc3ul_t t;
    modredc3ul_get_uls (t, r, m);
    ASSERT_EXPENSIVE (t[0] == s && t[1] == 0UL && t[2] == 0UL);
  }
#endif
}


/* Sets the residueredc3ul_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc3ul_set_ul_reduced (residueredc3ul_t r, const unsigned long s, 
			    const modulusredc3ul_t m MAYBE_UNUSED)
{
  modredc3ul_set_ul (r, s, m);
}


MAYBE_UNUSED
static inline void
modredc3ul_set_uls (residueredc3ul_t r, const modintredc3ul_t s, 
		     const modulusredc3ul_t m)
{
  if (!modredc3ul_intlt (s, m[0].m))
    modredc3ul_intmod (r, s, m[0].m);
  else
    modredc3ul_intset (r, s);

  modredc3ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc3ul_set_uls_reduced (residueredc3ul_t r, const modintredc3ul_t s, 
			     const modulusredc3ul_t m)
{
  ASSERT (modredc3ul_intlt (s, m[0].m));
  modredc3ul_intset (r, s);
  modredc3ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED 
static inline void 
modredc3ul_set0 (residueredc3ul_t r, const modulusredc3ul_t m MAYBE_UNUSED) 
{ 
  modredc3ul_intset_ul (r, 0UL);
}


MAYBE_UNUSED 
static inline void 
modredc3ul_set1 (residueredc3ul_t r, const modulusredc3ul_t m) 
{ 
  modredc3ul_intset (r, m[0].one);
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc3ul_swap (residueredc3ul_t a, residueredc3ul_t b, 
		  const modulusredc3ul_t m MAYBE_UNUSED)
{
  modintredc3ul_t t;
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc3ul_intlt (b, m[0].m));
  modredc3ul_intset (t, a);
  modredc3ul_intset (a, b);
  modredc3ul_intset (b, t);
}


/* Returns the least significant unsigned long of the residue. How to signal
   if the residue does not fit in one unsigned long? */

MAYBE_UNUSED
static inline unsigned long
modredc3ul_get_ul (const residueredc3ul_t s, 
		    const modulusredc3ul_t m MAYBE_UNUSED)
{
  unsigned long t[3];
  ASSERT_EXPENSIVE (modredc3ul_intlt (s, m[0].m));
  modredc3ul_frommontgomery (t, s, m);
  ASSERT (t[2] == 0UL && t[1] == 0UL);
  return t[0];
}


/* Returns the residue into an array of unsigned longs */

MAYBE_UNUSED
static inline void
modredc3ul_get_uls (modintredc3ul_t r, const residueredc3ul_t s, 
		     const modulusredc3ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (s, m[0].m));
  modredc3ul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc3ul_equal (const residueredc3ul_t a, const residueredc3ul_t b, 
		   const modulusredc3ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc3ul_intlt (b, m[0].m));
  return (modredc3ul_intequal(a, b));
}


MAYBE_UNUSED
static inline int
modredc3ul_is0 (const residueredc3ul_t a, const modulusredc3ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  return (modredc3ul_intequal_ul(a, 0UL));
}


MAYBE_UNUSED
static inline int
modredc3ul_is1 (const residueredc3ul_t a, const modulusredc3ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  return (modredc3ul_intequal(a, m[0].one));
}


MAYBE_UNUSED
static inline void
modredc3ul_add (residueredc3ul_t r, const residueredc3ul_t a, 
		 const residueredc3ul_t b, const modulusredc3ul_t m)
{
  /* r, a, and/or b may overlap */
  const unsigned long t0 = b[0], t1 = b[1], t2 = b[2];
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc3ul_intlt (b, m[0].m));

  modredc3ul_intset (r, a);
  ularith_add_3ul_3ul (&(r[0]), &(r[1]), &(r[2]), t0, t1, t2);
  ularith_sub_3ul_3ul_ge (&(r[0]), &(r[1]), &(r[2]),
          m[0].m[0], m[0].m[1], m[0].m[2]);
  ASSERT_EXPENSIVE (modredc3ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc3ul_sub (residueredc3ul_t r, const residueredc3ul_t a, 
		 const residueredc3ul_t b, const modulusredc3ul_t m)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc3ul_intlt (b, m[0].m));

#if 0
//#ifdef HAVE_GCC_STYLE_AMD64_ASM
  {
    unsigned long s1 = m[0].m[0], s2 = m[0].m[1], t1 = a[0], t2 = a[1];
    
    __asm__ (
	     "subq %4, %0\n\t"
	     "sbbq %5, %1\n\t"    /* r -= b */
	     "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
	     "cmovncq %6, %3\n"
	     : "+&r" (t1), "+&r" (t2), "+&r" (s1), "+r" (s2)
	     : "g" (b[0]), "g" (b[1]), "rm" (0UL)
	     : "cc"
	     );
    ularith_add_2ul_2ul (&t1, &t2, s1, s2);
    r[0] = t1;
    r[1] = t2;
  }
#else
  {
    unsigned long t1 = a[0], t2 = a[1], t3 = a[2];
    ularith_sub_3ul_3ul (&t1, &t2, &t3, b[0], b[1], b[2]);
    
    if (t3 > a[2] || (t3 == a[2] && t2 > a[1]) ||
            (t3 == a[2] && t2 == a[1] && t1 > a[0]))
      ularith_add_3ul_3ul (&t1, &t2, &t3, m[0].m[0], m[0].m[1], m[0].m[2]);

    r[0] = t1;
    r[1] = t2;
    r[2] = t3;
  }
#endif
}

MAYBE_UNUSED
static inline void
modredc3ul_add1 (residueredc3ul_t r, const residueredc3ul_t a, 
		  const modulusredc3ul_t m)
{
  modredc3ul_add(r, a, m[0].one, m);
}


MAYBE_UNUSED
static inline void
modredc3ul_add_ul (residueredc3ul_t r, const residueredc3ul_t a,
		    const unsigned long b, const modulusredc3ul_t m)
{
  residueredc3ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  modredc3ul_init_noset0 (t, m);
  modredc3ul_set_ul (t, b, m);
  modredc3ul_add (r, a, t, m);
  modredc3ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc3ul_sub_ul (residueredc3ul_t r, const residueredc3ul_t a,
		    const unsigned long b, const modulusredc3ul_t m)
{
  residueredc3ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  modredc3ul_init_noset0 (t, m);
  modredc3ul_set_ul (t, b, m);
  modredc3ul_sub (r, a, t, m);
  modredc3ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc3ul_neg (residueredc3ul_t r, const residueredc3ul_t a, 
		 const modulusredc3ul_t m)
{
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  if (a[0] == 0UL && a[1] == 0UL && a[2] == 0UL)
    modredc3ul_set (r, a, m);
  else
    {
      unsigned long t1 = m[0].m[0], t2 = m[0].m[1], t3 = m[0].m[2];
      ularith_sub_3ul_3ul (&t1, &t2, &t3, a[0], a[1], a[2]);
      r[0] = t1;
      r[1] = t2;
      r[2] = t3;
    }
}


MAYBE_UNUSED
static inline void
modredc3ul_div2 (residueredc3ul_t r, const residueredc3ul_t a, 
		  const modulusredc3ul_t m)
{
  ASSERT_EXPENSIVE (m[0].m[0] % 2UL != 0UL);
  ASSERT_EXPENSIVE (modredc3ul_intlt (a, m[0].m));
  modredc3ul_intset (r, a);
  if (r[0] % 2UL == 1UL)
    ularith_add_3ul_3ul (&(r[0]), &(r[1]), &(r[2]),
            m[0].m[0], m[0].m[1], m[0].m[2]);
  ularith_shrd (&(r[0]), r[1], 1);
  ularith_shrd (&(r[1]), r[2], 1);
  r[2] >>= 1;
}


#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif

MAYBE_UNUSED
static inline void
modredc3ul_mul (residueredc3ul_t r, const residueredc3ul_t a, 
               const residueredc3ul_t b, const modulusredc3ul_t m)
{
  unsigned long t[6];
  mul_3(t, a, b);
  redc_3(r, t, &m[0].invm, m[0].m);
}


MAYBE_UNUSED
static inline void
modredc3ul_sqr (residueredc3ul_t r, const residueredc3ul_t a, 
                 const modulusredc3ul_t m)
{
  unsigned long t[6];
  sqr_3(t, a);
  redc_3(r, t, &m[0].invm, m[0].m);
}


MAYBE_UNUSED
static inline int
modredc3ul_next (residueredc3ul_t r, const modulusredc3ul_t m)
{
  ularith_add_ul_3ul (&(r[0]), &(r[1]), &(r[2]), 1UL);
  return (modredc3ul_intequal (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc3ul_finished (const residueredc3ul_t r, const modulusredc3ul_t m)
{
  return (modredc3ul_intequal (r, m[0].m));
}

/* Division by small integer n, where (n-1)*m may overflow the most 
   significant word. Returns 1 if n is invertible modulo m, 0 if not. 
   
   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime 
   to n, for 0 <= i < n
   c = n^(-1) (mod word base)
*/

MAYBE_UNUSED
static inline int
modredc3ul_divn (residueredc3ul_t r, const residueredc3ul_t a, 
		  const unsigned long n, const unsigned long w_mod_n, 
		  const unsigned long *inv_n, MAYBE_UNUSED const unsigned long c,
		  const modulusredc3ul_t m)
{
  const unsigned long an = (((a[2] % n) * w_mod_n + (a[1] % n)) * w_mod_n
      + a[0] % n) % n;
  
  const unsigned long mn = (((m[0].m[2] % n) * w_mod_n + (m[0].m[1] % n))
          * w_mod_n + m[0].m[0] % n) % n;

  if (inv_n[mn] == 0)
    return 0;

  // build t of the form a + km such that t == 0 mod n
  unsigned long k = (inv_n[mn] * an) % n;
  unsigned long t[4];
  t[3] = mpn_mul_1(t, m[0].m, 3, k);
  t[3] += mpn_add_n(t, t, a, 3);

  // divide by n
  mpn_tdiv_qr(r, t, 0, t, 4, &n, 1);
  return 1;
}

MAYBE_UNUSED
static inline int
modredc3ul_inv (residueredc3ul_t r, const residueredc3ul_t s, 
        const modulusredc3ul_t m) {
  int ret;
  unsigned long t[3];
  modredc3ul_frommontgomery(t, s, m);
  ret = invmod_3(t, t, m[0].m);
  if (ret) {
    modredc3ul_tomontgomery(r, t, m);
    return 1;
  } else {
    r[0] = t[0]; r[1] = t[1]; r[2] = t[2];
    return 0;
  }
}


/* prototypes of non-inline functions */
int modredc3ul_div3 (residueredc3ul_t, const residueredc3ul_t, 
		      const modulusredc3ul_t);
int modredc3ul_div5 (residueredc3ul_t, const residueredc3ul_t, 
		       const modulusredc3ul_t);
int modredc3ul_div7 (residueredc3ul_t, const residueredc3ul_t, 
		      const modulusredc3ul_t);
int modredc3ul_div11 (residueredc3ul_t, const residueredc3ul_t, 
		       const modulusredc3ul_t);
int modredc3ul_div13 (residueredc3ul_t, const residueredc3ul_t, 
		       const modulusredc3ul_t);
void modredc3ul_gcd (modintredc3ul_t, const residueredc3ul_t, 
		      const modulusredc3ul_t);
void modredc3ul_pow_ul (residueredc3ul_t, const residueredc3ul_t, 
			 const unsigned long, const modulusredc3ul_t);
void modredc3ul_2pow_ul (residueredc3ul_t, const unsigned long, 
                          const modulusredc3ul_t);
void modredc3ul_pow_mp (residueredc3ul_t, const residueredc3ul_t, 
			 const unsigned long *, const int, 
			 const modulusredc3ul_t);
void modredc3ul_2pow_mp (residueredc3ul_t, const unsigned long *, const int, 
			  const modulusredc3ul_t);
void modredc3ul_V_ul (residueredc3ul_t, const residueredc3ul_t, 
		       const unsigned long, const modulusredc3ul_t);
void modredc3ul_V_mp (residueredc3ul_t, const residueredc3ul_t, 
		       const unsigned long *, const int, 
		       const modulusredc3ul_t);
int modredc3ul_sprp (const residueredc3ul_t, const modulusredc3ul_t);
int modredc3ul_sprp2 (const modulusredc3ul_t);
int modredc3ul_isprime (const modulusredc3ul_t);
int modredc3ul_jacobi (const residueredc3ul_t, const modulusredc3ul_t);
#endif  /* MODREDC_3UL3_H */
