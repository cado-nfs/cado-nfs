#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

/* The tests are based on assert() */
#ifdef NDEBUG
#  undef NDEBUG
#endif

#include <assert.h>
#include "mpfq_fixmp.h"

/* COMMON INPUT DATA
 * All test routines work on data areas whose name prescribes the length.
 * In order to provide fresh random inputs to all routines, we work with
 * read-only areas on input. These are all accessed via const pointers
 * which are global variables in the test file, and owned by main() only.
 * 
 * x, y are n-word long.
 * z is 2n-word long.
 * wx is a full limb.
 * 
 * For the hw routines, the following modifications apply:
 * 
 * the top word of x and y is a half-word.
 * z is only max(n+1, 2n-1)-word long.
 * in order to access a random half-limb, y[n-1] is recommended.
 * 
 * 
 * the output buffers are s, t, u, v, w. All have space allocated for 2n+1
 * words. s and t are generally used for the reference data computed by
 * gmp.
 */
const mp_limb_t * x, * y, * z;
mp_limb_t wx;
mp_limb_t * gx, * gy, * gz;



void test_fixmp_1() {
  mp_limb_t s[3];
  mp_limb_t t[3];
  mp_limb_t u[3];
  mp_limb_t v[3];
  mp_limb_t w[3];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[1] = { 2041087589UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[1] = {4929763703639915597UL};
#endif
  int j, k;
    // add
    u[1] = 0xdeadbeef;
    v[1] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 1);
    c2 = mpfq_fixmp_1_add(u, x, y);
    mpfq_fixmp_1_add_nc(v, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 1);
    c2 = mpfq_fixmp_1_add(u, x, x);
    mpfq_fixmp_1_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 1);
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c2 = mpfq_fixmp_1_add(v, v, y);
    mpfq_fixmp_1_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,1) == 0);
    assert (mpn_cmp(s,u,1) == 0);
    // add_ui
    u[1] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 1, wx);
    c2 = mpfq_fixmp_1_add_ui(u, x, wx);
    mpfq_fixmp_1_add_ui_nc(v, x, wx);
    assert (u[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c1 = mpn_add_1(s, x, 1, wx);
    c2 = mpfq_fixmp_1_add_ui(u, u, wx);
    mpfq_fixmp_1_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    // sub
    u[1] = 0xdeadbeef;
    v[1] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 1);
    c2 = mpfq_fixmp_1_sub(u, x, y);
    mpfq_fixmp_1_sub_nc(v, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 1);
    c2 = mpfq_fixmp_1_sub(u, x, x);
    mpfq_fixmp_1_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 1);
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c2 = mpfq_fixmp_1_sub(v, v, y);
    mpfq_fixmp_1_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,1) == 0);
    assert (mpn_cmp(s,u,1) == 0);
    // sub_ui
    u[1] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 1, wx);
    c2 = mpfq_fixmp_1_sub_ui(u, x, wx);
    mpfq_fixmp_1_sub_ui_nc(v, x, wx);
    assert (u[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c1 = mpn_sub_1(s, x, 1, wx);
    c2 = mpfq_fixmp_1_sub_ui(u, u, wx);
    mpfq_fixmp_1_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    // addmul1
    mpfq_copy(s, z, 2);
    mpfq_copy(u, z, 2);
    mpfq_copy(v, z, 2);
    mpfq_copy(w, z, 2);
    c1 = mpn_addmul_1(s, x, 1, wx);
    s[1] += c1;
    c3 = s[1] < c1;
    u[2]=0xdeadbeef;
    v[1]=0xdeadbeef;
    w[2]=0xdeadbeef;
    c4 = mpfq_fixmp_1_addmul1(u, x, wx);
    c2 = mpfq_fixmp_1_addmul1_shortz(v, x, wx);
    mpfq_fixmp_1_addmul1_nc(w, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (w[2] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 2) == 0);
    assert (mpn_cmp(s, v, 1) == 0);
    assert (mpn_cmp(s, w, 2) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    c1 = mpn_addmul_1(s, s, 1, wx);
    c2 = mpfq_fixmp_1_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 1) == 0);
    assert (c1 == c2);
    // mul1
    s[1] = mpn_mul_1(s, x, 1, wx);
    u[2] = 0xdeadbeef;
    mpfq_fixmp_1_mul1(u, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (mpn_cmp(s,u,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 1);
    mpfq_copy(v, x, 1);
    t[1] = mpn_mul_1(t, t, 1, wx);
    mpfq_fixmp_1_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 1 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 1);
    u[2]=0xdeadbeef;
    mpfq_fixmp_1_mul(u, x, y);
    assert (u[2]==0xdeadbeef);
    assert (mpn_cmp(s,u,2) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 1);
    u[1] = 0xdeadbeef;
    mpfq_fixmp_1_shortmul(u, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 1) == 0);
    // sqr
    mpn_mul_n(s, x, x, 1);
    u[2]=0xdeadbeef;
    mpfq_fixmp_1_sqr(u, x);
    assert (u[2]==0xdeadbeef);
    assert (mpn_cmp(s,u,2) == 0);
    // cmp
    j = mpn_cmp(x, y, 1);
    k = mpfq_fixmp_1_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 1);
    j = mpfq_fixmp_1_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 0);
    j = mpfq_fixmp_1_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_1_cmp_ui(u, ~x[0]);
    assert (j!=0);
    // mod
    mpfq_copy(v, y, 1);
    v[1-1] += !v[1-1];
    mpn_tdiv_qr(s+1, s, 0, z, 2, v, 1);
    mpfq_fixmp_1_mod(u, z, v);
    assert(mpn_cmp(s, u, 1) == 0);
    // inv
    mpfq_fixmp_1_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_1_invmod(v, u, P);
    mpfq_fixmp_1_invmod(v, v, P);
    assert(mpn_cmp(v, u, 1) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 1);
    mpfq_zero(v, 1);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_1_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_1_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 1);
    memset(u, ~0, 1 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 1);
    /* This creates an n-limb multiple of 257.  */
    v[1] = mpn_lshift(v, x, 1, 8);
    v[1] += mpfq_fixmp_1_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_1_add_ui(v, v, v[1]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_1_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 1, 257) == 0);
    //redc
    {
      mp_limb_t p[1], mip[1];
      mp_limb_t xe[1], ye[1], ze[1];
      mp_limb_t invR[1];
      
      // x[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 1);
      mpn_random2(p, 1);
      p[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[1-1] > p[1-1])
        p[1-1] = x[1-1];
      if (y[1-1] > p[1-1])
        p[1-1] = y[1-1];
      p[0] |= 1UL;
      p[1-1] += !p[1-1];
      mpfq_zero(w, 2*1);
      w[1]=1;
      mpfq_fixmp_1_mod(w, w, p);
      mpfq_fixmp_1_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 1);
        mpfq_fixmp_1_shortmul(t, mip, p);
	mpn_add_1(t, t, 1, 1);
	mpfq_fixmp_1_shortmul(u, t, mip);
	mpfq_fixmp_1_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 1));
      mpfq_fixmp_1_mgy_encode(xe, x, p);
      mpfq_fixmp_1_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_1_mul(s, x, y);
      mpfq_fixmp_1_mod(t, s, p);
      mpfq_fixmp_1_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_1_mul(s, xe, ye);
      mpfq_fixmp_1_mod(t, s, p);
      mpfq_fixmp_1_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 1) == 0);
      mpfq_fixmp_1_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 1) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_1_mul(t, xe, ye);
      s[2] = mpn_mul_1(s, t, 2, wx);
      mpfq_zero(t, 2*1+1);
      mpfq_fixmp_1_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_1_mod(s, t, p);
      mpfq_fixmp_1_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_1_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+1, v, 0, w, 1+1, p, 1);
      assert(mpn_cmp(v, s, 1) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 2 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 2, sat);
      mpfq_fixmp_1_redc(s, u, mip, p);
      mpfq_fixmp_1_mul1(w, s, sat);
      mpn_tdiv_qr(s+1, s, 0, w, 1+1, p, 1);
      mpfq_zero(w, 2*1+1);
      mpfq_fixmp_1_redc(w, v, mip, p);
      mpn_tdiv_qr(t+1, t, 0, w, 1+1, p, 1);
      assert(mpn_cmp(s, t, 1) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 2 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 2 + 1, sat);
      mpfq_fixmp_1_redc_ur(s, u, mip, p);
      mpfq_fixmp_1_mul1(w, s, sat);
      mpn_tdiv_qr(s+1, s, 0, w, 1+1, p, 1);
      mpfq_zero(w, 2*1+1);
      mpfq_fixmp_1_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+1, t, 0, w, 1+1, p, 1);
      assert(mpn_cmp(s, t, 1) == 0);
#ifdef  HAVE_native_1_mulredc
      mpfq_fixmp_1_mul(u, x, y);
      mpfq_fixmp_1_redc(s, u, mip, p);
      mpfq_fixmp_1_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 1) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 1, j);
    else
        mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_1_lshift(u, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 1, j);
    else
        mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_1_rshift(u, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // long_lshift
    j = wx % (1 * GMP_LIMB_BITS);
    mpfq_zero(s, 1);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 1 - k, j);
    else
        mpfq_copy(s + k, x, 1 - k);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_1_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // long_rshift
    j = wx % (1 * GMP_LIMB_BITS);
    mpfq_zero(s, 1);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 1 - k, j);
    else
        mpfq_copy(s, x + k, 1 - k);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_1_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 1) == 0);
}


void test_fixmp_2() {
  mp_limb_t s[5];
  mp_limb_t t[5];
  mp_limb_t u[5];
  mp_limb_t v[5];
  mp_limb_t w[5];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[2] = { 1737731653UL, 3654705850UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[2] = {14425133756266440979UL, 7028776506806380750UL};
#endif
  int j, k;
    // add
    u[2] = 0xdeadbeef;
    v[2] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 2);
    c2 = mpfq_fixmp_2_add(u, x, y);
    mpfq_fixmp_2_add_nc(v, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 2);
    c2 = mpfq_fixmp_2_add(u, x, x);
    mpfq_fixmp_2_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 2);
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c2 = mpfq_fixmp_2_add(v, v, y);
    mpfq_fixmp_2_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,2) == 0);
    assert (mpn_cmp(s,u,2) == 0);
    // add_ui
    u[2] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 2, wx);
    c2 = mpfq_fixmp_2_add_ui(u, x, wx);
    mpfq_fixmp_2_add_ui_nc(v, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c1 = mpn_add_1(s, x, 2, wx);
    c2 = mpfq_fixmp_2_add_ui(u, u, wx);
    mpfq_fixmp_2_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    // sub
    u[2] = 0xdeadbeef;
    v[2] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 2);
    c2 = mpfq_fixmp_2_sub(u, x, y);
    mpfq_fixmp_2_sub_nc(v, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 2);
    c2 = mpfq_fixmp_2_sub(u, x, x);
    mpfq_fixmp_2_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 2);
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c2 = mpfq_fixmp_2_sub(v, v, y);
    mpfq_fixmp_2_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,2) == 0);
    assert (mpn_cmp(s,u,2) == 0);
    // sub_ui
    u[2] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 2, wx);
    c2 = mpfq_fixmp_2_sub_ui(u, x, wx);
    mpfq_fixmp_2_sub_ui_nc(v, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c1 = mpn_sub_1(s, x, 2, wx);
    c2 = mpfq_fixmp_2_sub_ui(u, u, wx);
    mpfq_fixmp_2_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    // addmul1
    mpfq_copy(s, z, 3);
    mpfq_copy(u, z, 3);
    mpfq_copy(v, z, 3);
    mpfq_copy(w, z, 3);
    c1 = mpn_addmul_1(s, x, 2, wx);
    s[2] += c1;
    c3 = s[2] < c1;
    u[3]=0xdeadbeef;
    v[2]=0xdeadbeef;
    w[3]=0xdeadbeef;
    c4 = mpfq_fixmp_2_addmul1(u, x, wx);
    c2 = mpfq_fixmp_2_addmul1_shortz(v, x, wx);
    mpfq_fixmp_2_addmul1_nc(w, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (w[3] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 3) == 0);
    assert (mpn_cmp(s, v, 2) == 0);
    assert (mpn_cmp(s, w, 3) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    c1 = mpn_addmul_1(s, s, 2, wx);
    c2 = mpfq_fixmp_2_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 2) == 0);
    assert (c1 == c2);
    // mul1
    s[2] = mpn_mul_1(s, x, 2, wx);
    u[3] = 0xdeadbeef;
    mpfq_fixmp_2_mul1(u, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (mpn_cmp(s,u,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 2);
    mpfq_copy(v, x, 2);
    t[2] = mpn_mul_1(t, t, 2, wx);
    mpfq_fixmp_2_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 2 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 2);
    u[4]=0xdeadbeef;
    mpfq_fixmp_2_mul(u, x, y);
    assert (u[4]==0xdeadbeef);
    assert (mpn_cmp(s,u,4) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 2);
    u[2] = 0xdeadbeef;
    mpfq_fixmp_2_shortmul(u, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 2) == 0);
    // sqr
    mpn_mul_n(s, x, x, 2);
    u[4]=0xdeadbeef;
    mpfq_fixmp_2_sqr(u, x);
    assert (u[4]==0xdeadbeef);
    assert (mpn_cmp(s,u,4) == 0);
    // cmp
    j = mpn_cmp(x, y, 2);
    k = mpfq_fixmp_2_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 2);
    j = mpfq_fixmp_2_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 1);
    j = mpfq_fixmp_2_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_2_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_2_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 2);
    v[2-1] += !v[2-1];
    mpn_tdiv_qr(s+2, s, 0, z, 4, v, 2);
    mpfq_fixmp_2_mod(u, z, v);
    assert(mpn_cmp(s, u, 2) == 0);
    // inv
    mpfq_fixmp_2_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_2_invmod(v, u, P);
    mpfq_fixmp_2_invmod(v, v, P);
    assert(mpn_cmp(v, u, 2) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 2);
    mpfq_zero(v, 2);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_2_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_2_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 2);
    memset(u, ~0, 2 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 2);
    /* This creates an n-limb multiple of 257.  */
    v[2] = mpn_lshift(v, x, 2, 8);
    v[2] += mpfq_fixmp_2_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_2_add_ui(v, v, v[2]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_2_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 2, 257) == 0);
    //redc
    {
      mp_limb_t p[2], mip[2];
      mp_limb_t xe[2], ye[2], ze[2];
      mp_limb_t invR[2];
      
      // x[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 2);
      mpn_random2(p, 2);
      p[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[2-1] > p[2-1])
        p[2-1] = x[2-1];
      if (y[2-1] > p[2-1])
        p[2-1] = y[2-1];
      p[0] |= 1UL;
      p[2-1] += !p[2-1];
      mpfq_zero(w, 2*2);
      w[2]=1;
      mpfq_fixmp_2_mod(w, w, p);
      mpfq_fixmp_2_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 2);
        mpfq_fixmp_2_shortmul(t, mip, p);
	mpn_add_1(t, t, 2, 1);
	mpfq_fixmp_2_shortmul(u, t, mip);
	mpfq_fixmp_2_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 2));
      mpfq_fixmp_2_mgy_encode(xe, x, p);
      mpfq_fixmp_2_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_2_mul(s, x, y);
      mpfq_fixmp_2_mod(t, s, p);
      mpfq_fixmp_2_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_2_mul(s, xe, ye);
      mpfq_fixmp_2_mod(t, s, p);
      mpfq_fixmp_2_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 2) == 0);
      mpfq_fixmp_2_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 2) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_2_mul(t, xe, ye);
      s[4] = mpn_mul_1(s, t, 4, wx);
      mpfq_zero(t, 2*2+1);
      mpfq_fixmp_2_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_2_mod(s, t, p);
      mpfq_fixmp_2_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_2_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+2, v, 0, w, 2+1, p, 2);
      assert(mpn_cmp(v, s, 2) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 4 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 4, sat);
      mpfq_fixmp_2_redc(s, u, mip, p);
      mpfq_fixmp_2_mul1(w, s, sat);
      mpn_tdiv_qr(s+2, s, 0, w, 2+1, p, 2);
      mpfq_zero(w, 2*2+1);
      mpfq_fixmp_2_redc(w, v, mip, p);
      mpn_tdiv_qr(t+2, t, 0, w, 2+1, p, 2);
      assert(mpn_cmp(s, t, 2) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 4 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 4 + 1, sat);
      mpfq_fixmp_2_redc_ur(s, u, mip, p);
      mpfq_fixmp_2_mul1(w, s, sat);
      mpn_tdiv_qr(s+2, s, 0, w, 2+1, p, 2);
      mpfq_zero(w, 2*2+1);
      mpfq_fixmp_2_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+2, t, 0, w, 2+1, p, 2);
      assert(mpn_cmp(s, t, 2) == 0);
#ifdef  HAVE_native_2_mulredc
      mpfq_fixmp_2_mul(u, x, y);
      mpfq_fixmp_2_redc(s, u, mip, p);
      mpfq_fixmp_2_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 2) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 2, j);
    else
        mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_2_lshift(u, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 2, j);
    else
        mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_2_rshift(u, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // long_lshift
    j = wx % (2 * GMP_LIMB_BITS);
    mpfq_zero(s, 2);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 2 - k, j);
    else
        mpfq_copy(s + k, x, 2 - k);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_2_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // long_rshift
    j = wx % (2 * GMP_LIMB_BITS);
    mpfq_zero(s, 2);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 2 - k, j);
    else
        mpfq_copy(s, x + k, 2 - k);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_2_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 2) == 0);
}


void test_fixmp_3() {
  mp_limb_t s[7];
  mp_limb_t t[7];
  mp_limb_t u[7];
  mp_limb_t v[7];
  mp_limb_t w[7];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[3] = { 3826833745UL, 2279976717UL, 3984871455UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[3] = {9758664018554848775UL, 108797327114284110UL, 3855934483758865187UL};
#endif
  int j, k;
    // add
    u[3] = 0xdeadbeef;
    v[3] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 3);
    c2 = mpfq_fixmp_3_add(u, x, y);
    mpfq_fixmp_3_add_nc(v, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 3);
    c2 = mpfq_fixmp_3_add(u, x, x);
    mpfq_fixmp_3_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 3);
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c2 = mpfq_fixmp_3_add(v, v, y);
    mpfq_fixmp_3_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,3) == 0);
    assert (mpn_cmp(s,u,3) == 0);
    // add_ui
    u[3] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 3, wx);
    c2 = mpfq_fixmp_3_add_ui(u, x, wx);
    mpfq_fixmp_3_add_ui_nc(v, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c1 = mpn_add_1(s, x, 3, wx);
    c2 = mpfq_fixmp_3_add_ui(u, u, wx);
    mpfq_fixmp_3_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    // sub
    u[3] = 0xdeadbeef;
    v[3] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 3);
    c2 = mpfq_fixmp_3_sub(u, x, y);
    mpfq_fixmp_3_sub_nc(v, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 3);
    c2 = mpfq_fixmp_3_sub(u, x, x);
    mpfq_fixmp_3_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 3);
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c2 = mpfq_fixmp_3_sub(v, v, y);
    mpfq_fixmp_3_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,3) == 0);
    assert (mpn_cmp(s,u,3) == 0);
    // sub_ui
    u[3] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 3, wx);
    c2 = mpfq_fixmp_3_sub_ui(u, x, wx);
    mpfq_fixmp_3_sub_ui_nc(v, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c1 = mpn_sub_1(s, x, 3, wx);
    c2 = mpfq_fixmp_3_sub_ui(u, u, wx);
    mpfq_fixmp_3_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    // addmul1
    mpfq_copy(s, z, 4);
    mpfq_copy(u, z, 4);
    mpfq_copy(v, z, 4);
    mpfq_copy(w, z, 4);
    c1 = mpn_addmul_1(s, x, 3, wx);
    s[3] += c1;
    c3 = s[3] < c1;
    u[4]=0xdeadbeef;
    v[3]=0xdeadbeef;
    w[4]=0xdeadbeef;
    c4 = mpfq_fixmp_3_addmul1(u, x, wx);
    c2 = mpfq_fixmp_3_addmul1_shortz(v, x, wx);
    mpfq_fixmp_3_addmul1_nc(w, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (w[4] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 4) == 0);
    assert (mpn_cmp(s, v, 3) == 0);
    assert (mpn_cmp(s, w, 4) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    c1 = mpn_addmul_1(s, s, 3, wx);
    c2 = mpfq_fixmp_3_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 3) == 0);
    assert (c1 == c2);
    // mul1
    s[3] = mpn_mul_1(s, x, 3, wx);
    u[4] = 0xdeadbeef;
    mpfq_fixmp_3_mul1(u, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (mpn_cmp(s,u,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 3);
    mpfq_copy(v, x, 3);
    t[3] = mpn_mul_1(t, t, 3, wx);
    mpfq_fixmp_3_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 3 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 3);
    u[6]=0xdeadbeef;
    mpfq_fixmp_3_mul(u, x, y);
    assert (u[6]==0xdeadbeef);
    assert (mpn_cmp(s,u,6) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 3);
    u[3] = 0xdeadbeef;
    mpfq_fixmp_3_shortmul(u, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 3) == 0);
    // sqr
    mpn_mul_n(s, x, x, 3);
    u[6]=0xdeadbeef;
    mpfq_fixmp_3_sqr(u, x);
    assert (u[6]==0xdeadbeef);
    assert (mpn_cmp(s,u,6) == 0);
    // cmp
    j = mpn_cmp(x, y, 3);
    k = mpfq_fixmp_3_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 3);
    j = mpfq_fixmp_3_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 2);
    j = mpfq_fixmp_3_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_3_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_3_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 3);
    v[3-1] += !v[3-1];
    mpn_tdiv_qr(s+3, s, 0, z, 6, v, 3);
    mpfq_fixmp_3_mod(u, z, v);
    assert(mpn_cmp(s, u, 3) == 0);
    // inv
    mpfq_fixmp_3_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_3_invmod(v, u, P);
    mpfq_fixmp_3_invmod(v, v, P);
    assert(mpn_cmp(v, u, 3) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 3);
    mpfq_zero(v, 3);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_3_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_3_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 3);
    memset(u, ~0, 3 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 3);
    /* This creates an n-limb multiple of 257.  */
    v[3] = mpn_lshift(v, x, 3, 8);
    v[3] += mpfq_fixmp_3_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_3_add_ui(v, v, v[3]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_3_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 3, 257) == 0);
    //redc
    {
      mp_limb_t p[3], mip[3];
      mp_limb_t xe[3], ye[3], ze[3];
      mp_limb_t invR[3];
      
      // x[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 3);
      mpn_random2(p, 3);
      p[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[3-1] > p[3-1])
        p[3-1] = x[3-1];
      if (y[3-1] > p[3-1])
        p[3-1] = y[3-1];
      p[0] |= 1UL;
      p[3-1] += !p[3-1];
      mpfq_zero(w, 2*3);
      w[3]=1;
      mpfq_fixmp_3_mod(w, w, p);
      mpfq_fixmp_3_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 3);
        mpfq_fixmp_3_shortmul(t, mip, p);
	mpn_add_1(t, t, 3, 1);
	mpfq_fixmp_3_shortmul(u, t, mip);
	mpfq_fixmp_3_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 3));
      mpfq_fixmp_3_mgy_encode(xe, x, p);
      mpfq_fixmp_3_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_3_mul(s, x, y);
      mpfq_fixmp_3_mod(t, s, p);
      mpfq_fixmp_3_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_3_mul(s, xe, ye);
      mpfq_fixmp_3_mod(t, s, p);
      mpfq_fixmp_3_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 3) == 0);
      mpfq_fixmp_3_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 3) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_3_mul(t, xe, ye);
      s[6] = mpn_mul_1(s, t, 6, wx);
      mpfq_zero(t, 2*3+1);
      mpfq_fixmp_3_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_3_mod(s, t, p);
      mpfq_fixmp_3_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_3_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+3, v, 0, w, 3+1, p, 3);
      assert(mpn_cmp(v, s, 3) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 6 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 6, sat);
      mpfq_fixmp_3_redc(s, u, mip, p);
      mpfq_fixmp_3_mul1(w, s, sat);
      mpn_tdiv_qr(s+3, s, 0, w, 3+1, p, 3);
      mpfq_zero(w, 2*3+1);
      mpfq_fixmp_3_redc(w, v, mip, p);
      mpn_tdiv_qr(t+3, t, 0, w, 3+1, p, 3);
      assert(mpn_cmp(s, t, 3) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 6 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 6 + 1, sat);
      mpfq_fixmp_3_redc_ur(s, u, mip, p);
      mpfq_fixmp_3_mul1(w, s, sat);
      mpn_tdiv_qr(s+3, s, 0, w, 3+1, p, 3);
      mpfq_zero(w, 2*3+1);
      mpfq_fixmp_3_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+3, t, 0, w, 3+1, p, 3);
      assert(mpn_cmp(s, t, 3) == 0);
#ifdef  HAVE_native_3_mulredc
      mpfq_fixmp_3_mul(u, x, y);
      mpfq_fixmp_3_redc(s, u, mip, p);
      mpfq_fixmp_3_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 3) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 3, j);
    else
        mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_3_lshift(u, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 3, j);
    else
        mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_3_rshift(u, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // long_lshift
    j = wx % (3 * GMP_LIMB_BITS);
    mpfq_zero(s, 3);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 3 - k, j);
    else
        mpfq_copy(s + k, x, 3 - k);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_3_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // long_rshift
    j = wx % (3 * GMP_LIMB_BITS);
    mpfq_zero(s, 3);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 3 - k, j);
    else
        mpfq_copy(s, x + k, 3 - k);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_3_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 3) == 0);
}


void test_fixmp_4() {
  mp_limb_t s[9];
  mp_limb_t t[9];
  mp_limb_t u[9];
  mp_limb_t v[9];
  mp_limb_t w[9];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[4] = { 3662469475UL, 2096692762UL, 4151755841UL, 4009865730UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[4] = {12011675740079661751UL, 3294090837287775300UL, 9673935898323528142UL, 6244774036521541631UL };
#endif
  int j, k;
    // add
    u[4] = 0xdeadbeef;
    v[4] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 4);
    c2 = mpfq_fixmp_4_add(u, x, y);
    mpfq_fixmp_4_add_nc(v, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 4);
    c2 = mpfq_fixmp_4_add(u, x, x);
    mpfq_fixmp_4_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 4);
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c2 = mpfq_fixmp_4_add(v, v, y);
    mpfq_fixmp_4_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,4) == 0);
    assert (mpn_cmp(s,u,4) == 0);
    // add_ui
    u[4] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 4, wx);
    c2 = mpfq_fixmp_4_add_ui(u, x, wx);
    mpfq_fixmp_4_add_ui_nc(v, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c1 = mpn_add_1(s, x, 4, wx);
    c2 = mpfq_fixmp_4_add_ui(u, u, wx);
    mpfq_fixmp_4_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    // sub
    u[4] = 0xdeadbeef;
    v[4] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 4);
    c2 = mpfq_fixmp_4_sub(u, x, y);
    mpfq_fixmp_4_sub_nc(v, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 4);
    c2 = mpfq_fixmp_4_sub(u, x, x);
    mpfq_fixmp_4_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 4);
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c2 = mpfq_fixmp_4_sub(v, v, y);
    mpfq_fixmp_4_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,4) == 0);
    assert (mpn_cmp(s,u,4) == 0);
    // sub_ui
    u[4] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 4, wx);
    c2 = mpfq_fixmp_4_sub_ui(u, x, wx);
    mpfq_fixmp_4_sub_ui_nc(v, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c1 = mpn_sub_1(s, x, 4, wx);
    c2 = mpfq_fixmp_4_sub_ui(u, u, wx);
    mpfq_fixmp_4_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    // addmul1
    mpfq_copy(s, z, 5);
    mpfq_copy(u, z, 5);
    mpfq_copy(v, z, 5);
    mpfq_copy(w, z, 5);
    c1 = mpn_addmul_1(s, x, 4, wx);
    s[4] += c1;
    c3 = s[4] < c1;
    u[5]=0xdeadbeef;
    v[4]=0xdeadbeef;
    w[5]=0xdeadbeef;
    c4 = mpfq_fixmp_4_addmul1(u, x, wx);
    c2 = mpfq_fixmp_4_addmul1_shortz(v, x, wx);
    mpfq_fixmp_4_addmul1_nc(w, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (w[5] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 5) == 0);
    assert (mpn_cmp(s, v, 4) == 0);
    assert (mpn_cmp(s, w, 5) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    c1 = mpn_addmul_1(s, s, 4, wx);
    c2 = mpfq_fixmp_4_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 4) == 0);
    assert (c1 == c2);
    // mul1
    s[4] = mpn_mul_1(s, x, 4, wx);
    u[5] = 0xdeadbeef;
    mpfq_fixmp_4_mul1(u, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (mpn_cmp(s,u,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 4);
    mpfq_copy(v, x, 4);
    t[4] = mpn_mul_1(t, t, 4, wx);
    mpfq_fixmp_4_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 4 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 4);
    u[8]=0xdeadbeef;
    mpfq_fixmp_4_mul(u, x, y);
    assert (u[8]==0xdeadbeef);
    assert (mpn_cmp(s,u,8) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 4);
    u[4] = 0xdeadbeef;
    mpfq_fixmp_4_shortmul(u, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 4) == 0);
    // sqr
    mpn_mul_n(s, x, x, 4);
    u[8]=0xdeadbeef;
    mpfq_fixmp_4_sqr(u, x);
    assert (u[8]==0xdeadbeef);
    assert (mpn_cmp(s,u,8) == 0);
    // cmp
    j = mpn_cmp(x, y, 4);
    k = mpfq_fixmp_4_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 4);
    j = mpfq_fixmp_4_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 3);
    j = mpfq_fixmp_4_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_4_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_4_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 4);
    v[4-1] += !v[4-1];
    mpn_tdiv_qr(s+4, s, 0, z, 8, v, 4);
    mpfq_fixmp_4_mod(u, z, v);
    assert(mpn_cmp(s, u, 4) == 0);
    // inv
    mpfq_fixmp_4_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_4_invmod(v, u, P);
    mpfq_fixmp_4_invmod(v, v, P);
    assert(mpn_cmp(v, u, 4) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 4);
    mpfq_zero(v, 4);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_4_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_4_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 4);
    memset(u, ~0, 4 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 4);
    /* This creates an n-limb multiple of 257.  */
    v[4] = mpn_lshift(v, x, 4, 8);
    v[4] += mpfq_fixmp_4_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_4_add_ui(v, v, v[4]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_4_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 4, 257) == 0);
    //redc
    {
      mp_limb_t p[4], mip[4];
      mp_limb_t xe[4], ye[4], ze[4];
      mp_limb_t invR[4];
      
      // x[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 4);
      mpn_random2(p, 4);
      p[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[4-1] > p[4-1])
        p[4-1] = x[4-1];
      if (y[4-1] > p[4-1])
        p[4-1] = y[4-1];
      p[0] |= 1UL;
      p[4-1] += !p[4-1];
      mpfq_zero(w, 2*4);
      w[4]=1;
      mpfq_fixmp_4_mod(w, w, p);
      mpfq_fixmp_4_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 4);
        mpfq_fixmp_4_shortmul(t, mip, p);
	mpn_add_1(t, t, 4, 1);
	mpfq_fixmp_4_shortmul(u, t, mip);
	mpfq_fixmp_4_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 4));
      mpfq_fixmp_4_mgy_encode(xe, x, p);
      mpfq_fixmp_4_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_4_mul(s, x, y);
      mpfq_fixmp_4_mod(t, s, p);
      mpfq_fixmp_4_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_4_mul(s, xe, ye);
      mpfq_fixmp_4_mod(t, s, p);
      mpfq_fixmp_4_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 4) == 0);
      mpfq_fixmp_4_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 4) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_4_mul(t, xe, ye);
      s[8] = mpn_mul_1(s, t, 8, wx);
      mpfq_zero(t, 2*4+1);
      mpfq_fixmp_4_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_4_mod(s, t, p);
      mpfq_fixmp_4_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_4_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+4, v, 0, w, 4+1, p, 4);
      assert(mpn_cmp(v, s, 4) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 8 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 8, sat);
      mpfq_fixmp_4_redc(s, u, mip, p);
      mpfq_fixmp_4_mul1(w, s, sat);
      mpn_tdiv_qr(s+4, s, 0, w, 4+1, p, 4);
      mpfq_zero(w, 2*4+1);
      mpfq_fixmp_4_redc(w, v, mip, p);
      mpn_tdiv_qr(t+4, t, 0, w, 4+1, p, 4);
      assert(mpn_cmp(s, t, 4) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 8 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 8 + 1, sat);
      mpfq_fixmp_4_redc_ur(s, u, mip, p);
      mpfq_fixmp_4_mul1(w, s, sat);
      mpn_tdiv_qr(s+4, s, 0, w, 4+1, p, 4);
      mpfq_zero(w, 2*4+1);
      mpfq_fixmp_4_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+4, t, 0, w, 4+1, p, 4);
      assert(mpn_cmp(s, t, 4) == 0);
#ifdef  HAVE_native_4_mulredc
      mpfq_fixmp_4_mul(u, x, y);
      mpfq_fixmp_4_redc(s, u, mip, p);
      mpfq_fixmp_4_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 4) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 4, j);
    else
        mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_4_lshift(u, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 4, j);
    else
        mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_4_rshift(u, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // long_lshift
    j = wx % (4 * GMP_LIMB_BITS);
    mpfq_zero(s, 4);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 4 - k, j);
    else
        mpfq_copy(s + k, x, 4 - k);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_4_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // long_rshift
    j = wx % (4 * GMP_LIMB_BITS);
    mpfq_zero(s, 4);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 4 - k, j);
    else
        mpfq_copy(s, x + k, 4 - k);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_4_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 4) == 0);
}


void test_fixmp_5() {
  mp_limb_t s[11];
  mp_limb_t t[11];
  mp_limb_t u[11];
  mp_limb_t v[11];
  mp_limb_t w[11];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[5] = { 4034459419UL, 3797792253UL, 1419478273UL, 2675749510UL, 3664727098UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[5] = {11766063743839433833UL, 17023338808517031849UL, 6384879829007101141UL, 9814014250957810811UL,5856459693223253397UL };
#endif
  int j, k;
    // add
    u[5] = 0xdeadbeef;
    v[5] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 5);
    c2 = mpfq_fixmp_5_add(u, x, y);
    mpfq_fixmp_5_add_nc(v, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 5);
    c2 = mpfq_fixmp_5_add(u, x, x);
    mpfq_fixmp_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 5);
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c2 = mpfq_fixmp_5_add(v, v, y);
    mpfq_fixmp_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,5) == 0);
    assert (mpn_cmp(s,u,5) == 0);
    // add_ui
    u[5] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 5, wx);
    c2 = mpfq_fixmp_5_add_ui(u, x, wx);
    mpfq_fixmp_5_add_ui_nc(v, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c1 = mpn_add_1(s, x, 5, wx);
    c2 = mpfq_fixmp_5_add_ui(u, u, wx);
    mpfq_fixmp_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    // sub
    u[5] = 0xdeadbeef;
    v[5] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 5);
    c2 = mpfq_fixmp_5_sub(u, x, y);
    mpfq_fixmp_5_sub_nc(v, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 5);
    c2 = mpfq_fixmp_5_sub(u, x, x);
    mpfq_fixmp_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 5);
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c2 = mpfq_fixmp_5_sub(v, v, y);
    mpfq_fixmp_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,5) == 0);
    assert (mpn_cmp(s,u,5) == 0);
    // sub_ui
    u[5] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 5, wx);
    c2 = mpfq_fixmp_5_sub_ui(u, x, wx);
    mpfq_fixmp_5_sub_ui_nc(v, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c1 = mpn_sub_1(s, x, 5, wx);
    c2 = mpfq_fixmp_5_sub_ui(u, u, wx);
    mpfq_fixmp_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    // addmul1
    mpfq_copy(s, z, 6);
    mpfq_copy(u, z, 6);
    mpfq_copy(v, z, 6);
    mpfq_copy(w, z, 6);
    c1 = mpn_addmul_1(s, x, 5, wx);
    s[5] += c1;
    c3 = s[5] < c1;
    u[6]=0xdeadbeef;
    v[5]=0xdeadbeef;
    w[6]=0xdeadbeef;
    c4 = mpfq_fixmp_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_5_addmul1_nc(w, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (w[6] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 6) == 0);
    assert (mpn_cmp(s, v, 5) == 0);
    assert (mpn_cmp(s, w, 6) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    c1 = mpn_addmul_1(s, s, 5, wx);
    c2 = mpfq_fixmp_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 5) == 0);
    assert (c1 == c2);
    // mul1
    s[5] = mpn_mul_1(s, x, 5, wx);
    u[6] = 0xdeadbeef;
    mpfq_fixmp_5_mul1(u, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (mpn_cmp(s,u,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 5);
    mpfq_copy(v, x, 5);
    t[5] = mpn_mul_1(t, t, 5, wx);
    mpfq_fixmp_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 5 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 5);
    u[10]=0xdeadbeef;
    mpfq_fixmp_5_mul(u, x, y);
    assert (u[10]==0xdeadbeef);
    assert (mpn_cmp(s,u,10) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 5);
    u[5] = 0xdeadbeef;
    mpfq_fixmp_5_shortmul(u, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 5) == 0);
    // sqr
    mpn_mul_n(s, x, x, 5);
    u[10]=0xdeadbeef;
    mpfq_fixmp_5_sqr(u, x);
    assert (u[10]==0xdeadbeef);
    assert (mpn_cmp(s,u,10) == 0);
    // cmp
    j = mpn_cmp(x, y, 5);
    k = mpfq_fixmp_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 5);
    j = mpfq_fixmp_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 4);
    j = mpfq_fixmp_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 5);
    v[5-1] += !v[5-1];
    mpn_tdiv_qr(s+5, s, 0, z, 10, v, 5);
    mpfq_fixmp_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 5) == 0);
    // inv
    mpfq_fixmp_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_5_invmod(v, u, P);
    mpfq_fixmp_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 5) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 5);
    mpfq_zero(v, 5);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 5);
    memset(u, ~0, 5 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 5);
    /* This creates an n-limb multiple of 257.  */
    v[5] = mpn_lshift(v, x, 5, 8);
    v[5] += mpfq_fixmp_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_5_add_ui(v, v, v[5]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 5, 257) == 0);
    //redc
    {
      mp_limb_t p[5], mip[5];
      mp_limb_t xe[5], ye[5], ze[5];
      mp_limb_t invR[5];
      
      // x[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 5);
      mpn_random2(p, 5);
      p[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[5-1] > p[5-1])
        p[5-1] = x[5-1];
      if (y[5-1] > p[5-1])
        p[5-1] = y[5-1];
      p[0] |= 1UL;
      p[5-1] += !p[5-1];
      mpfq_zero(w, 2*5);
      w[5]=1;
      mpfq_fixmp_5_mod(w, w, p);
      mpfq_fixmp_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 5);
        mpfq_fixmp_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 5, 1);
	mpfq_fixmp_5_shortmul(u, t, mip);
	mpfq_fixmp_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 5));
      mpfq_fixmp_5_mgy_encode(xe, x, p);
      mpfq_fixmp_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_5_mul(s, x, y);
      mpfq_fixmp_5_mod(t, s, p);
      mpfq_fixmp_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_5_mul(s, xe, ye);
      mpfq_fixmp_5_mod(t, s, p);
      mpfq_fixmp_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 5) == 0);
      mpfq_fixmp_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 5) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_5_mul(t, xe, ye);
      s[10] = mpn_mul_1(s, t, 10, wx);
      mpfq_zero(t, 2*5+1);
      mpfq_fixmp_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_5_mod(s, t, p);
      mpfq_fixmp_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+5, v, 0, w, 5+1, p, 5);
      assert(mpn_cmp(v, s, 5) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 10 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 10, sat);
      mpfq_fixmp_5_redc(s, u, mip, p);
      mpfq_fixmp_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+5, s, 0, w, 5+1, p, 5);
      mpfq_zero(w, 2*5+1);
      mpfq_fixmp_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+5, t, 0, w, 5+1, p, 5);
      assert(mpn_cmp(s, t, 5) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 10 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 10 + 1, sat);
      mpfq_fixmp_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+5, s, 0, w, 5+1, p, 5);
      mpfq_zero(w, 2*5+1);
      mpfq_fixmp_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+5, t, 0, w, 5+1, p, 5);
      assert(mpn_cmp(s, t, 5) == 0);
#ifdef  HAVE_native_5_mulredc
      mpfq_fixmp_5_mul(u, x, y);
      mpfq_fixmp_5_redc(s, u, mip, p);
      mpfq_fixmp_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 5) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 5, j);
    else
        mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_5_lshift(u, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 5, j);
    else
        mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_5_rshift(u, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // long_lshift
    j = wx % (5 * GMP_LIMB_BITS);
    mpfq_zero(s, 5);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 5 - k, j);
    else
        mpfq_copy(s + k, x, 5 - k);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // long_rshift
    j = wx % (5 * GMP_LIMB_BITS);
    mpfq_zero(s, 5);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 5 - k, j);
    else
        mpfq_copy(s, x + k, 5 - k);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 5) == 0);
}


void test_fixmp_6() {
  mp_limb_t s[13];
  mp_limb_t t[13];
  mp_limb_t u[13];
  mp_limb_t v[13];
  mp_limb_t w[13];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[6] = { 2813719779UL, 3907769622UL, 704006380UL, 1485932037UL, 661860009UL, 2968664580UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[6] = {9332489496020345727UL, 13059118375404545793UL, 543826843599586942UL, 568657921352937073UL, 8714542686157595041UL, 8377129812810584371UL };
#endif
  int j, k;
    // add
    u[6] = 0xdeadbeef;
    v[6] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 6);
    c2 = mpfq_fixmp_6_add(u, x, y);
    mpfq_fixmp_6_add_nc(v, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 6);
    c2 = mpfq_fixmp_6_add(u, x, x);
    mpfq_fixmp_6_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 6);
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c2 = mpfq_fixmp_6_add(v, v, y);
    mpfq_fixmp_6_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,6) == 0);
    assert (mpn_cmp(s,u,6) == 0);
    // add_ui
    u[6] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 6, wx);
    c2 = mpfq_fixmp_6_add_ui(u, x, wx);
    mpfq_fixmp_6_add_ui_nc(v, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c1 = mpn_add_1(s, x, 6, wx);
    c2 = mpfq_fixmp_6_add_ui(u, u, wx);
    mpfq_fixmp_6_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    // sub
    u[6] = 0xdeadbeef;
    v[6] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 6);
    c2 = mpfq_fixmp_6_sub(u, x, y);
    mpfq_fixmp_6_sub_nc(v, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 6);
    c2 = mpfq_fixmp_6_sub(u, x, x);
    mpfq_fixmp_6_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 6);
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c2 = mpfq_fixmp_6_sub(v, v, y);
    mpfq_fixmp_6_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,6) == 0);
    assert (mpn_cmp(s,u,6) == 0);
    // sub_ui
    u[6] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 6, wx);
    c2 = mpfq_fixmp_6_sub_ui(u, x, wx);
    mpfq_fixmp_6_sub_ui_nc(v, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c1 = mpn_sub_1(s, x, 6, wx);
    c2 = mpfq_fixmp_6_sub_ui(u, u, wx);
    mpfq_fixmp_6_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    // addmul1
    mpfq_copy(s, z, 7);
    mpfq_copy(u, z, 7);
    mpfq_copy(v, z, 7);
    mpfq_copy(w, z, 7);
    c1 = mpn_addmul_1(s, x, 6, wx);
    s[6] += c1;
    c3 = s[6] < c1;
    u[7]=0xdeadbeef;
    v[6]=0xdeadbeef;
    w[7]=0xdeadbeef;
    c4 = mpfq_fixmp_6_addmul1(u, x, wx);
    c2 = mpfq_fixmp_6_addmul1_shortz(v, x, wx);
    mpfq_fixmp_6_addmul1_nc(w, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (w[7] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 7) == 0);
    assert (mpn_cmp(s, v, 6) == 0);
    assert (mpn_cmp(s, w, 7) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    c1 = mpn_addmul_1(s, s, 6, wx);
    c2 = mpfq_fixmp_6_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 6) == 0);
    assert (c1 == c2);
    // mul1
    s[6] = mpn_mul_1(s, x, 6, wx);
    u[7] = 0xdeadbeef;
    mpfq_fixmp_6_mul1(u, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (mpn_cmp(s,u,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 6);
    mpfq_copy(v, x, 6);
    t[6] = mpn_mul_1(t, t, 6, wx);
    mpfq_fixmp_6_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 6 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 6);
    u[12]=0xdeadbeef;
    mpfq_fixmp_6_mul(u, x, y);
    assert (u[12]==0xdeadbeef);
    assert (mpn_cmp(s,u,12) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 6);
    u[6] = 0xdeadbeef;
    mpfq_fixmp_6_shortmul(u, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 6) == 0);
    // sqr
    mpn_mul_n(s, x, x, 6);
    u[12]=0xdeadbeef;
    mpfq_fixmp_6_sqr(u, x);
    assert (u[12]==0xdeadbeef);
    assert (mpn_cmp(s,u,12) == 0);
    // cmp
    j = mpn_cmp(x, y, 6);
    k = mpfq_fixmp_6_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 6);
    j = mpfq_fixmp_6_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 5);
    j = mpfq_fixmp_6_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_6_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_6_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 6);
    v[6-1] += !v[6-1];
    mpn_tdiv_qr(s+6, s, 0, z, 12, v, 6);
    mpfq_fixmp_6_mod(u, z, v);
    assert(mpn_cmp(s, u, 6) == 0);
    // inv
    mpfq_fixmp_6_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_6_invmod(v, u, P);
    mpfq_fixmp_6_invmod(v, v, P);
    assert(mpn_cmp(v, u, 6) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 6);
    mpfq_zero(v, 6);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_6_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_6_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 6);
    memset(u, ~0, 6 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 6);
    /* This creates an n-limb multiple of 257.  */
    v[6] = mpn_lshift(v, x, 6, 8);
    v[6] += mpfq_fixmp_6_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_6_add_ui(v, v, v[6]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_6_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 6, 257) == 0);
    //redc
    {
      mp_limb_t p[6], mip[6];
      mp_limb_t xe[6], ye[6], ze[6];
      mp_limb_t invR[6];
      
      // x[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 6);
      mpn_random2(p, 6);
      p[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[6-1] > p[6-1])
        p[6-1] = x[6-1];
      if (y[6-1] > p[6-1])
        p[6-1] = y[6-1];
      p[0] |= 1UL;
      p[6-1] += !p[6-1];
      mpfq_zero(w, 2*6);
      w[6]=1;
      mpfq_fixmp_6_mod(w, w, p);
      mpfq_fixmp_6_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 6);
        mpfq_fixmp_6_shortmul(t, mip, p);
	mpn_add_1(t, t, 6, 1);
	mpfq_fixmp_6_shortmul(u, t, mip);
	mpfq_fixmp_6_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 6));
      mpfq_fixmp_6_mgy_encode(xe, x, p);
      mpfq_fixmp_6_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_6_mul(s, x, y);
      mpfq_fixmp_6_mod(t, s, p);
      mpfq_fixmp_6_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_6_mul(s, xe, ye);
      mpfq_fixmp_6_mod(t, s, p);
      mpfq_fixmp_6_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 6) == 0);
      mpfq_fixmp_6_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 6) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_6_mul(t, xe, ye);
      s[12] = mpn_mul_1(s, t, 12, wx);
      mpfq_zero(t, 2*6+1);
      mpfq_fixmp_6_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_6_mod(s, t, p);
      mpfq_fixmp_6_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_6_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+6, v, 0, w, 6+1, p, 6);
      assert(mpn_cmp(v, s, 6) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 12 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 12, sat);
      mpfq_fixmp_6_redc(s, u, mip, p);
      mpfq_fixmp_6_mul1(w, s, sat);
      mpn_tdiv_qr(s+6, s, 0, w, 6+1, p, 6);
      mpfq_zero(w, 2*6+1);
      mpfq_fixmp_6_redc(w, v, mip, p);
      mpn_tdiv_qr(t+6, t, 0, w, 6+1, p, 6);
      assert(mpn_cmp(s, t, 6) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 12 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 12 + 1, sat);
      mpfq_fixmp_6_redc_ur(s, u, mip, p);
      mpfq_fixmp_6_mul1(w, s, sat);
      mpn_tdiv_qr(s+6, s, 0, w, 6+1, p, 6);
      mpfq_zero(w, 2*6+1);
      mpfq_fixmp_6_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+6, t, 0, w, 6+1, p, 6);
      assert(mpn_cmp(s, t, 6) == 0);
#ifdef  HAVE_native_6_mulredc
      mpfq_fixmp_6_mul(u, x, y);
      mpfq_fixmp_6_redc(s, u, mip, p);
      mpfq_fixmp_6_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 6) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 6, j);
    else
        mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_6_lshift(u, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 6, j);
    else
        mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_6_rshift(u, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // long_lshift
    j = wx % (6 * GMP_LIMB_BITS);
    mpfq_zero(s, 6);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 6 - k, j);
    else
        mpfq_copy(s + k, x, 6 - k);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_6_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // long_rshift
    j = wx % (6 * GMP_LIMB_BITS);
    mpfq_zero(s, 6);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 6 - k, j);
    else
        mpfq_copy(s, x + k, 6 - k);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_6_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 6) == 0);
}


void test_fixmp_7() {
  mp_limb_t s[15];
  mp_limb_t t[15];
  mp_limb_t u[15];
  mp_limb_t v[15];
  mp_limb_t w[15];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[7] = { 3784709369UL, 269443326UL, 4028649229UL, 2906318846UL, 1307656400UL, 167308958UL, 3095675918UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[7] = {8305431681600953837UL, 8511912794376737076UL, 5827616680491403508UL, 11764963549898802560UL, 9952224619298044241UL, 2593919323804169004UL, 5707166315511930231UL };
#endif
  int j, k;
    // add
    u[7] = 0xdeadbeef;
    v[7] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 7);
    c2 = mpfq_fixmp_7_add(u, x, y);
    mpfq_fixmp_7_add_nc(v, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 7);
    c2 = mpfq_fixmp_7_add(u, x, x);
    mpfq_fixmp_7_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 7);
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c2 = mpfq_fixmp_7_add(v, v, y);
    mpfq_fixmp_7_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,7) == 0);
    assert (mpn_cmp(s,u,7) == 0);
    // add_ui
    u[7] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 7, wx);
    c2 = mpfq_fixmp_7_add_ui(u, x, wx);
    mpfq_fixmp_7_add_ui_nc(v, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c1 = mpn_add_1(s, x, 7, wx);
    c2 = mpfq_fixmp_7_add_ui(u, u, wx);
    mpfq_fixmp_7_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    // sub
    u[7] = 0xdeadbeef;
    v[7] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 7);
    c2 = mpfq_fixmp_7_sub(u, x, y);
    mpfq_fixmp_7_sub_nc(v, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 7);
    c2 = mpfq_fixmp_7_sub(u, x, x);
    mpfq_fixmp_7_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 7);
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c2 = mpfq_fixmp_7_sub(v, v, y);
    mpfq_fixmp_7_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,7) == 0);
    assert (mpn_cmp(s,u,7) == 0);
    // sub_ui
    u[7] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 7, wx);
    c2 = mpfq_fixmp_7_sub_ui(u, x, wx);
    mpfq_fixmp_7_sub_ui_nc(v, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c1 = mpn_sub_1(s, x, 7, wx);
    c2 = mpfq_fixmp_7_sub_ui(u, u, wx);
    mpfq_fixmp_7_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    // addmul1
    mpfq_copy(s, z, 8);
    mpfq_copy(u, z, 8);
    mpfq_copy(v, z, 8);
    mpfq_copy(w, z, 8);
    c1 = mpn_addmul_1(s, x, 7, wx);
    s[7] += c1;
    c3 = s[7] < c1;
    u[8]=0xdeadbeef;
    v[7]=0xdeadbeef;
    w[8]=0xdeadbeef;
    c4 = mpfq_fixmp_7_addmul1(u, x, wx);
    c2 = mpfq_fixmp_7_addmul1_shortz(v, x, wx);
    mpfq_fixmp_7_addmul1_nc(w, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (w[8] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 8) == 0);
    assert (mpn_cmp(s, v, 7) == 0);
    assert (mpn_cmp(s, w, 8) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    c1 = mpn_addmul_1(s, s, 7, wx);
    c2 = mpfq_fixmp_7_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 7) == 0);
    assert (c1 == c2);
    // mul1
    s[7] = mpn_mul_1(s, x, 7, wx);
    u[8] = 0xdeadbeef;
    mpfq_fixmp_7_mul1(u, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (mpn_cmp(s,u,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 7);
    mpfq_copy(v, x, 7);
    t[7] = mpn_mul_1(t, t, 7, wx);
    mpfq_fixmp_7_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 7 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 7);
    u[14]=0xdeadbeef;
    mpfq_fixmp_7_mul(u, x, y);
    assert (u[14]==0xdeadbeef);
    assert (mpn_cmp(s,u,14) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 7);
    u[7] = 0xdeadbeef;
    mpfq_fixmp_7_shortmul(u, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 7) == 0);
    // sqr
    mpn_mul_n(s, x, x, 7);
    u[14]=0xdeadbeef;
    mpfq_fixmp_7_sqr(u, x);
    assert (u[14]==0xdeadbeef);
    assert (mpn_cmp(s,u,14) == 0);
    // cmp
    j = mpn_cmp(x, y, 7);
    k = mpfq_fixmp_7_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 7);
    j = mpfq_fixmp_7_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 6);
    j = mpfq_fixmp_7_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_7_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_7_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 7);
    v[7-1] += !v[7-1];
    mpn_tdiv_qr(s+7, s, 0, z, 14, v, 7);
    mpfq_fixmp_7_mod(u, z, v);
    assert(mpn_cmp(s, u, 7) == 0);
    // inv
    mpfq_fixmp_7_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_7_invmod(v, u, P);
    mpfq_fixmp_7_invmod(v, v, P);
    assert(mpn_cmp(v, u, 7) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 7);
    mpfq_zero(v, 7);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_7_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_7_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 7);
    memset(u, ~0, 7 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 7);
    /* This creates an n-limb multiple of 257.  */
    v[7] = mpn_lshift(v, x, 7, 8);
    v[7] += mpfq_fixmp_7_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_7_add_ui(v, v, v[7]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_7_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 7, 257) == 0);
    //redc
    {
      mp_limb_t p[7], mip[7];
      mp_limb_t xe[7], ye[7], ze[7];
      mp_limb_t invR[7];
      
      // x[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 7);
      mpn_random2(p, 7);
      p[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[7-1] > p[7-1])
        p[7-1] = x[7-1];
      if (y[7-1] > p[7-1])
        p[7-1] = y[7-1];
      p[0] |= 1UL;
      p[7-1] += !p[7-1];
      mpfq_zero(w, 2*7);
      w[7]=1;
      mpfq_fixmp_7_mod(w, w, p);
      mpfq_fixmp_7_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 7);
        mpfq_fixmp_7_shortmul(t, mip, p);
	mpn_add_1(t, t, 7, 1);
	mpfq_fixmp_7_shortmul(u, t, mip);
	mpfq_fixmp_7_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 7));
      mpfq_fixmp_7_mgy_encode(xe, x, p);
      mpfq_fixmp_7_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_7_mul(s, x, y);
      mpfq_fixmp_7_mod(t, s, p);
      mpfq_fixmp_7_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_7_mul(s, xe, ye);
      mpfq_fixmp_7_mod(t, s, p);
      mpfq_fixmp_7_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 7) == 0);
      mpfq_fixmp_7_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 7) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_7_mul(t, xe, ye);
      s[14] = mpn_mul_1(s, t, 14, wx);
      mpfq_zero(t, 2*7+1);
      mpfq_fixmp_7_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_7_mod(s, t, p);
      mpfq_fixmp_7_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_7_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+7, v, 0, w, 7+1, p, 7);
      assert(mpn_cmp(v, s, 7) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 14 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 14, sat);
      mpfq_fixmp_7_redc(s, u, mip, p);
      mpfq_fixmp_7_mul1(w, s, sat);
      mpn_tdiv_qr(s+7, s, 0, w, 7+1, p, 7);
      mpfq_zero(w, 2*7+1);
      mpfq_fixmp_7_redc(w, v, mip, p);
      mpn_tdiv_qr(t+7, t, 0, w, 7+1, p, 7);
      assert(mpn_cmp(s, t, 7) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 14 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 14 + 1, sat);
      mpfq_fixmp_7_redc_ur(s, u, mip, p);
      mpfq_fixmp_7_mul1(w, s, sat);
      mpn_tdiv_qr(s+7, s, 0, w, 7+1, p, 7);
      mpfq_zero(w, 2*7+1);
      mpfq_fixmp_7_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+7, t, 0, w, 7+1, p, 7);
      assert(mpn_cmp(s, t, 7) == 0);
#ifdef  HAVE_native_7_mulredc
      mpfq_fixmp_7_mul(u, x, y);
      mpfq_fixmp_7_redc(s, u, mip, p);
      mpfq_fixmp_7_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 7) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 7, j);
    else
        mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_7_lshift(u, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 7, j);
    else
        mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_7_rshift(u, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // long_lshift
    j = wx % (7 * GMP_LIMB_BITS);
    mpfq_zero(s, 7);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 7 - k, j);
    else
        mpfq_copy(s + k, x, 7 - k);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_7_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // long_rshift
    j = wx % (7 * GMP_LIMB_BITS);
    mpfq_zero(s, 7);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 7 - k, j);
    else
        mpfq_copy(s, x + k, 7 - k);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_7_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 7) == 0);
}


void test_fixmp_8() {
  mp_limb_t s[17];
  mp_limb_t t[17];
  mp_limb_t u[17];
  mp_limb_t v[17];
  mp_limb_t w[17];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[8] = { 4093166397UL, 205402748UL, 1827875733UL, 2591432089UL, 498572719UL, 2575114975UL, 3040974997UL, 3977792999UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[8] = {16672777963903890445UL, 14321543724342516978UL, 5190009058579841038UL, 16894467406687282692UL, 5579682454395466331UL, 3120279582727612446UL, 2933066969036697885UL, 2125779597467003446UL };
#endif
  int j, k;
    // add
    u[8] = 0xdeadbeef;
    v[8] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 8);
    c2 = mpfq_fixmp_8_add(u, x, y);
    mpfq_fixmp_8_add_nc(v, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 8);
    c2 = mpfq_fixmp_8_add(u, x, x);
    mpfq_fixmp_8_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 8);
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c2 = mpfq_fixmp_8_add(v, v, y);
    mpfq_fixmp_8_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,8) == 0);
    assert (mpn_cmp(s,u,8) == 0);
    // add_ui
    u[8] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 8, wx);
    c2 = mpfq_fixmp_8_add_ui(u, x, wx);
    mpfq_fixmp_8_add_ui_nc(v, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c1 = mpn_add_1(s, x, 8, wx);
    c2 = mpfq_fixmp_8_add_ui(u, u, wx);
    mpfq_fixmp_8_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    // sub
    u[8] = 0xdeadbeef;
    v[8] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 8);
    c2 = mpfq_fixmp_8_sub(u, x, y);
    mpfq_fixmp_8_sub_nc(v, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 8);
    c2 = mpfq_fixmp_8_sub(u, x, x);
    mpfq_fixmp_8_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 8);
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c2 = mpfq_fixmp_8_sub(v, v, y);
    mpfq_fixmp_8_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,8) == 0);
    assert (mpn_cmp(s,u,8) == 0);
    // sub_ui
    u[8] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 8, wx);
    c2 = mpfq_fixmp_8_sub_ui(u, x, wx);
    mpfq_fixmp_8_sub_ui_nc(v, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c1 = mpn_sub_1(s, x, 8, wx);
    c2 = mpfq_fixmp_8_sub_ui(u, u, wx);
    mpfq_fixmp_8_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    // addmul1
    mpfq_copy(s, z, 9);
    mpfq_copy(u, z, 9);
    mpfq_copy(v, z, 9);
    mpfq_copy(w, z, 9);
    c1 = mpn_addmul_1(s, x, 8, wx);
    s[8] += c1;
    c3 = s[8] < c1;
    u[9]=0xdeadbeef;
    v[8]=0xdeadbeef;
    w[9]=0xdeadbeef;
    c4 = mpfq_fixmp_8_addmul1(u, x, wx);
    c2 = mpfq_fixmp_8_addmul1_shortz(v, x, wx);
    mpfq_fixmp_8_addmul1_nc(w, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (w[9] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 9) == 0);
    assert (mpn_cmp(s, v, 8) == 0);
    assert (mpn_cmp(s, w, 9) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    c1 = mpn_addmul_1(s, s, 8, wx);
    c2 = mpfq_fixmp_8_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 8) == 0);
    assert (c1 == c2);
    // mul1
    s[8] = mpn_mul_1(s, x, 8, wx);
    u[9] = 0xdeadbeef;
    mpfq_fixmp_8_mul1(u, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (mpn_cmp(s,u,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 8);
    mpfq_copy(v, x, 8);
    t[8] = mpn_mul_1(t, t, 8, wx);
    mpfq_fixmp_8_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 8 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 8);
    u[16]=0xdeadbeef;
    mpfq_fixmp_8_mul(u, x, y);
    assert (u[16]==0xdeadbeef);
    assert (mpn_cmp(s,u,16) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 8);
    u[8] = 0xdeadbeef;
    mpfq_fixmp_8_shortmul(u, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 8) == 0);
    // sqr
    mpn_mul_n(s, x, x, 8);
    u[16]=0xdeadbeef;
    mpfq_fixmp_8_sqr(u, x);
    assert (u[16]==0xdeadbeef);
    assert (mpn_cmp(s,u,16) == 0);
    // cmp
    j = mpn_cmp(x, y, 8);
    k = mpfq_fixmp_8_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 8);
    j = mpfq_fixmp_8_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 7);
    j = mpfq_fixmp_8_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_8_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_8_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 8);
    v[8-1] += !v[8-1];
    mpn_tdiv_qr(s+8, s, 0, z, 16, v, 8);
    mpfq_fixmp_8_mod(u, z, v);
    assert(mpn_cmp(s, u, 8) == 0);
    // inv
    mpfq_fixmp_8_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_8_invmod(v, u, P);
    mpfq_fixmp_8_invmod(v, v, P);
    assert(mpn_cmp(v, u, 8) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 8);
    mpfq_zero(v, 8);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_8_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_8_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 8);
    memset(u, ~0, 8 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 8);
    /* This creates an n-limb multiple of 257.  */
    v[8] = mpn_lshift(v, x, 8, 8);
    v[8] += mpfq_fixmp_8_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_8_add_ui(v, v, v[8]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_8_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 8, 257) == 0);
    //redc
    {
      mp_limb_t p[8], mip[8];
      mp_limb_t xe[8], ye[8], ze[8];
      mp_limb_t invR[8];
      
      // x[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 8);
      mpn_random2(p, 8);
      p[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[8-1] > p[8-1])
        p[8-1] = x[8-1];
      if (y[8-1] > p[8-1])
        p[8-1] = y[8-1];
      p[0] |= 1UL;
      p[8-1] += !p[8-1];
      mpfq_zero(w, 2*8);
      w[8]=1;
      mpfq_fixmp_8_mod(w, w, p);
      mpfq_fixmp_8_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 8);
        mpfq_fixmp_8_shortmul(t, mip, p);
	mpn_add_1(t, t, 8, 1);
	mpfq_fixmp_8_shortmul(u, t, mip);
	mpfq_fixmp_8_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 8));
      mpfq_fixmp_8_mgy_encode(xe, x, p);
      mpfq_fixmp_8_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_8_mul(s, x, y);
      mpfq_fixmp_8_mod(t, s, p);
      mpfq_fixmp_8_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_8_mul(s, xe, ye);
      mpfq_fixmp_8_mod(t, s, p);
      mpfq_fixmp_8_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 8) == 0);
      mpfq_fixmp_8_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 8) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_8_mul(t, xe, ye);
      s[16] = mpn_mul_1(s, t, 16, wx);
      mpfq_zero(t, 2*8+1);
      mpfq_fixmp_8_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_8_mod(s, t, p);
      mpfq_fixmp_8_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_8_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+8, v, 0, w, 8+1, p, 8);
      assert(mpn_cmp(v, s, 8) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 16 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 16, sat);
      mpfq_fixmp_8_redc(s, u, mip, p);
      mpfq_fixmp_8_mul1(w, s, sat);
      mpn_tdiv_qr(s+8, s, 0, w, 8+1, p, 8);
      mpfq_zero(w, 2*8+1);
      mpfq_fixmp_8_redc(w, v, mip, p);
      mpn_tdiv_qr(t+8, t, 0, w, 8+1, p, 8);
      assert(mpn_cmp(s, t, 8) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 16 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 16 + 1, sat);
      mpfq_fixmp_8_redc_ur(s, u, mip, p);
      mpfq_fixmp_8_mul1(w, s, sat);
      mpn_tdiv_qr(s+8, s, 0, w, 8+1, p, 8);
      mpfq_zero(w, 2*8+1);
      mpfq_fixmp_8_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+8, t, 0, w, 8+1, p, 8);
      assert(mpn_cmp(s, t, 8) == 0);
#ifdef  HAVE_native_8_mulredc
      mpfq_fixmp_8_mul(u, x, y);
      mpfq_fixmp_8_redc(s, u, mip, p);
      mpfq_fixmp_8_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 8) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 8, j);
    else
        mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_8_lshift(u, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 8, j);
    else
        mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_8_rshift(u, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // long_lshift
    j = wx % (8 * GMP_LIMB_BITS);
    mpfq_zero(s, 8);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 8 - k, j);
    else
        mpfq_copy(s + k, x, 8 - k);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_8_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // long_rshift
    j = wx % (8 * GMP_LIMB_BITS);
    mpfq_zero(s, 8);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 8 - k, j);
    else
        mpfq_copy(s, x + k, 8 - k);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_8_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 8) == 0);
}


void test_fixmp_9() {
  mp_limb_t s[19];
  mp_limb_t t[19];
  mp_limb_t u[19];
  mp_limb_t v[19];
  mp_limb_t w[19];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[9] = { 3897809411UL, 1993283498UL, 867915630UL, 886471665UL, 3987868346UL, 2967702854UL, 1194285669UL, 1588068146UL, 928806807UL, };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[9] = {6272071724723397689UL, 7097496403184731472UL, 6722451164852552420UL, 2557895735561628759UL, 11466998160538807963UL, 18232042263112599551UL, 4641538801156436724UL, 16426483130014462608UL, 7262099965674661736UL };
#endif
  int j, k;
    // add
    u[9] = 0xdeadbeef;
    v[9] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 9);
    c2 = mpfq_fixmp_9_add(u, x, y);
    mpfq_fixmp_9_add_nc(v, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 9);
    c2 = mpfq_fixmp_9_add(u, x, x);
    mpfq_fixmp_9_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 9);
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c2 = mpfq_fixmp_9_add(v, v, y);
    mpfq_fixmp_9_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,9) == 0);
    assert (mpn_cmp(s,u,9) == 0);
    // add_ui
    u[9] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 9, wx);
    c2 = mpfq_fixmp_9_add_ui(u, x, wx);
    mpfq_fixmp_9_add_ui_nc(v, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c1 = mpn_add_1(s, x, 9, wx);
    c2 = mpfq_fixmp_9_add_ui(u, u, wx);
    mpfq_fixmp_9_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    // sub
    u[9] = 0xdeadbeef;
    v[9] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 9);
    c2 = mpfq_fixmp_9_sub(u, x, y);
    mpfq_fixmp_9_sub_nc(v, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 9);
    c2 = mpfq_fixmp_9_sub(u, x, x);
    mpfq_fixmp_9_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 9);
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c2 = mpfq_fixmp_9_sub(v, v, y);
    mpfq_fixmp_9_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,9) == 0);
    assert (mpn_cmp(s,u,9) == 0);
    // sub_ui
    u[9] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 9, wx);
    c2 = mpfq_fixmp_9_sub_ui(u, x, wx);
    mpfq_fixmp_9_sub_ui_nc(v, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c1 = mpn_sub_1(s, x, 9, wx);
    c2 = mpfq_fixmp_9_sub_ui(u, u, wx);
    mpfq_fixmp_9_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    // addmul1
    mpfq_copy(s, z, 10);
    mpfq_copy(u, z, 10);
    mpfq_copy(v, z, 10);
    mpfq_copy(w, z, 10);
    c1 = mpn_addmul_1(s, x, 9, wx);
    s[9] += c1;
    c3 = s[9] < c1;
    u[10]=0xdeadbeef;
    v[9]=0xdeadbeef;
    w[10]=0xdeadbeef;
    c4 = mpfq_fixmp_9_addmul1(u, x, wx);
    c2 = mpfq_fixmp_9_addmul1_shortz(v, x, wx);
    mpfq_fixmp_9_addmul1_nc(w, x, wx);
    assert (u[10] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (w[10] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 10) == 0);
    assert (mpn_cmp(s, v, 9) == 0);
    assert (mpn_cmp(s, w, 10) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    c1 = mpn_addmul_1(s, s, 9, wx);
    c2 = mpfq_fixmp_9_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 9) == 0);
    assert (c1 == c2);
    // mul1
    s[9] = mpn_mul_1(s, x, 9, wx);
    u[10] = 0xdeadbeef;
    mpfq_fixmp_9_mul1(u, x, wx);
    assert (u[10] == 0xdeadbeef);
    assert (mpn_cmp(s,u,10) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 9);
    mpfq_copy(v, x, 9);
    t[9] = mpn_mul_1(t, t, 9, wx);
    mpfq_fixmp_9_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 9 + 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 9);
    u[18]=0xdeadbeef;
    mpfq_fixmp_9_mul(u, x, y);
    assert (u[18]==0xdeadbeef);
    assert (mpn_cmp(s,u,18) == 0);
    // shortmul
    mpn_mul_n(s, x, y, 9);
    u[9] = 0xdeadbeef;
    mpfq_fixmp_9_shortmul(u, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 9) == 0);
    // sqr
    mpn_mul_n(s, x, x, 9);
    u[18]=0xdeadbeef;
    mpfq_fixmp_9_sqr(u, x);
    assert (u[18]==0xdeadbeef);
    assert (mpn_cmp(s,u,18) == 0);
    // cmp
    j = mpn_cmp(x, y, 9);
    k = mpfq_fixmp_9_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 9);
    j = mpfq_fixmp_9_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 8);
    j = mpfq_fixmp_9_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_9_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_9_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 9);
    v[9-1] += !v[9-1];
    mpn_tdiv_qr(s+9, s, 0, z, 18, v, 9);
    mpfq_fixmp_9_mod(u, z, v);
    assert(mpn_cmp(s, u, 9) == 0);
    // inv
    mpfq_fixmp_9_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_9_invmod(v, u, P);
    mpfq_fixmp_9_invmod(v, v, P);
    assert(mpn_cmp(v, u, 9) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 9);
    mpfq_zero(v, 9);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_9_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_9_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 9);
    memset(u, ~0, 9 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 9);
    /* This creates an n-limb multiple of 257.  */
    v[9] = mpn_lshift(v, x, 9, 8);
    v[9] += mpfq_fixmp_9_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_9_add_ui(v, v, v[9]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_9_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 9, 257) == 0);
    //redc
    {
      mp_limb_t p[9], mip[9];
      mp_limb_t xe[9], ye[9], ze[9];
      mp_limb_t invR[9];
      
      // x[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 9);
      mpn_random2(p, 9);
      p[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[9-1] > p[9-1])
        p[9-1] = x[9-1];
      if (y[9-1] > p[9-1])
        p[9-1] = y[9-1];
      p[0] |= 1UL;
      p[9-1] += !p[9-1];
      mpfq_zero(w, 2*9);
      w[9]=1;
      mpfq_fixmp_9_mod(w, w, p);
      mpfq_fixmp_9_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 9);
        mpfq_fixmp_9_shortmul(t, mip, p);
	mpn_add_1(t, t, 9, 1);
	mpfq_fixmp_9_shortmul(u, t, mip);
	mpfq_fixmp_9_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 9));
      mpfq_fixmp_9_mgy_encode(xe, x, p);
      mpfq_fixmp_9_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_9_mul(s, x, y);
      mpfq_fixmp_9_mod(t, s, p);
      mpfq_fixmp_9_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_9_mul(s, xe, ye);
      mpfq_fixmp_9_mod(t, s, p);
      mpfq_fixmp_9_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 9) == 0);
      mpfq_fixmp_9_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 9) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_9_mul(t, xe, ye);
      s[18] = mpn_mul_1(s, t, 18, wx);
      mpfq_zero(t, 2*9+1);
      mpfq_fixmp_9_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_9_mod(s, t, p);
      mpfq_fixmp_9_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_9_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+9, v, 0, w, 9+1, p, 9);
      assert(mpn_cmp(v, s, 9) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 18 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 18, sat);
      mpfq_fixmp_9_redc(s, u, mip, p);
      mpfq_fixmp_9_mul1(w, s, sat);
      mpn_tdiv_qr(s+9, s, 0, w, 9+1, p, 9);
      mpfq_zero(w, 2*9+1);
      mpfq_fixmp_9_redc(w, v, mip, p);
      mpn_tdiv_qr(t+9, t, 0, w, 9+1, p, 9);
      assert(mpn_cmp(s, t, 9) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 18 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 18 + 1, sat);
      mpfq_fixmp_9_redc_ur(s, u, mip, p);
      mpfq_fixmp_9_mul1(w, s, sat);
      mpn_tdiv_qr(s+9, s, 0, w, 9+1, p, 9);
      mpfq_zero(w, 2*9+1);
      mpfq_fixmp_9_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+9, t, 0, w, 9+1, p, 9);
      assert(mpn_cmp(s, t, 9) == 0);
#ifdef  HAVE_native_9_mulredc
      mpfq_fixmp_9_mul(u, x, y);
      mpfq_fixmp_9_redc(s, u, mip, p);
      mpfq_fixmp_9_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 9) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 9, j);
    else
        mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_9_lshift(u, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 9, j);
    else
        mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_9_rshift(u, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // long_lshift
    j = wx % (9 * GMP_LIMB_BITS);
    mpfq_zero(s, 9);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 9 - k, j);
    else
        mpfq_copy(s + k, x, 9 - k);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_9_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // long_rshift
    j = wx % (9 * GMP_LIMB_BITS);
    mpfq_zero(s, 9);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 9 - k, j);
    else
        mpfq_copy(s, x + k, 9 - k);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_9_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 9) == 0);
}


void test_fixmp_0_5() {
  mp_limb_t s[3];
  mp_limb_t t[3];
  mp_limb_t u[3];
  mp_limb_t v[3];
  mp_limb_t w[3];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[1] = { 51157UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[1] = {665100797UL };
#endif
  int j, k;
    // add
    u[1] = 0xdeadbeef;
    v[1] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 1);
    c2 = mpfq_fixmp_0_5_add(u, x, y);
    mpfq_fixmp_0_5_add_nc(v, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 1);
    c2 = mpfq_fixmp_0_5_add(u, x, x);
    mpfq_fixmp_0_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 1);
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c2 = mpfq_fixmp_0_5_add(v, v, y);
    mpfq_fixmp_0_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,1) == 0);
    assert (mpn_cmp(s,u,1) == 0);
    // add_ui
    u[1] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 1, wx);
    c2 = mpfq_fixmp_0_5_add_ui(u, x, wx);
    mpfq_fixmp_0_5_add_ui_nc(v, x, wx);
    assert (u[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c1 = mpn_add_1(s, x, 1, wx);
    c2 = mpfq_fixmp_0_5_add_ui(u, u, wx);
    mpfq_fixmp_0_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    // sub
    u[1] = 0xdeadbeef;
    v[1] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 1);
    c2 = mpfq_fixmp_0_5_sub(u, x, y);
    mpfq_fixmp_0_5_sub_nc(v, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 1);
    c2 = mpfq_fixmp_0_5_sub(u, x, x);
    mpfq_fixmp_0_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 1);
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c2 = mpfq_fixmp_0_5_sub(v, v, y);
    mpfq_fixmp_0_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,1) == 0);
    assert (mpn_cmp(s,u,1) == 0);
    // sub_ui
    u[1] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 1, wx);
    c2 = mpfq_fixmp_0_5_sub_ui(u, x, wx);
    mpfq_fixmp_0_5_sub_ui_nc(v, x, wx);
    assert (u[1] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 1);
    mpfq_copy(v, x, 1);
    c1 = mpn_sub_1(s, x, 1, wx);
    c2 = mpfq_fixmp_0_5_sub_ui(u, u, wx);
    mpfq_fixmp_0_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,1) == 0);
    assert (mpn_cmp(s,v,1) == 0);
    // addmul1
    mpfq_copy(s, z, 2);
    mpfq_copy(u, z, 2);
    mpfq_copy(v, z, 2);
    mpfq_copy(w, z, 2);
    c1 = mpn_addmul_1(s, x, 1, wx);
    s[1] += c1;
    c3 = s[1] < c1;
    u[2]=0xdeadbeef;
    v[1]=0xdeadbeef;
    w[2]=0xdeadbeef;
    c4 = mpfq_fixmp_0_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_0_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_0_5_addmul1_nc(w, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (v[1] == 0xdeadbeef);
    assert (w[2] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 2) == 0);
    assert (mpn_cmp(s, v, 1) == 0);
    assert (mpn_cmp(s, w, 2) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    c1 = mpn_addmul_1(s, s, 1, wx);
    c2 = mpfq_fixmp_0_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 1) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 1);
        mpfq_copy(u, z, 1);
        mpfq_copy(w, z, 1);
        c1 = mpn_addmul_1(s, x, 1, y[0]);
        u[1]=0xdeadbeef;
        w[1]=0xdeadbeef;
        c2 = mpfq_fixmp_0_5_addmul05(u, x, y[0]);
        mpfq_fixmp_0_5_addmul05_nc(w, x, y[0]);
        assert (u[1] == 0xdeadbeef);
        assert (w[1] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 1) == 0);
        assert (mpn_cmp(s, w, 1) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 1);
        mpfq_copy(u, x, 1);
        c1 = mpn_addmul_1(s, s, 1, y[0]);
        c2 = mpfq_fixmp_0_5_addmul05(u, u, y[0]);
        assert (mpn_cmp(s, u, 1) == 0);
        assert (c1 == c2);
    // mul1
    s[1] = mpn_mul_1(s, x, 1, wx);
    u[2] = 0xdeadbeef;
    mpfq_fixmp_0_5_mul1(u, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (mpn_cmp(s,u,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 1);
    mpfq_copy(v, x, 1);
    t[1] = mpn_mul_1(t, t, 1, wx);
    mpfq_fixmp_0_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 1 + 1) == 0);
        // mul05
        s[1] = mpn_mul_1(s, x, 1, y[0]);
        u[1] = 0xdeadbeef;
        mpfq_fixmp_0_5_mul05(u, x, y[0]);
        assert (u[1] == 0xdeadbeef);
        assert (mpn_cmp(s,u,1) == 0);
        assert (s[1] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 1);
        mpfq_copy(v, x, 1);
        t[1] = mpn_mul_1(t, t, 1, y[0]);
        mpfq_fixmp_0_5_mul05(v, v, y[0]);
        assert (mpn_cmp(t, v, 1) == 0);
    // mul
    mpn_mul_n(s, x, y, 1);
    u[1]=0xdeadbeef;
    mpfq_fixmp_0_5_mul(u, x, y);
    assert (u[1]==0xdeadbeef);
    assert (mpn_cmp(s,u,1) == 0);
assert (s[1] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 1);
    u[1] = 0xdeadbeef;
    mpfq_fixmp_0_5_shortmul(u, x, y);
    assert (u[1] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 1) == 0);
    // sqr
    mpn_mul_n(s, x, x, 1);
    u[1]=0xdeadbeef;
    mpfq_fixmp_0_5_sqr(u, x);
    assert (u[1]==0xdeadbeef);
    assert (mpn_cmp(s,u,1) == 0);
assert (s[1] == 0);
    // cmp
    j = mpn_cmp(x, y, 1);
    k = mpfq_fixmp_0_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 1);
    j = mpfq_fixmp_0_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 0);
    j = mpfq_fixmp_0_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_0_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    // mod
    mpfq_copy(v, y, 1);
    v[1-1] += !v[1-1];
    mpn_tdiv_qr(s+1, s, 0, z, 1, v, 1);
    mpfq_fixmp_0_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 1) == 0);
    // inv
    mpfq_fixmp_0_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_0_5_invmod(v, u, P);
    mpfq_fixmp_0_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 1) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 1);
    mpfq_zero(v, 1);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_0_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_0_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 1);
    memset(u, ~0, 1 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 1);
    /* This creates an n-limb multiple of 257.  */
    v[1] = mpn_lshift(v, x, 1, 8);
    v[1] += mpfq_fixmp_0_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_0_5_add_ui(v, v, v[1]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_0_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 1, 257) == 0);
    //redc
    {
      mp_limb_t p[1], mip[1];
      mp_limb_t xe[1], ye[1], ze[1];
      mp_limb_t invR[1];
      
      // x[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 1);
      mpn_random2(p, 1);
      p[1-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[1-1] > p[1-1])
        p[1-1] = x[1-1];
      if (y[1-1] > p[1-1])
        p[1-1] = y[1-1];
      p[0] |= 1UL;
p[1-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[1-1] += !p[1-1];
        mpfq_zero(w, 2*1);
        w[0] = 1UL << (GMP_LIMB_BITS>>1);
        mpfq_fixmp_0_5_mod(w, w, p);
        mpfq_fixmp_0_5_mul(v, w, w);
        mpfq_fixmp_0_5_mod(w, v, p);
        mpfq_fixmp_0_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 1);
        mpfq_fixmp_0_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 1, 1);
	mpfq_fixmp_0_5_shortmul(u, t, mip);
	mpfq_fixmp_0_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 1));
      mpfq_fixmp_0_5_mgy_encode(xe, x, p);
      mpfq_fixmp_0_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_0_5_mul(s, x, y);
      mpfq_fixmp_0_5_mod(t, s, p);
      mpfq_fixmp_0_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_0_5_mul(s, xe, ye);
      mpfq_fixmp_0_5_mod(t, s, p);
      mpfq_fixmp_0_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 1) == 0);
      mpfq_fixmp_0_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 1) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_0_5_mul(t, xe, ye);
      s[1] = mpn_mul_1(s, t, 1, wx);
      mpfq_zero(t, 2*1+1);
      mpfq_fixmp_0_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_0_5_mod(s, t, p);
      mpfq_fixmp_0_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_0_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+1, v, 0, w, 1+1, p, 1);
      assert(mpn_cmp(v, s, 1) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 1, sat);
      mpfq_fixmp_0_5_redc(s, u, mip, p);
      mpfq_fixmp_0_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+1, s, 0, w, 1+1, p, 1);
      mpfq_zero(w, 2*1+1);
      mpfq_fixmp_0_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+1, t, 0, w, 1+1, p, 1);
      assert(mpn_cmp(s, t, 1) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 1 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 1 + 1, sat);
      mpfq_fixmp_0_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_0_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+1, s, 0, w, 1+1, p, 1);
      mpfq_zero(w, 2*1+1);
      mpfq_fixmp_0_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+1, t, 0, w, 1+1, p, 1);
      assert(mpn_cmp(s, t, 1) == 0);
#ifdef  HAVE_native_0_5_mulredc
      mpfq_fixmp_0_5_mul(u, x, y);
      mpfq_fixmp_0_5_redc(s, u, mip, p);
      mpfq_fixmp_0_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 1) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 1, j);
    else
        mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_0_5_lshift(u, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 1, j);
    else
        mpfq_copy(s, x, 1);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_0_5_rshift(u, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // long_lshift
    j = wx % (1 * GMP_LIMB_BITS);
    mpfq_zero(s, 1);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 1 - k, j);
    else
        mpfq_copy(s + k, x, 1 - k);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_0_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 1) == 0);
    // long_rshift
    j = wx % (1 * GMP_LIMB_BITS);
    mpfq_zero(s, 1);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 1 - k, j);
    else
        mpfq_copy(s, x + k, 1 - k);
    mpfq_copy(u, x, 1);
    mpfq_fixmp_0_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 1) == 0);
}


void test_fixmp_1_5() {
  mp_limb_t s[5];
  mp_limb_t t[5];
  mp_limb_t u[5];
  mp_limb_t v[5];
  mp_limb_t w[5];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[2] = { 3980356625UL, 58899UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[2] = {3107852032399716593UL, 561560781UL };
#endif
  int j, k;
    // add
    u[2] = 0xdeadbeef;
    v[2] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 2);
    c2 = mpfq_fixmp_1_5_add(u, x, y);
    mpfq_fixmp_1_5_add_nc(v, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 2);
    c2 = mpfq_fixmp_1_5_add(u, x, x);
    mpfq_fixmp_1_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 2);
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c2 = mpfq_fixmp_1_5_add(v, v, y);
    mpfq_fixmp_1_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,2) == 0);
    assert (mpn_cmp(s,u,2) == 0);
    // add_ui
    u[2] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 2, wx);
    c2 = mpfq_fixmp_1_5_add_ui(u, x, wx);
    mpfq_fixmp_1_5_add_ui_nc(v, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c1 = mpn_add_1(s, x, 2, wx);
    c2 = mpfq_fixmp_1_5_add_ui(u, u, wx);
    mpfq_fixmp_1_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    // sub
    u[2] = 0xdeadbeef;
    v[2] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 2);
    c2 = mpfq_fixmp_1_5_sub(u, x, y);
    mpfq_fixmp_1_5_sub_nc(v, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 2);
    c2 = mpfq_fixmp_1_5_sub(u, x, x);
    mpfq_fixmp_1_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 2);
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c2 = mpfq_fixmp_1_5_sub(v, v, y);
    mpfq_fixmp_1_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,2) == 0);
    assert (mpn_cmp(s,u,2) == 0);
    // sub_ui
    u[2] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 2, wx);
    c2 = mpfq_fixmp_1_5_sub_ui(u, x, wx);
    mpfq_fixmp_1_5_sub_ui_nc(v, x, wx);
    assert (u[2] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 2);
    mpfq_copy(v, x, 2);
    c1 = mpn_sub_1(s, x, 2, wx);
    c2 = mpfq_fixmp_1_5_sub_ui(u, u, wx);
    mpfq_fixmp_1_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,2) == 0);
    assert (mpn_cmp(s,v,2) == 0);
    // addmul1
    mpfq_copy(s, z, 3);
    mpfq_copy(u, z, 3);
    mpfq_copy(v, z, 3);
    mpfq_copy(w, z, 3);
    c1 = mpn_addmul_1(s, x, 2, wx);
    s[2] += c1;
    c3 = s[2] < c1;
    u[3]=0xdeadbeef;
    v[2]=0xdeadbeef;
    w[3]=0xdeadbeef;
    c4 = mpfq_fixmp_1_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_1_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_1_5_addmul1_nc(w, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (v[2] == 0xdeadbeef);
    assert (w[3] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 3) == 0);
    assert (mpn_cmp(s, v, 2) == 0);
    assert (mpn_cmp(s, w, 3) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    c1 = mpn_addmul_1(s, s, 2, wx);
    c2 = mpfq_fixmp_1_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 2) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 2);
        mpfq_copy(u, z, 2);
        mpfq_copy(w, z, 2);
        c1 = mpn_addmul_1(s, x, 2, y[1]);
        u[2]=0xdeadbeef;
        w[2]=0xdeadbeef;
        c2 = mpfq_fixmp_1_5_addmul05(u, x, y[1]);
        mpfq_fixmp_1_5_addmul05_nc(w, x, y[1]);
        assert (u[2] == 0xdeadbeef);
        assert (w[2] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 2) == 0);
        assert (mpn_cmp(s, w, 2) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 2);
        mpfq_copy(u, x, 2);
        c1 = mpn_addmul_1(s, s, 2, y[1]);
        c2 = mpfq_fixmp_1_5_addmul05(u, u, y[1]);
        assert (mpn_cmp(s, u, 2) == 0);
        assert (c1 == c2);
    // mul1
    s[2] = mpn_mul_1(s, x, 2, wx);
    u[3] = 0xdeadbeef;
    mpfq_fixmp_1_5_mul1(u, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (mpn_cmp(s,u,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 2);
    mpfq_copy(v, x, 2);
    t[2] = mpn_mul_1(t, t, 2, wx);
    mpfq_fixmp_1_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 2 + 1) == 0);
        // mul05
        s[2] = mpn_mul_1(s, x, 2, y[1]);
        u[2] = 0xdeadbeef;
        mpfq_fixmp_1_5_mul05(u, x, y[1]);
        assert (u[2] == 0xdeadbeef);
        assert (mpn_cmp(s,u,2) == 0);
        assert (s[2] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 2);
        mpfq_copy(v, x, 2);
        t[2] = mpn_mul_1(t, t, 2, y[1]);
        mpfq_fixmp_1_5_mul05(v, v, y[1]);
        assert (mpn_cmp(t, v, 2) == 0);
    // mul
    mpn_mul_n(s, x, y, 2);
    u[3]=0xdeadbeef;
    mpfq_fixmp_1_5_mul(u, x, y);
    assert (u[3]==0xdeadbeef);
    assert (mpn_cmp(s,u,3) == 0);
assert (s[3] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 2);
    u[2] = 0xdeadbeef;
    mpfq_fixmp_1_5_shortmul(u, x, y);
    assert (u[2] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 2) == 0);
    // sqr
    mpn_mul_n(s, x, x, 2);
    u[3]=0xdeadbeef;
    mpfq_fixmp_1_5_sqr(u, x);
    assert (u[3]==0xdeadbeef);
    assert (mpn_cmp(s,u,3) == 0);
assert (s[3] == 0);
    // cmp
    j = mpn_cmp(x, y, 2);
    k = mpfq_fixmp_1_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 2);
    j = mpfq_fixmp_1_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 1);
    j = mpfq_fixmp_1_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_1_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_1_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 2);
    v[2-1] += !v[2-1];
    mpn_tdiv_qr(s+2, s, 0, z, 3, v, 2);
    mpfq_fixmp_1_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 2) == 0);
    // inv
    mpfq_fixmp_1_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_1_5_invmod(v, u, P);
    mpfq_fixmp_1_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 2) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 2);
    mpfq_zero(v, 2);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_1_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_1_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 2);
    memset(u, ~0, 2 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 2);
    /* This creates an n-limb multiple of 257.  */
    v[2] = mpn_lshift(v, x, 2, 8);
    v[2] += mpfq_fixmp_1_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_1_5_add_ui(v, v, v[2]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_1_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 2, 257) == 0);
    //redc
    {
      mp_limb_t p[2], mip[2];
      mp_limb_t xe[2], ye[2], ze[2];
      mp_limb_t invR[2];
      
      // x[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 2);
      mpn_random2(p, 2);
      p[2-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[2-1] > p[2-1])
        p[2-1] = x[2-1];
      if (y[2-1] > p[2-1])
        p[2-1] = y[2-1];
      p[0] |= 1UL;
p[2-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[2-1] += !p[2-1];
      mpfq_zero(w, 2*2);
      w[2]=1;
      mpfq_fixmp_1_5_mod(w, w, p);
      mpfq_fixmp_1_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 2);
        mpfq_fixmp_1_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 2, 1);
	mpfq_fixmp_1_5_shortmul(u, t, mip);
	mpfq_fixmp_1_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 2));
      mpfq_fixmp_1_5_mgy_encode(xe, x, p);
      mpfq_fixmp_1_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_1_5_mul(s, x, y);
      mpfq_fixmp_1_5_mod(t, s, p);
      mpfq_fixmp_1_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_1_5_mul(s, xe, ye);
      mpfq_fixmp_1_5_mod(t, s, p);
      mpfq_fixmp_1_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 2) == 0);
      mpfq_fixmp_1_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 2) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_1_5_mul(t, xe, ye);
      s[3] = mpn_mul_1(s, t, 3, wx);
      mpfq_zero(t, 2*2+1);
      mpfq_fixmp_1_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_1_5_mod(s, t, p);
      mpfq_fixmp_1_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_1_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+2, v, 0, w, 2+1, p, 2);
      assert(mpn_cmp(v, s, 2) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 3 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 3, sat);
      mpfq_fixmp_1_5_redc(s, u, mip, p);
      mpfq_fixmp_1_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+2, s, 0, w, 2+1, p, 2);
      mpfq_zero(w, 2*2+1);
      mpfq_fixmp_1_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+2, t, 0, w, 2+1, p, 2);
      assert(mpn_cmp(s, t, 2) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 3 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 3 + 1, sat);
      mpfq_fixmp_1_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_1_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+2, s, 0, w, 2+1, p, 2);
      mpfq_zero(w, 2*2+1);
      mpfq_fixmp_1_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+2, t, 0, w, 2+1, p, 2);
      assert(mpn_cmp(s, t, 2) == 0);
#ifdef  HAVE_native_1_5_mulredc
      mpfq_fixmp_1_5_mul(u, x, y);
      mpfq_fixmp_1_5_redc(s, u, mip, p);
      mpfq_fixmp_1_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 2) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 2, j);
    else
        mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_1_5_lshift(u, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 2, j);
    else
        mpfq_copy(s, x, 2);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_1_5_rshift(u, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // long_lshift
    j = wx % (2 * GMP_LIMB_BITS);
    mpfq_zero(s, 2);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 2 - k, j);
    else
        mpfq_copy(s + k, x, 2 - k);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_1_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 2) == 0);
    // long_rshift
    j = wx % (2 * GMP_LIMB_BITS);
    mpfq_zero(s, 2);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 2 - k, j);
    else
        mpfq_copy(s, x + k, 2 - k);
    mpfq_copy(u, x, 2);
    mpfq_fixmp_1_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 2) == 0);
}


void test_fixmp_2_5() {
  mp_limb_t s[7];
  mp_limb_t t[7];
  mp_limb_t u[7];
  mp_limb_t v[7];
  mp_limb_t w[7];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[3] = { 2705370375UL, 3814976481UL, 12767UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[3] = {3022170862732092607UL, 2943310936735432501UL, 2389540284UL };
#endif
  int j, k;
    // add
    u[3] = 0xdeadbeef;
    v[3] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 3);
    c2 = mpfq_fixmp_2_5_add(u, x, y);
    mpfq_fixmp_2_5_add_nc(v, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 3);
    c2 = mpfq_fixmp_2_5_add(u, x, x);
    mpfq_fixmp_2_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 3);
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c2 = mpfq_fixmp_2_5_add(v, v, y);
    mpfq_fixmp_2_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,3) == 0);
    assert (mpn_cmp(s,u,3) == 0);
    // add_ui
    u[3] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 3, wx);
    c2 = mpfq_fixmp_2_5_add_ui(u, x, wx);
    mpfq_fixmp_2_5_add_ui_nc(v, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c1 = mpn_add_1(s, x, 3, wx);
    c2 = mpfq_fixmp_2_5_add_ui(u, u, wx);
    mpfq_fixmp_2_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    // sub
    u[3] = 0xdeadbeef;
    v[3] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 3);
    c2 = mpfq_fixmp_2_5_sub(u, x, y);
    mpfq_fixmp_2_5_sub_nc(v, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 3);
    c2 = mpfq_fixmp_2_5_sub(u, x, x);
    mpfq_fixmp_2_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 3);
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c2 = mpfq_fixmp_2_5_sub(v, v, y);
    mpfq_fixmp_2_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,3) == 0);
    assert (mpn_cmp(s,u,3) == 0);
    // sub_ui
    u[3] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 3, wx);
    c2 = mpfq_fixmp_2_5_sub_ui(u, x, wx);
    mpfq_fixmp_2_5_sub_ui_nc(v, x, wx);
    assert (u[3] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 3);
    mpfq_copy(v, x, 3);
    c1 = mpn_sub_1(s, x, 3, wx);
    c2 = mpfq_fixmp_2_5_sub_ui(u, u, wx);
    mpfq_fixmp_2_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,3) == 0);
    assert (mpn_cmp(s,v,3) == 0);
    // addmul1
    mpfq_copy(s, z, 4);
    mpfq_copy(u, z, 4);
    mpfq_copy(v, z, 4);
    mpfq_copy(w, z, 4);
    c1 = mpn_addmul_1(s, x, 3, wx);
    s[3] += c1;
    c3 = s[3] < c1;
    u[4]=0xdeadbeef;
    v[3]=0xdeadbeef;
    w[4]=0xdeadbeef;
    c4 = mpfq_fixmp_2_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_2_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_2_5_addmul1_nc(w, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (v[3] == 0xdeadbeef);
    assert (w[4] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 4) == 0);
    assert (mpn_cmp(s, v, 3) == 0);
    assert (mpn_cmp(s, w, 4) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    c1 = mpn_addmul_1(s, s, 3, wx);
    c2 = mpfq_fixmp_2_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 3) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 3);
        mpfq_copy(u, z, 3);
        mpfq_copy(w, z, 3);
        c1 = mpn_addmul_1(s, x, 3, y[2]);
        u[3]=0xdeadbeef;
        w[3]=0xdeadbeef;
        c2 = mpfq_fixmp_2_5_addmul05(u, x, y[2]);
        mpfq_fixmp_2_5_addmul05_nc(w, x, y[2]);
        assert (u[3] == 0xdeadbeef);
        assert (w[3] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 3) == 0);
        assert (mpn_cmp(s, w, 3) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 3);
        mpfq_copy(u, x, 3);
        c1 = mpn_addmul_1(s, s, 3, y[2]);
        c2 = mpfq_fixmp_2_5_addmul05(u, u, y[2]);
        assert (mpn_cmp(s, u, 3) == 0);
        assert (c1 == c2);
    // mul1
    s[3] = mpn_mul_1(s, x, 3, wx);
    u[4] = 0xdeadbeef;
    mpfq_fixmp_2_5_mul1(u, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (mpn_cmp(s,u,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 3);
    mpfq_copy(v, x, 3);
    t[3] = mpn_mul_1(t, t, 3, wx);
    mpfq_fixmp_2_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 3 + 1) == 0);
        // mul05
        s[3] = mpn_mul_1(s, x, 3, y[2]);
        u[3] = 0xdeadbeef;
        mpfq_fixmp_2_5_mul05(u, x, y[2]);
        assert (u[3] == 0xdeadbeef);
        assert (mpn_cmp(s,u,3) == 0);
        assert (s[3] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 3);
        mpfq_copy(v, x, 3);
        t[3] = mpn_mul_1(t, t, 3, y[2]);
        mpfq_fixmp_2_5_mul05(v, v, y[2]);
        assert (mpn_cmp(t, v, 3) == 0);
    // mul
    mpn_mul_n(s, x, y, 3);
    u[5]=0xdeadbeef;
    mpfq_fixmp_2_5_mul(u, x, y);
    assert (u[5]==0xdeadbeef);
    assert (mpn_cmp(s,u,5) == 0);
assert (s[5] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 3);
    u[3] = 0xdeadbeef;
    mpfq_fixmp_2_5_shortmul(u, x, y);
    assert (u[3] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 3) == 0);
    // sqr
    mpn_mul_n(s, x, x, 3);
    u[5]=0xdeadbeef;
    mpfq_fixmp_2_5_sqr(u, x);
    assert (u[5]==0xdeadbeef);
    assert (mpn_cmp(s,u,5) == 0);
assert (s[5] == 0);
    // cmp
    j = mpn_cmp(x, y, 3);
    k = mpfq_fixmp_2_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 3);
    j = mpfq_fixmp_2_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 2);
    j = mpfq_fixmp_2_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_2_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_2_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 3);
    v[3-1] += !v[3-1];
    mpn_tdiv_qr(s+3, s, 0, z, 5, v, 3);
    mpfq_fixmp_2_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 3) == 0);
    // inv
    mpfq_fixmp_2_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_2_5_invmod(v, u, P);
    mpfq_fixmp_2_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 3) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 3);
    mpfq_zero(v, 3);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_2_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_2_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 3);
    memset(u, ~0, 3 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 3);
    /* This creates an n-limb multiple of 257.  */
    v[3] = mpn_lshift(v, x, 3, 8);
    v[3] += mpfq_fixmp_2_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_2_5_add_ui(v, v, v[3]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_2_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 3, 257) == 0);
    //redc
    {
      mp_limb_t p[3], mip[3];
      mp_limb_t xe[3], ye[3], ze[3];
      mp_limb_t invR[3];
      
      // x[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 3);
      mpn_random2(p, 3);
      p[3-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[3-1] > p[3-1])
        p[3-1] = x[3-1];
      if (y[3-1] > p[3-1])
        p[3-1] = y[3-1];
      p[0] |= 1UL;
p[3-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[3-1] += !p[3-1];
      mpfq_zero(w, 2*3);
      w[3]=1;
      mpfq_fixmp_2_5_mod(w, w, p);
      mpfq_fixmp_2_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 3);
        mpfq_fixmp_2_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 3, 1);
	mpfq_fixmp_2_5_shortmul(u, t, mip);
	mpfq_fixmp_2_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 3));
      mpfq_fixmp_2_5_mgy_encode(xe, x, p);
      mpfq_fixmp_2_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_2_5_mul(s, x, y);
      mpfq_fixmp_2_5_mod(t, s, p);
      mpfq_fixmp_2_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_2_5_mul(s, xe, ye);
      mpfq_fixmp_2_5_mod(t, s, p);
      mpfq_fixmp_2_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 3) == 0);
      mpfq_fixmp_2_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 3) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_2_5_mul(t, xe, ye);
      s[5] = mpn_mul_1(s, t, 5, wx);
      mpfq_zero(t, 2*3+1);
      mpfq_fixmp_2_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_2_5_mod(s, t, p);
      mpfq_fixmp_2_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_2_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+3, v, 0, w, 3+1, p, 3);
      assert(mpn_cmp(v, s, 3) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 5 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 5, sat);
      mpfq_fixmp_2_5_redc(s, u, mip, p);
      mpfq_fixmp_2_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+3, s, 0, w, 3+1, p, 3);
      mpfq_zero(w, 2*3+1);
      mpfq_fixmp_2_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+3, t, 0, w, 3+1, p, 3);
      assert(mpn_cmp(s, t, 3) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 5 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 5 + 1, sat);
      mpfq_fixmp_2_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_2_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+3, s, 0, w, 3+1, p, 3);
      mpfq_zero(w, 2*3+1);
      mpfq_fixmp_2_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+3, t, 0, w, 3+1, p, 3);
      assert(mpn_cmp(s, t, 3) == 0);
#ifdef  HAVE_native_2_5_mulredc
      mpfq_fixmp_2_5_mul(u, x, y);
      mpfq_fixmp_2_5_redc(s, u, mip, p);
      mpfq_fixmp_2_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 3) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 3, j);
    else
        mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_2_5_lshift(u, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 3, j);
    else
        mpfq_copy(s, x, 3);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_2_5_rshift(u, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // long_lshift
    j = wx % (3 * GMP_LIMB_BITS);
    mpfq_zero(s, 3);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 3 - k, j);
    else
        mpfq_copy(s + k, x, 3 - k);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_2_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 3) == 0);
    // long_rshift
    j = wx % (3 * GMP_LIMB_BITS);
    mpfq_zero(s, 3);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 3 - k, j);
    else
        mpfq_copy(s, x + k, 3 - k);
    mpfq_copy(u, x, 3);
    mpfq_fixmp_2_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 3) == 0);
}


void test_fixmp_3_5() {
  mp_limb_t s[9];
  mp_limb_t t[9];
  mp_limb_t u[9];
  mp_limb_t v[9];
  mp_limb_t w[9];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[4] = { 1824034149UL, 3737757448UL, 3181502520UL, 10185UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[4] = {3365129921614404491UL, 3435152464474079436UL, 15222848675840898267UL, 757047498UL };
#endif
  int j, k;
    // add
    u[4] = 0xdeadbeef;
    v[4] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 4);
    c2 = mpfq_fixmp_3_5_add(u, x, y);
    mpfq_fixmp_3_5_add_nc(v, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 4);
    c2 = mpfq_fixmp_3_5_add(u, x, x);
    mpfq_fixmp_3_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 4);
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c2 = mpfq_fixmp_3_5_add(v, v, y);
    mpfq_fixmp_3_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,4) == 0);
    assert (mpn_cmp(s,u,4) == 0);
    // add_ui
    u[4] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 4, wx);
    c2 = mpfq_fixmp_3_5_add_ui(u, x, wx);
    mpfq_fixmp_3_5_add_ui_nc(v, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c1 = mpn_add_1(s, x, 4, wx);
    c2 = mpfq_fixmp_3_5_add_ui(u, u, wx);
    mpfq_fixmp_3_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    // sub
    u[4] = 0xdeadbeef;
    v[4] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 4);
    c2 = mpfq_fixmp_3_5_sub(u, x, y);
    mpfq_fixmp_3_5_sub_nc(v, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 4);
    c2 = mpfq_fixmp_3_5_sub(u, x, x);
    mpfq_fixmp_3_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 4);
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c2 = mpfq_fixmp_3_5_sub(v, v, y);
    mpfq_fixmp_3_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,4) == 0);
    assert (mpn_cmp(s,u,4) == 0);
    // sub_ui
    u[4] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 4, wx);
    c2 = mpfq_fixmp_3_5_sub_ui(u, x, wx);
    mpfq_fixmp_3_5_sub_ui_nc(v, x, wx);
    assert (u[4] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 4);
    mpfq_copy(v, x, 4);
    c1 = mpn_sub_1(s, x, 4, wx);
    c2 = mpfq_fixmp_3_5_sub_ui(u, u, wx);
    mpfq_fixmp_3_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,4) == 0);
    assert (mpn_cmp(s,v,4) == 0);
    // addmul1
    mpfq_copy(s, z, 5);
    mpfq_copy(u, z, 5);
    mpfq_copy(v, z, 5);
    mpfq_copy(w, z, 5);
    c1 = mpn_addmul_1(s, x, 4, wx);
    s[4] += c1;
    c3 = s[4] < c1;
    u[5]=0xdeadbeef;
    v[4]=0xdeadbeef;
    w[5]=0xdeadbeef;
    c4 = mpfq_fixmp_3_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_3_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_3_5_addmul1_nc(w, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (v[4] == 0xdeadbeef);
    assert (w[5] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 5) == 0);
    assert (mpn_cmp(s, v, 4) == 0);
    assert (mpn_cmp(s, w, 5) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    c1 = mpn_addmul_1(s, s, 4, wx);
    c2 = mpfq_fixmp_3_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 4) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 4);
        mpfq_copy(u, z, 4);
        mpfq_copy(w, z, 4);
        c1 = mpn_addmul_1(s, x, 4, y[3]);
        u[4]=0xdeadbeef;
        w[4]=0xdeadbeef;
        c2 = mpfq_fixmp_3_5_addmul05(u, x, y[3]);
        mpfq_fixmp_3_5_addmul05_nc(w, x, y[3]);
        assert (u[4] == 0xdeadbeef);
        assert (w[4] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 4) == 0);
        assert (mpn_cmp(s, w, 4) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 4);
        mpfq_copy(u, x, 4);
        c1 = mpn_addmul_1(s, s, 4, y[3]);
        c2 = mpfq_fixmp_3_5_addmul05(u, u, y[3]);
        assert (mpn_cmp(s, u, 4) == 0);
        assert (c1 == c2);
    // mul1
    s[4] = mpn_mul_1(s, x, 4, wx);
    u[5] = 0xdeadbeef;
    mpfq_fixmp_3_5_mul1(u, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (mpn_cmp(s,u,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 4);
    mpfq_copy(v, x, 4);
    t[4] = mpn_mul_1(t, t, 4, wx);
    mpfq_fixmp_3_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 4 + 1) == 0);
        // mul05
        s[4] = mpn_mul_1(s, x, 4, y[3]);
        u[4] = 0xdeadbeef;
        mpfq_fixmp_3_5_mul05(u, x, y[3]);
        assert (u[4] == 0xdeadbeef);
        assert (mpn_cmp(s,u,4) == 0);
        assert (s[4] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 4);
        mpfq_copy(v, x, 4);
        t[4] = mpn_mul_1(t, t, 4, y[3]);
        mpfq_fixmp_3_5_mul05(v, v, y[3]);
        assert (mpn_cmp(t, v, 4) == 0);
    // mul
    mpn_mul_n(s, x, y, 4);
    u[7]=0xdeadbeef;
    mpfq_fixmp_3_5_mul(u, x, y);
    assert (u[7]==0xdeadbeef);
    assert (mpn_cmp(s,u,7) == 0);
assert (s[7] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 4);
    u[4] = 0xdeadbeef;
    mpfq_fixmp_3_5_shortmul(u, x, y);
    assert (u[4] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 4) == 0);
    // sqr
    mpn_mul_n(s, x, x, 4);
    u[7]=0xdeadbeef;
    mpfq_fixmp_3_5_sqr(u, x);
    assert (u[7]==0xdeadbeef);
    assert (mpn_cmp(s,u,7) == 0);
assert (s[7] == 0);
    // cmp
    j = mpn_cmp(x, y, 4);
    k = mpfq_fixmp_3_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 4);
    j = mpfq_fixmp_3_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 3);
    j = mpfq_fixmp_3_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_3_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_3_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 4);
    v[4-1] += !v[4-1];
    mpn_tdiv_qr(s+4, s, 0, z, 7, v, 4);
    mpfq_fixmp_3_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 4) == 0);
    // inv
    mpfq_fixmp_3_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_3_5_invmod(v, u, P);
    mpfq_fixmp_3_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 4) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 4);
    mpfq_zero(v, 4);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_3_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_3_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 4);
    memset(u, ~0, 4 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 4);
    /* This creates an n-limb multiple of 257.  */
    v[4] = mpn_lshift(v, x, 4, 8);
    v[4] += mpfq_fixmp_3_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_3_5_add_ui(v, v, v[4]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_3_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 4, 257) == 0);
    //redc
    {
      mp_limb_t p[4], mip[4];
      mp_limb_t xe[4], ye[4], ze[4];
      mp_limb_t invR[4];
      
      // x[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 4);
      mpn_random2(p, 4);
      p[4-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[4-1] > p[4-1])
        p[4-1] = x[4-1];
      if (y[4-1] > p[4-1])
        p[4-1] = y[4-1];
      p[0] |= 1UL;
p[4-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[4-1] += !p[4-1];
      mpfq_zero(w, 2*4);
      w[4]=1;
      mpfq_fixmp_3_5_mod(w, w, p);
      mpfq_fixmp_3_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 4);
        mpfq_fixmp_3_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 4, 1);
	mpfq_fixmp_3_5_shortmul(u, t, mip);
	mpfq_fixmp_3_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 4));
      mpfq_fixmp_3_5_mgy_encode(xe, x, p);
      mpfq_fixmp_3_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_3_5_mul(s, x, y);
      mpfq_fixmp_3_5_mod(t, s, p);
      mpfq_fixmp_3_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_3_5_mul(s, xe, ye);
      mpfq_fixmp_3_5_mod(t, s, p);
      mpfq_fixmp_3_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 4) == 0);
      mpfq_fixmp_3_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 4) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_3_5_mul(t, xe, ye);
      s[7] = mpn_mul_1(s, t, 7, wx);
      mpfq_zero(t, 2*4+1);
      mpfq_fixmp_3_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_3_5_mod(s, t, p);
      mpfq_fixmp_3_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_3_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+4, v, 0, w, 4+1, p, 4);
      assert(mpn_cmp(v, s, 4) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 7 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 7, sat);
      mpfq_fixmp_3_5_redc(s, u, mip, p);
      mpfq_fixmp_3_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+4, s, 0, w, 4+1, p, 4);
      mpfq_zero(w, 2*4+1);
      mpfq_fixmp_3_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+4, t, 0, w, 4+1, p, 4);
      assert(mpn_cmp(s, t, 4) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 7 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 7 + 1, sat);
      mpfq_fixmp_3_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_3_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+4, s, 0, w, 4+1, p, 4);
      mpfq_zero(w, 2*4+1);
      mpfq_fixmp_3_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+4, t, 0, w, 4+1, p, 4);
      assert(mpn_cmp(s, t, 4) == 0);
#ifdef  HAVE_native_3_5_mulredc
      mpfq_fixmp_3_5_mul(u, x, y);
      mpfq_fixmp_3_5_redc(s, u, mip, p);
      mpfq_fixmp_3_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 4) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 4, j);
    else
        mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_3_5_lshift(u, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 4, j);
    else
        mpfq_copy(s, x, 4);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_3_5_rshift(u, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // long_lshift
    j = wx % (4 * GMP_LIMB_BITS);
    mpfq_zero(s, 4);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 4 - k, j);
    else
        mpfq_copy(s + k, x, 4 - k);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_3_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 4) == 0);
    // long_rshift
    j = wx % (4 * GMP_LIMB_BITS);
    mpfq_zero(s, 4);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 4 - k, j);
    else
        mpfq_copy(s, x + k, 4 - k);
    mpfq_copy(u, x, 4);
    mpfq_fixmp_3_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 4) == 0);
}


void test_fixmp_4_5() {
  mp_limb_t s[11];
  mp_limb_t t[11];
  mp_limb_t u[11];
  mp_limb_t v[11];
  mp_limb_t w[11];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[5] = { 855954195UL, 1867864209UL, 238425956UL, 1079706229UL, 46658UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[5] = {7237051979824649275UL, 4624810970717405228UL, 1884741607445120980UL, 13362986272152613141UL, 4289448930UL };
#endif
  int j, k;
    // add
    u[5] = 0xdeadbeef;
    v[5] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 5);
    c2 = mpfq_fixmp_4_5_add(u, x, y);
    mpfq_fixmp_4_5_add_nc(v, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 5);
    c2 = mpfq_fixmp_4_5_add(u, x, x);
    mpfq_fixmp_4_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 5);
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c2 = mpfq_fixmp_4_5_add(v, v, y);
    mpfq_fixmp_4_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,5) == 0);
    assert (mpn_cmp(s,u,5) == 0);
    // add_ui
    u[5] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 5, wx);
    c2 = mpfq_fixmp_4_5_add_ui(u, x, wx);
    mpfq_fixmp_4_5_add_ui_nc(v, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c1 = mpn_add_1(s, x, 5, wx);
    c2 = mpfq_fixmp_4_5_add_ui(u, u, wx);
    mpfq_fixmp_4_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    // sub
    u[5] = 0xdeadbeef;
    v[5] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 5);
    c2 = mpfq_fixmp_4_5_sub(u, x, y);
    mpfq_fixmp_4_5_sub_nc(v, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 5);
    c2 = mpfq_fixmp_4_5_sub(u, x, x);
    mpfq_fixmp_4_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 5);
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c2 = mpfq_fixmp_4_5_sub(v, v, y);
    mpfq_fixmp_4_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,5) == 0);
    assert (mpn_cmp(s,u,5) == 0);
    // sub_ui
    u[5] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 5, wx);
    c2 = mpfq_fixmp_4_5_sub_ui(u, x, wx);
    mpfq_fixmp_4_5_sub_ui_nc(v, x, wx);
    assert (u[5] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 5);
    mpfq_copy(v, x, 5);
    c1 = mpn_sub_1(s, x, 5, wx);
    c2 = mpfq_fixmp_4_5_sub_ui(u, u, wx);
    mpfq_fixmp_4_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,5) == 0);
    assert (mpn_cmp(s,v,5) == 0);
    // addmul1
    mpfq_copy(s, z, 6);
    mpfq_copy(u, z, 6);
    mpfq_copy(v, z, 6);
    mpfq_copy(w, z, 6);
    c1 = mpn_addmul_1(s, x, 5, wx);
    s[5] += c1;
    c3 = s[5] < c1;
    u[6]=0xdeadbeef;
    v[5]=0xdeadbeef;
    w[6]=0xdeadbeef;
    c4 = mpfq_fixmp_4_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_4_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_4_5_addmul1_nc(w, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (v[5] == 0xdeadbeef);
    assert (w[6] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 6) == 0);
    assert (mpn_cmp(s, v, 5) == 0);
    assert (mpn_cmp(s, w, 6) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    c1 = mpn_addmul_1(s, s, 5, wx);
    c2 = mpfq_fixmp_4_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 5) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 5);
        mpfq_copy(u, z, 5);
        mpfq_copy(w, z, 5);
        c1 = mpn_addmul_1(s, x, 5, y[4]);
        u[5]=0xdeadbeef;
        w[5]=0xdeadbeef;
        c2 = mpfq_fixmp_4_5_addmul05(u, x, y[4]);
        mpfq_fixmp_4_5_addmul05_nc(w, x, y[4]);
        assert (u[5] == 0xdeadbeef);
        assert (w[5] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 5) == 0);
        assert (mpn_cmp(s, w, 5) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 5);
        mpfq_copy(u, x, 5);
        c1 = mpn_addmul_1(s, s, 5, y[4]);
        c2 = mpfq_fixmp_4_5_addmul05(u, u, y[4]);
        assert (mpn_cmp(s, u, 5) == 0);
        assert (c1 == c2);
    // mul1
    s[5] = mpn_mul_1(s, x, 5, wx);
    u[6] = 0xdeadbeef;
    mpfq_fixmp_4_5_mul1(u, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (mpn_cmp(s,u,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 5);
    mpfq_copy(v, x, 5);
    t[5] = mpn_mul_1(t, t, 5, wx);
    mpfq_fixmp_4_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 5 + 1) == 0);
        // mul05
        s[5] = mpn_mul_1(s, x, 5, y[4]);
        u[5] = 0xdeadbeef;
        mpfq_fixmp_4_5_mul05(u, x, y[4]);
        assert (u[5] == 0xdeadbeef);
        assert (mpn_cmp(s,u,5) == 0);
        assert (s[5] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 5);
        mpfq_copy(v, x, 5);
        t[5] = mpn_mul_1(t, t, 5, y[4]);
        mpfq_fixmp_4_5_mul05(v, v, y[4]);
        assert (mpn_cmp(t, v, 5) == 0);
    // mul
    mpn_mul_n(s, x, y, 5);
    u[9]=0xdeadbeef;
    mpfq_fixmp_4_5_mul(u, x, y);
    assert (u[9]==0xdeadbeef);
    assert (mpn_cmp(s,u,9) == 0);
assert (s[9] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 5);
    u[5] = 0xdeadbeef;
    mpfq_fixmp_4_5_shortmul(u, x, y);
    assert (u[5] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 5) == 0);
    // sqr
    mpn_mul_n(s, x, x, 5);
    u[9]=0xdeadbeef;
    mpfq_fixmp_4_5_sqr(u, x);
    assert (u[9]==0xdeadbeef);
    assert (mpn_cmp(s,u,9) == 0);
assert (s[9] == 0);
    // cmp
    j = mpn_cmp(x, y, 5);
    k = mpfq_fixmp_4_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 5);
    j = mpfq_fixmp_4_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 4);
    j = mpfq_fixmp_4_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_4_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_4_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 5);
    v[5-1] += !v[5-1];
    mpn_tdiv_qr(s+5, s, 0, z, 9, v, 5);
    mpfq_fixmp_4_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 5) == 0);
    // inv
    mpfq_fixmp_4_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_4_5_invmod(v, u, P);
    mpfq_fixmp_4_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 5) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 5);
    mpfq_zero(v, 5);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_4_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_4_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 5);
    memset(u, ~0, 5 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 5);
    /* This creates an n-limb multiple of 257.  */
    v[5] = mpn_lshift(v, x, 5, 8);
    v[5] += mpfq_fixmp_4_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_4_5_add_ui(v, v, v[5]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_4_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 5, 257) == 0);
    //redc
    {
      mp_limb_t p[5], mip[5];
      mp_limb_t xe[5], ye[5], ze[5];
      mp_limb_t invR[5];
      
      // x[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 5);
      mpn_random2(p, 5);
      p[5-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[5-1] > p[5-1])
        p[5-1] = x[5-1];
      if (y[5-1] > p[5-1])
        p[5-1] = y[5-1];
      p[0] |= 1UL;
p[5-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[5-1] += !p[5-1];
      mpfq_zero(w, 2*5);
      w[5]=1;
      mpfq_fixmp_4_5_mod(w, w, p);
      mpfq_fixmp_4_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 5);
        mpfq_fixmp_4_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 5, 1);
	mpfq_fixmp_4_5_shortmul(u, t, mip);
	mpfq_fixmp_4_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 5));
      mpfq_fixmp_4_5_mgy_encode(xe, x, p);
      mpfq_fixmp_4_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_4_5_mul(s, x, y);
      mpfq_fixmp_4_5_mod(t, s, p);
      mpfq_fixmp_4_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_4_5_mul(s, xe, ye);
      mpfq_fixmp_4_5_mod(t, s, p);
      mpfq_fixmp_4_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 5) == 0);
      mpfq_fixmp_4_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 5) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_4_5_mul(t, xe, ye);
      s[9] = mpn_mul_1(s, t, 9, wx);
      mpfq_zero(t, 2*5+1);
      mpfq_fixmp_4_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_4_5_mod(s, t, p);
      mpfq_fixmp_4_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_4_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+5, v, 0, w, 5+1, p, 5);
      assert(mpn_cmp(v, s, 5) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 9 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 9, sat);
      mpfq_fixmp_4_5_redc(s, u, mip, p);
      mpfq_fixmp_4_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+5, s, 0, w, 5+1, p, 5);
      mpfq_zero(w, 2*5+1);
      mpfq_fixmp_4_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+5, t, 0, w, 5+1, p, 5);
      assert(mpn_cmp(s, t, 5) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 9 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 9 + 1, sat);
      mpfq_fixmp_4_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_4_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+5, s, 0, w, 5+1, p, 5);
      mpfq_zero(w, 2*5+1);
      mpfq_fixmp_4_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+5, t, 0, w, 5+1, p, 5);
      assert(mpn_cmp(s, t, 5) == 0);
#ifdef  HAVE_native_4_5_mulredc
      mpfq_fixmp_4_5_mul(u, x, y);
      mpfq_fixmp_4_5_redc(s, u, mip, p);
      mpfq_fixmp_4_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 5) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 5, j);
    else
        mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_4_5_lshift(u, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 5, j);
    else
        mpfq_copy(s, x, 5);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_4_5_rshift(u, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // long_lshift
    j = wx % (5 * GMP_LIMB_BITS);
    mpfq_zero(s, 5);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 5 - k, j);
    else
        mpfq_copy(s + k, x, 5 - k);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_4_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 5) == 0);
    // long_rshift
    j = wx % (5 * GMP_LIMB_BITS);
    mpfq_zero(s, 5);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 5 - k, j);
    else
        mpfq_copy(s, x + k, 5 - k);
    mpfq_copy(u, x, 5);
    mpfq_fixmp_4_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 5) == 0);
}


void test_fixmp_5_5() {
  mp_limb_t s[13];
  mp_limb_t t[13];
  mp_limb_t u[13];
  mp_limb_t v[13];
  mp_limb_t w[13];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[6] = { 1374420825UL, 2881636544UL, 355812594UL, 997066857UL, 3035701496UL, 32232UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[6] = {6517825047235477399UL, 10380803732985570388UL, 10617540373906134554UL, 9832882603074502134UL, 17768816977972787415UL, 3347008897UL };
#endif
  int j, k;
    // add
    u[6] = 0xdeadbeef;
    v[6] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 6);
    c2 = mpfq_fixmp_5_5_add(u, x, y);
    mpfq_fixmp_5_5_add_nc(v, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 6);
    c2 = mpfq_fixmp_5_5_add(u, x, x);
    mpfq_fixmp_5_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 6);
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c2 = mpfq_fixmp_5_5_add(v, v, y);
    mpfq_fixmp_5_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,6) == 0);
    assert (mpn_cmp(s,u,6) == 0);
    // add_ui
    u[6] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 6, wx);
    c2 = mpfq_fixmp_5_5_add_ui(u, x, wx);
    mpfq_fixmp_5_5_add_ui_nc(v, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c1 = mpn_add_1(s, x, 6, wx);
    c2 = mpfq_fixmp_5_5_add_ui(u, u, wx);
    mpfq_fixmp_5_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    // sub
    u[6] = 0xdeadbeef;
    v[6] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 6);
    c2 = mpfq_fixmp_5_5_sub(u, x, y);
    mpfq_fixmp_5_5_sub_nc(v, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 6);
    c2 = mpfq_fixmp_5_5_sub(u, x, x);
    mpfq_fixmp_5_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 6);
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c2 = mpfq_fixmp_5_5_sub(v, v, y);
    mpfq_fixmp_5_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,6) == 0);
    assert (mpn_cmp(s,u,6) == 0);
    // sub_ui
    u[6] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 6, wx);
    c2 = mpfq_fixmp_5_5_sub_ui(u, x, wx);
    mpfq_fixmp_5_5_sub_ui_nc(v, x, wx);
    assert (u[6] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 6);
    mpfq_copy(v, x, 6);
    c1 = mpn_sub_1(s, x, 6, wx);
    c2 = mpfq_fixmp_5_5_sub_ui(u, u, wx);
    mpfq_fixmp_5_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,6) == 0);
    assert (mpn_cmp(s,v,6) == 0);
    // addmul1
    mpfq_copy(s, z, 7);
    mpfq_copy(u, z, 7);
    mpfq_copy(v, z, 7);
    mpfq_copy(w, z, 7);
    c1 = mpn_addmul_1(s, x, 6, wx);
    s[6] += c1;
    c3 = s[6] < c1;
    u[7]=0xdeadbeef;
    v[6]=0xdeadbeef;
    w[7]=0xdeadbeef;
    c4 = mpfq_fixmp_5_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_5_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_5_5_addmul1_nc(w, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (v[6] == 0xdeadbeef);
    assert (w[7] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 7) == 0);
    assert (mpn_cmp(s, v, 6) == 0);
    assert (mpn_cmp(s, w, 7) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    c1 = mpn_addmul_1(s, s, 6, wx);
    c2 = mpfq_fixmp_5_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 6) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 6);
        mpfq_copy(u, z, 6);
        mpfq_copy(w, z, 6);
        c1 = mpn_addmul_1(s, x, 6, y[5]);
        u[6]=0xdeadbeef;
        w[6]=0xdeadbeef;
        c2 = mpfq_fixmp_5_5_addmul05(u, x, y[5]);
        mpfq_fixmp_5_5_addmul05_nc(w, x, y[5]);
        assert (u[6] == 0xdeadbeef);
        assert (w[6] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 6) == 0);
        assert (mpn_cmp(s, w, 6) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 6);
        mpfq_copy(u, x, 6);
        c1 = mpn_addmul_1(s, s, 6, y[5]);
        c2 = mpfq_fixmp_5_5_addmul05(u, u, y[5]);
        assert (mpn_cmp(s, u, 6) == 0);
        assert (c1 == c2);
    // mul1
    s[6] = mpn_mul_1(s, x, 6, wx);
    u[7] = 0xdeadbeef;
    mpfq_fixmp_5_5_mul1(u, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (mpn_cmp(s,u,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 6);
    mpfq_copy(v, x, 6);
    t[6] = mpn_mul_1(t, t, 6, wx);
    mpfq_fixmp_5_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 6 + 1) == 0);
        // mul05
        s[6] = mpn_mul_1(s, x, 6, y[5]);
        u[6] = 0xdeadbeef;
        mpfq_fixmp_5_5_mul05(u, x, y[5]);
        assert (u[6] == 0xdeadbeef);
        assert (mpn_cmp(s,u,6) == 0);
        assert (s[6] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 6);
        mpfq_copy(v, x, 6);
        t[6] = mpn_mul_1(t, t, 6, y[5]);
        mpfq_fixmp_5_5_mul05(v, v, y[5]);
        assert (mpn_cmp(t, v, 6) == 0);
    // mul
    mpn_mul_n(s, x, y, 6);
    u[11]=0xdeadbeef;
    mpfq_fixmp_5_5_mul(u, x, y);
    assert (u[11]==0xdeadbeef);
    assert (mpn_cmp(s,u,11) == 0);
assert (s[11] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 6);
    u[6] = 0xdeadbeef;
    mpfq_fixmp_5_5_shortmul(u, x, y);
    assert (u[6] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 6) == 0);
    // sqr
    mpn_mul_n(s, x, x, 6);
    u[11]=0xdeadbeef;
    mpfq_fixmp_5_5_sqr(u, x);
    assert (u[11]==0xdeadbeef);
    assert (mpn_cmp(s,u,11) == 0);
assert (s[11] == 0);
    // cmp
    j = mpn_cmp(x, y, 6);
    k = mpfq_fixmp_5_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 6);
    j = mpfq_fixmp_5_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 5);
    j = mpfq_fixmp_5_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_5_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_5_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 6);
    v[6-1] += !v[6-1];
    mpn_tdiv_qr(s+6, s, 0, z, 11, v, 6);
    mpfq_fixmp_5_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 6) == 0);
    // inv
    mpfq_fixmp_5_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_5_5_invmod(v, u, P);
    mpfq_fixmp_5_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 6) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 6);
    mpfq_zero(v, 6);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_5_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_5_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 6);
    memset(u, ~0, 6 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 6);
    /* This creates an n-limb multiple of 257.  */
    v[6] = mpn_lshift(v, x, 6, 8);
    v[6] += mpfq_fixmp_5_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_5_5_add_ui(v, v, v[6]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_5_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 6, 257) == 0);
    //redc
    {
      mp_limb_t p[6], mip[6];
      mp_limb_t xe[6], ye[6], ze[6];
      mp_limb_t invR[6];
      
      // x[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 6);
      mpn_random2(p, 6);
      p[6-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[6-1] > p[6-1])
        p[6-1] = x[6-1];
      if (y[6-1] > p[6-1])
        p[6-1] = y[6-1];
      p[0] |= 1UL;
p[6-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[6-1] += !p[6-1];
      mpfq_zero(w, 2*6);
      w[6]=1;
      mpfq_fixmp_5_5_mod(w, w, p);
      mpfq_fixmp_5_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 6);
        mpfq_fixmp_5_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 6, 1);
	mpfq_fixmp_5_5_shortmul(u, t, mip);
	mpfq_fixmp_5_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 6));
      mpfq_fixmp_5_5_mgy_encode(xe, x, p);
      mpfq_fixmp_5_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_5_5_mul(s, x, y);
      mpfq_fixmp_5_5_mod(t, s, p);
      mpfq_fixmp_5_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_5_5_mul(s, xe, ye);
      mpfq_fixmp_5_5_mod(t, s, p);
      mpfq_fixmp_5_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 6) == 0);
      mpfq_fixmp_5_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 6) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_5_5_mul(t, xe, ye);
      s[11] = mpn_mul_1(s, t, 11, wx);
      mpfq_zero(t, 2*6+1);
      mpfq_fixmp_5_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_5_5_mod(s, t, p);
      mpfq_fixmp_5_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_5_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+6, v, 0, w, 6+1, p, 6);
      assert(mpn_cmp(v, s, 6) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 11 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 11, sat);
      mpfq_fixmp_5_5_redc(s, u, mip, p);
      mpfq_fixmp_5_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+6, s, 0, w, 6+1, p, 6);
      mpfq_zero(w, 2*6+1);
      mpfq_fixmp_5_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+6, t, 0, w, 6+1, p, 6);
      assert(mpn_cmp(s, t, 6) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 11 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 11 + 1, sat);
      mpfq_fixmp_5_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_5_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+6, s, 0, w, 6+1, p, 6);
      mpfq_zero(w, 2*6+1);
      mpfq_fixmp_5_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+6, t, 0, w, 6+1, p, 6);
      assert(mpn_cmp(s, t, 6) == 0);
#ifdef  HAVE_native_5_5_mulredc
      mpfq_fixmp_5_5_mul(u, x, y);
      mpfq_fixmp_5_5_redc(s, u, mip, p);
      mpfq_fixmp_5_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 6) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 6, j);
    else
        mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_5_5_lshift(u, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 6, j);
    else
        mpfq_copy(s, x, 6);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_5_5_rshift(u, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // long_lshift
    j = wx % (6 * GMP_LIMB_BITS);
    mpfq_zero(s, 6);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 6 - k, j);
    else
        mpfq_copy(s + k, x, 6 - k);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_5_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 6) == 0);
    // long_rshift
    j = wx % (6 * GMP_LIMB_BITS);
    mpfq_zero(s, 6);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 6 - k, j);
    else
        mpfq_copy(s, x + k, 6 - k);
    mpfq_copy(u, x, 6);
    mpfq_fixmp_5_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 6) == 0);
}


void test_fixmp_6_5() {
  mp_limb_t s[15];
  mp_limb_t t[15];
  mp_limb_t u[15];
  mp_limb_t v[15];
  mp_limb_t w[15];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[7] = { 4269688555UL, 2942999559UL, 634598091UL, 3827993518UL, 782471155UL, 4124854894UL, 45865UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[7] = {14474939553991838547UL, 13588299506894413430UL, 2138988088321967398UL, 4885350755630930435UL, 11400929609603244885UL, 9580305905884609731UL, 3424080812UL };
#endif
  int j, k;
    // add
    u[7] = 0xdeadbeef;
    v[7] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 7);
    c2 = mpfq_fixmp_6_5_add(u, x, y);
    mpfq_fixmp_6_5_add_nc(v, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 7);
    c2 = mpfq_fixmp_6_5_add(u, x, x);
    mpfq_fixmp_6_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 7);
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c2 = mpfq_fixmp_6_5_add(v, v, y);
    mpfq_fixmp_6_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,7) == 0);
    assert (mpn_cmp(s,u,7) == 0);
    // add_ui
    u[7] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 7, wx);
    c2 = mpfq_fixmp_6_5_add_ui(u, x, wx);
    mpfq_fixmp_6_5_add_ui_nc(v, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c1 = mpn_add_1(s, x, 7, wx);
    c2 = mpfq_fixmp_6_5_add_ui(u, u, wx);
    mpfq_fixmp_6_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    // sub
    u[7] = 0xdeadbeef;
    v[7] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 7);
    c2 = mpfq_fixmp_6_5_sub(u, x, y);
    mpfq_fixmp_6_5_sub_nc(v, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 7);
    c2 = mpfq_fixmp_6_5_sub(u, x, x);
    mpfq_fixmp_6_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 7);
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c2 = mpfq_fixmp_6_5_sub(v, v, y);
    mpfq_fixmp_6_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,7) == 0);
    assert (mpn_cmp(s,u,7) == 0);
    // sub_ui
    u[7] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 7, wx);
    c2 = mpfq_fixmp_6_5_sub_ui(u, x, wx);
    mpfq_fixmp_6_5_sub_ui_nc(v, x, wx);
    assert (u[7] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 7);
    mpfq_copy(v, x, 7);
    c1 = mpn_sub_1(s, x, 7, wx);
    c2 = mpfq_fixmp_6_5_sub_ui(u, u, wx);
    mpfq_fixmp_6_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,7) == 0);
    assert (mpn_cmp(s,v,7) == 0);
    // addmul1
    mpfq_copy(s, z, 8);
    mpfq_copy(u, z, 8);
    mpfq_copy(v, z, 8);
    mpfq_copy(w, z, 8);
    c1 = mpn_addmul_1(s, x, 7, wx);
    s[7] += c1;
    c3 = s[7] < c1;
    u[8]=0xdeadbeef;
    v[7]=0xdeadbeef;
    w[8]=0xdeadbeef;
    c4 = mpfq_fixmp_6_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_6_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_6_5_addmul1_nc(w, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (v[7] == 0xdeadbeef);
    assert (w[8] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 8) == 0);
    assert (mpn_cmp(s, v, 7) == 0);
    assert (mpn_cmp(s, w, 8) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    c1 = mpn_addmul_1(s, s, 7, wx);
    c2 = mpfq_fixmp_6_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 7) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 7);
        mpfq_copy(u, z, 7);
        mpfq_copy(w, z, 7);
        c1 = mpn_addmul_1(s, x, 7, y[6]);
        u[7]=0xdeadbeef;
        w[7]=0xdeadbeef;
        c2 = mpfq_fixmp_6_5_addmul05(u, x, y[6]);
        mpfq_fixmp_6_5_addmul05_nc(w, x, y[6]);
        assert (u[7] == 0xdeadbeef);
        assert (w[7] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 7) == 0);
        assert (mpn_cmp(s, w, 7) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 7);
        mpfq_copy(u, x, 7);
        c1 = mpn_addmul_1(s, s, 7, y[6]);
        c2 = mpfq_fixmp_6_5_addmul05(u, u, y[6]);
        assert (mpn_cmp(s, u, 7) == 0);
        assert (c1 == c2);
    // mul1
    s[7] = mpn_mul_1(s, x, 7, wx);
    u[8] = 0xdeadbeef;
    mpfq_fixmp_6_5_mul1(u, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (mpn_cmp(s,u,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 7);
    mpfq_copy(v, x, 7);
    t[7] = mpn_mul_1(t, t, 7, wx);
    mpfq_fixmp_6_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 7 + 1) == 0);
        // mul05
        s[7] = mpn_mul_1(s, x, 7, y[6]);
        u[7] = 0xdeadbeef;
        mpfq_fixmp_6_5_mul05(u, x, y[6]);
        assert (u[7] == 0xdeadbeef);
        assert (mpn_cmp(s,u,7) == 0);
        assert (s[7] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 7);
        mpfq_copy(v, x, 7);
        t[7] = mpn_mul_1(t, t, 7, y[6]);
        mpfq_fixmp_6_5_mul05(v, v, y[6]);
        assert (mpn_cmp(t, v, 7) == 0);
    // mul
    mpn_mul_n(s, x, y, 7);
    u[13]=0xdeadbeef;
    mpfq_fixmp_6_5_mul(u, x, y);
    assert (u[13]==0xdeadbeef);
    assert (mpn_cmp(s,u,13) == 0);
assert (s[13] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 7);
    u[7] = 0xdeadbeef;
    mpfq_fixmp_6_5_shortmul(u, x, y);
    assert (u[7] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 7) == 0);
    // sqr
    mpn_mul_n(s, x, x, 7);
    u[13]=0xdeadbeef;
    mpfq_fixmp_6_5_sqr(u, x);
    assert (u[13]==0xdeadbeef);
    assert (mpn_cmp(s,u,13) == 0);
assert (s[13] == 0);
    // cmp
    j = mpn_cmp(x, y, 7);
    k = mpfq_fixmp_6_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 7);
    j = mpfq_fixmp_6_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 6);
    j = mpfq_fixmp_6_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_6_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_6_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 7);
    v[7-1] += !v[7-1];
    mpn_tdiv_qr(s+7, s, 0, z, 13, v, 7);
    mpfq_fixmp_6_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 7) == 0);
    // inv
    mpfq_fixmp_6_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_6_5_invmod(v, u, P);
    mpfq_fixmp_6_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 7) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 7);
    mpfq_zero(v, 7);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_6_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_6_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 7);
    memset(u, ~0, 7 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 7);
    /* This creates an n-limb multiple of 257.  */
    v[7] = mpn_lshift(v, x, 7, 8);
    v[7] += mpfq_fixmp_6_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_6_5_add_ui(v, v, v[7]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_6_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 7, 257) == 0);
    //redc
    {
      mp_limb_t p[7], mip[7];
      mp_limb_t xe[7], ye[7], ze[7];
      mp_limb_t invR[7];
      
      // x[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 7);
      mpn_random2(p, 7);
      p[7-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[7-1] > p[7-1])
        p[7-1] = x[7-1];
      if (y[7-1] > p[7-1])
        p[7-1] = y[7-1];
      p[0] |= 1UL;
p[7-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[7-1] += !p[7-1];
      mpfq_zero(w, 2*7);
      w[7]=1;
      mpfq_fixmp_6_5_mod(w, w, p);
      mpfq_fixmp_6_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 7);
        mpfq_fixmp_6_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 7, 1);
	mpfq_fixmp_6_5_shortmul(u, t, mip);
	mpfq_fixmp_6_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 7));
      mpfq_fixmp_6_5_mgy_encode(xe, x, p);
      mpfq_fixmp_6_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_6_5_mul(s, x, y);
      mpfq_fixmp_6_5_mod(t, s, p);
      mpfq_fixmp_6_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_6_5_mul(s, xe, ye);
      mpfq_fixmp_6_5_mod(t, s, p);
      mpfq_fixmp_6_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 7) == 0);
      mpfq_fixmp_6_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 7) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_6_5_mul(t, xe, ye);
      s[13] = mpn_mul_1(s, t, 13, wx);
      mpfq_zero(t, 2*7+1);
      mpfq_fixmp_6_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_6_5_mod(s, t, p);
      mpfq_fixmp_6_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_6_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+7, v, 0, w, 7+1, p, 7);
      assert(mpn_cmp(v, s, 7) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 13 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 13, sat);
      mpfq_fixmp_6_5_redc(s, u, mip, p);
      mpfq_fixmp_6_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+7, s, 0, w, 7+1, p, 7);
      mpfq_zero(w, 2*7+1);
      mpfq_fixmp_6_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+7, t, 0, w, 7+1, p, 7);
      assert(mpn_cmp(s, t, 7) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 13 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 13 + 1, sat);
      mpfq_fixmp_6_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_6_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+7, s, 0, w, 7+1, p, 7);
      mpfq_zero(w, 2*7+1);
      mpfq_fixmp_6_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+7, t, 0, w, 7+1, p, 7);
      assert(mpn_cmp(s, t, 7) == 0);
#ifdef  HAVE_native_6_5_mulredc
      mpfq_fixmp_6_5_mul(u, x, y);
      mpfq_fixmp_6_5_redc(s, u, mip, p);
      mpfq_fixmp_6_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 7) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 7, j);
    else
        mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_6_5_lshift(u, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 7, j);
    else
        mpfq_copy(s, x, 7);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_6_5_rshift(u, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // long_lshift
    j = wx % (7 * GMP_LIMB_BITS);
    mpfq_zero(s, 7);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 7 - k, j);
    else
        mpfq_copy(s + k, x, 7 - k);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_6_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 7) == 0);
    // long_rshift
    j = wx % (7 * GMP_LIMB_BITS);
    mpfq_zero(s, 7);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 7 - k, j);
    else
        mpfq_copy(s, x + k, 7 - k);
    mpfq_copy(u, x, 7);
    mpfq_fixmp_6_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 7) == 0);
}


void test_fixmp_7_5() {
  mp_limb_t s[17];
  mp_limb_t t[17];
  mp_limb_t u[17];
  mp_limb_t v[17];
  mp_limb_t w[17];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[8] = { 1034369339UL, 1441735391UL, 87723275UL, 1402427089UL, 2414785066UL, 3785364497UL, 3077063242UL, 5378UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[8] = {16041218152472538879UL, 7399024743779162188UL, 13475004769486012038UL, 11599672397031013245UL, 8334617143536304840UL, 914620958498863738UL, 7941356306452725792UL, 2100729329UL };
#endif
  int j, k;
    // add
    u[8] = 0xdeadbeef;
    v[8] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 8);
    c2 = mpfq_fixmp_7_5_add(u, x, y);
    mpfq_fixmp_7_5_add_nc(v, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 8);
    c2 = mpfq_fixmp_7_5_add(u, x, x);
    mpfq_fixmp_7_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 8);
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c2 = mpfq_fixmp_7_5_add(v, v, y);
    mpfq_fixmp_7_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,8) == 0);
    assert (mpn_cmp(s,u,8) == 0);
    // add_ui
    u[8] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 8, wx);
    c2 = mpfq_fixmp_7_5_add_ui(u, x, wx);
    mpfq_fixmp_7_5_add_ui_nc(v, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c1 = mpn_add_1(s, x, 8, wx);
    c2 = mpfq_fixmp_7_5_add_ui(u, u, wx);
    mpfq_fixmp_7_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    // sub
    u[8] = 0xdeadbeef;
    v[8] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 8);
    c2 = mpfq_fixmp_7_5_sub(u, x, y);
    mpfq_fixmp_7_5_sub_nc(v, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 8);
    c2 = mpfq_fixmp_7_5_sub(u, x, x);
    mpfq_fixmp_7_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 8);
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c2 = mpfq_fixmp_7_5_sub(v, v, y);
    mpfq_fixmp_7_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,8) == 0);
    assert (mpn_cmp(s,u,8) == 0);
    // sub_ui
    u[8] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 8, wx);
    c2 = mpfq_fixmp_7_5_sub_ui(u, x, wx);
    mpfq_fixmp_7_5_sub_ui_nc(v, x, wx);
    assert (u[8] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 8);
    mpfq_copy(v, x, 8);
    c1 = mpn_sub_1(s, x, 8, wx);
    c2 = mpfq_fixmp_7_5_sub_ui(u, u, wx);
    mpfq_fixmp_7_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,8) == 0);
    assert (mpn_cmp(s,v,8) == 0);
    // addmul1
    mpfq_copy(s, z, 9);
    mpfq_copy(u, z, 9);
    mpfq_copy(v, z, 9);
    mpfq_copy(w, z, 9);
    c1 = mpn_addmul_1(s, x, 8, wx);
    s[8] += c1;
    c3 = s[8] < c1;
    u[9]=0xdeadbeef;
    v[8]=0xdeadbeef;
    w[9]=0xdeadbeef;
    c4 = mpfq_fixmp_7_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_7_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_7_5_addmul1_nc(w, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (v[8] == 0xdeadbeef);
    assert (w[9] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 9) == 0);
    assert (mpn_cmp(s, v, 8) == 0);
    assert (mpn_cmp(s, w, 9) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    c1 = mpn_addmul_1(s, s, 8, wx);
    c2 = mpfq_fixmp_7_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 8) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 8);
        mpfq_copy(u, z, 8);
        mpfq_copy(w, z, 8);
        c1 = mpn_addmul_1(s, x, 8, y[7]);
        u[8]=0xdeadbeef;
        w[8]=0xdeadbeef;
        c2 = mpfq_fixmp_7_5_addmul05(u, x, y[7]);
        mpfq_fixmp_7_5_addmul05_nc(w, x, y[7]);
        assert (u[8] == 0xdeadbeef);
        assert (w[8] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 8) == 0);
        assert (mpn_cmp(s, w, 8) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 8);
        mpfq_copy(u, x, 8);
        c1 = mpn_addmul_1(s, s, 8, y[7]);
        c2 = mpfq_fixmp_7_5_addmul05(u, u, y[7]);
        assert (mpn_cmp(s, u, 8) == 0);
        assert (c1 == c2);
    // mul1
    s[8] = mpn_mul_1(s, x, 8, wx);
    u[9] = 0xdeadbeef;
    mpfq_fixmp_7_5_mul1(u, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (mpn_cmp(s,u,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 8);
    mpfq_copy(v, x, 8);
    t[8] = mpn_mul_1(t, t, 8, wx);
    mpfq_fixmp_7_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 8 + 1) == 0);
        // mul05
        s[8] = mpn_mul_1(s, x, 8, y[7]);
        u[8] = 0xdeadbeef;
        mpfq_fixmp_7_5_mul05(u, x, y[7]);
        assert (u[8] == 0xdeadbeef);
        assert (mpn_cmp(s,u,8) == 0);
        assert (s[8] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 8);
        mpfq_copy(v, x, 8);
        t[8] = mpn_mul_1(t, t, 8, y[7]);
        mpfq_fixmp_7_5_mul05(v, v, y[7]);
        assert (mpn_cmp(t, v, 8) == 0);
    // mul
    mpn_mul_n(s, x, y, 8);
    u[15]=0xdeadbeef;
    mpfq_fixmp_7_5_mul(u, x, y);
    assert (u[15]==0xdeadbeef);
    assert (mpn_cmp(s,u,15) == 0);
assert (s[15] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 8);
    u[8] = 0xdeadbeef;
    mpfq_fixmp_7_5_shortmul(u, x, y);
    assert (u[8] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 8) == 0);
    // sqr
    mpn_mul_n(s, x, x, 8);
    u[15]=0xdeadbeef;
    mpfq_fixmp_7_5_sqr(u, x);
    assert (u[15]==0xdeadbeef);
    assert (mpn_cmp(s,u,15) == 0);
assert (s[15] == 0);
    // cmp
    j = mpn_cmp(x, y, 8);
    k = mpfq_fixmp_7_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 8);
    j = mpfq_fixmp_7_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 7);
    j = mpfq_fixmp_7_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_7_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_7_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 8);
    v[8-1] += !v[8-1];
    mpn_tdiv_qr(s+8, s, 0, z, 15, v, 8);
    mpfq_fixmp_7_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 8) == 0);
    // inv
    mpfq_fixmp_7_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_7_5_invmod(v, u, P);
    mpfq_fixmp_7_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 8) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 8);
    mpfq_zero(v, 8);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_7_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_7_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 8);
    memset(u, ~0, 8 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 8);
    /* This creates an n-limb multiple of 257.  */
    v[8] = mpn_lshift(v, x, 8, 8);
    v[8] += mpfq_fixmp_7_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_7_5_add_ui(v, v, v[8]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_7_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 8, 257) == 0);
    //redc
    {
      mp_limb_t p[8], mip[8];
      mp_limb_t xe[8], ye[8], ze[8];
      mp_limb_t invR[8];
      
      // x[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 8);
      mpn_random2(p, 8);
      p[8-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[8-1] > p[8-1])
        p[8-1] = x[8-1];
      if (y[8-1] > p[8-1])
        p[8-1] = y[8-1];
      p[0] |= 1UL;
p[8-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[8-1] += !p[8-1];
      mpfq_zero(w, 2*8);
      w[8]=1;
      mpfq_fixmp_7_5_mod(w, w, p);
      mpfq_fixmp_7_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 8);
        mpfq_fixmp_7_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 8, 1);
	mpfq_fixmp_7_5_shortmul(u, t, mip);
	mpfq_fixmp_7_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 8));
      mpfq_fixmp_7_5_mgy_encode(xe, x, p);
      mpfq_fixmp_7_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_7_5_mul(s, x, y);
      mpfq_fixmp_7_5_mod(t, s, p);
      mpfq_fixmp_7_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_7_5_mul(s, xe, ye);
      mpfq_fixmp_7_5_mod(t, s, p);
      mpfq_fixmp_7_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 8) == 0);
      mpfq_fixmp_7_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 8) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_7_5_mul(t, xe, ye);
      s[15] = mpn_mul_1(s, t, 15, wx);
      mpfq_zero(t, 2*8+1);
      mpfq_fixmp_7_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_7_5_mod(s, t, p);
      mpfq_fixmp_7_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_7_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+8, v, 0, w, 8+1, p, 8);
      assert(mpn_cmp(v, s, 8) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 15 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 15, sat);
      mpfq_fixmp_7_5_redc(s, u, mip, p);
      mpfq_fixmp_7_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+8, s, 0, w, 8+1, p, 8);
      mpfq_zero(w, 2*8+1);
      mpfq_fixmp_7_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+8, t, 0, w, 8+1, p, 8);
      assert(mpn_cmp(s, t, 8) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 15 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 15 + 1, sat);
      mpfq_fixmp_7_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_7_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+8, s, 0, w, 8+1, p, 8);
      mpfq_zero(w, 2*8+1);
      mpfq_fixmp_7_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+8, t, 0, w, 8+1, p, 8);
      assert(mpn_cmp(s, t, 8) == 0);
#ifdef  HAVE_native_7_5_mulredc
      mpfq_fixmp_7_5_mul(u, x, y);
      mpfq_fixmp_7_5_redc(s, u, mip, p);
      mpfq_fixmp_7_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 8) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 8, j);
    else
        mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_7_5_lshift(u, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 8, j);
    else
        mpfq_copy(s, x, 8);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_7_5_rshift(u, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // long_lshift
    j = wx % (8 * GMP_LIMB_BITS);
    mpfq_zero(s, 8);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 8 - k, j);
    else
        mpfq_copy(s + k, x, 8 - k);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_7_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 8) == 0);
    // long_rshift
    j = wx % (8 * GMP_LIMB_BITS);
    mpfq_zero(s, 8);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 8 - k, j);
    else
        mpfq_copy(s, x + k, 8 - k);
    mpfq_copy(u, x, 8);
    mpfq_fixmp_7_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 8) == 0);
}


void test_fixmp_8_5() {
  mp_limb_t s[19];
  mp_limb_t t[19];
  mp_limb_t u[19];
  mp_limb_t v[19];
  mp_limb_t w[19];
  mp_limb_t c1, c2, c3, c4;
#if GMP_LIMB_BITS == 32
  mp_limb_t P[9] = { 1576205429UL, 1707228548UL, 3046795374UL, 846023214UL, 3621694389UL, 53080496UL, 1314150003UL, 322039988UL, 42563UL };
#elif GMP_LIMB_BITS == 64
  mp_limb_t P[9] = {9735718323941594257UL, 17406359535604746461UL, 12157468420192567345UL, 1042532625931608455UL, 13551854600515349299UL, 5602118475269555012UL, 9091328115459418722UL, 12828508824304788877UL, 3205924263UL };
#endif
  int j, k;
    // add
    u[9] = 0xdeadbeef;
    v[9] = 0xdeadbeef;
    c1 = mpn_add_n(s, x, y, 9);
    c2 = mpfq_fixmp_8_5_add(u, x, y);
    mpfq_fixmp_8_5_add_nc(v, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing x==y */
    c1 = mpn_add_n(s, x, x, 9);
    c2 = mpfq_fixmp_8_5_add(u, x, x);
    mpfq_fixmp_8_5_add_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    c1 = mpn_add_n(s, x, y, 9);
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c2 = mpfq_fixmp_8_5_add(v, v, y);
    mpfq_fixmp_8_5_add_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,9) == 0);
    assert (mpn_cmp(s,u,9) == 0);
    // add_ui
    u[9] = 0xdeadbeef;
    c1 = mpn_add_1(s, x, 9, wx);
    c2 = mpfq_fixmp_8_5_add_ui(u, x, wx);
    mpfq_fixmp_8_5_add_ui_nc(v, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c1 = mpn_add_1(s, x, 9, wx);
    c2 = mpfq_fixmp_8_5_add_ui(u, u, wx);
    mpfq_fixmp_8_5_add_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    // sub
    u[9] = 0xdeadbeef;
    v[9] = 0xdeadbeef;
    c1 = mpn_sub_n(s, x, y, 9);
    c2 = mpfq_fixmp_8_5_sub(u, x, y);
    mpfq_fixmp_8_5_sub_nc(v, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing x==y */
    c1 = mpn_sub_n(s, x, x, 9);
    c2 = mpfq_fixmp_8_5_sub(u, x, x);
    mpfq_fixmp_8_5_sub_nc(v, x, x);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    c1 = mpn_sub_n(s, x, y, 9);
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c2 = mpfq_fixmp_8_5_sub(v, v, y);
    mpfq_fixmp_8_5_sub_nc(u, u, y);
    assert (c1 == c2);
    assert (mpn_cmp(s,v,9) == 0);
    assert (mpn_cmp(s,u,9) == 0);
    // sub_ui
    u[9] = 0xdeadbeef;
    c1 = mpn_sub_1(s, x, 9, wx);
    c2 = mpfq_fixmp_8_5_sub_ui(u, x, wx);
    mpfq_fixmp_8_5_sub_ui_nc(v, x, wx);
    assert (u[9] == 0xdeadbeef);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    /* check aliasing z==x */
    mpfq_copy(u, x, 9);
    mpfq_copy(v, x, 9);
    c1 = mpn_sub_1(s, x, 9, wx);
    c2 = mpfq_fixmp_8_5_sub_ui(u, u, wx);
    mpfq_fixmp_8_5_sub_ui_nc(v, v, wx);
    assert (c1 == c2);
    assert (mpn_cmp(s,u,9) == 0);
    assert (mpn_cmp(s,v,9) == 0);
    // addmul1
    mpfq_copy(s, z, 10);
    mpfq_copy(u, z, 10);
    mpfq_copy(v, z, 10);
    mpfq_copy(w, z, 10);
    c1 = mpn_addmul_1(s, x, 9, wx);
    s[9] += c1;
    c3 = s[9] < c1;
    u[10]=0xdeadbeef;
    v[9]=0xdeadbeef;
    w[10]=0xdeadbeef;
    c4 = mpfq_fixmp_8_5_addmul1(u, x, wx);
    c2 = mpfq_fixmp_8_5_addmul1_shortz(v, x, wx);
    mpfq_fixmp_8_5_addmul1_nc(w, x, wx);
    assert (u[10] == 0xdeadbeef);
    assert (v[9] == 0xdeadbeef);
    assert (w[10] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 10) == 0);
    assert (mpn_cmp(s, v, 9) == 0);
    assert (mpn_cmp(s, w, 10) == 0);
    assert (c1 == c2);
    assert (c3 == c4);
    /* check aliasing z==x */
    mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    c1 = mpn_addmul_1(s, s, 9, wx);
    c2 = mpfq_fixmp_8_5_addmul1_shortz(u, u, wx);
    assert (mpn_cmp(s, u, 9) == 0);
    assert (c1 == c2);
        // addmul05
        mpfq_copy(s, z, 9);
        mpfq_copy(u, z, 9);
        mpfq_copy(w, z, 9);
        c1 = mpn_addmul_1(s, x, 9, y[8]);
        u[9]=0xdeadbeef;
        w[9]=0xdeadbeef;
        c2 = mpfq_fixmp_8_5_addmul05(u, x, y[8]);
        mpfq_fixmp_8_5_addmul05_nc(w, x, y[8]);
        assert (u[9] == 0xdeadbeef);
        assert (w[9] == 0xdeadbeef);
        assert (mpn_cmp(s, u, 9) == 0);
        assert (mpn_cmp(s, w, 9) == 0);
        assert (c1 == c2);
        /* check aliasing z==x */
        mpfq_copy(s, x, 9);
        mpfq_copy(u, x, 9);
        c1 = mpn_addmul_1(s, s, 9, y[8]);
        c2 = mpfq_fixmp_8_5_addmul05(u, u, y[8]);
        assert (mpn_cmp(s, u, 9) == 0);
        assert (c1 == c2);
    // mul1
    s[9] = mpn_mul_1(s, x, 9, wx);
    u[10] = 0xdeadbeef;
    mpfq_fixmp_8_5_mul1(u, x, wx);
    assert (u[10] == 0xdeadbeef);
    assert (mpn_cmp(s,u,10) == 0);
    /* check aliasing z==x */
    mpfq_copy(t, x, 9);
    mpfq_copy(v, x, 9);
    t[9] = mpn_mul_1(t, t, 9, wx);
    mpfq_fixmp_8_5_mul1(v, v, wx);
    assert (mpn_cmp(t, v, 9 + 1) == 0);
        // mul05
        s[9] = mpn_mul_1(s, x, 9, y[8]);
        u[9] = 0xdeadbeef;
        mpfq_fixmp_8_5_mul05(u, x, y[8]);
        assert (u[9] == 0xdeadbeef);
        assert (mpn_cmp(s,u,9) == 0);
        assert (s[9] == 0);
        /* check aliasing z==x */
        mpfq_copy(t, x, 9);
        mpfq_copy(v, x, 9);
        t[9] = mpn_mul_1(t, t, 9, y[8]);
        mpfq_fixmp_8_5_mul05(v, v, y[8]);
        assert (mpn_cmp(t, v, 9) == 0);
    // mul
    mpn_mul_n(s, x, y, 9);
    u[17]=0xdeadbeef;
    mpfq_fixmp_8_5_mul(u, x, y);
    assert (u[17]==0xdeadbeef);
    assert (mpn_cmp(s,u,17) == 0);
assert (s[17] == 0);
    // shortmul
    mpn_mul_n(s, x, y, 9);
    u[9] = 0xdeadbeef;
    mpfq_fixmp_8_5_shortmul(u, x, y);
    assert (u[9] == 0xdeadbeef);
    assert (mpn_cmp(s, u, 9) == 0);
    // sqr
    mpn_mul_n(s, x, x, 9);
    u[17]=0xdeadbeef;
    mpfq_fixmp_8_5_sqr(u, x);
    assert (u[17]==0xdeadbeef);
    assert (mpn_cmp(s,u,17) == 0);
assert (s[17] == 0);
    // cmp
    j = mpn_cmp(x, y, 9);
    k = mpfq_fixmp_8_5_cmp(x, y);
    assert (j==k);
    mpfq_copy(u, x, 9);
    j = mpfq_fixmp_8_5_cmp(x, u);
    assert (j==0);
    mpfq_zero(u+1, 8);
    j = mpfq_fixmp_8_5_cmp_ui(u, x[0]);
    assert (j==0);
    j = mpfq_fixmp_8_5_cmp_ui(u, ~x[0]);
    assert (j!=0);
    u[1]=1;
    j = mpfq_fixmp_8_5_cmp_ui(u, x[0]);
    assert (j>0);
    // mod
    mpfq_copy(v, y, 9);
    v[9-1] += !v[9-1];
    mpn_tdiv_qr(s+9, s, 0, z, 17, v, 9);
    mpfq_fixmp_8_5_mod(u, z, v);
    assert(mpn_cmp(s, u, 9) == 0);
    // inv
    mpfq_fixmp_8_5_mod(u, z, P);
    u[0] += !u[0];
    mpfq_fixmp_8_5_invmod(v, u, P);
    mpfq_fixmp_8_5_invmod(v, v, P);
    assert(mpn_cmp(v, u, 9) == 0);
    /* also test the degenerate case invmod(0) or invmod(non invertible) */
    mpfq_zero(u, 9);
    mpfq_zero(v, 9);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_8_5_invmod(v, u, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    v[0] = ~0UL;
    c1 = mpfq_fixmp_8_5_invmod(v, P, P);
    assert(c1 == 0);
    assert(v[0] == 0);
    /* non invertible ; We'll make it slightly artificial */
    mpfq_zero(u, 9);
    memset(u, ~0, 9 * sizeof(mp_limb_t)); /* multiple of 257 */
    mpfq_copy(v, x, 9);
    /* This creates an n-limb multiple of 257.  */
    v[9] = mpn_lshift(v, x, 9, 8);
    v[9] += mpfq_fixmp_8_5_add(v, v, x);
    /* the add_ui below can't overflow since it would imply (for k=257)
     * k*x + cy >= k*2^(n*w).
     * Since cy <= k-1 and x <= 2^(n*w)-1, that can't be. */
    /* note though that v == 2^(n*w)-1 is possible... */
    mpfq_fixmp_8_5_add_ui(v, v, v[9]);
    w[0] = 1UL;
    c1 = mpfq_fixmp_8_5_invmod(w, v, u);
    assert(c1 == 0);
    assert(mpn_mod_1(w, 9, 257) == 0);
    //redc
    {
      mp_limb_t p[9], mip[9];
      mp_limb_t xe[9], ye[9], ze[9];
      mp_limb_t invR[9];
      
      // x[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      // y[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      mpfq_zero(mip, 9);
      mpn_random2(p, 9);
      p[9-1] &= (1UL<<(GMP_LIMB_BITS-1)) -1;
      if (x[9-1] > p[9-1])
        p[9-1] = x[9-1];
      if (y[9-1] > p[9-1])
        p[9-1] = y[9-1];
      p[0] |= 1UL;
p[9-1] &= (1UL<<(GMP_LIMB_BITS>>1))-1;
      p[9-1] += !p[9-1];
      mpfq_zero(w, 2*9);
      w[9]=1;
      mpfq_fixmp_8_5_mod(w, w, p);
      mpfq_fixmp_8_5_invmod(invR, w, p);
      // compute inverse of opposite of p mod R, 
      // with iterated  x <- x + x*(p*x+1) mod R
      mip[0] = 1UL;
      do{ 
        mpfq_copy(v, mip, 9);
        mpfq_fixmp_8_5_shortmul(t, mip, p);
	mpn_add_1(t, t, 9, 1);
	mpfq_fixmp_8_5_shortmul(u, t, mip);
	mpfq_fixmp_8_5_add(mip, u, mip);
      } while (mpn_cmp(v, mip, 9));
      mpfq_fixmp_8_5_mgy_encode(xe, x, p);
      mpfq_fixmp_8_5_mgy_encode(ye, y, p);
      // encode x*y mod p
      mpfq_fixmp_8_5_mul(s, x, y);
      mpfq_fixmp_8_5_mod(t, s, p);
      mpfq_fixmp_8_5_mgy_encode(ze, t, p);
      // do the product in Mgy form
      mpfq_fixmp_8_5_mul(s, xe, ye);
      mpfq_fixmp_8_5_mod(t, s, p);
      mpfq_fixmp_8_5_mgy_decode(u, t, invR, p);
      assert(mpn_cmp(ze, u, 9) == 0);
      mpfq_fixmp_8_5_redc(t, s, mip, p);
      assert(mpn_cmp(ze, t, 9) == 0);
      /* test redc_ur too.
       * We should be able to accumulate many additions here.
       * Let's cheat a bit, and accumulate a number of additions which
       * is simulated. */
      mpfq_fixmp_8_5_mul(t, xe, ye);
      s[17] = mpn_mul_1(s, t, 17, wx);
      mpfq_zero(t, 2*9+1);
      mpfq_fixmp_8_5_redc_ur(t, s, mip, p); /* input is 2n+1-hw words */
      mpfq_fixmp_8_5_mod(s, t, p);
      mpfq_fixmp_8_5_mul1(w, ze, wx);       /* output is n+1-hw/2 words == n+1 */
      /* we can't do mpfq_fixmp_8_5_mod(v, w, p) because that would only eat 2n-hw
       * input limbs, and we want n+1 (which is larger for n=1 and hw=1) */
      mpn_tdiv_qr(v+9, v, 0, w, 9+1, p, 9);
      assert(mpn_cmp(v, s, 9) == 0);
      /* Do the same for something which is carry-friendly */
      mp_limb_t sat = wx;
      /* Saturate the high bits of the random word. This is a very strong
       * degradation of randomness, but its intent is to pinpoint corner
       * cases (we could as well set sat=~0UL, that is). */
      for(j = GMP_LIMB_BITS >> 1 ; j >= 2 ; j >>= 1) {
        sat |= sat << j;
      }
      /* first for plain redc */
      for(int i = 0 ; i < 17 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 17, sat);
      mpfq_fixmp_8_5_redc(s, u, mip, p);
      mpfq_fixmp_8_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+9, s, 0, w, 9+1, p, 9);
      mpfq_zero(w, 2*9+1);
      mpfq_fixmp_8_5_redc(w, v, mip, p);
      mpn_tdiv_qr(t+9, t, 0, w, 9+1, p, 9);
      assert(mpn_cmp(s, t, 9) == 0);
      /* then for redc_ur */
      for(int i = 0 ; i < 17 + 1 ; u[i++] = 1) ;
      mpn_mul_1(v, u, 17 + 1, sat);
      mpfq_fixmp_8_5_redc_ur(s, u, mip, p);
      mpfq_fixmp_8_5_mul1(w, s, sat);
      mpn_tdiv_qr(s+9, s, 0, w, 9+1, p, 9);
      mpfq_zero(w, 2*9+1);
      mpfq_fixmp_8_5_redc_ur(w, v, mip, p);
      mpn_tdiv_qr(t+9, t, 0, w, 9+1, p, 9);
      assert(mpn_cmp(s, t, 9) == 0);
#ifdef  HAVE_native_8_5_mulredc
      mpfq_fixmp_8_5_mul(u, x, y);
      mpfq_fixmp_8_5_redc(s, u, mip, p);
      mpfq_fixmp_8_5_mulredc(t, x, y, mip, p);
      assert(mpn_cmp(s, t, 9) == 0);
#endif
    }
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j) 
        mpn_lshift(s, x, 9, j);
    else
        mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_8_5_lshift(u, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // lshift
    j = wx % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x, 9, j);
    else
        mpfq_copy(s, x, 9);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_8_5_rshift(u, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // long_lshift
    j = wx % (9 * GMP_LIMB_BITS);
    mpfq_zero(s, 9);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_lshift(s + k, x, 9 - k, j);
    else
        mpfq_copy(s + k, x, 9 - k);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_8_5_long_lshift(u, k, j);
    assert (mpn_cmp(s, u, 9) == 0);
    // long_rshift
    j = wx % (9 * GMP_LIMB_BITS);
    mpfq_zero(s, 9);
    k = j / GMP_LIMB_BITS;
    j = j % GMP_LIMB_BITS;
    if (j)
        mpn_rshift(s, x + k, 9 - k, j);
    else
        mpfq_copy(s, x + k, 9 - k);
    mpfq_copy(u, x, 9);
    mpfq_fixmp_8_5_long_rshift(u, k, j);
    assert (mpn_cmp(s, u, 9) == 0);
}





void do_test(int N, int k, int hw, void (*func)())
{
    int kk = 2*k-hw;
    if (kk < k+1) kk = k + 1;

    gx = (mp_limb_t *) malloc(k * sizeof(mp_limb_t));
    gy = (mp_limb_t *) malloc(k * sizeof(mp_limb_t));
    gz = (mp_limb_t *) malloc(2*k * sizeof(mp_limb_t));

    for(int i = 0 ; i < N ; i++) {
        mpn_random2(gx, k);
        mpn_random2(gy, k);
        mpn_random2(gz, kk);

        if (hw) {
            mp_limb_t hwmask = (1UL<<(GMP_LIMB_BITS>>1))-1;
            gx[k-1] &= hwmask; gy[k-1] &= hwmask;
            if (k>1) { gz[2*k-1] = 0xbada55; }
        }

        mpn_random2(&wx, 1);
        x=gx; y=gy; z=gz;
        (*func)();

        if (hw && k>1) { assert(gz[2*k-1] == 0xbada55); }
    }

    free(gx);
    free(gy);
    free(gz);
    printf("."); fflush(stdout);
}
int main(int argc, char **argv) {
  int k=100;

  if (argc==2)
    k = atoi(argv[1]);

  do_test(k, 1, 0, &test_fixmp_1);
  do_test(k, 2, 0, &test_fixmp_2);
  do_test(k, 3, 0, &test_fixmp_3);
  do_test(k, 4, 0, &test_fixmp_4);
  do_test(k, 5, 0, &test_fixmp_5);
  do_test(k, 6, 0, &test_fixmp_6);
  do_test(k, 7, 0, &test_fixmp_7);
  do_test(k, 8, 0, &test_fixmp_8);
  do_test(k, 9, 0, &test_fixmp_9);
  do_test(k, 1, 1, &test_fixmp_0_5);
  do_test(k, 2, 1, &test_fixmp_1_5);
  do_test(k, 3, 1, &test_fixmp_2_5);
  do_test(k, 4, 1, &test_fixmp_3_5);
  do_test(k, 5, 1, &test_fixmp_4_5);
  do_test(k, 6, 1, &test_fixmp_5_5);
  do_test(k, 7, 1, &test_fixmp_6_5);
  do_test(k, 8, 1, &test_fixmp_7_5);
  do_test(k, 9, 1, &test_fixmp_8_5);
  printf("\n");
  return 0;
}
