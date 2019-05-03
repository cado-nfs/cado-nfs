#include "cado.h"
/* This file defines some functions that work more or less the same
   with mod_ul.h and modredc_ul.h. I.e. mod_div3() and mod_gcd() work
   unchanged with plain and Montgomery representation (so we can work on
   the stored residue directly, whatever its representation is);
   mod_jacobi() converts to plain "unsigned long" first, the others use
   only mod_*() inline functions.
   Speed-critical functions need to be rewritten in assembly for REDC,
   but this is a start.
*/

// #include "mod_64.hpp"
#include "mod_common.cpp"

int
Modulus::div3 (Residue &r, const Residue a) const
{
  const uint64_t a3 = a.r % 3;
  uint64_t ml, m3;
  Residue t;

  assertValid(a);

  ml = getmod_u64 ();
  m3 = ml % 3;
  if (m3 == 0)
    return 0;

  init_noset0 (t);

  set (t, a);
  if (a3 != 0) {
    if (a3 + m3 == 3) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
      t.r += ml;
    else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
      t.r += 2 * ml;
  }

  /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 3.
     (a+k*m)/3 < 2^w, so doing a division (mod 2^w) produces the
     correct result. */

#ifdef WANT_ASSERT_EXPENSIVE
  sub (r, a, t);
  sub (r, r, t);
  sub (r, r, t);
  ASSERT_EXPENSIVE (is0 (r));
#endif

  set (r, t);
  clear (t);

  return 1;
}


int
Modulus::div5 (Residue &r, const Residue a) const
{
  uint64_t ml, m5, k;
  Residue t;
  const uint64_t a5 = a.r % 5;
  const uint64_t inv5[5] = {0,4,2,3,1}; /* inv5[i] = -1/i (mod 5) */

  assertValid(a);

  ml = getmod_u64 ();
  m5 = ml % 5;
  if (m5 == 0)
    return 0;

  init_noset0 (t);
  set (t, a);
  if (a5 != 0)
    {
      /* We want a+km == 0 (mod 5), so k = -a*m^{-1} (mod 5) */
      k = (a5 * inv5[m5]) % 5;
      ASSERT_EXPENSIVE ((k*m5 + a5) % 5 == 0);
      t.r = a.r + k * ml;
    }

  /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 5.
     (a+k*m)/5 < 2^w, so doing a division (mod 2^w) produces the
     correct result. */

    t.r *= 0xcccccccccccccccdULL; /* 1/5 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
  ASSERT_EXPENSIVE (t.r < mod_getmod_ul (m));
  assertValid(t);
  mod_sub (r, a, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  ASSERT_EXPENSIVE (mod_is0 (r, m));
#endif

  set (r, t);
  clear (t);

  return 1;
}


int
Modulus::div7 (Residue &r, const Residue a) const
{
  uint64_t ml, m7, k;
  Residue t;
  const uint64_t a7 = a.r % 7;
  const uint64_t inv7[7] = {0,6,3,2,5,4,1}; /* inv7[i] = -1/i (mod 7) */

  assertValid(a);

  ml = getmod_u64 ();
  m7 = ml % 7;
  if (m7 == 0)
    return 0;

  init_noset0 (t);
  set (t, a);
  if (a7 != 0)
    {
      /* We want a+km == 0 (mod 7), so k = -a*m^{-1} (mod 7) */
      k = (a7 * inv7[m7]) % 7;
      ASSERT_EXPENSIVE ((k*m7 + a7) % 7 == 0);
      t.r = a.r + k * ml;
    }

  /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 7.
     (a+k*m)/7 < 2^w, so doing a division (mod 2^w) produces the
     correct result. */

    t.r *= 0x6db6db6db6db6db7ULL; /* 1/7 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
  assertValid(t);
  mod_sub (r, a, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  mod_sub (r, r, t, m);
  ASSERT_EXPENSIVE (is0 (r));
#endif

  set (r, t);
  clear (t);

  return 1;
}


int
Modulus::div11 (Residue &r, const Residue a) const
{
  uint64_t ml, m11, k;
  Residue t;
  const uint64_t a11 = a.r % 11;
  /* inv11[i] = -1/i (mod 11) */
  const uint64_t inv11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1};

  assertValid(a);

  ml = getmod_u64();
  m11 = ml % 11;
  if (m11 == 0)
    return 0;

  init_noset0 (t);
  set (t, a);
  if (a11 != 0)
    {
      /* We want a+km == 0 (mod 11), so k = -a*m^{-1} (mod 11) */
      k = (a11 * inv11[m11]) % 11;
      ASSERT_EXPENSIVE ((k*m11 + a11) % 11 == 0);
      t.r = a.r + k * ml;
    }

  /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 11.
     (a+k*m)/11 < 2^w, so doing a division (mod 2^w) produces the
     correct result. */

    t.r *= 0x2e8ba2e8ba2e8ba3ULL; /* 1/11 (mod 2^64) */

  set (r, t);
  clear (t);

  return 1;
}


int
Modulus::div13 (Residue &r, const Residue a) const
{
  uint64_t ml, m13, k;
  Residue t;
  const uint64_t a13 = a.r % 13UL;
  /* inv13[i] = -1/i (mod 13) */
  const uint64_t inv13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1};

  assertValid(a);

  ml = getmod_u64 ();
  m13 = ml % 13;
  if (m13 == 0)
    return 0;

  init_noset0 (t);
  set (t, a);
  if (a13 != 0)
    {
      /* We want a+km == 0 (mod 13), so k = -a*m^{-1} (mod 13) */
      k = (a13 * inv13[m13]) % 13;
      ASSERT_EXPENSIVE ((k*m13 + a13) % 13 == 0);
      t.r = a.r + k * ml;
    }

  /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 13.
     (a+k*m)/13 < 2^w, so doing a division (mod 2^w) produces the
     correct result. */

    t.r *= 0x4ec4ec4ec4ec4ec5ULL; /* 1/13 (mod 2^64) */

  set (r, t);
  clear (t);

  return 1;
}


void
Modulus::gcd (Integer &g, const Residue r) const
{
  uint64_t a, b, t;

  a = r.r; /* This works the same for "a" in plain or Montgomery
               representation */
  b = getmod_u64 ();
  /* ASSERT (a < b); Should we require this? */
  ASSERT (b > 0);

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a;
      b = a;
      a = t;
    }

  g = b;
}


/* Compute r = 2^e. Here, e is a uint64_t */
void
Modulus::pow2_u64 (Residue &r, const uint64_t e) const
{
  uint64_t mask;
  Residue t, u;

  if (e == 0)
    {
      set1 (r);
      return;
    }

  mask = (1ULL << 63) >> ularith_clz (e);

  init_noset0 (t);
  init_noset0 (u);
  set1 (t);
  add (t, t, t);
  mask >>= 1;

  while (mask > 0)
    {
      sqr (t, t);
      add (u, t, t);
      if (e & mask)
        set (t, u);
      mask >>= 1;
    }
  set (r, t);
  clear (t);
  clear (u);
}


/* Compute r = 3^e. Here, e is a uint64_t */
void
Modulus::pow3_u64 (Residue &r, const uint64_t e) const
{
  uint64_t mask;
  Residue t, u;

  if (e == 0)
    {
      set1 (r);
      return;
    }

  mask = (1ULL << 63) >> ularith_clz (e);

  init_noset0 (t);
  init_noset0 (u);
  set1 (u);
  add (t, u, u);
  add (t, t, u);
  mask >>= 1;

  while (mask > 0)
    {
      sqr (t, t);
      add (u, t, t);
      add (u, u, t);
      if (e & mask)
        set (t, u);
      mask >>= 1;
    }
  set (r, t);
  clear (t);
  clear (u);
}


/* Computes 2^e (mod m), where e is a multiple precision integer.
   Requires e != 0. The value of 2 in Montgomery representation
   (i.e. 2*2^w (mod m) must be passed. */

void
Modulus::pow2_mp (Residue &r, const uint64_t *e, const int e_nrwords) const
{
  Residue t, u;
  uint64_t mask;
  int i = e_nrwords - 1;

  ASSERT (e_nrwords != 0 && e[i] != 0);

  mask = (1ULL << 63) >> ularith_clz (e[i]);

  init_noset0 (t);
  init_noset0 (u);
  set1 (t);
  add (t, t, t);
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          sqr (t, t);
          add (u, t, t);
          if (e[i] & mask)
            set (t, u);
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }

  set (r, t);
  clear (t);
  clear (u);
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
int
Modulus::sprp (const Residue b) const
{
  Residue r1, minusone;
  int i = 0, po2 = 1;
  uint64_t mm1;

  mm1 = getmod_u64 ();

  /* Set mm1 to the odd part of m-1 */
  mm1 = (mm1 - 1) >> 1;
  while (mm1 % 2 == 0)
    {
      po2++;
      mm1 >>= 1;
    }
  /* Hence, m-1 = mm1 * 2^po2 */

  init_noset0 (r1);
  init_noset0 (minusone);
  set1 (minusone);
  neg (minusone, minusone);

  /* Exponentiate */
  pow_u64 (r1, b, mm1);

  /* Now r1 == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(%lu,%lu) ^ %lu) == %lu /* PARI */\n",
	  get_u64 (b, m), getmod_u64 (m), mm1, get_u64 (r1, m));
#endif

  i = find_minus1 (r1, minusone, po2);

  clear (r1);
  clear (minusone);
  return i;
}


/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise.
   Assumes m > 1 and is odd.
 */
int
Modulus::sprp2 () const
{
  Residue r, minusone;
  int i = 0, po2 = 1;
  uint64_t mm1;

  mm1 = getmod_u64 ();

  /* If m == 1,7 (mod 8), then 2 is a quadratic residue, and we must find
     -1 with one less squaring. This does not reduce the number of
     pseudo-primes because strong pseudo-primes are also Euler pseudo-primes,
     but makes identifying composites a little faster on average. */
  if (mm1 % 8 == 1 || mm1 % 8 == 7)
    po2--;

  /* Set mm1 to the odd part of m-1 */
  mm1 = (mm1 - 1) >> 1;
  while (mm1 % 2UL == 0UL)
    {
      po2++;
      mm1 >>= 1;
    }
  /* Hence, m-1 = mm1 * 2^po2 */

  init_noset0 (r);
  init_noset0 (minusone);
  set1 (minusone);
  neg (minusone, minusone);

  /* Exponentiate */
  pow2_u64 (r, mm1);

  /* Now r == b^mm1 (mod m) */
#if defined(PARI)
  printf ("(Mod(2,%lu) ^ %lu) == %lu /* PARI */\n",
	  getmod_u64 (m), mm1, get_u64 (r, m));
#endif

  i = find_minus1 (r, minusone, po2);

  clear (r);
  clear (minusone);
  return i;
}


int
Modulus::isprime () const
{
  Residue b, minusone, r1;
  const uint64_t n = getmod_u64 ();
  uint64_t mm1;
  int r = 0, po2;

  if (n == 1)
    return 0;

  if (n % 2 == 0)
    {
      r = (n == 2);
#if defined(PARI)
      printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
      return r;
    }

  /* Set mm1 to the odd part of m-1 */
  mm1 = n - 1;
  po2 = u64arith_ctz (mm1);
  mm1 >>= po2;

  init_noset0 (b);
  init_noset0 (minusone);
  init_noset0 (r1);
  set1 (minusone);
  neg (minusone, minusone);

  /* Do base 2 SPRP test */
  pow2_u64 (r1, mm1);   /* r = 2^mm1 mod m */
  /* If n is prime and 1 or 7 (mod 8), then 2 is a square (mod n)
     and one less squaring must suffice. This does not strengthen the
     test but saves one squaring for composite input */
  if (n % 8 == 7)
    {
      if (!is1 (r1))
        goto end;
    }
  else if (!find_minus1 (r1, minusone, po2 - ((n % 8 == 1) ? 1 : 0)))
    goto end; /* Not prime */

  if (n < 2047UL)
    {
      r = 1;
      goto end;
    }

  /* Base 3 is poor at identifying composites == 1 (mod 3), but good at
     identifying composites == 2 (mod 3). Thus we use it only for 2 (mod 3) */
  if (n % 3 == 1)
    {
      set1 (b);
      add (b, b, b);
      add (b, b, b);
      add (b, b, b);
      add (b, b, minusone);  /* b = 7 */
      pow_u64 (r1, b, mm1);   /* r = 7^mm1 mod m */
      if (!find_minus1 (r1, minusone, po2))
        goto end; /* Not prime */

      if (n < 2269093UL) {
        r = (n != 314821UL);
        goto end;
      }

      /* b is still 7 here */
      add (b, b, b); /* 14 */
      sub (b, b, minusone); /* 15 */
      add (b, b, b); /* 30 */
      add (b, b, b); /* 60 */
      sub (b, b, minusone); /* 61 */
      pow_u64 (r1, b, mm1);   /* r = 61^mm1 mod m */
      if (!find_minus1 (r1, minusone, po2))
	goto end; /* Not prime */

#if (ULONG_MAX > 4294967295UL)
      if (n != 4759123141UL && n != 8411807377UL && n < 11207066041UL)
        {
            r = 1;
            goto end;
        }

      set1 (b);
      add (b, b, b);
      add (b, b, b);
      sub (b, b, minusone);    /* b = 5 */
      pow_u64 (r1, b, mm1);   /* r = 5^mm1 mod m */
      if (!find_minus1 (r1, minusone, po2))
	goto end; /* Not prime */
	
          /* These are the base 5,7,61 SPSP < 10^13 and n == 1 (mod 3) */
	  r = (n != 30926647201UL && n != 45821738881UL &&
	       n != 74359744201UL && n != 90528271681UL &&
	       n != 110330267041UL && n != 373303331521UL &&
	       n != 440478111067UL && n != 1436309367751UL &&
	       n != 1437328758421UL && n != 1858903385041UL &&
	       n != 4897239482521UL && n != 5026103290981UL &&
	       n != 5219055617887UL && n != 5660137043641UL &&
	       n != 6385803726241UL);
#else
	      r = 1;
#endif
    }
  else
    {
      /* Case n % 3 == 0, 2 */

      pow3_u64 (r1, mm1);   /* r = 3^mm1 mod m */
      if (!find_minus1 (r1, minusone, po2))
	goto end; /* Not prime */

      if (n < 102690677UL && n != 5173601UL && n != 16070429UL &&
          n != 54029741UL)
        {
	  r = 1;
	  goto end;
	}

      set1 (b);
      add (b, b, b);
      add (b, b, b);
      sub (b, b, minusone);    /* b = 5 */
      pow_u64 (r1, b, mm1);   /* r = 5^mm1 mod m */
      if (!find_minus1 (r1, minusone, po2))
	goto end; /* Not prime */

#if (ULONG_MAX > 4294967295UL)
      /* These are the base 3,5 SPSP < 10^13 with n == 2 (mod 3) */
      r = (n != 244970876021UL && n != 405439595861UL &&
	   n != 1566655993781UL && n != 3857382025841UL &&
	   n != 4074652846961UL && n != 5783688565841UL);
#else
      r = 1;
#endif
    }

 end:
#if defined(PARI)
  printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
  clear (b);
  clear (minusone);
  clear (r1);
  return r;
}

#if 1

int
Modulus::jacobi (const Residue a_par) const
{
  uint64_t x;
  uint64_t mm;
  unsigned int s, j;

  /* Get residue in Montgomery form directly without converting */
  x = a_par.r;
  mm = getmod_u64 ();
  ASSERT (x < mm);
  ASSERT(mm % 2 == 1);

  j = ularith_ctz(x);
  x = x >> j;
  /* If we divide by an odd power of 2, and 2 is a QNR, flip sign */
  /* 2 is a QNR (mod mm) iff mm = 3,5 (mod 8)
     mm = 1 = 001b:   1
     mm = 3 = 011b:  -1
     mm = 5 = 101b:  -1
     mm = 7 = 111b:   1
     Hence we can store in s the exponent of -1, i.e., s=0 for jacobi()=1
     and s=1 for jacobi()=-1, and update s ^= (mm>>1) & (mm>>2) & 1.
     We can do the &1 at the very end.
     In fact, we store the exponent of -1 in the second bit of s.
     The s ^= ((j<<1) & (mm ^ (mm>>1))) still needs 2 shift but one of them can
     be done with LEA, and f = s ^ (x&mm) needs no shift */

  s = ((j<<1) & (mm ^ (mm>>1)));

  while (x > 1) {
    /* Here, x < mm, x and mm are odd */

    /* Implicitly swap by reversing roles of x and mm in next loop */
    /* Flip sign if both are 3 (mod 4) */
    s = s ^ (x&mm);

    /* Make mm smaller by subtracting and shifting */
    do {
      mm -= x; /* Difference is even */
      if (mm == 0)
        break;
      /* Make odd again */
      j = ularith_ctz(mm);
      s ^= ((j<<1) & (x ^ (x>>1)));
      mm >>= j;
    } while (mm >= x);

    if (mm <= 1) {
      x = mm;
      break;
    }

    /* Flip sign if both are 3 (mod 4) */
    /* Implicitly swap again */
    s = s ^ (x&mm);

    /* Make x<mm by subtracting and shifting */
    do {
      x -= mm; /* Difference is even */
      if (x == 0)
        break;
      /* Make odd again */
      j = ularith_ctz(x);
      s ^= ((j<<1) & (mm ^ (mm>>1)));
      x >>= j;
    } while (x >= mm);
  }

  if (x == 0)
    return 0;
  return ((s & 2) == 0) ? 1 : -1;
}

#else

/*
#!/usr/bin/env python3
# Python program to create mult.h

def rate(k,b):
  # The 0.25 magic constant here tries to estimate the ratio m/x, 
  # to minimize (x+c*m)/2^b
  r = (abs(k)*0.25 + 1.)/2**b
  # print ("rate(" + str(k) + ", " + str(b) + ") = " + str(r))
  return(r)

def bestb(k):
  best_b = 0
  best_r = 1
  best_c = 0
  for b in range(1, 8):
    c = k % (2**b)
    r = rate(c, b) 
    if r < best_r: 
      best_r=r 
      best_b = b 
      best_c = c
    c = - (2**b - c)
    r = rate(c, b)
    if r < best_r: 
      best_r=r 
      best_b = b 
      best_c = c
  # print ("bestb(" + str(k) + ") = " + str(best_b))
  return([k, best_b, best_c, 1/best_r])


r = [str(-bestb(2*i+1)[2]) for i in range(0, 128) ]
print("static char mult[128] = {" + ", ".join(r) + "};")
*/

#include "mult.h"
static unsigned char invmod8[256] = {
0, 1, 0, 171, 0, 205, 0, 183, 0, 57, 0, 163, 0, 197, 0, 239, 0, 241, 0, 27, 0, 61, 0, 167, 0, 41, 0, 19, 0, 53, 0, 223, 0, 225, 0, 139, 0, 173, 0, 151, 0, 25, 0, 131, 0, 165, 0, 207, 0, 209, 0, 251, 0, 29, 0, 135, 0, 9, 0, 243, 0, 21, 0, 191, 0, 193, 0, 107, 0, 141, 0, 119, 0, 249, 0, 99, 0, 133, 0, 175, 0, 177, 0, 219, 0, 253, 0, 103, 0, 233, 0, 211, 0, 245, 0, 159, 0, 161, 0, 75, 0, 109, 0, 87, 0, 217, 0, 67, 0, 101, 0, 143, 0, 145, 0, 187, 0, 221, 0, 71, 0, 201, 0, 179, 0, 213, 0, 127, 0, 129, 0, 43, 0, 77, 0, 55, 0, 185, 0, 35, 0, 69, 0, 111, 0, 113, 0, 155, 0, 189, 0, 39, 0, 169, 0, 147, 0, 181, 0, 95, 0, 97, 0, 11, 0, 45, 0, 23, 0, 153, 0, 3, 0, 37, 0, 79, 0, 81, 0, 123, 0, 157, 0, 7, 0, 137, 0, 115, 0, 149, 0, 63, 0, 65, 0, 235, 0, 13, 0, 247, 0, 121, 0, 227, 0, 5, 0, 47, 0, 49, 0, 91, 0, 125, 0, 231, 0, 105, 0, 83, 0, 117, 0, 31, 0, 33, 0, 203, 0, 237, 0, 215, 0, 89, 0, 195, 0, 229, 0, 15, 0, 17, 0, 59, 0, 93, 0, 199, 0, 73, 0, 51, 0, 85, 0, 255
};

static inline int
s_val(unsigned int s) {
  return ((s & 2) == 0) ? 1 : -1;
}

static int
Modulus::jacobi1 (const Residue a_par, const Modulus m_par)
{
  unsigned long x, m;
  unsigned int s, j;

  x = mod_get_ul (a_par, m_par);
  m = mod_getmod_ul (m_par);
  ASSERT (x < m);
  ASSERT(m % 2 == 1);

  j = ularith_ctz(x);
  x = x >> j;

  s = ((j<<1) & (m ^ (m>>1)));

  while (x > 1) {
    unsigned long t;
    unsigned char inv;

    // printf ("kronecker(%lu, %lu) == %d * kronecker(%lu, %lu)\n",
    //        mod_get_ul(a_par, m_par), mod_getmod_ul (m_par), s_val(s), x, m);
    /* Here x < m. Swap to make x > m */
    t = m;
    m = x; 
    x = t;
    s = s ^ (x&m);

    /* Now x > m */
    inv = invmod8[(unsigned char)m];
    do {
      /* Do a REDC-like step. We want a multiplier k such that the low 
         8 bits of x+k*m are all zero. 
         That is, we want k = -x/m (mod 2^8). */
      unsigned char k;
      unsigned long t1;
      long int c, t2;
      // const unsigned long old_x = x;
      
      k = inv * (unsigned char)x;
      // ASSERT_ALWAYS((k & 1) == 1);
      c = mult[k / 2];
      
      /* Compute x+cm */
      long tmp = c >> 63;
      ularith_mul_ul_ul_2ul (&t1, (unsigned long *)&t2, (c ^ tmp) - tmp, m);
      t2 ^= tmp;
      t1 ^= tmp;
      ularith_add_ul_2ul (&t1, (unsigned long *)&t2, x-tmp);
      tmp = ((long) t2) >> 63;
      
      t2 ^= tmp;
      t1 ^= tmp;
      s ^= m & tmp;
      ularith_add_ul_2ul (&t1, (unsigned long *)&t2, -tmp);
      // ASSERT_ALWAYS(t2 >= 0);

      if (t1 == 0) {
        if (t2 == 0) {
          x = 0;
          break;
        }
        t1 = t2;
        /* Divided by 2^64 which is square, so no adjustment to s */
        t2 = 0;
      }

      j = ularith_ctz(t1);
      ularith_shrd (&t1, t2, j);
      // ASSERT_ALWAYS((t2 >> j) == 0);
      x = t1;
      s ^= ((j<<1) & (m ^ (m>>1)));
      // ASSERT_ALWAYS(x < old_x);
      // printf ("%f\n", (double)old_x / (double)x);
    } while (x >= m);
  }

  if (x == 0)
    return 0;
  return s_val(s);
}

int
Modulus::jacobi (const Residue a_par, const Modulus m_par)
{
  unsigned long x, m;
  unsigned int s, j;

  x = mod_get_ul (a_par, m_par);
  m = mod_getmod_ul (m_par);
  ASSERT (x < m);
  ASSERT(m % 2 == 1);

  if ((LONG_MAX - x) / 50 < m)
    return mod_jacobi1 (a_par, m_par);

  j = ularith_ctz(x);
  x = x >> j;

  s = ((j<<1) & (m ^ (m>>1)));

  while (x > 1) {
    unsigned long t;
    unsigned char inv;

    // printf ("kronecker(%lu, %lu) == %d * kronecker(%lu, %lu)\n",
    //        mod_get_ul(a_par, m_par), mod_getmod_ul (m_par), s_val(s), x, m);
    /* Here x < m. Swap to make x > m */
    t = m;
    m = x; 
    x = t;
    s = s ^ (x&m);

    /* Now x > m */
    inv = invmod8[(unsigned char)m];
    do {
      /* Do a REDC-like step. We want a multiplier k such that the low 
         8 bits of x+k*m are all zero. 
         That is, we want k = -x/m (mod 2^8). */
      unsigned char k;
      long int c;
      // const unsigned long old_x = x;
      
      k = inv * x;
      // ASSERT_ALWAYS((k & 1) == 1);
      c = mult[k / 2];
      
      c = x + c*m;
      x = c;
      c >>= 63;
      x = (x ^ c) - c;
      
      if (x == 0) {
        break;
      }
      s ^= m & c;

      j = ularith_ctz(x);
      x >>= j;
      s ^= ((j<<1) & (m ^ (m>>1)));
      // printf ("%f\n", (double)old_x / (double)x);
    } while (x >= m);
  }

  if (x == 0)
    return 0;
  return s_val(s);
}

#endif
