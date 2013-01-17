/* Routines that are the same for modredc_15ul.c and modredc_2ul2.c */


/* Divide residue by 3. Returns 1 if division is possible, 0 otherwise.
   Assumes that a+3m does not overflow */

int
mod_div3 (residue_t r, const residue_t a, const modulus_t m)
{
  residue_t t;
  unsigned long a3 = (a[1] % 256UL + a[1] / 256UL + 
		      a[0] % 256UL + a[0] / 256UL) % 3UL;
  const unsigned long m3 = (m[0].m[0] % 256UL + m[0].m[0] / 256UL +
			    m[0].m[1] % 256UL + m[0].m[1] / 256UL) % 3UL;
#ifdef WANT_ASSERT_EXPENSIVE
  residue_t a_backup;
#endif

  if (m3 == 0)
    return 0;

#ifdef WANT_ASSERT_EXPENSIVE
  mod_init_noset0 (a_backup, m);
  mod_set (a_backup, a, m);
#endif

  mod_init_noset0 (t, m);
  t[1] = a[1];
  t[0] = a[0];

  /* Make t[1]:t[0] divisible by 3 */
  if (a3 != 0UL)
    {
      if (a3 + m3 == 3UL) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
	{
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	}
      else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
	{
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	  ularith_add_2ul_2ul (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	}
      
      /* Now t[1]:t[0] is divisible by 3 */
      ASSERT_EXPENSIVE ((t[0] % 3UL + t[1] % 3UL) % 3UL == 0UL);
    }
  
  /* a = a1 * 2^w + a0, 3|a
     Let a = a' * 3 * 2^w + a'', a'' < 3 * 2^w. 
     3 | a'', a'' / 3 < 2^w
     So a / 3 = a' * w + a'' / 3
     a' = trunc(a1 / 3)
     a'' = a0 * 3^{-1} (mod 2^w)
     Hence we get the correct result with one one-word multiplication
     and one one-word truncating division by a small constant.
  */
  
  r[1] = t[1] / 3UL;
#if LONG_BIT == 32
    r[0] = t[0] * 0xaaaaaaabUL; /* 1/3 (mod 2^32) */
#elif LONG_BIT == 64
    r[0] = t[0] * 0xaaaaaaaaaaaaaaabUL; /* 1/3 (mod 2^64) */
#else
#error Unknown word size
#endif

#ifdef WANT_ASSERT_EXPENSIVE
  mod_add (t, r, r, m);
  mod_add (t, t, r, m);
  ASSERT_EXPENSIVE (mod_equal (a_backup, t, m));
  mod_clear (a_backup, m);
#endif
  mod_clear (t, m);

  return 1;
}


/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */

int
mod_div5 (residue_t r, const residue_t a, const modulus_t m)
{
  /* inv_5[i] = -1/i (mod 5) */
  const unsigned long inv_5[5] = {0,4,2,3,1};
#if LONG_BIT == 32
  const unsigned long c = 0xcccccccdUL; /* 1/5 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0xcccccccccccccccdUL; /* 1/5 (mod 2^64) */
#else
#error Unknown word size
#endif
  
  return mod_divn (r, a, 5UL, 1UL, inv_5, c, m);
}


/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */

int
mod_div7 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_7 = (sizeof (unsigned long) == 4) ? 4UL : 2UL;
  /* inv_7[i] = -1/i (mod 7) */
  const unsigned long inv_7[7] = {0,6,3,2,5,4,1};
#if LONG_BIT == 32
  const unsigned long c = 0xb6db6db7UL; /* 1/7 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x6db6db6db6db6db7UL; /* 1/7 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 7UL, w_mod_7, inv_7, c, m);
}


/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */

int
mod_div11 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_11 = (sizeof (unsigned long) == 4) ? 4UL : 5UL;
  /* inv_11[i] = -1/i (mod 11) */
  const unsigned long inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1}; 
#if LONG_BIT == 32
  const unsigned long c = 0xba2e8ba3UL; /* 1/11 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x2e8ba2e8ba2e8ba3UL; /* 1/11 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 11UL, w_mod_11, inv_11, c, m);
}


/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */

int
mod_div13 (residue_t r, const residue_t a, const modulus_t m)
{
  const unsigned long w_mod_13 = (sizeof (unsigned long) == 4) ? 9UL : 3UL;
  /* inv_13[i] = -1/i (mod 13) */
  const unsigned long inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1}; 
#if LONG_BIT == 32
  const unsigned long c = 0xc4ec4ec5UL; /* 1/13 (mod 2^32) */
#elif LONG_BIT == 64
  const unsigned long c = 0x4ec4ec4ec4ec4ec5UL; /* 1/13 (mod 2^64) */
#else
#error Unknown word size
#endif

  return mod_divn (r, a, 13UL, w_mod_13, inv_13, c, m);
}


void
mod_gcd (unsigned long *r, const residue_t A, const modulus_t m)
{
  unsigned long a[2], b[2];
  int sh;

  /* Since we do REDC arithmetic, we must have m odd */
  ASSERT_EXPENSIVE (m[0].m[0] & 1UL);

  if (A[0] == 0UL && A[1] == 0UL)
    {
      r[0] = m[0].m[0];
      r[1] = m[0].m[1];
      return;
    }

  a[0] = A[0];
  a[1] = A[1];
  b[0] = m[0].m[0];
  b[1] = m[0].m[1];

  while (a[1] != 0UL || a[0] != 0UL)
    {
      /* Make a odd */
#if LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) a[0]];
	ularith_shrd (&(a[0]), a[1], sh);
	*(long *) &(a[1]) >>= sh;
      } while (sh == 8);
#else
      if (a[0] == 0UL) /* ctzl does not like zero input */
	{
	  a[0] = a[1];
	  a[1] = ((long)a[1] < 0L) ? (unsigned long) (-1L) : 0UL;
	}
      sh = ularith_ctz (a[0]);
      ularith_shrd (&(a[0]), a[1], sh);
      *(long *) &(a[1]) >>= sh;
#endif
      
      /* Try to make the low two bits of b[0] zero */
      ASSERT_EXPENSIVE (a[0] % 2UL == 1UL);
      ASSERT_EXPENSIVE (b[0] % 2UL == 1UL);
      if ((a[0] ^ b[0]) & 2UL)
	ularith_add_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);
      else
	ularith_sub_2ul_2ul (&(b[0]), &(b[1]), a[0], a[1]);

      if (b[0] == 0UL && b[1] == 0UL)
	{
	  if ((long) a[1] < 0L)
	    {
	      a[1] = -a[1];
	      if (a[0] != 0UL)
		a[1]--;
	      a[0] = -a[0];
	    }
	  r[0] = a[0];
	  r[1] = a[1];
	  return;
	}

      /* Make b odd */
#if LOOKUP_TRAILING_ZEROS
      do {
	sh = trailing_zeros [(unsigned char) b[0]];
	ularith_shrd (&(b[0]), b[1], sh);
	*(long *) &(b[1]) >>= sh;
      } while (sh == 8);
#else
      if (b[0] == 0UL) /* ctzl does not like zero input */
	{
	  b[0] = b[1];
	  b[1] = ((long)b[1] < 0) ? (unsigned long) (-1L) : 0UL;
	}
      sh = ularith_ctz (b[0]);
      ularith_shrd (&(b[0]), b[1], sh);
      *(long *) &(b[1]) >>= sh;
#endif
      ASSERT_EXPENSIVE (a[0] % 2UL == 1UL);
      ASSERT_EXPENSIVE (b[0] % 2UL == 1UL);

      if ((a[0] ^ b[0]) & 2)
	ularith_add_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
      else
	ularith_sub_2ul_2ul (&(a[0]), &(a[1]), b[0], b[1]);
    }

  if ((long) b[1] < 0)
    {
      b[1] = -b[1];
      if (b[0] != 0UL)
	b[1]--;
      b[0] = -b[0];
    }
  r[0] = b[0];
  r[1] = b[1];

  return;
}


/* Compute r = b^e. Here, e is an unsigned long */
void
mod_pow_ul (residue_t r, const residue_t b, const unsigned long e, 
	    const modulus_t m)
{
  unsigned long mask;
  residue_t t;
#ifndef NDEBUG
  unsigned long e1 = e;
#endif
  
  if (e == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  /* Assume t = 1 here for the invariant.
     r^mask * b^e is invariant, and is the result we want */

  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */

  mod_init (t, m);
  mod_set (t, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
#ifndef NDEBUG
  e1 -= mask;
#endif

  while (mask > 1UL)
    {
      mod_sqr (t, t, m);
      mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
      if (e & mask)
        {
	  mod_mul (t, t, b, m);
#ifndef NDEBUG
          e1 -= mask;
#endif
        }
    }
  ASSERT (e1 == 0UL && mask == 1UL);
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Simple addition chains for small multipliers, used in the powering 
   functions with small bases */
static inline void
simple_mul (residue_t r, const residue_t a, const unsigned long b, 
	    const modulus_t m)
{
  if (b == 2UL) {
    mod_add (r, a, a, m);
  } else if (b == 3UL) {
    mod_add (r, a, a, m);
    mod_add (r, r, a, m);
  } else if (b == 5UL) {
    mod_add (r, a, a, m);
    mod_add (r, r, r, m);
    mod_add (r, r, a, m);
  } else {
    ASSERT (b == 7UL);
    mod_add (r, a, a, m);
    mod_add (r, r, r, m);
    mod_add (r, r, r, m);
    mod_sub (r, r, a, m);
  }
}

/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are 
   implemented. Here, e is an unsigned long */
static inline void
mod_npow_ul (residue_t r, const unsigned long b, const unsigned long e, 
	     const modulus_t m)
{
  unsigned long mask;
  residue_t t, u;

  if (e == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  mod_init_noset0 (t, m);
  
  if (b == 2UL) {
    mod_set1 (t, m);
    mod_add (t, t, t, m); /* t = 2 */  
  } else {
    mod_init_noset0 (u, m);
    mod_set1 (u, m);
    simple_mul (t, u, b, m); /* t = b */
  }

  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  ASSERT (e & mask);
  mask >>= 1;

  while (mask > 0UL)
    {
      mod_sqr (t, t, m);
      if (b == 2UL) {
	mod_intshl (t, t, (e & mask) ? 1 : 0);
	ularith_sub_2ul_2ul_ge (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
      } else {
	simple_mul (u, t, b, m);
	if (e & mask)
	  mod_set (t, u, m);
      }
      mask >>= 1;
    }
  mod_set (r, t, m);
  mod_clear (t, m);
  if (b != 2UL)
    mod_clear (u, m);
}


/* Compute r = 2^e mod m. Here, e is an unsigned long */
void
mod_2pow_ul (residue_t r, const unsigned long e, const modulus_t m)
{
  mod_npow_ul (r, 2UL, e, m);
}


/* Compute r = b^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
void
mod_pow_mp (residue_t r, const residue_t b, const unsigned long *e, 
	    const int e_nrwords, const modulus_t m)
{
  unsigned long mask;
  residue_t t;
  int i = e_nrwords - 1;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  /* Find highest set bit in e[i]. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */

  /* Exponentiate */

  mod_init (t, m);
  mod_set (t, b, m);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          mod_sqr (t, t, m);
          if (e[i] & mask)
            mod_mul (t, t, b, m);
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
  mod_set (r, t, m);
  mod_clear (t, m);
}


/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are 
   implemented.  Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
static inline void
mod_npow_mp (residue_t r, const unsigned long b, const unsigned long *e, 
	     const int e_nrwords, const modulus_t m)
{
  residue_t t, u;
  int i = e_nrwords - 1;
  unsigned long mask, ei;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      return;
    }

  mod_init_noset0 (t, m);
  if (b == 2UL) {
    mod_set1 (t, m);
    mod_add (t, t, t, m); /* t = 2 */  
  } else {
    mod_init_noset0 (u, m);
    mod_set1 (u, m);
    simple_mul (t, u, b, m); /* t = b */
  }

  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  mask >>= 1;

  for ( ; i >= 0; i--)
    {
      ei = e[i];
      while (mask > 0UL)
        {
          mod_sqr (t, t, m);
	  if (b == 2UL) {
	    mod_intshl (t, t, (ei & mask) ? 1 : 0);
	    ularith_sub_2ul_2ul_ge (&(t[0]), &(t[1]), m[0].m[0], m[0].m[1]);
	  } else {
	    simple_mul (u, t, b, m);
	    if (ei & mask)
	      mod_set (t, u, m);
	  }
          mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
      mask = ~0UL - (~0UL >> 1);
    }
  mod_set (r, t, m);
  mod_clear (t, m);
  if (b != 2UL)
    mod_clear (u, m);
}


/* Compute r = 2^e mod m.  Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
void
mod_2pow_mp (residue_t r, const unsigned long *e, const int e_nrwords, 
	     const modulus_t m)
{
  mod_npow_mp (r, 2UL, e, e_nrwords, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is an unsigned long. */

void
mod_V_ul (residue_t r, const residue_t b, 
		  const unsigned long e, const modulus_t m)
{
  unsigned long mask;
  residue_t t, t1, two;

  if (e == 0UL)
    {
      mod_set1 (r, m);
      mod_add (r, r, r, m);
      return;
    }
  
  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e);
  /* r = 1, so r^(mask/2) * b^e = r^mask * b^e  */

  /* Exponentiate */
  
  mod_init_noset0 (t, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);
  mod_set (t, b, m);        /* t = b = V_1 (b) */
  mod_sqr (t1, b, m);
  mod_sub (t1, t1, two, m);  /* r1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here t = V_j (b) and t1 = V_{j+1} (b) for j = 1 */

  while (mask > 0UL)
    {
      if (e & mask)
        {
          /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
          mod_mul (t, t, t1, m);
          mod_sub (t, t, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
          mod_sqr (t1, t1, m);
          mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
        }
      else
        {
          /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
          mod_mul (t1, t1, t, m);
          mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
          mod_sqr (t, t, m);
          mod_sub (t, t, two, m);
        }
      mask >>= 1;
    }

  mod_set (r, t, m);
  mod_clear (t, m);
  mod_clear (t1, m);
  mod_clear (two, m);
}


/* Compute r = V_e(b), where V_e(x) is the Chebyshev polynomial defined by
   V_e(x + 1/x) = x^e + 1/x^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords} e[i] * (machine word base)^i */

void
mod_V_mp (residue_t r, const residue_t b, 
		  const unsigned long *e, const int e_nrwords, 
		  const modulus_t m)
{
  unsigned long mask;
  int i = e_nrwords - 1;
  residue_t t, t1, two;

  if (e_nrwords == 0 || e[i] == 0UL)
    {
      mod_set1 (r, m);
      mod_add (r, r, r, m);
      return;
    }

  /* Find highest set bit in e. */
  mask = (1UL << (LONG_BIT - 1)) >> ularith_clz (e[i]);
  /* t = 1, so r^(mask/2) * b^e = t^mask * b^e  */

  /* Exponentiate */

  mod_init_noset0 (t, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (two, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);
  mod_set (t, b, m);         /* t = b = V_1 (b) */
  mod_sqr (t1, b, m);
  mod_sub (t1, t1, two, m);  /* t1 = b^2 - 2 = V_2 (b) */
  mask >>= 1;

  /* Here t = V_j (b) and t1 = V_{j+1} (b) for j = 1 */

  for ( ; i >= 0; i--)
    {
      while (mask > 0UL)
        {
          if (e[i] & mask)
	    {
	      /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
	      mod_mul (t, t, t1, m);
	      mod_sub (t, t, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1} */
	      mod_sqr (t1, t1, m);
	      mod_sub (t1, t1, two, m); /* (V_{j+1})^2 - 2 = V_{2j+2} */
	    }
	  else
	    {
	      /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
	      mod_mul (t1, t1, t, m);
	      mod_sub (t1, t1, b, m); /* V_j * V_{j+1} - V_1 = V_{2j+1}*/
	      mod_sqr (t, t, m);
	      mod_sub (t, t, two, m);
	    }
          mask >>= 1;
        }
      mask = ~0UL - (~0UL >> 1);
    }
  
  mod_set (r, t, m);
  mod_clear (t, m);
  mod_clear (t1, m);
  mod_clear (two, m);
}


/* Returns 1 if r1 == 1 (mod m) or if r1 == -1 (mod m) or if
   one of r1^(2^1), r1^(2^2), ..., r1^(2^(po2-1)) == -1 (mod m),
   zero otherwise. Requires -1 (mod m) in minusone. */

static inline int
find_minus1 (residue_t r1, const residue_t minusone, const int po2, 
             const modulus_t m)
{
  int i;

  if (mod_is1 (r1, m) || mod_equal (r1, minusone, m))
    return 1;

  for (i = 1 ; i < po2; i++)
    {
      mod_sqr (r1, r1, m);
      if (mod_equal (r1, minusone, m))
        break;
    }

  return (i < po2) ? 1 : 0;
}


/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
int
mod_sprp (const residue_t b, const modulus_t m)
{
  residue_t r, minusone;
  int i = 0, po2 = 0;
  modint_t mm1;

  mod_getmod_uls (mm1, m);

  if (mod_intequal_ul (mm1, 1UL))
    return 0;

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1[0]--; /* No borrow since m is odd */
  if (mm1[0] == 0UL)
    {
      mm1[0] = mm1[1];
      mm1[1] = 0UL;
      po2 += LONG_BIT;
    }
  ASSERT (mm1[0] != 0UL);
  i = ularith_clz (mm1[0]);
  mod_intshr (mm1, mm1, i);
  po2 += i;

  mod_init_noset0 (r, m);
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  if (mm1[1] > 0UL)
    mod_pow_mp (r, b, mm1, 2UL, m);
  else
    mod_pow_ul (r, b, mm1[0], m);
  
  i = find_minus1 (r, minusone, po2, m);

  mod_clear (r, m);
  mod_clear (minusone, m);
  return i;
}


/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise. 
   We assume m is odd. */
int
mod_sprp2 (const modulus_t m)
{
  residue_t r, minusone;
  int i = 0, po2 = 0;
  modint_t mm1;

  mod_getmod_uls (mm1, m);

  /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
  mm1[0]--; /* No borrow since m is odd */

  if (mm1[0] == 0UL)
    {
      mm1[0] = mm1[1];
      mm1[1] = 0UL;
      po2 += LONG_BIT;
    }
  ASSERT (mm1[0] != 0UL);
  i = ularith_clz (mm1[0]);
  mod_intshr (mm1, mm1, i);
  po2 += i;

  mod_init_noset0 (r, m);
  mod_init_noset0 (minusone, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Exponentiate */
  if (mm1[1] > 0UL)
    mod_2pow_mp (r, mm1, 2UL, m);
  else
    mod_2pow_ul (r, mm1[0], m);

  i = find_minus1 (r, minusone, po2, m);

  mod_clear (r, m);
  mod_clear (minusone, m);
  return i;
}


int
mod_isprime (const modulus_t m)
{
  residue_t b, minusone, r1;
  modint_t n, mm1;
  int r = 0, po2 = 0, i;
  
  mod_getmod_uls (n, m);

  if (mod_intcmp_ul (n, 1UL) == 0)
    return 0;

  if (n[0] % 2UL == 0UL)
    {
      r = (mod_intcmp_ul(n, 2UL) == 0);
#if defined(PARI)
      printf ("isprime(%lu) == %d /* PARI */\n", n[0], r);
#endif
      return r;
    }

  /* Set mm1 to the odd part of m-1 */
  mod_intset (mm1, n);
  mm1[0]--;
  if (mm1[0] == 0UL)
    {
      mm1[0] = mm1[1];
      mm1[1] = 0UL;
      po2 += LONG_BIT;
    }
  ASSERT (mm1[0] != 0UL);
  i = ularith_ctz (mm1[0]);
  mod_intshr (mm1, mm1, i);
  po2 += i;
  
  mod_init_noset0 (b, m);
  mod_init_noset0 (minusone, m);
  mod_init_noset0 (r1, m);
  mod_set1 (minusone, m);
  mod_neg (minusone, minusone, m);

  /* Do base 2 SPRP test */
  if (mm1[1] != 0UL)
    mod_2pow_mp (r1, mm1, 2UL, m);
  else
    mod_2pow_ul (r1, mm1[0], m);
  /* If n is prime and 1 or 7 (mod 8), then 2 is a square (mod n)
     and one less squaring must suffice. This does not strengthen the
     test but saves one squaring for composite input */
  if (n[0] % 8 == 7)
    { 
      if (!mod_is1 (r1, m))
        goto end;
    }
  else if (!find_minus1 (r1, minusone, po2 - ((n[0] % 8 == 7) ? 1 : 0), m))
    goto end; /* Not prime */

  /* Base 3 is poor at identifying composites == 1 (mod 3), but good at
     identifying composites == 2 (mod 3). Thus we use it only for 2 (mod 3) */
  i = n[0] % 3UL + n[1] % 3;
  if (i == 1 || i == 4)
    {
      if (mm1[1] != 0UL)
	mod_npow_mp (r1, 7UL, mm1, 2UL, m); /* r = 7^mm1 mod m */
      else
	mod_npow_ul (r1, 7UL, mm1[0], m);
      if (!find_minus1 (r1, minusone, po2, m))
	goto end; /* Not prime */

      mod_set_ul_reduced (b, 61UL, m); /* Use addition chain? */
      if (mm1[1] != 0UL)
	mod_pow_mp (r1, b, mm1, 2UL, m); /* r = 61^mm1 mod m */
      else
	mod_pow_ul (r1, b, mm1[0], m);
      if (!find_minus1 (r1, minusone, po2, m))
	goto end; /* Not prime */

#if LONG_BIT == 32
      {
	/* These are the two base 2,7,61 SPSP below 11207066041 */
	const modint_t 
	  c4759123141 = {464155845UL, 1UL},
	  c8411807377 = {4116840081UL, 1UL},
	  c11207066041 = {2617131449UL, 2UL};
	if (mod_intcmp (n, c11207066041) < 0)
	  {
	    r = mod_intcmp (n, c4759123141) != 0 &&
	      mod_intcmp (n, c8411807377) != 0;
	    goto end;
	  }
      }
#endif
      
      if (mm1[1] != 0UL)
	mod_npow_mp (r1, 5UL, mm1, 2UL, m); /* r = 5^mm1 mod m */
      else
	mod_npow_ul (r1, 5UL, mm1[0], m);
      if (!find_minus1 (r1, minusone, po2, m))
	goto end; /* Not prime */
	  
#if LONG_BIT == 32
      {
	/* These are the base 2,5,7,61 SPSP < 10^13 and n == 1 (mod 3) */
	const modint_t 
	  c30926647201 = {861876129UL,7UL},
	  c45821738881 = {2872065921UL,10UL},
	  c74359744201 = {1345300169UL,17UL},
	  c90528271681 = {333958465UL,21UL},
	  c110330267041 = {2956084641UL,25UL},
	  c373303331521 = {3936144065UL,86UL},
	  c440478111067 = {2391446875UL,102UL},
	  c1436309367751 = {1790290887UL,334UL},
	  c1437328758421 = {2809681557UL,334UL},
	  c1858903385041 = {3477513169UL,432UL},
	  c4897239482521 = {976765081UL,1140UL},
	  c5026103290981 = {991554661UL,1170UL},
	  c5219055617887 = {670353247UL,1215UL},
	  c5660137043641 = {3665114809UL,1317UL},
	  c6385803726241 = {3482324385UL,1486UL};
				    
	  r = mod_intcmp (n, c30926647201) != 0 &&
	    mod_intcmp (n, c45821738881) != 0 &&
	    mod_intcmp (n, c74359744201) != 0 &&
	    mod_intcmp (n, c90528271681) != 0 &&
	    mod_intcmp (n, c110330267041) != 0 &&
	    mod_intcmp (n, c373303331521) != 0 &&
	    mod_intcmp (n, c440478111067) != 0 &&
	    mod_intcmp (n, c1436309367751) != 0 &&
	    mod_intcmp (n, c1437328758421) != 0 &&
	    mod_intcmp (n, c1858903385041) != 0 &&
	    mod_intcmp (n, c4897239482521) != 0 &&
	    mod_intcmp (n, c5026103290981) != 0 &&
	    mod_intcmp (n, c5219055617887) != 0 &&
	    mod_intcmp (n, c5660137043641) != 0 &&
	    mod_intcmp (n, c6385803726241) != 0;
      }
#else
      /* For 64 bit arithmetic, a two-word modulus is neccessarily too 
	 large for any deterministic test (with the lists of SPSP 
	 currently available). A modulus >2^64 and == 1 (mod 3) that 
	 survived base 2,5,7,61 is assumed to be prime. */
      r = 1;
#endif
    }
  else
    {
      /* Case n % 3 == 0, 2 */
      
      if (mm1[1] != 0UL)
	mod_npow_mp (r1, 3UL, mm1, 2UL, m); /* r = 3^mm1 mod m */
      else
	mod_npow_ul (r1, 3UL, mm1[0], m);
      if (!find_minus1 (r1, minusone, po2, m))
	goto end; /* Not prime */
      
      if (mm1[1] != 0UL)
	mod_npow_mp (r1, 5UL, mm1, 2UL, m); /* r = 5^mm1 mod m */
      else
	mod_npow_ul (r1, 5UL, mm1[0], m);
      if (!find_minus1 (r1, minusone, po2, m))
	goto end; /* Not prime */

#if LONG_BIT == 32
      {
	/* These are the base 2,3,5 SPSP < 10^13 and n == 2 (mod 3) */
	const modint_t 
	  c244970876021 = {157740149UL,57UL},
	  c405439595861 = {1712670037UL,94UL},
	  c1566655993781 = {3287898037UL,364UL},
	  c3857382025841 = {501394033UL,898UL},
	  c4074652846961 = {3023850353UL,948UL},
	  c5783688565841 = {2662585425UL,1346UL};

	r = mod_intcmp (n, c244970876021) != 0 &&
	  mod_intcmp (n, c405439595861) != 0 &&
	  mod_intcmp (n, c1566655993781) != 0 &&
	  mod_intcmp (n, c3857382025841) != 0 &&
	  mod_intcmp (n, c4074652846961) != 0 &&
	  mod_intcmp (n, c5783688565841) != 0;
      }
#else
      r = 1;
#endif
    }
 
 end:
#if defined(PARI)
  printf ("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
  mod_clear (b, m);
  mod_clear (minusone, m);
  mod_clear (r1, m);
  return r;
}



int
mod_jacobi (const residue_t a_par, const modulus_t m_par)
{
  modint_t a, m, s;
  int t = 1;

  mod_get_uls (a, a_par, m_par);
  mod_getmod_uls (m, m_par);
  
  while (!mod_intequal_ul (a, 0UL))
  {
    while (a[0] % 2UL == 0UL) /* TODO speedup */
    {
      mod_intshr (a, a, 1);
      if (m[0] % 8UL == 3UL || m[0] % 8UL == 5UL)
        t = -t;
    }
    mod_intset (s, a); /* swap a and m */
    mod_intset (a, m);
    mod_intset (m, s);
    if (a[0] % 4UL == 3UL && m[0] % 4UL == 3UL)
      t = -t;
    
    /* m is odd here */
    if (mod_intcmp (a, m) >= 0)
      mod_intmod (a, a, m);
  }
  if (m[1] != 0UL || m[0] != 1UL)
    t = 0;
  
#ifdef MODTRACE
  printf ("kronecker(%lu, %lu) == %d\n", 
          mod_get_ul (a_par, m_par), mod_getmod_ul (m_par), t);
#endif
  return t;
}
