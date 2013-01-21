#include "cado.h"
#include <stdio.h>
#include "modredc_2ul2_default.h"
#include "modredc_2ul_common.c"

#define PARI 0
#if PARI
#define MODINV_PRINT_PARI_M \
    printf ("m = (%lu << %d) + %lu; /* PARI %d */\n", m[0].m[1], LONG_BIT, m[0].m[0], __LINE__)
#define MODINV_PRINT_PARI_x \
    printf ("x = (%lu << %d) + %lu; /* PARI %d */\n", a[1], LONG_BIT, a[0], __LINE__);
#define MODINV_PRINT_PARI_X \
    printf ("X = (%lu << %d) + %lu; /* PARI %d */\n", a[1], LONG_BIT, a[0], __LINE__);
#define MODINV_PRINT_PARI_INVARIANT_A \
    printf ("a = %lu *2^%d + %lu; u = %lu *2^%d + %lu; Mod(u, m) * X == a << %d /* PARIC %d */\n", a[1], LONG_BIT, a[0], u[1], LONG_BIT, u[0], t, __LINE__)
#define MODINV_PRINT_PARI_INVARIANT_B \
    printf ("b = %lu *2^%d + %lu; v = %lu *2^%d + %lu; -Mod(v, m) * X == b << %d /* PARIC %d */\n", b[1], LONG_BIT, b[0], v[1], LONG_BIT, v[0], t, __LINE__)
#else
#define MODINV_PRINT_PARI_M
#define MODINV_PRINT_PARI_x
#define MODINV_PRINT_PARI_X
#define MODINV_PRINT_PARI_INVARIANT_A
#define MODINV_PRINT_PARI_INVARIANT_B
#endif

int
modredc2ul2_inv (residueredc2ul2_t r, const residueredc2ul2_t A, 
		 const modulusredc2ul2_t m) 
{
  modintredc2ul2_t a, b, u, v;
  int t, lsh;
#ifdef WANT_ASSERT_EXPENSIVE
  residueredc2ul2_t tmp;
  
  modredc2ul2_init_noset0 (tmp, m);
  modredc2ul2_set (tmp, A, m);
#endif

  ASSERT_EXPENSIVE (modredc2ul2_intcmp (A, m[0].m) < 0);
  ASSERT_EXPENSIVE (m[0].m[0] & 1UL);

  MODINV_PRINT_PARI_M;

  if (A[0] == 0UL && A[1] == 0UL)
    return 0;

  modredc2ul2_getmod_uls (b, m);

  /* Let A = x*2^{2w}, so we want the Montgomery representation of 1/x, 
     which is 2^{2w}/x. We start by getting a = x */ 
  modredc2ul2_get_uls (a, A, m);
  MODINV_PRINT_PARI_x;

  /* We simply set a = x/2^{2w} and t=0. The result before correction 
     will be 2^(2w+t)/x so we have to divide by t, which may be >64, 
     so we may have to do one or more full and a variable width REDC. */
  /* TODO: If b[1] > 1, we could skip one of the two REDC */
  modredc2ul2_redc1 (a, a, m);
  /* Now a = x/2^w */
  MODINV_PRINT_PARI_X;
  t = -LONG_BIT;

  modredc2ul2_intset_ul (u, 1UL);
  modredc2ul2_intset_ul (v, 0UL);

  MODINV_PRINT_PARI_INVARIANT_A;
  MODINV_PRINT_PARI_INVARIANT_B;

  /* make a odd */
  if (a[0] == 0UL)
    {
      /* x86 bsf gives undefined result for zero input */
      a[0] = a[1];
      a[1] = 0UL;
      t += LONG_BIT;
    }
  ASSERT_EXPENSIVE (a[0] != 0UL);
  lsh = ularith_ctz (a[0]);
  modredc2ul2_intshr (a, a, lsh);
  t += lsh;

  // Here a and b are odd, and a < b
  do {
    /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
    ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, b) < 0);
    ASSERT_EXPENSIVE ((a[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((b[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((u[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((v[0] & 1UL) == 0UL);
    
    MODINV_PRINT_PARI_INVARIANT_A;
    MODINV_PRINT_PARI_INVARIANT_B;

    do {
      modredc2ul2_intsub (b, b, a);
      modredc2ul2_intadd (v, v, u);

      MODINV_PRINT_PARI_INVARIANT_A;
      MODINV_PRINT_PARI_INVARIANT_B;

      if (b[0] == 0UL)
	{
	  b[0] = b[1]; /* b[0] can be odd now, so lsh might be 0 below! */
	  b[1] = 0UL;
	  ASSERT_EXPENSIVE (u[1] == 0UL);
	  u[1] = u[0]; /* Shift left u by LONG_BIT */
	  u[0] = 0UL;
	  t += LONG_BIT;
	}
      else
        {
	  ASSERT_EXPENSIVE (ularith_ctz (b[0]) > 0);
        }
      lsh = ularith_ctz (b[0]);
      ASSERT_EXPENSIVE ((b[0] & ((1UL << lsh) - 1UL)) == 0UL);
      modredc2ul2_intshr (b, b, lsh);
      t += lsh;
      modredc2ul2_intshl (u, u, lsh);
      MODINV_PRINT_PARI_INVARIANT_A;
      MODINV_PRINT_PARI_INVARIANT_B;
    } while (modredc2ul2_intlt (a, b)); /* ~50% branch taken :( */
    
    /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
    ASSERT_EXPENSIVE ((a[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((b[0] & 1UL) == 1UL);
    ASSERT_EXPENSIVE ((u[0] & 1UL) == 0UL);
    ASSERT_EXPENSIVE ((v[0] & 1UL) == 1UL);
    
    if (modredc2ul2_intequal (a, b))
      break;
    ASSERT_EXPENSIVE (modredc2ul2_intcmp (a, b) > 0);
    
    /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
    do {
      modredc2ul2_intsub (a, a, b);
      modredc2ul2_intadd (u, u, v);
      MODINV_PRINT_PARI_INVARIANT_A;
      MODINV_PRINT_PARI_INVARIANT_B;
      
      if (a[0] == 0UL)
	{
	  a[0] = a[1];
	  a[1] = 0UL;
	  v[1] = v[0]; /* Shift left v by LONG_BIT */
	  v[0] = 0UL;
	  t += LONG_BIT;
	}
      else
        {
	  ASSERT_EXPENSIVE (ularith_ctz (a[0]) > 0);
        }
	lsh = ularith_ctz (a[0]);
        ASSERT_EXPENSIVE ((a[0] & ((1UL << lsh) - 1UL)) == 0UL);
	modredc2ul2_intshr (a, a, lsh);
	t += lsh;
	modredc2ul2_intshl (v, v, lsh);
	MODINV_PRINT_PARI_INVARIANT_A;
	MODINV_PRINT_PARI_INVARIANT_B;
    } while (modredc2ul2_intlt (b, a)); /* about 50% branch taken :( */
    /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
  } while (!modredc2ul2_intequal (a, b));
  
  if (modredc2ul2_intcmp_ul (a, 1UL) != 0) /* Non-trivial GCD */
    return 0;

  ASSERT (t >= 0);

  /* Here, the inverse of a is u/2^t mod m. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t
     with impunity. */
  while (t >= LONG_BIT)
    {
      modredc2ul2_redc1 (u, u, m);
      t -= LONG_BIT;
    }

  if (t > 0)
    {
      unsigned long s[5], k;
      k = ((u[0] * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      ularith_mul_ul_ul_2ul (&(s[0]), &(s[1]), k, m[0].m[0]); 
      /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
      ularith_add_2ul_2ul (&(s[0]), &(s[1]), u[0], u[1]);
      /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
      /* s[0] == 0 (mod 2^t) */
      ASSERT_EXPENSIVE ((s[0] & ((1UL << t) - 1UL)) == 0UL);
      s[2] = 0;
      ularith_mul_ul_ul_2ul (&(s[3]), &(s[4]), k, m[0].m[1]);
      ularith_add_2ul_2ul (&(s[1]), &(s[2]), s[3], s[4]);

      /* Now shift s[2]:s[1]:s[0] right by t */
      ularith_shrd (&(s[0]), s[1], t);
      ularith_shrd (&(s[1]), s[2], t);

      u[0] = s[0];
      u[1] = s[1];
      t = 0;
      MODINV_PRINT_PARI_INVARIANT_A;
    }

#ifdef WANT_ASSERT_EXPENSIVE
  modredc2ul2_mul (tmp, tmp, u, m);
  if (!modredc2ul2_is1 (tmp, m))
    {
      modintredc2ul2_t tmpi;
      modredc2ul2_get_uls (tmpi, tmp, m);
      fprintf (stderr, "Error, Mod(1/(%lu + 2^%d * %lu), %lu + 2^%d * %lu) == "
               "%lu + 2^%d * %lu\n",
               A[0], LONG_BIT, A[1], m[0].m[0], LONG_BIT, m[0].m[1],
               tmpi[0], LONG_BIT, tmpi[1]);
      ASSERT_EXPENSIVE (modredc2ul2_intcmp_ul (tmpi, 1UL) == 0);
    }
  modredc2ul2_clear (tmp, m);
#endif

  r[0] = u[0];
  r[1] = u[1];
  return 1;
}


