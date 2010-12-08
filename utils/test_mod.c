#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include <time.h>
#include ARITHMETIC


/* Mask the (possibly multi-word) integer m so that only the
   low maxbits can be non-zero */

static void 
limit_integer (modint_t m, const int maxbits)
{
  int i, j;
  for (i = 0, j = maxbits; i < MOD_SIZE; i++)
    {
      if (j < LONG_BIT)
	{
	  m[i] &= (1UL << j) - 1UL;
	  j = 0;
	}
      else
	j -= LONG_BIT;
    }
}


static void 
limit_modulus (modint_t m)
{
  limit_integer (m, MOD_MAXBITS);
}


static void
random_modulus (modint_t m)
{
  int i;

#ifdef MOD_MINBITS
  modint_t minmod;

  mod_intset_ul (minmod, 0UL);
  minmod[MOD_MINBITS / LONG_BIT] = 1UL << (MOD_MINBITS % LONG_BIT);
  do {
#endif

  mod_intset_ul (m, 0UL);
  for (i = 0; i < MOD_SIZE; i++)
    {
      m[i] = (unsigned long) random() + 
	(RAND_MAX + 1UL) * (unsigned long) random();
    }
  m[0] |= 1UL;

  limit_modulus (m);

#ifdef MOD_MINBITS
  } while (mod_intcmp (m, minmod) < 0);
#endif
  

#if 0
  printf ("Random modulus is ");
  for (i = MOD_SIZE - 1; i >= 0; i--)
    printf ("%lx%s", m[i], (i>0) ? ":" : "\n");
#endif
}


static void
random_integer (modint_t z)
{
  int i;

  for (i = 0; i < MOD_SIZE; i++)
    {
      z[i] = (unsigned long) random() + 
	(RAND_MAX + 1UL) * (unsigned long) random();
    }
}


static void
mod_setint_mpz (mpz_t r, const modint_t s)
{
  mpz_import (r, MOD_SIZE, -1, sizeof (unsigned long), 0, 0, s);
}



static void
test_mod_intmod (const modint_t la, const modint_t lm)
{
  modint_t s;
  mpz_t mm, ma, mr, ms;
  
  mod_intmod (s, la, lm);
  
  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (ms);

  mod_setint_mpz (mm, lm);
  mod_setint_mpz (ma, la);
  mpz_mod (mr, ma, mm);
  mod_setint_mpz (ms, s);

  if (mpz_cmp (mr, ms) != 0)
    {
      gmp_printf ("mod_intmod(%Zd, %Zd) wrong (%Zd), GMP has %Zd\n", 
	          ma, mm, ms, mr);
      abort ();
    }

  mpz_clear (mm);
  mpz_clear (ma);
  mpz_clear (mr);
  mpz_clear (ms);
}

void
tests_mod_intmod (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  /* Test large moduli */
  for (i = 0; i < MOD_SIZE; i++)
    tm[i] = ULONG_MAX;

  for (i = 1; i < iter; i++) /* Test i/2 and (m-i)/2 for i=0, ..., 99 */
    {
      mod_intset_ul (tr, (unsigned long) i);
      test_mod_intmod (tr, tm);
      mod_intset (tr, tm);
      tr[0] -= i;
      test_mod_intmod (tr, tm);
    }
  
  /* Test random moduli of different size */
  for (i = 2; i < MOD_MAXBITS; i++)
    {
      do {
        random_integer (tm);
        limit_integer (tm, i);
      } while (mod_intcmp_ul (tm, 0UL) == 0);
      
      /* Test random integers of different size */
      for (j = 1; j < MOD_MAXBITS; j++)
	{
	  random_integer (tr);
	  limit_integer (tr, j);
	  test_mod_intmod (tr, tm);
	}
    }
}


static void
test_mod_set_uls(const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t r;
  modint_t s;
  mpz_t mm, ma, mr, ms;

  mod_initmod_uls (m, lm);
  mod_init (r, m);
  mod_set_uls (r, la, m);
  mod_get_uls (s, r, m);

  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (ms);
  mod_setint_mpz (mm, lm);
  mod_setint_mpz (ma, la);
  mpz_mod (mr, ma, mm);
  mod_setint_mpz (ms, s);

  if (mpz_cmp (mr, ms) != 0)
    {
      gmp_printf ("mod_set(%Zd, %Zd) wrong (%Zd), GMP has %Zd\n", 
	          ma, mm, ms, mr);
      abort ();
    }
  mpz_clear (mm);
  mpz_clear (ma);
  mpz_clear (mr);
  mpz_clear (ms);
  mod_clear (r, m);
  mod_clearmod (m);
}


void
tests_mod_set_uls(int iter)
{
  modint_t tm, tr;
  int i, j;
  
  /* Test large moduli */
  for (i = 0; i < MOD_SIZE; i++)
    tm[i] = ULONG_MAX;
  limit_modulus (tm);

  for (i = 1; i < 100; i++) /* Test i/2 and (m-i)/2 for i=0, ..., 99 */
    {
      mod_intset_ul (tr, (unsigned long) i);
      test_mod_set_uls (tr, tm);
      mod_intset (tr, tm);
      tr[0] -= i;
      test_mod_set_uls (tr, tm);
    }
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_set_uls (tr, tm);
      mod_intset_ul (tr, 1UL);
      test_mod_set_uls (tr, tm);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_set_uls (tr, tm);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_set_uls (tr, tm);
	}
    }
}


static void
test_mod_divn(const modint_t la, const modint_t lm, const unsigned long n)
{
  modulus_t m;
  residue_t r, s, t;

  mod_initmod_uls (m, lm);
  mod_init (r, m);
  mod_init (s, m);
  mod_init (t, m);
  mod_set_uls (r, la, m);
  if (n == 2UL)
    mod_div2 (s, r, m);   /* s = r/n */
  else if (n == 3UL)
    mod_div3 (s, r, m);
  else if (n == 5UL)
    mod_div5 (s, r, m);
  else if (n == 7UL)
    mod_div7 (s, r, m);
  else if (n == 11UL)
    mod_div11 (s, r, m);
  else if (n == 13UL)
    mod_div13 (s, r, m);
  
  mod_set_ul (t, n, m);
  mod_mul (t, t, s, m); /* t = s*n = r */
  
  if (!mod_equal (r, t, m))
    {
      printf ("mod_div%lu(%lu, %lu) wrong\n", 
	      n, la[0], lm[0]);
      abort ();
    }
  mod_clear (r, m);
  mod_clear (s, m);
  mod_clearmod (m);
}


void
tests_mod_divn(const int iter, const unsigned long n)
{
  modint_t tm, tr;
  int i, j;
  mpz_t mod;

  mpz_init (mod);
  
  /* Test large moduli */
  for (i = 0; i < MOD_SIZE; i++)
    tm[i] = ULONG_MAX;
  limit_modulus (tm);
  mod_setint_mpz (mod, tm);
  while (mpz_divisible_ui_p (mod, n))
    {
      tm[0] = ULONG_MAX - 2UL;
      mod_setint_mpz (mod, tm);
    }

  for (i = 0; i < 100; i++) /* Test i/n and (m-i)/n for i=0, ..., 99 */
    {
      mod_intset_ul (tr, (unsigned long) i);
      test_mod_divn (tr, tm, n);
      mod_intset (tr, tm);
      tr[0] -= i;
      test_mod_divn (tr, tm, n);
    }
  
  for (i = 0; i < iter; i++)
    {
      do {
	random_modulus (tm);
	mod_setint_mpz (mod, tm);
      } while (mpz_divisible_ui_p (mod, n));

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_divn (tr, tm, n);
      mod_intset_ul (tr, 1UL);
      test_mod_divn (tr, tm, n);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_divn (tr, tm, n);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_divn (tr, tm, n);
	}
    }
  mpz_clear (mod);
}


static void
test_mod_gcd (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a;
  modint_t lt;
  mpz_t mt, ma, mr, mm;

  mod_initmod_uls (m, lm);
  mod_init (a, m);
  mod_set_uls (a, la, m);
  
  mod_gcd (lt, a, m);
  
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mm);
  mpz_init (mt);
  mod_setint_mpz (mm, lm); /* modulus */
  mod_setint_mpz (ma, la); /* input */
  mod_setint_mpz (mt, lt); /* gcd */
  mpz_gcd (mr, ma, mm);    /* re-compute gcd */
  if (mpz_cmp (mr, mt) != 0)
    {
      gmp_printf ("mod_gcd(%Zd, %Zd) wrong (%Zd), GMP has %Zd\n", 
                  ma, mm, mt, mr);
      abort ();
    }
  mpz_clear (mt);
  mpz_clear (mm);
  mpz_clear (mr);
  mpz_clear (ma);
  
  mod_clear (a, m);
  mod_clearmod (m);
}


void
tests_mod_gcd (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_gcd (tr, tm);
      mod_intset_ul (tr, 1UL);
      test_mod_gcd (tr, tm);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_gcd (tr, tm);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_gcd (tr, tm);
	}
    }
}


static void
test_mod_pow_ul (const modint_t la, const modint_t lm, 
		 const unsigned long e)
{
  modulus_t m;
  residue_t a, r;
  modint_t lt;
  mpz_t mt, ma, mr, mm;

  mod_initmod_uls (m, lm);
  mod_init (a, m);
  mod_set_uls (a, la, m);
  mod_init (r, m);
  
  mod_pow_ul (r, a, e, m); /* s = a^e % m */
  mod_get_uls (lt, r, m);   /* result as integer in t */
  
  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mt);
  mod_setint_mpz (mm, lm); /* modulus */
  mod_setint_mpz (ma, la); /* input */
  mod_setint_mpz (mt, lt);  /* result of modul function */
  
  mpz_powm_ui (mr, ma, e, mm);
  if (mpz_cmp (mr, mt) != 0)
    {
      gmp_printf ("mod_pow_ui(%Zd, %lu, %Zd) wrong (%Zd)\n", 
	          ma, e, mm, mt);
      abort ();
    }
  mpz_clear (mt);
  mpz_clear (mr);
  mpz_clear (ma);
  mpz_clear (mm);
  
  mod_clear (a, m);
  mod_clear (r, m);
  mod_clearmod (m);
}


void
tests_mod_pow_ul (int iter)
{
  modint_t tm, tr;
  unsigned long e;
  int i, j;
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);
      e = (unsigned long) random();

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_pow_ul (tr, tm, e);
      mod_intset_ul (tr, 1UL);
      test_mod_pow_ul (tr, tm, e);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_pow_ul (tr, tm, e);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_pow_ul (tr, tm, e);
	}
    }
}


static void
test_mod_2pow_ul (const modint_t lm, const unsigned long e)
{
  modulus_t m;
  residue_t r;
  modint_t lt;
  mpz_t mt, mr, mm, m2;

  mod_initmod_uls (m, lm);
  mod_init (r, m);
  
  mod_2pow_ul (r, e, m);  /* s = 2^e % m */
  mod_get_uls (lt, r, m); /* result as integer in t */
  
  mpz_init (mm);
  mpz_init (mr);
  mpz_init (mt);
  mpz_init (m2);
  mod_setint_mpz (mm, lm); /* modulus */
  mpz_set_ui (m2, 2UL);
  mod_setint_mpz (mt, lt);  /* result of mod_2pow_ul() */
  
  mpz_powm_ui (mr, m2, e, mm);
  if (mpz_cmp (mr, mt) != 0)
    {
      gmp_printf ("mod_pow_ui(2, %lu, %Zd) wrong (%Zd)\n", 
	          e, mm, mt);
      abort ();
    }
  mpz_clear (mt);
  mpz_clear (mr);
  mpz_clear (m2);
  mpz_clear (mm);
  
  mod_clear (r, m);
  mod_clearmod (m);
}


void
tests_mod_2pow_ul (int iter)
{
  modint_t tm;
  unsigned long e;
  int i;
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);
      test_mod_2pow_ul (tm, 0UL);
      test_mod_2pow_ul (tm, 1UL);
      test_mod_2pow_ul (tm, 2UL);
      test_mod_2pow_ul (tm, 3UL);
      test_mod_2pow_ul (tm, 4UL);
      test_mod_2pow_ul (tm, ~0UL);
      e = (unsigned long) random();
      test_mod_2pow_ul (tm, e);
    }
}


static void 
test_mod_inv (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a, r;
  modint_t lt;
  mpz_t mr, ma, mm, mt;
  int ok1, ok2;

  mod_initmod_uls (m, lm);
  mod_init (a, m);
  mod_set_uls (a, la, m);
  mod_init (r, m);
  
  ok1 = mod_inv (r, a, m) ? 1 : 0; /* r = 1/a (mod lm) */
  mod_get_uls (lt, r, m);
  
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mm);
  mpz_init (mt);
  mod_setint_mpz (ma, la);
  mod_setint_mpz (mm, lm);
  mod_setint_mpz (mt, lt);
  ok2 = mpz_invert (mr, ma, mm) ? 1 : 0;

  if (ok1 != ok2)
    {
      gmp_printf ("mod_inv(%Zd, %Zd) wrong, return code %d, GMP has code %d\n", 
		  ma, mm, ok1, ok2);
      abort ();
    }
  else if (ok1 == 1)
    {
      if (mpz_cmp (mt, mr) != 0)
	{
	  gmp_printf ("mod_inv(%Zd, %Zd) wrong (%Zd), GMP has %Zd\n", 
		      ma, mm, mt, mr);
	  abort ();
	}
    }
  
  mpz_clear (mr);
  mpz_clear (ma);
  mpz_clear (mm);
  mpz_clear (mt);
  
  mod_clear (a, m);
  mod_clear (r, m);
  mod_clearmod (m); 
}


void
tests_mod_inv (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_inv (tr, tm);
      mod_intset_ul (tr, 1UL);
      test_mod_inv (tr, tm);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_inv (tr, tm);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_inv (tr, tm);
	}
    }
}


static void 
test_mod_jacobi (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a;
  mpz_t ma, mm;
  int j1, j2;

  mod_initmod_uls (m, lm);
  mod_init (a, m);
  mod_set_uls (a, la, m);
  
  j1 = mod_jacobi (a, m);
  
  mpz_init (ma);
  mpz_init (mm);
  mod_setint_mpz (ma, la);
  mod_setint_mpz (mm, lm);
  j2 = mpz_jacobi (ma, mm);

  if (j1 != j2)
    {
      gmp_printf ("mod_jacobi(%Zd, %Zd) wrong (%d), GMP has %d\n", 
		  ma, mm, j1, j2);
      abort ();
    }
  
  mpz_clear (ma);
  mpz_clear (mm);
  
  mod_clear (a, m);
  mod_clearmod (m); 
}


void
tests_mod_jacobi (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_jacobi (tr, tm);
      mod_intset_ul (tr, 1UL);
      test_mod_jacobi (tr, tm);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_jacobi (tr, tm);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_jacobi (tr, tm);
	}
    }
}


int main(int argc, char **argv)
{
  int iter = 10000;

  if (argc > 1 && argv[1][0] == '-' && argv[1][1] == 'r')
    srandom (time (NULL)); /* Seed RNG with nr of seconds since epoch */
  else if (argc > 1)
    iter = atoi (argv[1]);

  printf ("Testing mod_intmod()\n");
  tests_mod_intmod (iter);
  printf ("Testing mod_set_uls()\n");
  tests_mod_set_uls (iter);
  printf ("Testing mod_div2()\n");
  tests_mod_divn (iter, 2);
  printf ("Testing mod_div3()\n");
  tests_mod_divn (iter, 3);
  printf ("Testing mod_div5()\n");
  tests_mod_divn (iter, 5);
  printf ("Testing mod_div7()\n");
  tests_mod_divn (iter, 7);
  printf ("Testing mod_div11()\n");
  tests_mod_divn (iter, 11);
  printf ("Testing mod_div13()\n");
  tests_mod_divn (iter, 13);
  printf ("Testing mod_gcd()\n");
  tests_mod_gcd (iter);
  printf ("Testing mod_pow_ul()\n");
  tests_mod_pow_ul (iter);
  printf ("Testing mod_2pow_ul()\n");
  tests_mod_2pow_ul (iter);
  printf ("Testing mod_inv()\n");
  tests_mod_inv (iter);
  printf ("Testing mod_jacobi()\n");
  tests_mod_jacobi (iter);
  return 0;
}
