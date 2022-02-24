#include "cado.h" // IWYU pragma: keep
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include <time.h>
#include "tests_common.h"
#include ARITHMETIC
#include "macros.h"
#include "misc.h"


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

  mod_intinit (minmod);
  mod_intset_ul (minmod, 1);
  mod_intshl (minmod, minmod, MOD_MINBITS);
  do {
#endif

  mod_intset_ul (m, 0UL);
  for (i = 0; i < MOD_SIZE; i++)
    {
      m[i] = u64_random(state);
    }
  m[0] |= 1UL;

  limit_modulus (m);

#ifdef MOD_MINBITS
  } while (mod_intcmp (m, minmod) < 0);
  mod_intclear (minmod);
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
      z[i] = u64_random(state);
    }
}


static void
mpz_set_modint (mpz_t r, const modint_t s)
{
  mpz_import (r, MOD_SIZE, -1, sizeof (unsigned long), 0, 0, s);
}


static inline void
mod_intset_mpz (modint_t r, const mpz_t s)
{
  unsigned long t[MOD_SIZE];
  size_t written;
  mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, s);
  ASSERT_ALWAYS(written <= MOD_SIZE);
  mod_intset_uls(r, t, written);
}

static void
test_mod_intmod (const modint_t la, const modint_t lm)
{
  modint_t s;
  mpz_t mm, ma, mr, ms;
  
  mod_intinit (s);
  mod_intmod (s, la, lm);
  
  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (ms);

  mpz_set_modint (mm, lm);
  mpz_set_modint (ma, la);
  mpz_mod (mr, ma, mm);
  mpz_set_modint (ms, s);

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
  mod_intclear (s);
}

void
tests_mod_intmod (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  mod_intinit (tm);
  mod_intinit (tr);
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
  mod_intclear (tm);
  mod_intclear (tr);
}


static void
test_mod_set_int(const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t r;
  modint_t s;
  mpz_t mm, ma, mr, ms;

  mod_intinit (s);
  mod_initmod_int (m, lm);
  mod_init (r, m);
  mod_set_int (r, la, m);
  mod_get_int (s, r, m);

  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (ms);
  mpz_set_modint (mm, lm);
  mpz_set_modint (ma, la);
  mpz_mod (mr, ma, mm);
  mpz_set_modint (ms, s);

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
  mod_intclear (s);
}


void
tests_mod_set_int(int iter)
{
  modint_t tm, tr;
  int i, j;
  
  mod_intinit (tm);
  mod_intinit (tr);
  /* Test large moduli */
  for (i = 0; i < MOD_SIZE; i++)
    tm[i] = ULONG_MAX;
  limit_modulus (tm);

  for (i = 1; i < 100; i++) /* Test i/2 and (m-i)/2 for i=0, ..., 99 */
    {
      mod_intset_ul (tr, (unsigned long) i);
      test_mod_set_int (tr, tm);
      mod_intset (tr, tm);
      tr[0] -= i;
      test_mod_set_int (tr, tm);
    }
  
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);

      /* Test 0, 1 and -1 residue */
      mod_intset_ul (tr, 0UL);
      test_mod_set_int (tr, tm);
      mod_intset_ul (tr, 1UL);
      test_mod_set_int (tr, tm);
      mod_intset (tr, tm);
      tr[0]--;
      test_mod_set_int (tr, tm);

      /* Test 10 random residues */
      for (j = 0; j < 10; j++)
	{
	  random_integer (tr);
	  test_mod_set_int (tr, tm);
	}
    }
  mod_intclear (tm);
  mod_intclear (tr);
}


static void
test_mod_divn(const modint_t la, const modint_t lm, const unsigned long n)
{
  modulus_t m;
  residue_t r, s, t;

  mod_initmod_int (m, lm);
  mod_init (r, m);
  mod_init (s, m);
  mod_init (t, m);
  mod_set_int (r, la, m);
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

  mod_intinit (tm);
  mod_intinit (tr);
  mpz_init (mod);
  
  /* Test large moduli */
  for (i = 0; i < MOD_SIZE; i++)
    tm[i] = ULONG_MAX;
  limit_modulus (tm);
  mpz_set_modint (mod, tm);
  while (mpz_divisible_ui_p (mod, n))
    {
      tm[0] = ULONG_MAX - 2UL;
      mpz_set_modint (mod, tm);
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
	mpz_set_modint (mod, tm);
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
  mod_intclear (tm);
  mod_intclear (tr);
  mpz_clear (mod);
}


static void
test_mod_gcd (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a;
  modint_t lt;
  mpz_t mt, ma, mr, mm;

  mod_intinit (lt);
  mod_initmod_int (m, lm);
  mod_init (a, m);
  mod_set_int (a, la, m);
  
  mod_gcd (lt, a, m);
  
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mm);
  mpz_init (mt);
  mpz_set_modint (mm, lm); /* modulus */
  mpz_set_modint (ma, la); /* input */
  mpz_set_modint (mt, lt); /* gcd */
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
  mod_intclear (lt);
}


void
tests_mod_gcd (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  mod_intinit (tm);
  mod_intinit (tr);
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
  mod_intclear (tm);
  mod_intclear (tr);
}


static void
test_mod_pow_ul (const modint_t la, const modint_t lm, 
		 const unsigned long e)
{
  modulus_t m;
  residue_t a, r;
  modint_t lt;
  mpz_t mt, ma, mr, mm;

  mod_intinit (lt);
  mod_initmod_int (m, lm);
  mod_init (a, m);
  mod_set_int (a, la, m);
  mod_init (r, m);
  
  mod_pow_ul (r, a, e, m); /* s = a^e % m */
  mod_get_int (lt, r, m);   /* result as integer in t */
  
  mpz_init (mm);
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mt);
  mpz_set_modint (mm, lm); /* modulus */
  mpz_set_modint (ma, la); /* input */
  mpz_set_modint (mt, lt);  /* result of modul function */
  
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
  mod_intclear (lt);
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
      e = (unsigned long) rand ();

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

  mod_intinit (lt);
  mod_initmod_int (m, lm);
  mod_init (r, m);
  
  mod_2pow_ul (r, e, m);  /* s = 2^e % m */
  mod_get_int (lt, r, m); /* result as integer in t */
  
  mpz_init (mm);
  mpz_init (mr);
  mpz_init (mt);
  mpz_init (m2);
  mpz_set_modint (mm, lm); /* modulus */
  mpz_set_ui (m2, 2UL);
  mpz_set_modint (mt, lt);  /* result of mod_2pow_ul() */
  
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
  mod_intclear (lt);
}


void
tests_mod_2pow_ul (int iter)
{
  modint_t tm;
  unsigned long e;
  int i;
  
  mod_intinit (tm);
  for (i = 0; i < iter; i++)
    {
      random_modulus (tm);
      test_mod_2pow_ul (tm, 0UL);
      test_mod_2pow_ul (tm, 1UL);
      test_mod_2pow_ul (tm, 2UL);
      test_mod_2pow_ul (tm, 3UL);
      test_mod_2pow_ul (tm, 4UL);
      test_mod_2pow_ul (tm, ~0UL);
      e = (unsigned long) rand ();
      test_mod_2pow_ul (tm, e);
    }
  mod_intclear (tm);
}


static void 
test_mod_inv (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a, r;
  modint_t lt;
  mpz_t mr, ma, mm, mt;
  int ok1, ok2;

  mod_intinit (lt);
  mod_initmod_int (m, lm);
  mod_init (a, m);
  mod_set_int (a, la, m);
  mod_init (r, m);
  
  ok1 = mod_inv (r, a, m) ? 1 : 0; /* r = 1/a (mod lm) */
  mod_get_int (lt, r, m);
  
  mpz_init (ma);
  mpz_init (mr);
  mpz_init (mm);
  mpz_init (mt);
  mpz_set_modint (ma, la);
  mpz_set_modint (mm, lm);
  mpz_set_modint (mt, lt);
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
  mod_intclear (lt);
}


void
tests_mod_inv (int iter)
{
  modint_t tm, tr;
  int i, j;
  
  mod_intinit (tm);
  mod_intinit (tr);
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
  mod_intclear (tm);
  mod_intclear (tr);
}


static void 
test_mod_jacobi (const modint_t la, const modint_t lm)
{
  modulus_t m;
  residue_t a;
  mpz_t ma, mm;
  int j1, j2;

  mod_initmod_int (m, lm);
  mod_init (a, m);
  mod_set_int (a, la, m);
  
  j1 = mod_jacobi (a, m);
  
  mpz_init (ma);
  mpz_init (mm);
  mpz_set_modint (ma, la);
  mpz_set_modint (mm, lm);
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
  
  mod_intinit (tm);
  mod_intinit (tr);
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
  mod_intclear (tm);
  mod_intclear (tr);
}

#if MOD_MINBITS <= 65 &&  MOD_MAXBITS >= 65 && LONG_BIT == 64
void test_sprp(const mpz_t n, const int is_prime)
{
  const char *prime_str[2] = {"composite", "prime"};
  modint_t ni;
  modulus_t m;

  mod_intinit(ni);
  mod_intset_mpz(ni, n);
  mod_initmod_int(m, ni);
  mod_intclear(ni);

  if (mod_sprp2(m) != is_prime) {
    gmp_fprintf (stderr, "%Zd incorrectly declared %s by mod_sprp2()\n",
                 n, prime_str[!is_prime]);
    abort();
  }

  for (unsigned long b = 2; b < 10; b++) {
    residue_t r;
    mod_init(r, m);
    mod_set_ul(r, b, m);
    if (mod_sprp(r, m) != is_prime) {
      gmp_fprintf (stderr, "%Zd incorrectly declared %s by mod_sprp(%lu)\n",
              n, prime_str[!is_prime], b);
      abort();
    }
    mod_clear(r, m);
  }
  mod_clearmod(m);
}
#endif

void tests_sprp()
{
#if MOD_MINBITS <= 65 &&  MOD_MAXBITS >= 65 && LONG_BIT == 64
  mpz_t n;
  mpz_init(n);
  mpz_set_str(n, "22626675434590779179", 10); /* a prime */
  test_sprp(n, 1);
  mpz_set_str(n, "22626675434590779197", 10); /* a composite */
  test_sprp(n, 0);
  mpz_clear(n);
#endif
}


// coverity[root_function]
int main(int argc, const char **argv)
{
  unsigned long iter = 10000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  printf ("Testing mod_intmod()\n");
  tests_mod_intmod (iter);
  printf ("Testing mod_set_int()\n");
  tests_mod_set_int (iter);
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
  printf ("Testing mod_sprp2()\n");
  tests_sprp();
  tests_common_clear();
  return 0;
}
