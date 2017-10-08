#include "cado.h"
#include "utils.h"
#include "utils_int64.h"

uint64_t invmod_uint64(uint64_t xx, uint64_t mm)
{
  uint64_t yy;

#if LONG_BIT == 64
  modulus_t m;
  residue_t x, y;
  modul_initmod_ul(m, mm);        //we work mod mm
  modul_init_noset0(x, m);
  modul_init_noset0(y, m);
  modul_set_ul_reduced(x, xx, m);  // we know that xx < m

  modul_inv(y, x, m);

  yy = modul_get_ul(y, m);
  modul_clear(x, m);
  modul_clear(y, m);
  modul_clearmod(m);
#else
  /* otherwise use a stupid placeholder */
  mpz_t xz, mz;
  mpz_init(xz);
  mpz_init(mz);
  mpz_set_uint64(xz, xx);
  mpz_set_uint64(mz, mm);
  mpz_invert(xz, xz, mz);
  yy = mpz_get_uint64(xz);
  mpz_clear(mz);
  mpz_clear(xz);
#endif

  return yy;
}


void int64_fdiv_qr(int64_t * q, int64_t * r, int64_t n, int64_t d)
{
  ASSERT(d != 0);

  if ((n < 0) ^ (d < 0)) {
#ifndef NDEBUG
    if (n < 0) {
      ASSERT(d > 0);
    } else {
      ASSERT(d < 0);
    }
#endif // NDEBUG

    * q = n / d;
    if (* q * d != n) {
      * q = * q - 1;
    }

  } else {
#ifndef NDEBUG
    if (n < 0) {
      ASSERT(d < 0);
    } else {
      ASSERT(d > 0);
    }
#endif // NDEBUG

    * q = n / d;

  }

  * r = n - d * * q;

  ASSERT(ABS(* r) < ABS(d));

#ifndef NDEBUG
  mpz_t n_Z, d_Z, q_Z, r_Z;
  mpz_init(n_Z);
  mpz_init(d_Z);
  mpz_init(q_Z);
  mpz_init(r_Z);
  mpz_set_si(n_Z, n);
  mpz_set_si(d_Z, d);
  mpz_fdiv_qr(q_Z, r_Z, n_Z, d_Z);

  ASSERT(mpz_cmp_si(q_Z, * q) == 0);
  ASSERT(mpz_cmp_si(r_Z, * r) == 0);

  mpz_clear(n_Z);
  mpz_clear(d_Z);
  mpz_clear(q_Z);
  mpz_clear(r_Z);
#endif // NDEBUG
}

void int64_xgcd(int64_t * g, int64_t * u, int64_t * v, int64_t a, int64_t b)
{
  if (!a) {
    ASSERT(a == 0);

    * g = b;
    * u = 0;
    * v = 1;
  } else {
    int64_t x = 0;
    int64_t y = 0;
    int64_xgcd(g, &y, &x, b % a, a);
    * u = x - (b / a) * y;
    * v = y;
  }

  ASSERT(* g == * u * a + * v * b);
  ASSERT(ABS(* g) == ABS(gcd_int64(a, b)));
}

//TODO: change that!
void int64_gcdext(int64_t * e, int64_t * s, int64_t * t, int64_t a, int64_t b)
{
  mpz_t e_Z, s_Z, t_Z, a_Z, b_Z;

  mpz_init(e_Z);
  mpz_init(s_Z);
  mpz_init(t_Z);
  mpz_init(a_Z);
  mpz_init(b_Z);

  mpz_set_si(a_Z, a);
  mpz_set_si(b_Z, b);

  mpz_gcdext(e_Z, s_Z, t_Z, a_Z, b_Z);

  ASSERT(mpz_cmp_ui(e_Z, 0) > 0);
  ASSERT(mpz_sizeinbase(e_Z, 2) < 63);
  ASSERT(mpz_sizeinbase(s_Z, 2) < 63);
  ASSERT(mpz_sizeinbase(t_Z, 2) < 63);

  * e = mpz_get_si(e_Z);
  * s = mpz_get_si(s_Z);
  * t = mpz_get_si(t_Z);

  mpz_clear(e_Z);
  mpz_clear(s_Z);
  mpz_clear(t_Z);
  mpz_clear(a_Z);
  mpz_clear(b_Z);
}


