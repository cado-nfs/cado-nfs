#include <stdio.h>
#include <stdlib.h>
#ifndef EC_ARITH_WEIERSTRASS_H_
#define EC_ARITH_WEIERSTRASS_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "ec_arith_common.h"
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Short Weierstrass elliptic curves
 *
 * Affine coordinates, with equation
 *                   y^2 = x^3 + a*x + b
 * Projective coordinates, with equation:
 *                  Y^2*Z = X^3 + a*X*Z^Z + b*Z^3
 *
 * Curve coefficient needed in computation: a
 */


#define weierstrass_aff_curve_fprintf MOD_APPEND_TYPE(weierstrass_aff_curve_fprintf)
static inline void
weierstrass_aff_curve_fprintf (FILE *out, const char *prefix,
                                const residue_t a, ec_point_t P,
                                const modulus_t m)
{
  const char *pre = (prefix == NULL) ? "" : prefix;

  modint_t cc;
  mod_intinit (cc);

  mod_get_int (cc, a, m);

  mod_fprintf (out, "%sWeierstrass curve: y^2 = x^3 + a*x + b\n"
                    "%sa = 0x%" PRIMODx "\n", pre, pre, MOD_PRINT_INT (cc));


  mod_intclear (cc);

  if (P)
  {
    fprintf (out, "%swith point (x, y) = ", pre);
    ec_point_fprintf (out, P, SHORT_WEIERSTRASS_aff, m);
    fputc ('\n', out);
  }
}

/* Computes R=2P, with 1 inv, 4 muls (2 muls and 2 squares) and 8 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in R->x.
 *
 * It is permissible to let P and Q use the same memory.
 */
#define weierstrass_aff_dbl MOD_APPEND_TYPE(weierstrass_aff_dbl)
static int
weierstrass_aff_dbl (ec_point_t R, const ec_point_t P, const residue_t a,
                     const modulus_t m)
{
  int ret = 1;
  residue_t lambda, u, v;

  mod_init_noset0 (lambda, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_sqr (u, P->x, m);
  mod_add (v, u, u, m);
  mod_add (v, v, u, m);
  mod_add (v, v, a, m); /* 3x^2 + a */
  mod_add (u, P->y, P->y, m);
  ret = mod_inv (lambda, u, m);    /* 1/(2*y) */
  if (ret != 0)
  {
    mod_mul (lambda, lambda, v, m);
    mod_sqr (u, lambda, m);
    mod_sub (u, u, P->x, m);
    mod_sub (u, u, P->x, m);    /* x3 = u = lambda^2 - 2*x */
    mod_sub (v, P->x, u, m);
    mod_mul (v, v, lambda, m);
    mod_sub (R->y, v, P->y, m);
    mod_set (R->x, u, m);
  }
  else
    mod_set (R->x, u, m);

  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (lambda, m);
  return ret;
}

/* Computes R=P+Q, with 1 inv, 3 muls (2 muls and 1 square) and 6 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in R->x.
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
#define weierstrass_aff_add MOD_APPEND_TYPE(weierstrass_aff_add)
static int
weierstrass_aff_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
                     const residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;
  int ret;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (lambda, m);

  mod_sub (u, Q->x, P->x, m);
  ret = mod_inv (v, u, m);
  if (ret != 0)
  {
    mod_sub (u, Q->y, P->y, m);
    mod_mul (lambda, u, v, m);
    mod_sqr (u, lambda, m);
    mod_sub (u, u, P->x, m);
    mod_sub (u, u, Q->x, m);    /* x3 = u = lambda^2 - P->x - Q->x */
    mod_sub (v, P->x, u, m);
    mod_mul (v, v, lambda, m);
    mod_sub (R->y, v, P->y, m);
    mod_set (R->x, u, m);
  }
  else if (mod_equal (P->x, Q->x, m) && mod_equal (P->y, Q->y, m))
    ret = weierstrass_aff_dbl (R, P, a, m);
  else
    mod_set (R->x, u, m);

  mod_clear (lambda, m);
  mod_clear (v, m);
  mod_clear (u, m);
  return ret;
}


/* Convert the Montgomery curve B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2 with a valid
 * point Pm into a affine Weierstrass curve y^2 = x^3 + a*x + b with a
 * valid affine point Pw.
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in Pw->x.
 *
 * The curve coefficient b of the short Weierstrass curve will not be computed.
 * a and A can be the same variable.
 * Pm and Pw can be the same variable.
 */
#define weierstrass_aff_from_montgomery MOD_APPEND_TYPE(weierstrass_aff_from_montgomery)
MAYBE_UNUSED
static int
weierstrass_aff_from_montgomery (residue_t a, ec_point_t Pw, const residue_t A,
                                 ec_point_t Pm, const modulus_t m)
{
  residue_t B, one, t, x;
  int ret;

  mod_init_noset0 (B, m);
  mod_init_noset0 (one, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (x, m);
  mod_set1 (one, m);

  ret = mod_inv (t, Pm->z, m);
  if (ret == 0)
  {
    fprintf (stderr, "%s: could not invert Z\n", __func__);
    mod_set (Pw->x, Pm->z, m);
  }
  else
  {
    mod_mul (x, Pm->x, t, m); /* x = X/Z */
    mod_add (B, x, A, m);
    mod_mul (B, B, x, m);
    mod_add (B, B, one, m);
    mod_mul (B, B, x, m); /* B = x^3 + A*x^2 + x */

    /* Now (x,1) is on the curve B*y^2 = x^3 + A*x^2 + x. */
    ret = mod_inv (Pw->y, B, m);    /* y = 1/B */
    if (ret == 0)
    {
      fprintf (stderr, "%s: could not invert B\n", __func__);
      mod_set (Pw->x, B, m);
    }
    else
    {
      mod_div3 (t, A, m);
      mod_add (Pw->x, x, t, m);
      mod_mul (Pw->x, Pw->x, Pw->y, m); /* x = (X + A/3)/B */
      mod_mul (a, t, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, Pw->y, m);
      mod_mul (a, a, Pw->y, m);         /* a = (1 - (A^2)/3)/B^2 */
    }
  }

  mod_clear (one, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (x, m);

  return ret;
}


/* Computes P<-eP, with double-and-add algorithm.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 0 if e*P is the point at infinity, else return nonzero.
 * If the point at infinity is due to a failed inversion, the non-invertible
 * value is returned in P->x.
 */
#define weierstrass_aff_smul_ui MOD_APPEND_TYPE(weierstrass_aff_smul_ui)
static int
weierstrass_aff_smul_ui (ec_point_t P, const unsigned long e, const residue_t a,
                         const modulus_t m)
{
  unsigned long i;
  ec_point_t T;
  int tfinite; /* Nonzero iff T is NOT point at infinity */

  if (e == 0)
    return 0; /* signal point at infinity */

  ec_point_init (T, m);

  i = ~(0UL);
  i -= i/2;   /* Now the most significant bit of i is set */
  while ((i & e) == 0)
    i >>= 1;

  ec_point_set (T, P, m, SHORT_WEIERSTRASS_aff);
  tfinite = 1;
  i >>= 1;

  while (i > 0)
  {
    if (tfinite)
      tfinite = weierstrass_aff_dbl (T, T, a, m);
    if (e & i)
    {
      if (tfinite)
        tfinite = weierstrass_aff_add (T, T, P, a, m);
      else
      {
        ec_point_set (T, P, m, SHORT_WEIERSTRASS_aff);
        tfinite = 1;
      }
    }
    i >>= 1;
  }

  if (tfinite)
    ec_point_set (P, T, m, SHORT_WEIERSTRASS_aff);
  else
    mod_set (P->x, T->x, m);

  ec_point_clear (T, m);

  return tfinite;
}

/* Return the order of the point P on the Weierstrass curve defined by the curve
 * coefficient a modulo the prime m.
 * Looks for i in Hasse interval so that i*P = O, has complexity O(m^(1/4)).
 * If the group order is known to be == r (mod m), this can be supplied in
 * the variables "known_r" and" known_m".
 * XXX For now, this function only works for _ul arithmetic, so only the _ul
 * version is defined.
 */
#if defined(MOD_SIZE) && MOD_SIZE == 1
#define weierstrass_aff_point_order MOD_APPEND_TYPE(weierstrass_aff_point_order)
unsigned long
weierstrass_aff_point_order (const residue_t a, ec_point_t P,
                             const unsigned long known_m,
                             const unsigned long known_r, const modulus_t m,
                             const int verbose)
{
  ec_point_t Pi, Pg, Q, *baby;
  residue_t x, d;
  unsigned long min, max, i, j, order, cof, p;
  unsigned long giant_step, giant_min, baby_len;
  modint_t tm;

  ASSERT (known_r < known_m);

  mod_intinit (tm);
  mod_getmod_int (tm, m);

  if (verbose >= 2)
  {
    printf ("%s: ", __func__);
    weierstrass_aff_curve_fprintf (stdout, NULL, a, P, m);
  }

  mod_init (x, m);
  mod_init (d, m);
  ec_point_init (Pi, m);
  ec_point_init (Pg, m);

  /* XXX Here we assume m fits in an unsigned long */
  i = (unsigned long) (2. * sqrt((double) mod_intget_ul(tm)));
  min = mod_intget_ul(tm) - i;
  max = mod_intget_ul(tm) + i;

  /* Giant steps visit values == r (mod m), baby steps values == 0 (mod m) */
  giant_step = ceil(sqrt(2.*(double)i / (double) known_m));
  /* Round up to multiple of m */
  giant_step = ((giant_step - 1) / known_m + 1) * known_m;

  /* We test Pi +- Pj, where Pi = (giant_min + i*giant_step), i >= 0,
     and Pj = j*P, 0 <= j <= giant_step / 2.
     To ensure we can find all values >= min, ensure
     giant_min <= min + giant_step / 2.
     We also want giant_min == r (mod m) */
  giant_min = ((min + giant_step / 2) / known_m) * known_m + known_r;
  if (giant_min > min + giant_step / 2)
    giant_min -= known_m;
  if (verbose >= 2)
    printf ("known_m = %lu, known_r = %lu, giant_step = %lu, "
            "giant_min = %lu\n", known_m, known_r, giant_step, giant_min);

  baby_len = giant_step / known_m / 2 + 1;
  baby = (ec_point_t *) malloc (baby_len * sizeof (ec_point_t));
  for (i = 0; i < baby_len; i++)
    ec_point_init (baby[i], m);

  ec_point_set (Pg, P, m, SHORT_WEIERSTRASS_aff);
  i = known_m;
  if (weierstrass_aff_smul_ui (Pg, i, a, m) == 0) /* Pg = m*P for now */
    goto found_inf;

  if (1 < baby_len)
    ec_point_set (baby[1], Pg, m, SHORT_WEIERSTRASS_aff);

  if (2 < baby_len)
    {
      if (weierstrass_aff_dbl (Pi, Pg, a, m) == 0)
        {
          i = 2 * known_m;
          goto found_inf;
        }
      ec_point_set (baby[2], Pi, m, SHORT_WEIERSTRASS_aff);
    }

  for (i = 3; i < baby_len; i++)
    {
      if (weierstrass_aff_add (Pi, Pi, Pg, a, m) == 0)
        {
          i *= known_m;
          goto found_inf;
        }
      ec_point_set (baby[i], Pi, m, SHORT_WEIERSTRASS_aff);
    }

  /* Now compute the giant steps in [giant_min, giant_max] */
  i = giant_step;
  ec_point_set (Pg, P, m, SHORT_WEIERSTRASS_aff);
  if (weierstrass_aff_smul_ui (Pg, i, a, m) == 0)
    goto found_inf;

  i = giant_min;
  ec_point_set (Pi, P, m, SHORT_WEIERSTRASS_aff);
  if (weierstrass_aff_smul_ui (Pi, i, a, m) == 0)
    goto found_inf;

  while (i <= max + giant_step - 1)
    {
      /* Compare x-coordinate with stored baby steps. This makes it
         O(sqrt(p)) complexity, strictly speaking. */
      for (j = 1; j < baby_len; j++)
        if (mod_equal (Pi[0].x, baby[j]->x, m))
          {
            if (mod_equal (Pi[0].y, baby[j]->y, m))
              i -= j * known_m; /* Equal, so iP = jP and (i-j)P = 0 */
            else
              {
                mod_neg (Pi[0].y, Pi[0].y, m);
                if (mod_equal (Pi[0].y, baby[j]->y, m))
                  i += j * known_m; /* Negatives, so iP = -jP and (i+j)P = 0 */
                else
                  {
                    fprintf (stderr, "Matching x-coordinates, but y neither "
                             "equal nor negatives\n");
                    abort();
                  }
              }
            goto found_inf;
          }

      i += giant_step;
      if (!weierstrass_aff_add (Pi, Pi, Pg, a, m))
        goto found_inf;
    }

  if (i > max)
  {
      fprintf (stderr, "ec_order: Error, no match found for p = %lu, "
               "min = %lu, max = %lu, giant_step = %lu, giant_min = %lu\n",
               mod_intget_ul(tm), min, max, giant_step, giant_min);
      abort ();
  }

found_inf:
  /* Check that i is a multiple of the order */
  ec_point_set (Pi, P, m, SHORT_WEIERSTRASS_aff);
  if (weierstrass_aff_smul_ui (Pi, i, a, m) != 0)
    {
      modint_t tx1, ty1;
      for (i = 0; i < giant_step; i++)
          ec_point_clear (baby[i], m);
      free (baby);
      baby = NULL;
      mod_intinit (tx1);
      mod_intinit (ty1);
      mod_get_int (tx1, P[0].x, m);
      mod_get_int (ty1, P[0].y, m);
#ifndef MODMPZ_MAXBITS
      fprintf (stderr, "ec_order: Error, %ld*(%ld, %ld) (mod %ld) is "
               "not the point at infinity\n",
               i, tx1[0], ty1[0], tm[0]);
#endif
      mod_intclear (tx1);
      mod_intclear (ty1);
      return 0UL;
    }

  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  /* For each prime factor of the order, reduce the exponent of
     that prime factor as far as possible */

  cof = order = i;
  for (p = 2; p * p <= cof; p += 1 + p%2)
    if (cof % p == 0)
      {
        ASSERT (order % p == 0);
        /* Remove all factors of p */
        for (order /= p, cof /= p; order % p == 0; order /= p)
          {
            ASSERT(cof % p == 0);
            cof /= p;
          }
        ASSERT (cof % p != 0);

        /* Add factors of p again one by one, stopping when we hit
           point at infinity */
        ec_point_set (Pi, P, m, SHORT_WEIERSTRASS_aff);
        if (weierstrass_aff_smul_ui (Pi, order, a, m) != 0)
          {
            order *= p;
            while (weierstrass_aff_smul_ui (Pi, p, a, m) != 0)
              order *= p;
          }
      }
  /* Now cof is 1 or a prime */
  if (cof > 1)
    {
      ec_point_set (Pi, P, m, SHORT_WEIERSTRASS_aff);
      ASSERT (order % cof == 0);
      if (weierstrass_aff_smul_ui (Pi, order / cof, a, m) == 0)
        order /= cof;
    }


  /* One last check that order divides real order */
  ec_point_set (Pi, P, m, SHORT_WEIERSTRASS_aff);
  if (weierstrass_aff_smul_ui (Pi, order, a, m) != 0)
    {
      modint_t tx1, ty1;
      mod_intinit (tx1);
      mod_intinit (ty1);
      mod_get_int (tx1, P[0].x, m);
      mod_get_int (ty1, P[0].y, m);
      fprintf (stderr, "ec_order: Error, final order %ld is wrong\n",
               order);
      mod_intclear (tx1);
      mod_intclear (ty1);
      abort ();
    }

  for (i = 0; i < giant_step; i++)
    ec_point_clear (baby[i], m);
  free (baby);
  baby = NULL;
  mod_clear (x, m);
  mod_clear (d, m);
  mod_intclear (tm);
  ec_point_clear (Pi, m);
  ec_point_clear (Pg, m);
  ec_point_clear (Q, m);

  return order;
}
#endif /* defined(MOD_SIZE) && MOD_SIZE == 1 */

/***************************** projective version *****************************/

/* Set P to zero (the neutral point): (0:1:0) */
#define weierstrass_proj_point_set_zero MOD_APPEND_TYPE(weierstrass_proj_point_set_zero)
static inline void
weierstrass_proj_point_set_zero (ec_point_t P, const modulus_t m)
{
  mod_set0 (P->x, m);
  mod_set1 (P->y, m);
  mod_set0 (P->z, m);
}

/* Computes R=2P, with ? muls (? muls and ? squares) and ? add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * It is permissible to let P and Q use the same memory.
 */
#define weierstrass_proj_dbl MOD_APPEND_TYPE(weierstrass_proj_dbl)
static void
weierstrass_proj_dbl (ec_point_t R, const ec_point_t P, const residue_t a,
                      const modulus_t m)
{
  residue_t xx, zz, w, u, s, ss, sss, r, rr, B, t, h;
  /* TODO reduce number of var */

  mod_init_noset0 (xx, m);
  mod_init_noset0 (zz, m);
  mod_init_noset0 (w, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (s, m);
  mod_init_noset0 (ss, m);
  mod_init_noset0 (sss, m);
  mod_init_noset0 (r, m);
  mod_init_noset0 (rr, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (h, m);

  mod_sqr (xx, P->x, m);
  mod_sqr (zz, P->z, m);

  mod_mul (w, a, zz, m);
  mod_add (t, xx, xx, m);
  mod_add (t, t, xx, m);
  mod_add (w, w, t, m);

  mod_mul (s, P->y, P->z, m);
  mod_add (s, s, s, m);
  mod_sqr (ss, s, m);
  mod_mul (sss, ss, s, m);
  mod_mul (r, P->y, s, m);
  mod_sqr (rr, r, m);

  mod_add (B, P->x, r, m);
  mod_sqr (B, B, m);
  mod_sub (B, B, xx, m);
  mod_sub (B, B, rr, m);
  mod_sqr (h, w, m);
  mod_sub (h, h, B, m);
  mod_sub (h, h, B, m);

  mod_mul (R->x, h, s, m);
  mod_sub (R->y, B, h, m);
  mod_mul (R->y, R->y, w, m);
  mod_sub (R->y, R->y, rr, m);
  mod_sub (R->y, R->y, rr, m);
  mod_set (R->z, sss, m);

  mod_clear (xx, m);
  mod_clear (zz, m);
  mod_clear (w, m);
  mod_clear (u, m);
  mod_clear (s, m);
  mod_clear (ss, m);
  mod_clear (sss, m);
  mod_clear (r, m);
  mod_clear (rr, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (h, m);
}

/* Computes R=P+Q, with 14 muls (12 muls and 2 squares) and 7 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
#define weierstrass_proj_add MOD_APPEND_TYPE(weierstrass_proj_add)
static void
weierstrass_proj_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
                      const modulus_t m)
{
  residue_t t1, t2, t3, u, uu, v, vv, vvv, r, A; /* TODO reduce number of var */

  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (uu, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (vv, m);
  mod_init_noset0 (vvv, m);
  mod_init_noset0 (r, m);
  mod_init_noset0 (A, m);

  mod_mul (t1, P->x, Q->z, m);
  mod_mul (t2, P->y, Q->z, m);
  mod_mul (t3, P->z, Q->z, m);

  mod_mul (u, Q->y, P->z, m);
  mod_sub (u, u, t2, m);
  mod_sqr (uu, u, m);

  mod_mul (v, Q->x, P->z, m);
  mod_sub (v, v, t1, m);
  mod_sqr (vv, v, m);
  mod_mul (vvv, vv, v, m);

  mod_mul (r, vv, t1, m);

  mod_mul (A, uu, t3, m);
  mod_sub (A, A, vvv, m);
  mod_sub (A, A, r, m);
  mod_sub (A, A, r, m);

  mod_mul (R->x, v, A, m);

  mod_sub (R->y, r, A, m);
  mod_mul (R->y, R->y, u, m);
  mod_mul (t2, t2, vvv, m);
  mod_sub (R->y, R->y, t2, m);

  mod_mul (R->z, vvv, t3, m);

  mod_clear (t1, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
  mod_clear (u, m);
  mod_clear (uu, m);
  mod_clear (v, m);
  mod_clear (vv, m);
  mod_clear (vvv, m);
  mod_clear (r, m);
  mod_clear (A, m);
}

/* Computes P<-eP, with double-and-add algorithm.
 *    - m : modulus
 *    - a : curve coefficient
 */
#define weierstrass_proj_smul_ui MOD_APPEND_TYPE(weierstrass_proj_smul_ui)
static void
weierstrass_proj_smul_ui (ec_point_t P, const unsigned long e,
                          const residue_t a, const modulus_t m)
{
  if (e == 0)
    weierstrass_proj_point_set_zero (P, m);
  else if (e > 1)
  {
    unsigned long i;
    ec_point_t T;

    ec_point_init (T, m);

    i = ~(0UL);
    i -= i/2;   /* Now the most significant bit of i is set */
    while ((i & e) == 0)
      i >>= 1;

    ec_point_set (T, P, m, SHORT_WEIERSTRASS_proj);
    i >>= 1; /* skip most significant bit of e */

    for (; i > 0; i >>= 1)
    {
      weierstrass_proj_dbl (T, T, a, m);
      if (e & i)
        weierstrass_proj_add (T, T, P, m);
    }

    ec_point_set (P, T, m, SHORT_WEIERSTRASS_proj);
    ec_point_clear (T, m);
  }
  /* else do nothing for e == 1 */
}

#ifdef __cplusplus
}
#endif


#endif /* EC_ARITH_WEIERSTRASS_H_ */
