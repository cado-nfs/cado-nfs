#include <stdio.h>
#ifndef EC_ARITH_MONTGOMERY_H_
#define EC_ARITH_MONTGOMERY_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "ec_arith_common.h"

#ifdef ECM_COUNT_OPS
#include "ec_arith_cost.h"
#include "macros.h"
static unsigned int _count_montgomery_dadd, _count_montgomery_dbl;
#define MONTGOMERY_COUNT_OPS_M _count_montgomery_dadd * MONTGOMERY_dADD \
                             + _count_montgomery_dbl * MONTGOMERY_DBL
#define MONTGOMERY_COUNT_OPS_RESET() do {                   \
      _count_montgomery_dadd = _count_montgomery_dbl = 0;   \
    } while (0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Montgomery elliptic curves
 *
 * XZ-only coordinates, with equation:
 *    B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2
 *
 * Constant needed in computation: b = (A+2)/4
 *
 */

/* Compute A = 4*b-2. A and b can be the same variable. */
#define montgomery_A_from_b MOD_APPEND_TYPE(montgomery_A_from_b)
static inline void
montgomery_A_from_b (residue_t A, const residue_t b, const modulus_t m)
{
  mod_add (A, b, b, m);    /* A <- b+b = 2b */
  mod_add (A, A, A, m);    /* A <- 4b */
  mod_sub_ul (A, A, 2, m); /* A <- 4b-2 */
}

#define montgomery_curve_fprintf MOD_APPEND_TYPE(montgomery_curve_fprintf)
static inline void
montgomery_curve_fprintf (FILE *out, const char *prefix, residue_t A,
                          ec_point_t P, const modulus_t m)
{
  modint_t cc;
  const char *pre = (prefix == NULL) ? "" : prefix;

  mod_intinit (cc);
  mod_get_int (cc, A, m);

  mod_fprintf (out, "%sMontgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2\n"
                    "%sA = 0x%" PRIMODx "\n", pre, pre, MOD_PRINT_INT (cc));

  mod_intclear (cc);

  if (P)
  {
    fprintf (out, "%swith point (X::Z) = ", pre);
    ec_point_fprintf (out, P, MONTGOMERY_xz, m);
    fputc ('\n', out);
  }
}

/* Set P to zero (the neutral point): (0::0) */
static inline void
montgomery_point_set_zero (ec_point_t P, const modulus_t m)
{
  mod_set0 (P->x, m);
  mod_set0 (P->z, m);
}

/* Set Q to the same point as P but with z = 1.
 * Return 1 if it worked, 0 if the computation of the modular inverse of P->z
 * failed.
 * P and Q can be the same variables
 */
#define montgomery_point_to_affine MOD_APPEND_TYPE(montgomery_point_to_affine)
static inline int
montgomery_point_to_affine (ec_point_t Q, ec_point_t P, const modulus_t m)
{
  residue_t t;
  mod_init_noset0 (t, m);

  int ret = mod_inv (t, P->z, m);
  if (ret)
  {
    mod_mul (Q->x, P->x, t, m);
    mod_set1 (Q->z, m);
  }

  mod_clear (t, m);
  return ret;
}

/* Convert the point P in affine before printing it (mostly used for debug) */
#define montgomery_point_fprintf_affine MOD_APPEND_TYPE(montgomery_point_fprintf_affine)
static inline void
montgomery_point_fprintf_affine (FILE *out, ec_point_t P, const modulus_t m)
{
  ec_point_t Paff;
  ec_point_init (Paff, m);
  int ret = montgomery_point_to_affine (Paff, P, m);
  if (ret)
    ec_point_fprintf (out, Paff, MONTGOMERY_xz, m);
  else
    fprintf (out, "(not a affine point)");
  ec_point_clear (Paff, m);
}

#define montgomery_point_from_edwards_point MOD_APPEND_TYPE(montgomery_point_from_edwards_point)
static inline void
montgomery_point_from_edwards_point (ec_point_t PM, ec_point_t PE, int out_aff,
                                     const modulus_t m)
{
  mod_add (PM->x, PE->z, PE->y, m);
  mod_sub (PM->z, PE->z, PE->y, m);
  if (out_aff)
    montgomery_point_to_affine (PM, PM, m);
}

/* montgomery_dbl (Q, P)
 *     Q <- 2*P
 * It is permissible to let P and Q use the same memory.
 * Cost:
 *    3M + 2S + 4add/sub
 */
#define montgomery_dbl MOD_APPEND_TYPE(montgomery_dbl)
static inline void
montgomery_dbl (ec_point_t Q, const ec_point_t P, const modulus_t m,
             const residue_t b)
{
  residue_t u, v, w;

#if ELLM_SAFE_ADD
  if (mod_is0 (P->z, m))
    {
      ASSERT (mod_is0 (P->x, m));
    }
#endif

#ifdef ECM_COUNT_OPS
  _count_montgomery_dbl++;
#endif

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_add (u, P->x, P->z, m);
  mod_sqr (u, u, m);          /* u = (x + z)^2 */
  mod_sub (v, P->x, P->z, m);
  mod_sqr (v, v, m);          /* v = (x - z)^2 */
  mod_mul (Q->x, u, v, m);    /* x2 = (x^2 - z^2)^2 */
  mod_sub (w, u, v, m);       /* w = 4 * x * z */
  mod_mul (u, w, b, m);       /* u = x * z * (A + 2) */
  mod_add (u, u, v, m);       /* u = x^2 + x * z * A + z^2 */
  mod_mul (Q->z, w, u, m);    /* Q_z = (4xz) * (x^2 + xzA + z^2) */

#if ELLM_SAFE_ADD
  if (mod_is0 (Q->z, m))
    mod_set0 (Q->x, m);
#endif

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* montgomery_dadd (R, P, Q, D)
 *     R <- P+Q with D=P-Q or D=Q-P
 * R may be identical to P, Q and/or D.
 * This function assumes that P !~= Q, i.e. that there is
 * no t!=0 so that P->x = t*Q->x and P->z = t*Q->z, for otherwise the result
 * is (0:0) although it shouldn't be (which actually is good for factoring!).
 * Cost:
 *    4M + 2S + 6add/sub
 */
#define montgomery_dadd MOD_APPEND_TYPE(montgomery_dadd)
static inline void
montgomery_dadd (ec_point_t R, const ec_point_t P, const ec_point_t Q,
          const ec_point_t D, MAYBE_UNUSED const residue_t b,
          const modulus_t m)
{
  residue_t u, v, w;

#ifdef ECM_COUNT_OPS
  _count_montgomery_dadd++;
#endif

#if ELLM_SAFE_ADD
  /* Handle case where at least one input point is point at infinity */
  if (mod_is0 (P->z, m))
    {
      ASSERT (mod_is0 (P->x, m));
      ellM_set (R, Q, m);
      return;
    }
  if (mod_is0 (Q->z, m))
    {
      ASSERT (mod_is0 (Q->x, m));
      ellM_set (R, P, m);
      return;
    }
#endif

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_sub (u, P->x, P->z, m);
  mod_add (v, Q->x, Q->z, m);
  mod_mul (u, u, v, m);      /* u = (Px-Pz)*(Qx+Qz) */
  mod_add (w, P->x, P->z, m);
  mod_sub (v, Q->x, Q->z, m);
  mod_mul (v, w, v, m);      /* v = (Px+Pz)*(Qx-Qz) */
  mod_add (w, u, v, m);      /* w = 2*(Qx*Px - Qz*Pz)*/
  mod_sub (v, u, v, m);      /* v = 2*(Qz*Px - Qx*Pz) */
#if ELLM_SAFE_ADD
  /* Check if v == 0, which happens if P=Q or P=-Q.
     If P=-Q, set result to point at infinity.
     If P=Q, use ellM_double() instead.
     This test only works if P=Q on the pseudo-curve modulo N, i.e.,
     if N has several prime factors p, q, ... and P=Q or P=-Q on E_p but
     not on E_q, this test won't notice it. */
  if (mod_is0 (v, m))
    {
      mod_clear (w, m);
      mod_clear (v, m);
      mod_clear (u, m);
      /* Test if difference is point at infinity */
      if (mod_is0 (D->z, m))
        {
          ASSERT (mod_is0 (D->x, m));
          montgomery_dbl (R, P, m, b); /* Yes, points are identical, use doubling */
        }
      else
        montgomery_point_set_zero (R, m); /* Set result to point at infinity */
      return;
    }
#endif
  mod_sqr (w, w, m);          /* w = 4*(Qx*Px - Qz*Pz)^2 */
  mod_sqr (v, v, m);          /* v = 4*(Qz*Px - Qx*Pz)^2 */
  mod_set (u, D->x, m);       /* save D->x */
  mod_mul (R->x, w, D->z, m); /* may overwrite D->x */
  mod_mul (R->z, u, v, m);

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* montgomery_smul_ul (R, Rp1, P, k)
 *     R <- k*P
 *     Rp1 <- (k+1)*P       if Rp1 != NULL
 * R or Rp1 can be the same variable as P.
 * Cost:
 *    log2(k) DBL and log2(k)-1 dADD to compute R and Rp1
 *    computing only R save 1 DBL or 1 dADD depending on the parity of k
 *    If the second most significant bit of k is 0, we save 1 DBL
 */
#define montgomery_smul_ul MOD_APPEND_TYPE(montgomery_smul_ul)
static inline void
montgomery_smul_ul (ec_point_t R, ec_point_t Rp1, const ec_point_t P,
                    const unsigned long k, const modulus_t m, const residue_t b)
{
  ec_point_t T0, T1;

  ec_point_init (T0, m);
  ec_point_init (T1, m);

  if (k == 0UL)
  {
    montgomery_point_set_zero (R, m);
    if (Rp1)
      ec_point_set (Rp1, P, m, MONTGOMERY_xz);
  }
  else if (k == 1UL)
  {
    ec_point_set (R, P, m, MONTGOMERY_xz);
    if (Rp1)
      montgomery_dbl (Rp1, P, m, b);
  }
  else if (k == 2UL)
  {
    if (Rp1)
    {
      montgomery_dbl (T1, P, m, b);
      montgomery_dadd (Rp1, T1, P, P, b, m);
      ec_point_set (R, T1, m, MONTGOMERY_xz);
    }
    else
      montgomery_dbl (R, P, m, b);
  }
  else /* k >= 3 */
  {
    /* Montgomery Ladder */
    unsigned long mask;

    mask = ~(0UL);
    mask -= mask/2;   /* Now the most significant bit of i is set */
    while ((mask & k) == 0)
      mask >>= 1;

    /* Most significant bit of k is 1, do it outside the loop */
    ec_point_set (T0, P, m, MONTGOMERY_xz); /* starting value T0 = P */
    montgomery_dbl (T1, P, m, b);           /* starting value T1 = 2P */
    mask >>= 1;

    /* If the second most significant bit of k is 0, then we do the iteration
     * manually (to avoid to compute again 2P)
     * As k >= 3, we know that in this case k has at least 3 bits.
     */
    if (!(k & mask)) /* (T0,T1) <- (2*P, 3*P) */
    {
      ec_point_set (T0, T1, m, MONTGOMERY_xz);
      montgomery_dadd (T1, T1, P, P, b, m);
      mask >>= 1;
    }

    for ( ; mask > 1; mask >>= 1) /* loop invariant: T1-T0 = P */
    {
      if (k & mask) /* (T0,T1) <- (T1+T0, 2*T1) */
      {
        montgomery_dadd (T0, T1, T0, P, b, m);
        montgomery_dbl (T1, T1, m, b);
      }
      else /* (T0,T1) <- (2*T0, T1+T0) */
      {
        montgomery_dadd (T1, T1, T0, P, b, m);
        montgomery_dbl (T0, T0, m, b);
      }
    }

    /* Deal with least significant bit outside the loop */
    if (k & mask)
    {
      montgomery_dadd (R, T1, T0, P, b, m);
      if (Rp1)
        montgomery_dbl (Rp1, T1, m, b);
    }
    else
    {
      montgomery_dbl (R, T0, m, b);
      if (Rp1)
        montgomery_dadd (Rp1, T1, T0, P, b, m);
    }
  }

  ec_point_clear (T0, m);
  ec_point_clear (T1, m);
}

/* Return the order of a Mongomery curve, using the Jacobi symbol.
 * The curve is described by the curve coefficient A and a point P (which gives
 * the Jacobi symbol of the curve coefficient B).
 * P must not be a point of order 2 and P->z must be invertible modulo m.
 * This has complexity O(m).
 * This function only works for _ul arithmetic, so only the _ul version is
 * defined. Do we need the others versions (_15ul, _2ul2, _mpz) ? Not really,
 * because due to the complexity of this function, it will be too expensive to
 * use it with large value of m.
 */
#if defined(MOD_SIZE) && MOD_SIZE == 1
#define montgomery_curve_order MOD_APPEND_TYPE(montgomery_curve_order)
static inline unsigned long
montgomery_curve_order (residue_t A, ec_point_t P, const modulus_t m)
{
  residue_t x, t, one;
  unsigned long order, i;
  modint_t tm;
  int bchar;

  /* Compute x = P->x/P->z mod m */
  mod_init_noset0 (x, m);

  int ret = mod_inv (x, P->z, m);
  if (ret)
    mod_mul (x, x, P->x, m);
  else
  {
    mod_clear (x, m);
    return 0UL;
  }

  mod_init_noset0 (t, m);
  mod_init_noset0 (one, m);
  mod_set1 (one, m);
  mod_intinit (tm);


  /* Compute x^3 + A*x^2 + x and see if it is a square */
  mod_set (t, x, m);
  mod_add (t, t, A, m);
  mod_mul (t, t, x, m);
  mod_add (t, t, one, m);
  mod_mul (t, t, x, m);
  bchar = mod_jacobi (t, m);
  ASSERT (bchar != 0);

  order = 2; /* One for (0, 0, 1), one for the point at infinity */
  mod_getmod_int (tm, m);
  /* XXX Here we assume m fits in an unsigned long (or the loop never stops) */
  for (i = 1; mod_intcmp_ul (tm, i) > 0; i++)
  {
    mod_set_ul (x, i, m);
    mod_set (t, x, m);
    mod_add (t, t, A, m);
    mod_mul (t, t, x, m);
    mod_add (t, t, one, m);
    mod_mul (t, t, x, m);
    if (bchar == 1)
      order = order + (unsigned long) (1L + (long) mod_jacobi (t, m));
    else
      order = order + (unsigned long) (1L - (long) mod_jacobi (t, m));
      /* Brackets put like this to avoid signedness warning */
  }

  mod_intclear (tm);
  mod_clear (one, m);
  mod_clear (t, m);
  mod_clear (x, m);

  return order;
}
#endif /* defined(MOD_SIZE) && MOD_SIZE == 1 */


#ifdef __cplusplus
}
#endif

#endif /* EC_ARITH_MONTGOMERY_H_ */
