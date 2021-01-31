#include <stdio.h>
#ifndef EC_ARITH_EDWARDS_H_
#define EC_ARITH_EDWARDS_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "ec_arith_common.h"

/* a=-1 Twisted Edwards elliptic curves
 *
 * Extended coordinates, with equations
 *                  -X^2 + Y^2 = Z^2 + d*T^2
 *                  X*Y = Z*T
 * Projective coordinates, with equation:
 *                  -X^2*Z^2 + Y^2*Z^2 = Z^4+d*X^2*Y^2
 *
 * Constant needed in computation: none
 */

#ifdef ECM_COUNT_OPS
#include "ec_arith_cost.h"
static unsigned int _count_edwards_dbl, _count_edwards_dblext,
                    _count_edwards_tpl, _count_edwards_tplext,
                    _count_edwards_add, _count_edwards_addext,
                    _count_edwards_addmont;
#define EDWARDS_COUNT_OPS_M _count_edwards_dbl * EDWARDS_DBL \
                          + _count_edwards_dblext * EDWARDS_DBLext \
                          + _count_edwards_tpl * EDWARDS_TPL \
                          + _count_edwards_tplext * EDWARDS_TPLext \
                          + _count_edwards_add * EDWARDS_ADD \
                          + _count_edwards_addext * EDWARDS_ADDext \
                          + _count_edwards_addmont * EDWARDS_ADDmontgomery
#define EDWARDS_COUNT_OPS_RESET() do { \
      _count_edwards_dbl = _count_edwards_dblext = _count_edwards_tpl = 0;    \
      _count_edwards_tplext = _count_edwards_add = _count_edwards_addext = 0; \
      _count_edwards_addmont = 0;                                             \
    } while (0)
#endif

/* #define SAFE_TWISTED_EDWARDS_TO_MONTGOMERY */

#ifdef __cplusplus
extern "C" {
#endif

/* Compute d = -(A-2)/(A+2). A and d can be the same variable. */
#define edwards_d_from_montgomery_A MOD_APPEND_TYPE(edwards_d_from_montgomery_A)
static inline void
edwards_d_from_montgomery_A (residue_t d, const residue_t A, const modulus_t m)
{
  residue_t (t);
  mod_init_noset0 (t, m);

  mod_add_ul (t, A, 2UL, m);
  mod_inv (d, t, m);
  mod_sub_ul (t, t, 4UL, m); /* t = A-2 */
  mod_mul (d, d, t, m);
  mod_neg (d, d, m);

  mod_clear (t, m);
}

#define edwards_ext_curve_fprintf MOD_APPEND_TYPE(edwards_ext_curve_fprintf)
static inline void
edwards_ext_curve_fprintf (FILE *out, const char *prefix, residue_t d,
                           ec_point_t P, const modulus_t m)
{
  const char *pre = (prefix == NULL) ? "" : prefix;

  modint_t cc;
  mod_intinit (cc);

  mod_get_int (cc, d, m);

  mod_fprintf (out, "%sTwisted Edwards curve: -X^2 + Y^2 = Z^2 + d*T^2\n"
                    "%sXY = ZT (extended coordinates)\n%sd = 0x%" PRIMODx "\n",
                    pre, pre, pre, MOD_PRINT_INT (cc));

  mod_intclear (cc);

  if (P)
  {
    fprintf (out, "%swith point (X:Y:Z:T) = ", pre);
    ec_point_fprintf (out, P, TWISTED_EDWARDS_ext, m);
    fputc ('\n', out);
  }
}

/* Set P to zero (the neutral point):
 *    (0:1:1:0)     if coord is TWISTED_EDWARDS_ext
 *    (0:1:1)       if coord is TWISTED_EDWARDS_proj
 */
#define edwards_point_set_zero MOD_APPEND_TYPE(edwards_point_set_zero)
static inline void
edwards_point_set_zero (ec_point_t P, const modulus_t m,
                        const ec_point_coord_type_t coord)
{
  ASSERT_EXPENSIVE (coord==TWISTED_EDWARDS_ext || coord==TWISTED_EDWARDS_proj);
  mod_set0 (P->x, m);
  mod_set1 (P->y, m);
  mod_set1 (P->z, m);
  if (coord == TWISTED_EDWARDS_ext)
    mod_set0 (P->t, m);
}

#define edwards_neg MOD_APPEND_TYPE(edwards_neg)
static inline void
edwards_neg (ec_point_t Q, const ec_point_t P, const modulus_t m)
{
  mod_neg (Q->x, P->x, m);
  mod_set (Q->y, P->y, m);
  mod_set (Q->z, P->z, m);
  mod_neg (Q->t, P->t, m);
}

/* edwards_addsub (R:output_type, P:edwards_ext, Q:edwards_ext, sub,output_type)
 *     R <- P+Q      if sub == 0
 *     R <- P-Q      if sub != 0
 *     output_type can be edwards_proj, edwards_ext or montgomery
 * R can be the same variable as P or Q.
 * All coordinates of the output point R can be modified (because they may be
 * used as temporary variables).
 * Cost:
 *    7M + 8add + 2*2         if output_type == TWISTED_EDWARDS_proj
 *    8M + 8add + 2*2         if output_type == TWISTED_EDWARDS_ext
 *    4M + 10add + 2*2        if output_type == TWISTED_EDWARDS_ext
 * Notations in the comments come from:
 *    https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-add-2008-hwcd-4
 * Source: Hisil–Wong–Carter–Dawson, 2008, section 3.2 of
 *    http://eprint.iacr.org/2008/522
 */
#define edwards_addsub MOD_APPEND_TYPE(edwards_addsub)
static inline void
edwards_addsub (ec_point_t R, const ec_point_t P, const ec_point_t Q, int sub,
                const modulus_t m, const ec_point_coord_type_t output_type)
{
  ASSERT_EXPENSIVE (output_type == TWISTED_EDWARDS_ext ||
                    output_type == TWISTED_EDWARDS_proj ||
                    output_type == MONTGOMERY_xz);

#ifdef ECM_COUNT_OPS
  if (output_type == TWISTED_EDWARDS_proj)
    _count_edwards_add++;
  else if (output_type == TWISTED_EDWARDS_ext)
    _count_edwards_addext++;
  else /* if (output_type == MONTGOMERY_xz) */
    _count_edwards_addmont++;
#endif

  residue_t u0, u1, u2;

  mod_init_noset0 (u0, m);
  mod_init_noset0 (u1, m);
  mod_init_noset0 (u2, m);

  if (sub)
  {
    mod_add (u1, Q->y, Q->x, m);    /* u1 <-      (Y2+X2) */
    mod_sub (u0, Q->y, Q->x, m);    /* u0 <-      (Y2-X2) */
  }
  else
  {
    mod_add (u0, Q->y, Q->x, m);    /* u0 <-      (Y2+X2) */
    mod_sub (u1, Q->y, Q->x, m);    /* u1 <-      (Y2-X2) */
  }

  mod_sub (u2, P->y, P->x, m);      /* u2 <-      (Y1-X1) */
  mod_mul (u0, u0, u2, m);          /* u0 <- A := (Y1-X1)*(Y2+/-X2) */
  mod_add (u2, P->y, P->x, m);      /* u2 <-      (X1+Y1) */
  mod_mul (u1, u1, u2, m);          /* u1 <- B := (Y1+X1)*(Y2-/+X2) */
  mod_sub (R->x, u1, u0, m);        /* Rx <- F := B-A */
  mod_add (R->y, u1, u0, m);        /* Ry <- G := B+A */

  mod_add (u1, P->z, P->z, m);      /* u1 <-      2*Z1 */
  mod_mul (u1, u1, Q->t, m);        /* u1 <- C := 2*Z1*T2 */
  mod_add (u2, Q->z, Q->z, m);      /* u2 <-      2*Z2 */
  mod_mul (u2, u2, P->t, m);        /* u2 <- D := 2*T1*Z2 */
  if (sub)
  {
    mod_add (u0, u2, u1, m);        /* u0 <- H := D+C */
    mod_sub (u1, u2, u1, m);        /* u1 <- E := D-C */
  }
  else
  {
    mod_sub (u0, u2, u1, m);        /* u0 <- H := D-C */
    mod_add (u1, u2, u1, m);        /* u1 <- E := D+C */
  }


  if (output_type == TWISTED_EDWARDS_ext || output_type == TWISTED_EDWARDS_proj)
  {
    mod_mul (R->z, R->x, R->y, m);  /* Rz <- Z3 := F*G */
    mod_mul (R->x, R->x, u1, m);    /* Rx <- X3 := E*F */
    mod_mul (R->y, R->y, u0, m);    /* Ry <- Y3 := G*H */
    if (output_type == TWISTED_EDWARDS_ext)
      mod_mul (R->t, u0, u1, m);    /* Rt <- T3 := E*H */
  }
  else /* if (output_type == MONTGOMERY_xz) */
  {
#ifdef SAFE_TWISTED_EDWARDS_TO_MONTGOMERY
    mod_sub (R->z, R->x, u0, m);    /* Rz <-       F-H */
    mod_mul (R->z, R->z, R->y, m);  /* Rz <-       G*(F-H) */
    mod_add (R->x, R->x, u0, m);    /* Rx <-       F+H */
    mod_mul (R->x, R->x, R->y, m);  /* Rx <-       G*(F+H) */
#else
    /* Conversion: Edwards completed -> Montgomery XZ-only
     *                 ((E:G),(H,F)) -> (F+H :: F-H)
     *
     * Note: the forgotten Y-coordinate is G(F+H)/E.
     *
     * The above map is correct for every points except ((0:1),(1:1)) which is
     * sent to (2 :: 0) instead of the point at infinity (0 :: 0).
     * Nevertheless, in our context, computing with the point (2 :: 0)
     * in place of the point at infinity (0 :: 0) produces the same result
     * because both points have Z = 0 and this property is preserved after a
     * doubling dDBL or a differential addition dADD of two points having Z = 0.
     */
    mod_sub (R->z, R->x, u0, m);    /* Rz <-       F-H */
    mod_add (R->x, R->x, u0, m);    /* Rx <-       F+H */
#endif
  }

  mod_clear (u0, m);
  mod_clear (u1, m);
  mod_clear (u2, m);
}

#define edwards_add MOD_APPEND_TYPE(edwards_add)
static inline void
edwards_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
             const modulus_t m, const ec_point_coord_type_t output_type)
{
  edwards_addsub (R, P, Q, 0, m, output_type);
}

#define edwards_sub MOD_APPEND_TYPE(edwards_sub)
static inline void
edwards_sub (ec_point_t R, const ec_point_t P, const ec_point_t Q,
             const modulus_t m, const ec_point_coord_type_t output_type)
{
  edwards_addsub (R, P, Q, 1, m, output_type);
}

/* edwards_dbl (R:output_type, P:edwards_proj, output_type)
 *     R <- 2*P
 *     output_type can be edwards_proj or edwards_ext
 * R can be the same variable as P.
 * All coordinates of the output point R can be modified (because they may be
 * used as temporary variables).
 * Cost:
 *    3M + 4S + 5add + 1*2         if output_type == TWISTED_EDWARDS_proj
 *    4M + 4S + 5add + 1*2         if output_type == TWISTED_EDWARDS_ext
 * Notations in the comments come from:
 *    https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp
 * Source: Bernstein–Birkner–Joye–Lange–Peters, 2008
 *    http://eprint.iacr.org/2008/013.
 */
#define edwards_dbl MOD_APPEND_TYPE(edwards_dbl)
static inline void
edwards_dbl (ec_point_t R, const ec_point_t P,
             const modulus_t m, const ec_point_coord_type_t output_type)
{
  ASSERT_EXPENSIVE (output_type == TWISTED_EDWARDS_ext ||
                    output_type == TWISTED_EDWARDS_proj);

#ifdef ECM_COUNT_OPS
  if (output_type == TWISTED_EDWARDS_proj)
    _count_edwards_dbl++;
  else /* if (output_type == TWISTED_EDWARDS_ext) */
    _count_edwards_dblext++;
#endif

  residue_t u0, u1;

  mod_init_noset0 (u0, m);
  mod_init_noset0 (u1, m);

  mod_sqr (u0, P->x, m);            /* u0 <-  C := X1^2 */
  mod_sqr (u1, P->y, m);            /* u1 <-  D := Y1^2 */
  mod_add (R->x, P->x, P->y, m);    /* Rx <-       X1+Y1 */
  mod_sqr (R->x, R->x, m);          /* Rx <-  B := (X1+Y1)^2 */
  mod_add (R->y, u0, u1, m);        /* Ry <-       C+D */
  mod_sub (u0, u0, u1, m);          /* u0 <- -F := C-D */
  mod_sub (R->x, R->y, R->x, m);    /* Rx <-    := C+D-B */

  mod_sqr (u1, P->z, m);            /* u1 <-  H := Z1^2 */
  mod_add (u1, u1, u1, m);          /* u1 <-    := 2*H  */
  mod_add (u1, u0, u1, m);          /* u1 <- -J := -F + 2*H */

  if (output_type == TWISTED_EDWARDS_ext)
    mod_mul(R->t, R->x, R->y, m);   /* Rt <- T3 := (B-C-D)*(-C-D) */
  mod_mul (R->x, R->x, u1, m);      /* Rx <- X3 := (B-C-D) * J */
  mod_mul (R->y, R->y, u0, m);      /* Ry <- Y3 := F*(-C-D) */
  mod_mul (R->z, u0, u1, m);        /* Rz <- Z3 := F * J */

  mod_clear (u0, m);
  mod_clear (u1, m);
}


/* edwards_tpl (R:output_type, P:edwards_proj, output_type)
 *     R <- 3*P
 *     output_type can be edwards_proj or edwards_ext
 * R can be the same variable as P.
 * All coordinates of the output point R can be modified (because they may be
 * used as temporary variables).
 * Cost:
 *     9M + 3S + 7add + 2*2 + 1*-1      if output_type == TWISTED_EDWARDS_proj
 *    11M + 3S + 7add + 2*2 + 1*-1      if output_type == TWISTED_EDWARDS_ext
 * Notations in the comments come from:
 *    https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#tripling-tpl-2015-c
 * Source: 2015 Chuengsatiansup.
 */
#define edwards_tpl MOD_APPEND_TYPE(edwards_tpl)
static inline void
edwards_tpl (ec_point_t R, const ec_point_t P,
             const modulus_t m, const ec_point_coord_type_t output_type)
{
  ASSERT_EXPENSIVE (output_type == TWISTED_EDWARDS_ext ||
                    output_type == TWISTED_EDWARDS_proj);

#ifdef ECM_COUNT_OPS
  if (output_type == TWISTED_EDWARDS_proj)
    _count_edwards_tpl++;
  else /* if (output_type == TWISTED_EDWARDS_ext) */
    _count_edwards_tplext++;
#endif

  residue_t u0, u1, u2, u3;

  mod_init_noset0 (u0, m);
  mod_init_noset0 (u1, m);
  mod_init_noset0 (u2, m);
  mod_init_noset0 (u3, m);

  mod_sqr (u0, P->x, m);            /* u0 <-     := X1^2 */
  mod_neg (u0, u0, m);              /* u0 <- aXX := -X1^2 */
  mod_sqr (u1, P->y, m);            /* u1 <-  YY := Y1^2 */
  mod_add (u2, u1, u0, m);          /* u2 <-  Ap := YY+aXX */
  mod_sub (u3, u1, u0, m);          /* u3 <-        YY-aXX */
  mod_sqr (R->t, P->z, m);          /* Rt <-        Z1^2 */
  mod_add (R->t, R->t, R->t, m);    /* Rt <-        2*Z1^2 */
  mod_sub (R->t, R->t, u2, m);      /* Rt <-        2*Z1^2-Ap */
  mod_add (R->t, R->t, R->t, m);    /* Rt <-   B := 2*(2*Z1^2-Ap) */
  mod_mul (u2, u2, u3, m);          /* u2 <-  AA := Ap*(YY-aXX) */
  mod_mul (u0, u0, R->t, m);        /* u0 <-  xB := aXX*B */
  mod_mul (u1, u1, R->t, m);        /* u1 <-  yB := YY*B */
  mod_sub (R->t, u2, u1, m);        /* Rt <-   F := AA-yB */
  mod_add (u1, u1, u2, m);          /* u1 <-        yB+AA */
  mod_add (u3, u2, u0, m);          /* u3 <-   G := AA+xB */
  mod_sub (u2, u0, u2, m);          /* u2 <-        xB-AA */
  mod_mul (u1, P->x, u1, m);        /* u1 <-  xE := X1*(yB+AA) */
  mod_mul (u2, P->y, u2, m);        /* u2 <-  yH := Y1*(xB-AA) */
  mod_mul (u0, P->z, R->t, m);      /* u0 <-  zF := Z1*F */

  if (output_type == TWISTED_EDWARDS_proj)
  {
    mod_mul (R->x, u1, R->t, m);    /* Rx <-  X3 := xE*F */
    mod_mul (R->y, u2, u3, m);      /* Ry <-  Y3 := yH*G */
    mod_mul (R->z, u0, u3, m);      /* Rz <-  Z3 := zF*G */
  }
  else /* if (output_type == TWISTED_EDWARDS_ext) */
  {
    mod_mul (u3, P->z, u3, m);      /* u3 <-  zG := Z1*G */
    mod_mul (R->x, u1, u0, m);      /* Rx <-  X3 := xE*zF */
    mod_mul (R->y, u2, u3, m);      /* Ry <-  Y3 := yH*zG */
    mod_mul (R->z, u0, u3, m);      /* Rz <-  Z3 := zF*zG */
    mod_mul (R->t, u1, u2, m);      /* Rt <- T3 := xE*yH */
  }

  mod_clear (u0, m);
  mod_clear (u1, m);
  mod_clear (u2, m);
  mod_clear (u3, m);
}


/* edwards_smul_ui (R:edwards_ext, P:edwards_ext, k: unsigned long)
 *     R <- k*P
 * R can be the same variable as P.
 */
#define edwards_smul_ui MOD_APPEND_TYPE(edwards_smul_ui)
MAYBE_UNUSED
static inline void
edwards_smul_ui (ec_point_t R, const ec_point_t P, const unsigned long k,
                 const modulus_t m)
{
  if (k == 0UL)
    edwards_point_set_zero (R, m, TWISTED_EDWARDS_ext);
  else if (k == 1UL)
    ec_point_set (R, P, m, TWISTED_EDWARDS_ext);
  else if (k == 2UL)
    edwards_dbl (R, P, m, TWISTED_EDWARDS_ext);
  else if (k == 3UL)
    edwards_tpl (R, P, m, TWISTED_EDWARDS_ext);
  else
  {
    /* basic double-and-add */
    unsigned long mask;
    ec_point_t T;

    ec_point_init (T, m);

    mask = ~(0UL);
    mask -= mask/2;   /* Now the most significant bit of i is set */
    while ((mask & k) == 0)
      mask >>= 1;

    /* Most significant bit of k is 1, do it outside the loop */
    ec_point_set (T, P, m, TWISTED_EDWARDS_ext);
    mask >>= 1;

    for ( ; mask > 0; mask >>= 1)
    {
      edwards_dbl (T, T, m, TWISTED_EDWARDS_ext);
      if (k & mask) /* output in ext only on the last iteration */
        edwards_add (T, T, P, m, (mask > 1) ? TWISTED_EDWARDS_proj :
                                                           TWISTED_EDWARDS_ext);
    }

    ec_point_set (R, T, m, TWISTED_EDWARDS_ext);
    ec_point_clear (T, m);
  }
}

#ifdef __cplusplus
}
#endif


#endif /* EC_ARITH_EDWARDS_H_ */
