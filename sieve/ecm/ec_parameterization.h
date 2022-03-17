#include <stdlib.h>
#ifndef EC_PARAMETERIZATION_H_
#define EC_PARAMETERIZATION_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "ec_arith_common.h"
#include "ec_arith_Weierstrass.h"
#include "macros.h"


/******************************************************************************/
/*********************** Brent--Suyama parameterization ***********************/
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/* Produces curve in Montgomery form.
 *
 * Rational parameterization (parameter is called sigma).
 *
 * Valid parameter: sigma in Q \ { 0, -1, 1, -3, 3, -5, 5, -5/3, 5/3 }
 *
 * Note: sigma and -sigma give isomorphic curves.
 *
 * Parameterization:
      u = sigma^2-5
      v = 4*sigma
      x0 = u^3
      y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)
      z0 = v^3
      A = (v-u)^3*(3*u+v)/(4*u^3*v) - 2
      B = u/z0
      b = (v-u)^3*(3*u+v)/(16*u^3*v)                          # [ b = (A+2)/4 ]
 * The point of order 3 is given by:
      x3 = u
      y3 = 2*(sigma^2-25)*(sigma^2-1)*sigma/u
      z3 = v
 * To check with Sage:
      QQsigma.<sigma> = QQ[]
      # copy paste all the above formula
      (b-(A+2)/4).is_zero() # b is correct ?
      (B*y0^2*z0 - (x0^3+A*x0^2*z0+x0*z0^2)).is_zero() # P0 is on the curve ?
      (B*y3^2*z3 - (x3^3+A*x3^2*z3+x3*z3^2)).is_zero() # P3 is on the curve ?
      (4*A*(x3/z3)^3+3*(x3/z3)^4+6*(x3/z3)^2-1).is_zero() # P3 of order 3 ?
 */

static inline int
ec_parameterization_Brent_Suyama_is_valid (const unsigned long sigma)
{
  if (sigma == 0 || sigma == 1 || sigma == 3 || sigma == 5)
    return 0; /* invalid values */
  else
    return 1;
}

static inline unsigned long
ec_parameterization_Brent_Suyama_valid_parameter_from_sequence(const unsigned long sequence_value)
{
    switch(sequence_value) {
        case 0: return 11;
        case 1: return 2;
        case 3: return 4;
    }
    if (sequence_value <= 8)
        return sequence_value + 2;
    return sequence_value + 3;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
#define ec_parameterization_Brent_Suyama MOD_APPEND_TYPE(ec_parameterization_Brent_Suyama)
static int
ec_parameterization_Brent_Suyama (residue_t b, ec_point_t P0,
                                  const unsigned long sigma, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Brent_Suyama_is_valid (sigma));

  residue_t s, u, v, t1, t2;
  int ret;

  mod_init_noset0 (s, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);

  mod_set_ul (s, sigma, m);

  mod_add (v, s, s, m);
  mod_add (v, v, v, m);             /* v = 4*sigma */
  mod_sqr (u, s, m);
  mod_set_ul (t1, 5, m);
  mod_sub (u, u, t1, m);            /* u = sigma^2 - 5 */
  mod_sqr (t1, u, m);
  mod_mul (P0->x, t1, u, m);        /* x0 = u^3 */
  mod_sqr (t1, v, m);
  mod_mul (P0->z, t1, v, m);        /* z0 = v^3 */
  mod_mul (t1, P0->x, v, m);
  mod_2pow_ul (t2, 4, m);
  mod_mul (t1, t1, t2, m);          /* t1 = 16*x0*v = 16*u^3*v */
  mod_add (b, u, u, m);
  mod_add (b, b, u, m);
  mod_add (b, b, v, m);             /* b = 3*u+v (so far) */
  mod_sub (t2, v, u, m);
  mod_mul (b, b, t2, m);            /* b = (v-u)*(3*u+v) (so far) */
  mod_sqr (t2, t2, m);
  mod_mul (b, b, t2, m);            /* b = (v-u)^3*(3*u+v) (so far) */

  ret = mod_inv (t2, t1, m);         /* t2 = 1/t1 */
  if (ret == 0) /* non-trivial gcd */
    mod_set (P0->x, t1, m);
  else
    mod_mul (b, b, t2, m);          /* b = (v-u)^3*(3*u+v)/(16*u^3*v) */

  mod_clear (s, m);
  mod_clear (u, m);
  mod_clear (v, m);
  mod_clear (t1, m);
  mod_clear (t2, m);

  return ret;
}

/******************************************************************************/
/********************* Montgomery Z/12Z parameterization **********************/
/******************************************************************************/

/* Produces curve in Montgomery form.
 *
 * Elliptic parameterization (parameter is called k).
 *    E: y^2 = x^3 - 12*x
 *      rank = 1
 *      P = (-2, 4) is a point of infinite order
 *      Torsion points: P2 = (0:0:1) of order 2
 *
 * Reference: Montgomery's thesis (Section 6.1)
 *
 * Valid parameter: k in Z \ {-1, 0, 1}
 *
 * Note: Using P2 does not produce a valid curve.
 *       A point and its opposite produce the same curve
 *       Two points that differ by P2 produce the same curve (but not the same
 *        starting point)
 *
 * Parameterization (where (u,v) := k*P):
      d1 = v^2+12*u^2
      d2 = v^2-4*u^2
      x0 = 2*(u*(u^2+12))^2
      y0 = d1*u*(u^2+12)
      z0 = 2*d1*d2
      A = (4*(32*v*u^3)^2-2*d1*d2^3)/(d1*d2^3)
      B = (16*u^2*(v^2+4*u^2))^2/(d1*d2^3)
      b = (32*v*u^3)^2/(d1*d2^3)
 * The point of order 12 is given by:
      x12 = 2*u*(2*u + v)^2
      y12 = v*(2*u + v)^2
      z12 = 2*u*(-2*u + v)^2
 * To check with Sage:
      QQuv.<u,v> = QQ[]
      I = QQuv.ideal (v^2 - (u^3-12*u))
      # copy paste all the above formula
      (b-(A+2)/4).is_zero() # b is correct ?
      # P0 is on the curve ?
      eq_P0 = (B*y0^2*z0 - (x0^3+A*x0^2*z0+x0*z0^2))
      eq_P0.numerator() in I and not eq_P0.denominator() in I
      # P12 is on the curve ?
      eq_P12 = (B*y12^2*z12 - (x12^3+A*x12^2*z12+x12*z12^2))
      eq_P12.numerator() in I and not eq_P12.denominator() in I
      # TODO check that 4*P12 is of order 3 and 3*P12 is of order 4
 */

static inline int
ec_parameterization_Montgomery12_is_valid (const unsigned long k)
{
  if (k <= 1)
    return 0; /* invalid values */
  else
    return 1;
}
static inline unsigned long
ec_parameterization_Montgomery12_valid_parameter_from_sequence (const unsigned long sequence_value)
{
    return sequence_value + 2;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
#define ec_parameterization_Montgomery12 MOD_APPEND_TYPE(ec_parameterization_Montgomery12)
static int
ec_parameterization_Montgomery12 (residue_t b, ec_point_t P0,
                                  const unsigned long k, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Montgomery12_is_valid (k));

  residue_t u, v, a, uu, vv, d1, d2, t1, t2, one;
  ec_point_t T;
  int ret = 0;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (a, m);
  mod_init_noset0 (uu, m);
  mod_init_noset0 (vv, m);
  mod_init_noset0 (d1, m);
  mod_init_noset0 (d2, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (one, m);

  ec_point_init (T, m);

  mod_set1 (one, m);

  /* Compute u and v */
  mod_add (T->y, one, one, m);
  mod_neg (T->x, T->y, m);              /* x = -2 */
  mod_add (T->y, T->y, T->y, m);        /* y = 4 */
  mod_add (a, T->y, T->y, m);
  mod_add (a, a, T->y, m);
  mod_neg (a, a, m);                    /* a = -12 */

  ret = weierstrass_aff_smul_ui (T, k, a, m);
  if (ret == 0) /* non-trivial gcd */
    mod_set (P0->x, T->x, m);
  else
  {
    mod_set (u, T->x, m);
    mod_set (v, T->y, m);

    mod_sqr (uu, u, m);                   /* uu = u^2 */
    mod_sqr (vv, v, m);                   /* vv = v^2 */

    mod_add (t1, uu, uu, m);
    mod_add (t1, t1, t1, m);
    mod_sub (d2, vv, t1, m);              /* d2 = v^2-4*u^2 */
    mod_add (t2, t1, t1, m);
    mod_add (t2, t2, t1, m);
    mod_add (d1, vv, t2, m);              /* d1 = v^2+12*u^2 */

    mod_set_ul (t1, 12, m);               /* t1 = 12 */
    mod_add (P0->x, uu, t1, m);           /* x0 = (u^2+12) [so far] */
    mod_sqr (P0->x, P0->x, m);            /* x0 = (u^2+12)^2 [so far] */
    mod_mul (P0->x, P0->x, uu, m);        /* x0 = u^2*(u^2+12)^2 [so far] */
    mod_add (P0->x, P0->x, P0->x, m);     /* x0 = 2*u^2*(u^2+12)^2 */

    mod_mul (P0->z, d1, d2, m);           /* z0 = d1*d2 [so far] */
    mod_add (P0->z, P0->z, P0->z, m);     /* z0 = 2*d1*d2 */

    mod_sqr (t1, d2, m);
    mod_mul (t1, t1, d2, m);
    mod_mul (t1, t1, d1, m);              /* t1 = d1*d2^3 */

    ret = mod_inv (b, t1, m);             /* 1/t1 = 1/(d1*d2^3) (stored in b) */
    if (ret == 0) /* non-trivial gcd */
      mod_set (P0->x, t1, m);
    else
    {
      mod_sqr (t1, uu, m);
      mod_mul (t1, t1, uu, m);
      mod_mul (b, b, t1, m);              /* b = u^6/(d1*d2^3) [so far] */
      mod_mul (b, b, vv, m);              /* b = v^2*u^6/(d1*d2^3)  [so far] */
      mod_2pow_ul (t1, 10, m);            /* t1 = 2^10 = 1024 */
      mod_mul (b, b, t1, m);              /* b = 1024*v^2*u^6/(d1*d2^3) */
    }
  }

  mod_clear (u, m);
  mod_clear (v, m);
  mod_clear (a, m);
  mod_clear (uu, m);
  mod_clear (vv, m);
  mod_clear (d1, m);
  mod_clear (d2, m);
  mod_clear (t1, m);
  mod_clear (t2, m);
  mod_clear (one, m);

  ec_point_clear (T, m);

  return ret;
}

/******************************************************************************/
/********************* Montgomery Z/16Z parameterization **********************/
/******************************************************************************/

/* Produces curve in Montgomery form.
 *
 * Finite number of curves.
 *
 * Reference: Montgomery's thesis (Section 6.2)
 *
 * Valid parameter: k in { 1 }
 */


static inline int
ec_parameterization_Montgomery16_is_valid (const unsigned long k)
{
  if (k != 1)
    return 0; /* invalid values */
  else
    return 1;
}
static inline unsigned long
ec_parameterization_Montgomery16_valid_parameter_from_sequence (const unsigned long k MAYBE_UNUSED)
{
    /* there's only one curve. (well, it's more complicated, it seems)
     * Is it best to always return the same curve, or to error out
     * when k>0 ? */
    return 1;
}


#define ec_parameterization_Montgomery16 MOD_APPEND_TYPE(ec_parameterization_Montgomery16)
/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 */
static int
ec_parameterization_Montgomery16 (residue_t b, ec_point_t P0,
                                  const unsigned long k, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Montgomery16_is_valid (k));

  /* Make curve corresponding to (a,b,c) = (8, 15, 17) in Montgomery's thesis,
   * equation (6.2.2).
   *    A = 54721/14400
   *    b = (A+2)/4 = (289/240)^2
   *    x0 = 8
   *    z0 = 15
   * [ see table 6.2.1 ]
   * We only need b and the coordinates of the starting point (x0::z0)
   * This curve is cheap to initialise: four div2, one div3, one div5,
   * three add, one mul.
   */
  residue_t t, one;

  mod_init (t, m);
  mod_init (one, m);

  mod_set1 (one, m);

  mod_set_ul (P0->x, 8UL, m);     /* x0 = 8 */
  mod_set_ul (P0->z, 15UL, m);    /* z0 = 15 */

  mod_div3 (t, one, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);              /* t = 1/48 */
  mod_add (t, t, one, m);         /* t = 49/48 */
  mod_div5 (t, t, m);             /* t = 49/240 */
  mod_add (t, t, one, m);         /* t = 289/240 */
  mod_sqr (b, t, m);              /* b = 83521/57600 */

#ifdef WANT_ASSERT_EXPENSIVE
  residue_t t2;
  mod_init (t2, m);
  mod_set_ul (t, 57600UL, m);
  mod_mul (t2, b, t, m);
  mod_set_ul (t, 83521UL, m);
  ASSERT (mod_equal (t, t2, m));
  mod_clear (t2, m);
#endif

  mod_clear (t, m);
  mod_clear (one, m);

  return 1; /* we assume gcd (m, 2*3*5) = 1 */
}

/******************************************************************************/
/******************* Z/6Z-rational-torsion parameterization *******************/
/******************************************************************************/


static inline int
ec_parameterization_Z6_is_valid (const unsigned long k)
{
  if (k == 0)
    return 0; /* invalid values */
  else
    return 1;
}

static inline unsigned long
ec_parameterization_Z6_valid_parameter_from_sequence (const unsigned long k)
{
    return k + 1;
}

/* Produces curve
 *     - in "a=-1" Twisted Edwards form (in extended or projective coordinates)
 *     - or in Montgomery from
 * 
 * Elliptic parameterization (parameter is called k).
 *    E: y^2 = x^3 - 9747*x + 285714
 *       rank = 1
 *       P = (15, 378) is a point of infinite order
 *       Torsion points: P2 = (33:0:1) of order 2
 *                       Q2 = (78:0:1) of order 2
 *                       P2+Q2 = (-111:0:1) of order 2
 *
 * Reference: based on Theorem 5.4 of the article "Starfish on strike".
 *
 * Valid parameter: k in Z \ { 0 }
 *
 * Note: P2, Q2, P2+Q2 does not produce a valid curve.
 *       P+Q2, P+P2+Q2 does not produce a valid curve.
 *       Two points that differ by P2 produce the same curve.
 *       A point and its opposite produce the same curve (with opposite starting
 *        points)
 *       It generates curves isomorphic to Brent--Suyama curves with
 *       sigma = -(p-213)/(p+3)
 *       The value k=1 produces the same curve as the Brent-Suyama
 *       parameterization with sigma=11 (up to isomorphism, i.e. the values of
 *       B differ by a square factor).
 *

 * Parameterization (where (U, V, W) := [k]*P)

 * The following code can be pasted into Sage
 * -------------------------------------------------------------------------

 # Curve from Strafish on strike (theorem 5.4)
 E1 = EllipticCurve([0,-1/2304,0,-5/221184,1/28311552])

 # A curve is short Weierstrass form, isomorphic to E1 with integer coefficients
 E2 = EllipticCurve([-9747, 285714])

 # QQpqr.<p,q,r> = QQ[]
 # I = QQpqr.ideal(E2.defining_polynomial()(z=r,x=p,y=q))

 QQxyz.<x,y,z> = QQ[]
 I = QQxyz.ideal(E2.defining_polynomial())

 phi1 = E1.isomorphism_to(E2)
 phi2 = E2.isomorphism_to(E1)

 # Apply the Weierstrass transformation (u,r,s,t) in projective coordinate
 # (x,y,z)  |-->  (U,V,W) = (u^2x + r*z , u^3y + su^2x + t*z, z).

 U = (phi1).u^2*x+(phi1).r*z
 V = (phi1).u^3*y+(phi1).s*(phi1).u^2*x+(phi1).t*z
 W = z

 l = lcm([c.denominator() for c in U.coefficients()\
 + V.coefficients() + W.coefficients()])

 U *= l     # U = 144*(x + 3*z)
 V *= l     # V = y
 W *= l     # W = 2985984*z

 # -----------------------------------------------------------------------------

 u0 = 96*U                       # u0 = 96*U
 u1 = (W - u0)                   # u1 = (W - sigma_d)
 u2 = u1^2                       # u2 = u1^2
 u3 = u0^2                       # u3 = u0^2
 u4 = u2 - 5*u3                  # u4 = u2 - 5*u3
 u5 = u4^3                       # u5 = u4^3
 u6 = 4*u1                       # u6 = 4*u1

 # Compute M0 = (xM0::zE0) and A: Montgomery curve parameters
 # following 20 years of ECM by Zimmermann and Dodson

 xM0 = u5*u3*u0
 zM0 = (u6*u3)^3

 #A = (((beta-alpha)^3)*(3*alpha+beta) / (4*alpha3*beta)) - 2
 A = ((u6*u3-u4*u0)^3 * (3*u4*u0 + u6*u3)) / (4*u3*u5*u6*u2) - 2

 # Compute E0 = (xE0:yE0:zE0:tE0) base point on twisted Edwards curve
 # following starfish on strike

 Tx = (u1 - u0)*(u1 + 5*u0)*(u2 + 5*u3)
 Ty = (u2 - 5*u3)^3
 Tt = (u6*u0)^3

 u1 = 2*u1*u2*V*W
 u2 = u0^3
 u7 = Tx*U^2

 zE0 = (Ty + Tt)*u7
 xE0 = (Ty + Tt)*u1
 yE0 = (Ty - Tt)*u7
 tE0 = (Ty - Tt)*u1

 # -----------------------------------------------------------------------------

 # Check if PE0 = (xE0:yE0:zE0:tE0) is on the extended twisted
 # Edwards curve: (-1)*xE0^2 + yE0^2 - zE0^2 - d*tE0^2
 a = -1
 # Compute r = v/u^2 following Starfish on strike (Theorem 5.4)
 # with v = V/W and u = U/W
 r = (V*W)/U^2
 # Compute d following Starfish on strike (Theorem 5.4)
 alpha = u4 / u3
 beta = u6 / u0
 d = ((beta+alpha)^3*(beta-3*alpha)) / (r*(beta-alpha))^2

 eq_Ed = a*xE0^2 + yE0^2 - zE0^2 - d*tE0^2
 print "P0 is on curve:       ",
 print eq_Ed.numerator() in I\
 and not(eq_Ed.denominator() in I)\
 and (xE0*yE0 - zE0*tE0 in I)


 # Applying the morphism from the Montgomery curve By2 = x^3 + Ax^2 + x
 # to its equivalent Edwards form ax^2 + y^2 = 1 + dx^2y^2
 # sets a = (A+2)/B.
 # We define B so that the corresponding Edwards curve has a = -1.
 B = - (A+2)

 # Check if PM0 = (xM0:yM0:zM0) is on the Montgomery curve equivalent to
 # the Edwards curve eq_Ed
 alpha3 = alpha^3
 beta3 = beta^3
 sigma = u1 / u0
 _b = alpha/beta3
 y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)*(u3*u6)^3

 # Ad-hoc code for checking that b/B is a square in I
 e2 = E2.defining_polynomial()
 e3 = (-e2 + y^2*z)*z
 f = (_b/(B*e3)).factor()
 print "b/B is a square in I: ",\
 not(any([m[1] % 2 for m in list(f)])) and QQ(f.unit()).is_square()
 yM02 = (_b/B)*y0^2
 eq_M0 = B*yM02*zM0-xM0*(xM0^2+A*xM0*zM0+zM0^2)
 print "M0 in on curve:       ",
 print eq_M0.numerator() in I and not(eq_M0.denominator() in I)

 # Define P3 = (xE3, yE3) a point of order 3 on eq_Ed
 # Check that P3 in on the curve

 xE3 = -r/((sigma-1)*(sigma+5))
 yE3 = (alpha - beta)/(alpha + beta)
 eq_P3 = a*xE3^2+yE3^2 - (1 + d*xE3^2*yE3^2)
 print "P3 is on curve:       ",
 print eq_P3.numerator() in I\
 and not(eq_P3.denominator() in I)

 * -------------------------------------------------------------------------

 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * We only need the starting point (xE0:yE0:zE0:tE0) (for Edwards extented) or
 * (xE0:yE0:zE0) (for Edwards projective) or (xM0::zM0) (for Montgomery) and
 * the curve coefficient b = (A+2)/4 (if b is not NULL).
 */
#define ec_parameterization_Z6 MOD_APPEND_TYPE(ec_parameterization_Z6)
static int
ec_parameterization_Z6 (residue_t b, ec_point_t P0, const unsigned long k,
                        ec_point_coord_type_t coord, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Z6_is_valid (k));
  ASSERT_ALWAYS (coord == MONTGOMERY_xz || coord == TWISTED_EDWARDS_ext
                                        || coord == TWISTED_EDWARDS_proj);

  residue_t U, V, W;
  residue_t u0, u1, u2, u3, u4, u5, u6, u7;
  
  ec_point_t T;
  int ret = 0;

  ec_point_init (T, m);

  mod_init_noset0 (U, m);
  mod_init_noset0 (V, m);
  mod_init_noset0 (W, m);
  mod_init_noset0 (u0, m);
  mod_init_noset0 (u1, m);
  mod_init_noset0 (u2, m);
  mod_init_noset0 (u3, m);
  mod_init_noset0 (u4, m);
  mod_init_noset0 (u5, m);
  mod_init_noset0 (u6, m);
  mod_init_noset0 (u7, m);

  mod_set_ul (u0, 9747UL, m);
  mod_neg (u0, u0, m);
  mod_set_ul (T->x, 15UL, m);
  mod_set_ul (T->y, 378UL, m);
  mod_set1 (T->z, m);

  /* Inversion free scalar multiplication */
  /* T <- [k]T */
  weierstrass_proj_smul_ui (T, k, u0, m);

  mod_set_ul (u0, 144UL, m);
  mod_add (U, T->x, T->z, m);
  mod_add (U, U, T->z, m);
  mod_add (U, U, T->z, m);
  mod_mul (U, U, u0, m);                /* U = 144*(T->x + 3*T->z) */  
  mod_set (V, T->y, m);             	/* V = T->y */		   
  mod_set_ul (W, 2985984UL, m);		                               
  mod_mul (W, W, T->z, m);          	/* W = 12^6 * T->z*/	   
					                               
  mod_set_ul (u0, 96UL, m);                                            
  mod_mul (u0, u0, U, m);           	/* u0 = 96*U */		   
  mod_sub (u1, W, u0, m);           	/* u1 = (W - u0) */		   
  mod_sqr (u2, u1, m);              	/* u2 = u1^2 */		   
  mod_sqr (u3, u0, m);              	/* u3 = u0^2 */		   
  mod_set_ul (u4, 5UL, m);		                               
  mod_mul (u4, u4, u3, m);  		                               
  mod_sub (u4, u2, u4, m);          	/* u4 = u2 - 5*u3 */	   
  mod_sqr (u5, u4, m);     		                               
  mod_mul (u5, u5, u4, m);          	/* u5 = u4^3 */		   

  mod_add (u7, u1, u1, m);		                               
  mod_add (u7, u7, u7, m);          	/* u7 = 4*u1 */		   
  
  if (coord == MONTGOMERY_xz)
  {
    mod_mul (P0->x, u5, u3, m);
    mod_mul (P0->x, P0->x, u0, m);      /* xM0 = u5*u3*u0 */
    mod_mul (u6, u7, u3, m);
    mod_sqr (P0->z, u6, m);
    mod_mul (P0->z, P0->z, u6, m);      /* zM0 = (u7*u3)^3 */
    mod_mul (u2, u3, u0, m);            /* u2 = u0^3 */
  }
  else
  {
    /* At this point, T is not needed anymore.  */
    /* We use its coordinates for temporary variables.  */
  
    mod_sub (T->x, u1, u0, m);          /* T->x = u1 - u0 */
    mod_set_ul (u6, 5UL, m);                 
    mod_mul (u6, u6, u0, m);            /* u6 = 5*u0 */
    mod_add (u6, u1, u6, m);            /* u6 = u1 + 5*u0 */
    mod_mul (T->x, T->x, u6, m);        /* T->x = (u1 - u0)(u1 + 5*u0) */
    mod_set_ul (u6, 5UL, m);                 
    mod_mul (u6, u6, u3, m);             /* u6 = 5*u3 */
    mod_sub (T->y, u2, u6, m);           /* T->y = u2 - 5*u3 */
    mod_add (u6, u2, u6, m);             /* u6 = u2 + 5*u6 */
    mod_mul (T->x, T->x, u6, m);         /* T->x = (u1-u0)(u1+5*u0)(u2+5*u6) */
    mod_sqr (u6, T->y, m);
    mod_mul (T->y, T->y, u6, m);         /* T->y = u6^3 */

    mod_mul (u6, u7, u0, m);
    mod_sqr (T->t, u6, m);
    mod_mul (T->t, T->t, u6, m);         /* T->t = (u7*u0)^3 */
	
    mod_sqr (u6, U, m);
    mod_mul (u6, u6, T->x, m);           /* u6 = (T->x)*U^2 */
    
    mod_mul (u2, u3, u0, m);             /* u2 = u0^3 */
    mod_mul (u1, u1, u2, m);
    mod_mul (u1, u1, V, m);
    mod_mul (u1, u1, W, m);
    mod_add (u1, u1, u1, m);             /* u1 = 2*u1*u2*V*W */
    
    mod_add (P0->z, T->y, T->t, m);      /* P0->z = T->y + T->t */
    mod_set (T->z, P0->z, m);      
    mod_mul (P0->z, P0->z, u6, m);       /* P0->z = (T->y + T->t)*u6 */
    mod_mul (P0->x, T->z, u1, m);        /* P0->x = (T->y + T->t)*u1 */
    mod_sub (P0->y, T->y, T->t, m); 
    
    if (coord == TWISTED_EDWARDS_ext)
      mod_set (T->z, P0->y, m);

    mod_mul (P0->y, P0->y, u6, m);       /* P0->y = (T->y - T->z)*u6 */
    
    if (coord == TWISTED_EDWARDS_ext)
	mod_mul (P0->t, T->z, u1, m);    /* P0->t = (T->y - T->z)*u1 */
  }
  
  if (b != NULL)
  {

    // A = ((u7*u3 - u4*u0)^3 * (3*u4*u0 + u7*u3)) / (16*u3*u5*u7*u2)

    // u0, u2, u3, u4, u5, u7
    
    mod_mul (T->t, u3, u5, m);
    mod_mul (T->t, T->t, u7, m);
    mod_mul (T->t, T->t, u2, m);
    mod_add (T->t, T->t, T->t, m);
    mod_add (T->t, T->t, T->t, m);
    mod_add (T->t, T->t, T->t, m);
    mod_add (T->t, T->t, T->t, m);       /* 16*u3*u5*u7*u2 */
    ret = mod_inv (T->t, T->t, m);
    if (ret == 0)
      mod_set (P0->x, T->t, m);
    else
      {
	mod_mul (T->y, u7, u3, m);       /* T->y = u7*u3 */
	mod_mul (T->z, u4, u0, m);       /* T->z = u4*u0 */
	mod_add (T->x, T->z, T->z, m);
	mod_add (T->x, T->x, T->z, m);   /* T->x = 3*T->z */
	mod_add (T->x, T->x, T->y, m);   /* T->x = 3*T->z + T->y */
	mod_sub (T->y, T->y, T->z, m);   /* T->y = T->y - T->z */
	mod_sqr (b, T->y, m);      
	mod_mul (b, b, T->y, m);      
	mod_mul (b, b, T->x, m);
	mod_mul (b, b, T->t, m);         /* b = (T->y)^3*(T->x)*(T->t) */
      }
  }
  
  mod_clear (U, m);
  mod_clear (V, m);
  mod_clear (W, m);
  mod_clear (u0, m);
  mod_clear (u1, m);
  mod_clear (u2, m);
  mod_clear (u3, m);
  mod_clear (u4, m);
  mod_clear (u5, m);
  mod_clear (u6, m);
  mod_clear (u7, m);

  ec_point_clear (T, m);

  return ret;
}

#ifdef __cplusplus
}
#endif

#endif /* EC_PARAMETERIZATION_H_ */
