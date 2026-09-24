#ifndef CADO_SIEVE_ECM_EC_PARAMETERIZATION_HPP
#define CADO_SIEVE_ECM_EC_PARAMETERIZATION_HPP

#include <cstdint>
#include <cstdio>

#include "ec_arith_common.hpp"
#include "ec_arith_Weierstrass_new.hpp"
#include "macros.h"

/******************************************************************************/
/*********************** Brent--Suyama parameterization ***********************/
/******************************************************************************/

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
inline int ec_parameterization_Brent_Suyama_is_valid(unsigned long const sigma)
{
    return !(sigma == 0 || sigma == 1 || sigma == 3 || sigma == 5);
}

inline unsigned long
ec_parameterization_Brent_Suyama_valid_parameter_from_sequence(unsigned long const sequence_value)
{
    switch (sequence_value) {
        case 0: return 11;
        case 1: return 2;
        case 3: return 4;
    }
    if (sequence_value <= 8)
        return sequence_value + 2;
    return sequence_value + 3;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0.x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
template<typename layer>
int ec_parameterization_Brent_Suyama(typename layer::Residue & b,
        ec_point<layer> & P0, unsigned long const sigma,
        typename layer::Modulus const & m)
{
    ASSERT_ALWAYS(ec_parameterization_Brent_Suyama_is_valid(sigma));

    typename layer::Residue s(m), u(m), v(m), t1(m), t2(m);

    m.set(s, uint64_t(sigma));
    m.add(v, s, s);
    m.add(v, v, v);             /* v = 4*sigma */
    m.sqr(u, s);
    m.set(t1, uint64_t(5));
    m.sub(u, u, t1);            /* u = sigma^2 - 5 */
    m.sqr(t1, u);
    m.mul(P0.x, t1, u);         /* x0 = u^3 */
    m.sqr(t1, v);
    m.mul(P0.z, t1, v);         /* z0 = v^3 */

    m.mul(t1, P0.x, v);
    m.pow2(t2, uint64_t(4));
    m.mul(t1, t1, t2);          /* t1 = 16*x0*v = 16*u^3*v */
    m.add(b, u, u);
    m.add(b, b, u);
    m.add(b, b, v);             /* b = 3*u+v (so far) */
    m.sub(t2, v, u);
    m.mul(b, b, t2);            /* b = (v-u)*(3*u+v) (so far) */
    m.sqr(t2, t2);
    m.mul(b, b, t2);            /* b = (v-u)^3*(3*u+v) (so far) */
    int const ret = m.inv(t2, t1); /* t2 = 1/t1 */
    if (ret == 0) /* non-trivial gcd */
        m.set(P0.x, t1);
    else
        m.mul(b, b, t2);        /* b = (v-u)^3*(3*u+v)/(16*u^3*v) */

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
inline int ec_parameterization_Montgomery12_is_valid(unsigned long const k)
{
    return k > 1;
}

inline unsigned long
ec_parameterization_Montgomery12_valid_parameter_from_sequence(unsigned long const sequence_value)
{
    return sequence_value + 2;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0.x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
template<typename layer>
int ec_parameterization_Montgomery12(typename layer::Residue & b,
        ec_point<layer> & P0, unsigned long const k,
        typename layer::Modulus const & m)
{
    ASSERT_ALWAYS(ec_parameterization_Montgomery12_is_valid(k));

    using Residue = typename layer::Residue;
    Residue u(m), v(m), a(m), uu(m), vv(m), d1(m), d2(m), t1(m), t2(m), one(m);
    Residue x(m), y(m);

    m.set1(one);

    /* Compute u and v */
    m.add(y, one, one);
    m.neg(x, y);                /* x = -2 */
    m.add(y, y, y);             /* y = 4 */
    m.add(a, y, y);
    m.add(a, a, y);
    m.neg(a, a);                /* a = -12 */

    ECWeierstrass<layer> const E(m, a);
    typename ECWeierstrass<layer>::AffinePoint T(E, x, y);
    T.smul(T, k);
    if (T.is0()) { /* non-trivial gcd */
        m.set(P0.x, T.getX());
        return 0;
    }

    m.set(u, T.getX());
    m.set(v, T.getY());
    m.sqr(uu, u);               /* uu = u^2 */
    m.sqr(vv, v);               /* vv = v^2 */

    m.add(t1, uu, uu);
    m.add(t1, t1, t1);
    m.sub(d2, vv, t1);          /* d2 = v^2-4*u^2 */
    m.add(t2, t1, t1);
    m.add(t2, t2, t1);
    m.add(d1, vv, t2);          /* d1 = v^2+12*u^2 */

    m.set(t1, uint64_t(12));    /* t1 = 12 */
    m.add(P0.x, uu, t1);        /* x0 = (u^2+12) [so far] */
    m.sqr(P0.x, P0.x);          /* x0 = (u^2+12)^2 [so far] */
    m.mul(P0.x, P0.x, uu);      /* x0 = u^2*(u^2+12)^2 [so far] */
    m.add(P0.x, P0.x, P0.x);    /* x0 = 2*u^2*(u^2+12)^2 */

    m.mul(P0.z, d1, d2);        /* z0 = d1*d2 [so far] */
    m.add(P0.z, P0.z, P0.z);    /* z0 = 2*d1*d2 */

    m.sqr(t1, d2);
    m.mul(t1, t1, d2);
    m.mul(t1, t1, d1);          /* t1 = d1*d2^3 */
    int const ret = m.inv(b, t1); /* 1/t1 = 1/(d1*d2^3) (stored in b) */
    if (ret == 0) { /* non-trivial gcd */
        m.set(P0.x, t1);
    } else {
        m.sqr(t1, uu);
        m.mul(t1, t1, uu);
        m.mul(b, b, t1);        /* b = u^6/(d1*d2^3) [so far] */
        m.mul(b, b, vv);        /* b = v^2*u^6/(d1*d2^3)  [so far] */
        m.pow2(t1, uint64_t(10)); /* t1 = 2^10 = 1024 */
        m.mul(b, b, t1);        /* b = 1024*v^2*u^6/(d1*d2^3) */
    }
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
inline int ec_parameterization_Montgomery16_is_valid(unsigned long const k)
{
    return k == 1;
}

inline unsigned long
ec_parameterization_Montgomery16_valid_parameter_from_sequence(unsigned long const k MAYBE_UNUSED)
{
    /* there's only one curve. (well, it's more complicated, it seems)
     * Is it best to always return the same curve, or to error out
     * when k>0 ? */
    return 1;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0.x.
 */
template<typename layer>
int ec_parameterization_Montgomery16(typename layer::Residue & b,
        ec_point<layer> & P0, unsigned long const k MAYBE_UNUSED,
        typename layer::Modulus const & m)
{
    ASSERT_ALWAYS(ec_parameterization_Montgomery16_is_valid(k));

    /* Make curve corresponding to (a,b,c) = (8, 15, 17) in Montgomery's
     * thesis, equation (6.2.2).
     *    A = 54721/14400
     *    b = (A+2)/4 = (289/240)^2
     *    x0 = 8
     *    z0 = 15
     * [ see table 6.2.1 ]
     * We only need b and the coordinates of the starting point (x0::z0)
     * This curve is cheap to initialise: four div2, one div3, one div5,
     * three add, one mul.
     */
    typename layer::Residue t(m), one(m);

    m.set1(one);
    m.set(P0.x, uint64_t(8));   /* x0 = 8 */
    m.set(P0.z, uint64_t(15));  /* z0 = 15 */

    m.div3(t, one);
    m.div2(t, t);
    m.div2(t, t);
    m.div2(t, t);
    m.div2(t, t);               /* t = 1/48 */
    m.add(t, t, one);           /* t = 49/48 */
    m.div5(t, t);               /* t = 49/240 */
    m.add(t, t, one);           /* t = 289/240 */
    m.sqr(b, t);                /* b = 83521/57600 */

#ifdef WANT_ASSERT_EXPENSIVE
    typename layer::Residue t2(m);
    m.set(t, uint64_t(57600));
    m.mul(t2, b, t);
    m.set(t, uint64_t(83521));
    ASSERT(m.equal(t, t2));
#endif

    return 1; /* we assume gcd (m, 2*3*5) = 1 */
}

/******************************************************************************/
/******************* Z/6Z-rational-torsion parameterization *******************/
/******************************************************************************/

inline int ec_parameterization_Z6_is_valid(unsigned long const k)
{
    return k != 0;
}

inline unsigned long
ec_parameterization_Z6_valid_parameter_from_sequence(unsigned long const k)
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
 * If modular inverse failed, return non-invertible value in P0.x.
 *
 * We only need the starting point (xE0:yE0:zE0:tE0) (for Edwards extented) or
 * (xE0:yE0:zE0) (for Edwards projective) or (xM0::zM0) (for Montgomery) and
 * the curve coefficient b = (A+2)/4 (if b is not nullptr).
 */
template<typename layer>
int ec_parameterization_Z6(typename layer::Residue * b,
        ec_point<layer> & P0, unsigned long const k,
        ec_point_coord_type_t const coord,
        typename layer::Modulus const & m)
{
    ASSERT_ALWAYS(ec_parameterization_Z6_is_valid(k));
    ASSERT_ALWAYS(coord == MONTGOMERY_xz || coord == TWISTED_EDWARDS_ext
                                         || coord == TWISTED_EDWARDS_proj);

    using Residue = typename layer::Residue;
    Residue U(m), V(m), W(m);
    Residue u0(m), u1(m), u2(m), u3(m), u4(m), u5(m), u6(m), u7(m);
    Residue Tx(m), Ty(m), Tz(m), Tt(m);
    int ret = 0;

    {
        /* Inversion free scalar multiplication: T <- [k](15:378:1) */
        m.set(u0, uint64_t(9747));
        m.neg(u0, u0);
        m.set(Tx, uint64_t(15));
        m.set(Ty, uint64_t(378));
        m.set1(Tz);
        ECWeierstrass<layer> const E(m, u0);
        typename ECWeierstrass<layer>::ProjectivePoint T(E, Tx, Ty, Tz);
        T.smul(T, k);
        m.set(Tx, T.getX());
        m.set(Ty, T.getY());
        m.set(Tz, T.getZ());
    }

    m.set(u0, uint64_t(144));
    m.add(U, Tx, Tz);
    m.add(U, U, Tz);
    m.add(U, U, Tz);
    m.mul(U, U, u0);            /* U = 144*(T.x + 3*T.z) */
    m.set(V, Ty);               /* V = T.y */
    m.set(W, uint64_t(2985984));
    m.mul(W, W, Tz);            /* W = 12^6 * T.z*/

    m.set(u0, uint64_t(96));
    m.mul(u0, u0, U);           /* u0 = 96*U */
    m.sub(u1, W, u0);           /* u1 = (W - u0) */
    m.sqr(u2, u1);              /* u2 = u1^2 */
    m.sqr(u3, u0);              /* u3 = u0^2 */
    m.set(u4, uint64_t(5));
    m.mul(u4, u4, u3);
    m.sub(u4, u2, u4);          /* u4 = u2 - 5*u3 */
    m.sqr(u5, u4);
    m.mul(u5, u5, u4);          /* u5 = u4^3 */
    m.add(u7, u1, u1);
    m.add(u7, u7, u7);          /* u7 = 4*u1 */

    if (coord == MONTGOMERY_xz) {
        m.mul(P0.x, u5, u3);
        m.mul(P0.x, P0.x, u0);  /* xM0 = u5*u3*u0 */
        m.mul(u6, u7, u3);
        m.sqr(P0.z, u6);
        m.mul(P0.z, P0.z, u6);  /* zM0 = (u7*u3)^3 */
        m.mul(u2, u3, u0);      /* u2 = u0^3 */
    } else {
        /* At this point, T is not needed anymore.  */
        /* We use its coordinates for temporary variables.  */
        m.sub(Tx, u1, u0);      /* T.x = u1 - u0 */
        m.set(u6, uint64_t(5));
        m.mul(u6, u6, u0);      /* u6 = 5*u0 */
        m.add(u6, u1, u6);      /* u6 = u1 + 5*u0 */
        m.mul(Tx, Tx, u6);      /* T.x = (u1 - u0)(u1 + 5*u0) */
        m.set(u6, uint64_t(5));
        m.mul(u6, u6, u3);      /* u6 = 5*u3 */
        m.sub(Ty, u2, u6);      /* T.y = u2 - 5*u3 */
        m.add(u6, u2, u6);      /* u6 = u2 + 5*u6 */
        m.mul(Tx, Tx, u6);      /* T.x = (u1-u0)(u1+5*u0)(u2+5*u6) */
        m.sqr(u6, Ty);
        m.mul(Ty, Ty, u6);      /* T.y = u6^3 */
        m.mul(u6, u7, u0);
        m.sqr(Tt, u6);
        m.mul(Tt, Tt, u6);      /* T.t = (u7*u0)^3 */
        m.sqr(u6, U);
        m.mul(u6, u6, Tx);      /* u6 = (T.x)*U^2 */
        m.mul(u2, u3, u0);      /* u2 = u0^3 */
        m.mul(u1, u1, u2);
        m.mul(u1, u1, V);
        m.mul(u1, u1, W);
        m.add(u1, u1, u1);      /* u1 = 2*u1*u2*V*W */

        m.add(P0.z, Ty, Tt);    /* P0.z = T.y + T.t */
        m.set(Tz, P0.z);
        m.mul(P0.z, P0.z, u6);  /* P0.z = (T.y + T.t)*u6 */
        m.mul(P0.x, Tz, u1);    /* P0.x = (T.y + T.t)*u1 */
        m.sub(P0.y, Ty, Tt);
        if (coord == TWISTED_EDWARDS_ext)
            m.set(Tz, P0.y);
        m.mul(P0.y, P0.y, u6);  /* P0.y = (T.y - T.z)*u6 */
        if (coord == TWISTED_EDWARDS_ext)
            m.mul(P0.t, Tz, u1); /* P0.t = (T.y - T.z)*u1 */
    }

    if (b != nullptr) {
        // A = ((u7*u3 - u4*u0)^3 * (3*u4*u0 + u7*u3)) / (16*u3*u5*u7*u2)
        // u0, u2, u3, u4, u5, u7
        m.mul(Tt, u3, u5);
        m.mul(Tt, Tt, u7);
        m.mul(Tt, Tt, u2);
        m.add(Tt, Tt, Tt);
        m.add(Tt, Tt, Tt);
        m.add(Tt, Tt, Tt);
        m.add(Tt, Tt, Tt);      /* 16*u3*u5*u7*u2 */
        ret = m.inv(Tt, Tt);
        if (ret == 0) {
            m.set(P0.x, Tt);
        } else {
            m.mul(Ty, u7, u3);  /* T.y = u7*u3 */
            m.mul(Tz, u4, u0);  /* T.z = u4*u0 */
            m.add(Tx, Tz, Tz);
            m.add(Tx, Tx, Tz);  /* T.x = 3*T.z */
            m.add(Tx, Tx, Ty);  /* T.x = 3*T.z + T.y */
            m.sub(Ty, Ty, Tz);  /* T.y = T.y - T.z */
            m.sqr(*b, Ty);
            m.mul(*b, *b, Ty);
            m.mul(*b, *b, Tx);
            m.mul(*b, *b, Tt);  /* b = (T.y)^3*(T.x)*(T.t) */
        }
    }

    return ret;
}

/******************************************************************************/
/*************** from a Montgomery curve to a Weierstrass curve ***************/
/******************************************************************************/

/* Convert the Montgomery curve B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2 with a valid
 * point Pm into a affine Weierstrass curve y^2 = x^3 + a*x + b with a
 * valid affine point (xw, yw).
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in xw.
 *
 * The curve coefficient b of the short Weierstrass curve will not be
 * computed.
 */
template<typename layer>
int weierstrass_aff_from_montgomery(typename layer::Residue & a,
        typename layer::Residue & xw, typename layer::Residue & yw,
        typename layer::Residue const & A, ec_point<layer> const & Pm,
        typename layer::Modulus const & m)
{
    typename layer::Residue B(m), one(m), t(m), x(m);

    m.set1(one);

    if (!m.inv(t, Pm.z)) {
        fprintf(stderr, "%s: could not invert Z\n", __func__);
        m.set(xw, Pm.z);
        return 0;
    }

    m.mul(x, Pm.x, t); /* x = X/Z */
    m.add(B, x, A);
    m.mul(B, B, x);
    m.add(B, B, one);
    m.mul(B, B, x); /* B = x^3 + A*x^2 + x */

    /* Now (x,1) is on the curve B*y^2 = x^3 + A*x^2 + x. */
    if (!m.inv(yw, B)) {    /* y = 1/B */
        fprintf(stderr, "%s: could not invert B\n", __func__);
        m.set(xw, B);
        return 0;
    }

    m.div3(t, A);
    m.add(xw, x, t);
    m.mul(xw, xw, yw);      /* x = (X + A/3)/B */
    m.mul(a, t, A);
    m.sub(a, one, a);
    m.mul(a, a, yw);
    m.mul(a, a, yw);        /* a = (1 - (A^2)/3)/B^2 */

    return 1;
}

#endif /* CADO_SIEVE_ECM_EC_PARAMETERIZATION_HPP */
