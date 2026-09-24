#ifndef CADO_SIEVE_ECM_EC_ARITH_EDWARDS_HPP
#define CADO_SIEVE_ECM_EC_ARITH_EDWARDS_HPP

#include <cstdio>

#include "ec_arith_common.hpp"
#include "macros.h"
#ifdef ECM_COUNT_OPS
#include "ec_arith_cost.h"
#endif

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
/* not thread-safe */
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

/* Compute d = -(A-2)/(A+2). A and d can be the same variable. */
template<typename layer>
void edwards_d_from_montgomery_A(typename layer::Residue & d,
        typename layer::Residue const & A,
        typename layer::Modulus const & m)
{
    typename layer::Residue t(m), two(m);

    m.set1(two);
    m.add(two, two, two);
    m.add(t, A, two);           /* t = A+2 */
    m.inv(d, t);
    m.sub(t, t, two);
    m.sub(t, t, two);           /* t = A-2 */
    m.mul(d, d, t);
    m.neg(d, d);
}

template<typename layer>
void edwards_ext_curve_fprintf(FILE * out, char const * prefix,
        typename layer::Residue const & d, ec_point<layer> const * P,
        typename layer::Modulus const & m)
{
    char const * pre = (prefix == nullptr) ? "" : prefix;

    fmt::print(out, "{}Twisted Edwards curve: -X^2 + Y^2 = Z^2 + d*T^2\n"
                    "{}XY = ZT (extended coordinates)\n{}d = {:#x}\n",
                    pre, pre, pre, ec_residue_get_mpz<layer>(d, m));

    if (P) {
        fmt::print(out, "{}with point (X:Y:Z:T) = ", pre);
        ec_point_fprintf(out, *P, TWISTED_EDWARDS_ext, m);
        fputc('\n', out);
    }
}

/* Set P to zero (the neutral point):
 *    (0:1:1:0)     if coord is TWISTED_EDWARDS_ext
 *    (0:1:1)       if coord is TWISTED_EDWARDS_proj
 */
template<typename layer>
void edwards_point_set_zero(ec_point<layer> & P,
        typename layer::Modulus const & m,
        ec_point_coord_type_t const coord)
{
    ASSERT_EXPENSIVE(coord == TWISTED_EDWARDS_ext ||
                     coord == TWISTED_EDWARDS_proj);
    m.set0(P.x);
    m.set1(P.y);
    m.set1(P.z);
    if (coord == TWISTED_EDWARDS_ext)
        m.set0(P.t);
}

template<typename layer>
void edwards_neg(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m)
{
    m.neg(Q.x, P.x);
    m.set(Q.y, P.y);
    m.set(Q.z, P.z);
    m.neg(Q.t, P.t);
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
 *    4M + 10add + 2*2        if output_type == MONTGOMERY_xz
 * Notations in the comments come from:
 *    https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-add-2008-hwcd-4
 * Source: Hisil–Wong–Carter–Dawson, 2008, section 3.2 of
 *    http://eprint.iacr.org/2008/522
 */
template<typename layer>
void edwards_addsub(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, int const sub,
        typename layer::Modulus const & m,
        ec_point_coord_type_t const output_type)
{
    ASSERT_EXPENSIVE(output_type == TWISTED_EDWARDS_ext ||
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

    typename layer::Residue u0(m), u1(m), u2(m);

    if (sub) {
        m.add(u1, Q.y, Q.x);    /* u1 <-      (Y2+X2) */
        m.sub(u0, Q.y, Q.x);    /* u0 <-      (Y2-X2) */
    } else {
        m.add(u0, Q.y, Q.x);    /* u0 <-      (Y2+X2) */
        m.sub(u1, Q.y, Q.x);    /* u1 <-      (Y2-X2) */
    }

    m.sub(u2, P.y, P.x);        /* u2 <-      (Y1-X1) */
    m.mul(u0, u0, u2);          /* u0 <- A := (Y1-X1)*(Y2+/-X2) */
    m.add(u2, P.y, P.x);        /* u2 <-      (X1+Y1) */
    m.mul(u1, u1, u2);          /* u1 <- B := (Y1+X1)*(Y2-/+X2) */

    m.sub(R.x, u1, u0);         /* Rx <- F := B-A */
    m.add(R.y, u1, u0);         /* Ry <- G := B+A */

    m.add(u1, P.z, P.z);        /* u1 <-      2*Z1 */
    m.mul(u1, u1, Q.t);         /* u1 <- C := 2*Z1*T2 */
    m.add(u2, Q.z, Q.z);        /* u2 <-      2*Z2 */
    m.mul(u2, u2, P.t);         /* u2 <- D := 2*T1*Z2 */

    if (sub) {
        m.add(u0, u2, u1);      /* u0 <- H := D+C */
        m.sub(u1, u2, u1);      /* u1 <- E := D-C */
    } else {
        m.sub(u0, u2, u1);      /* u0 <- H := D-C */
        m.add(u1, u2, u1);      /* u1 <- E := D+C */
    }

    if (output_type == TWISTED_EDWARDS_ext || output_type == TWISTED_EDWARDS_proj) {
        m.mul(R.z, R.x, R.y);   /* Rz <- Z3 := F*G */
        m.mul(R.x, R.x, u1);    /* Rx <- X3 := E*F */
        m.mul(R.y, R.y, u0);    /* Ry <- Y3 := G*H */
        if (output_type == TWISTED_EDWARDS_ext)
            m.mul(R.t, u0, u1); /* Rt <- T3 := E*H */
    } else { /* if (output_type == MONTGOMERY_xz) */
#ifdef SAFE_TWISTED_EDWARDS_TO_MONTGOMERY
        m.sub(R.z, R.x, u0);    /* Rz <-       F-H */
        m.mul(R.z, R.z, R.y);   /* Rz <-       G*(F-H) */
        m.add(R.x, R.x, u0);    /* Rx <-       F+H */
        m.mul(R.x, R.x, R.y);   /* Rx <-       G*(F+H) */
#else
        /* Conversion: Edwards completed -> Montgomery XZ-only
         *                 ((E:G),(H,F)) -> (F+H :: F-H)
         *
         * Note: the forgotten Y-coordinate is G(F+H)/E.
         *
         * The above map is correct for every points except ((0:1),(1:1))
         * which is sent to (2 :: 0) instead of the point at infinity
         * (0 :: 0). Nevertheless, in our context, computing with the
         * point (2 :: 0) in place of the point at infinity (0 :: 0)
         * produces the same result because both points have Z = 0 and
         * this property is preserved after a doubling dDBL or a
         * differential addition dADD of two points having Z = 0.
         */
        m.sub(R.z, R.x, u0);    /* Rz <-       F-H */
        m.add(R.x, R.x, u0);    /* Rx <-       F+H */
#endif
    }
}

template<typename layer>
void edwards_add(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, typename layer::Modulus const & m,
        ec_point_coord_type_t const output_type)
{
    edwards_addsub(R, P, Q, 0, m, output_type);
}

template<typename layer>
void edwards_sub(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, typename layer::Modulus const & m,
        ec_point_coord_type_t const output_type)
{
    edwards_addsub(R, P, Q, 1, m, output_type);
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
template<typename layer>
void edwards_dbl(ec_point<layer> & R, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        ec_point_coord_type_t const output_type)
{
    ASSERT_EXPENSIVE(output_type == TWISTED_EDWARDS_ext ||
                     output_type == TWISTED_EDWARDS_proj);

#ifdef ECM_COUNT_OPS
    if (output_type == TWISTED_EDWARDS_proj)
        _count_edwards_dbl++;
    else /* if (output_type == TWISTED_EDWARDS_ext) */
        _count_edwards_dblext++;
#endif

    typename layer::Residue u0(m), u1(m);

    m.sqr(u0, P.x);             /* u0 <-  C := X1^2 */
    m.sqr(u1, P.y);             /* u1 <-  D := Y1^2 */
    m.add(R.x, P.x, P.y);       /* Rx <-       X1+Y1 */
    m.sqr(R.x, R.x);            /* Rx <-  B := (X1+Y1)^2 */
    m.add(R.y, u0, u1);         /* Ry <-       C+D */
    m.sub(u0, u0, u1);          /* u0 <- -F := C-D */
    m.sub(R.x, R.y, R.x);       /* Rx <-    := C+D-B */
    m.sqr(u1, P.z);             /* u1 <-  H := Z1^2 */
    m.add(u1, u1, u1);          /* u1 <-    := 2*H  */
    m.add(u1, u0, u1);          /* u1 <- -J := -F + 2*H */

    if (output_type == TWISTED_EDWARDS_ext)
        m.mul(R.t, R.x, R.y);   /* Rt <- T3 := (B-C-D)*(-C-D) */
    m.mul(R.x, R.x, u1);        /* Rx <- X3 := (B-C-D) * J */
    m.mul(R.y, R.y, u0);        /* Ry <- Y3 := F*(-C-D) */
    m.mul(R.z, u0, u1);         /* Rz <- Z3 := F * J */
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
template<typename layer>
void edwards_tpl(ec_point<layer> & R, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        ec_point_coord_type_t const output_type)
{
    ASSERT_EXPENSIVE(output_type == TWISTED_EDWARDS_ext ||
                     output_type == TWISTED_EDWARDS_proj);

#ifdef ECM_COUNT_OPS
    if (output_type == TWISTED_EDWARDS_proj)
        _count_edwards_tpl++;
    else /* if (output_type == TWISTED_EDWARDS_ext) */
        _count_edwards_tplext++;
#endif

    typename layer::Residue u0(m), u1(m), u2(m), u3(m);

    m.sqr(u0, P.x);             /* u0 <-     := X1^2 */
    m.neg(u0, u0);              /* u0 <- aXX := -X1^2 */
    m.sqr(u1, P.y);             /* u1 <-  YY := Y1^2 */
    m.add(u2, u1, u0);          /* u2 <-  Ap := YY+aXX */
    m.sub(u3, u1, u0);          /* u3 <-        YY-aXX */
    m.sqr(R.t, P.z);            /* Rt <-        Z1^2 */
    m.add(R.t, R.t, R.t);       /* Rt <-        2*Z1^2 */
    m.sub(R.t, R.t, u2);        /* Rt <-        2*Z1^2-Ap */
    m.add(R.t, R.t, R.t);       /* Rt <-   B := 2*(2*Z1^2-Ap) */
    m.mul(u2, u2, u3);          /* u2 <-  AA := Ap*(YY-aXX) */
    m.mul(u0, u0, R.t);         /* u0 <-  xB := aXX*B */
    m.mul(u1, u1, R.t);         /* u1 <-  yB := YY*B */
    m.sub(R.t, u2, u1);         /* Rt <-   F := AA-yB */
    m.add(u1, u1, u2);          /* u1 <-        yB+AA */
    m.add(u3, u2, u0);          /* u3 <-   G := AA+xB */
    m.sub(u2, u0, u2);          /* u2 <-        xB-AA */
    m.mul(u1, P.x, u1);         /* u1 <-  xE := X1*(yB+AA) */
    m.mul(u2, P.y, u2);         /* u2 <-  yH := Y1*(xB-AA) */
    m.mul(u0, P.z, R.t);        /* u0 <-  zF := Z1*F */

    if (output_type == TWISTED_EDWARDS_proj) {
        m.mul(R.x, u1, R.t);    /* Rx <-  X3 := xE*F */
        m.mul(R.y, u2, u3);     /* Ry <-  Y3 := yH*G */
        m.mul(R.z, u0, u3);     /* Rz <-  Z3 := zF*G */
    } else { /* if (output_type == TWISTED_EDWARDS_ext) */
        m.mul(u3, P.z, u3);     /* u3 <-  zG := Z1*G */
        m.mul(R.x, u1, u0);     /* Rx <-  X3 := xE*zF */
        m.mul(R.y, u2, u3);     /* Ry <-  Y3 := yH*zG */
        m.mul(R.z, u0, u3);     /* Rz <-  Z3 := zF*zG */
        m.mul(R.t, u1, u2);     /* Rt <- T3 := xE*yH */
    }
}

/* edwards_smul_ui (R:edwards_ext, P:edwards_ext, k: unsigned long)
 *     R <- k*P
 * R can be the same variable as P.
 */
template<typename layer>
void edwards_smul_ui(ec_point<layer> & R, ec_point<layer> const & P,
        unsigned long const k, typename layer::Modulus const & m)
{
    if (k == 0UL) {
        edwards_point_set_zero(R, m, TWISTED_EDWARDS_ext);
    } else if (k == 1UL) {
        ec_point_set(R, P, m, TWISTED_EDWARDS_ext);
    } else if (k == 2UL) {
        edwards_dbl(R, P, m, TWISTED_EDWARDS_ext);
    } else if (k == 3UL) {
        edwards_tpl(R, P, m, TWISTED_EDWARDS_ext);
    } else {
        /* basic double-and-add */
        ec_point<layer> T(m);

        unsigned long mask = ~(0UL);
        mask -= mask / 2; /* Now the most significant bit of i is set */
        while ((mask & k) == 0)
            mask >>= 1;

        /* Most significant bit of k is 1, do it outside the loop */
        ec_point_set(T, P, m, TWISTED_EDWARDS_ext);
        mask >>= 1;

        for (; mask > 0; mask >>= 1) {
            edwards_dbl(T, T, m, TWISTED_EDWARDS_ext);
            if (k & mask) /* output in ext only on the last iteration */
                edwards_add(T, T, P, m, (mask > 1) ? TWISTED_EDWARDS_proj
                                                   : TWISTED_EDWARDS_ext);
        }

        ec_point_set(R, T, m, TWISTED_EDWARDS_ext);
    }
}

#endif /* CADO_SIEVE_ECM_EC_ARITH_EDWARDS_HPP */
