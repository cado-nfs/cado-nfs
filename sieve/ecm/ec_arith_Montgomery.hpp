#ifndef CADO_SIEVE_ECM_EC_ARITH_MONTGOMERY_HPP
#define CADO_SIEVE_ECM_EC_ARITH_MONTGOMERY_HPP

#include <cstdint>
#include <cstdio>

#include <type_traits>

#include "cxx_mpz.hpp"

#include "ec_arith_common.hpp"
#include "macros.h"
#ifdef ECM_COUNT_OPS
#include "ec_arith_cost.h"
#endif

/* Define to 1 to make montgomery_dadd() test if the two points are
   identical, and call montgomery_dbl() if they are */
#ifndef ELLM_SAFE_ADD
#define ELLM_SAFE_ADD 0
#endif

#ifdef ECM_COUNT_OPS
/* not thread-safe */
static unsigned int _count_montgomery_dadd, _count_montgomery_dbl;
#define MONTGOMERY_COUNT_OPS_M _count_montgomery_dadd * MONTGOMERY_dADD \
                             + _count_montgomery_dbl * MONTGOMERY_DBL
#define MONTGOMERY_COUNT_OPS_RESET() do {                   \
      _count_montgomery_dadd = _count_montgomery_dbl = 0;   \
    } while (0)
#endif

/* Montgomery elliptic curves
 *
 * XZ-only coordinates, with equation:
 *    B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2
 *
 * Constant needed in computation: b = (A+2)/4
 */

/* Compute A = 4*b-2. A and b can be the same variable. */
template<typename layer>
void montgomery_A_from_b(typename layer::Residue & A,
        typename layer::Residue const & b,
        typename layer::Modulus const & m)
{
    m.add(A, b, b);    /* A <- b+b = 2b */
    m.add(A, A, A);    /* A <- 4b */
    m.sub1(A, A);
    m.sub1(A, A);      /* A <- 4b-2 */
}

template<typename layer>
void montgomery_curve_fprintf(FILE * out, char const * prefix,
        typename layer::Residue const & A, ec_point<layer> const * P,
        typename layer::Modulus const & m)
{
    char const * pre = (prefix == nullptr) ? "" : prefix;

    fmt::print(out, "{}Montgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2\n"
                    "{}A = {:#x}\n", pre, pre, ec_residue_get_mpz<layer>(A, m));

    if (P) {
        fmt::print(out, "{}with point (X::Z) = ", pre);
        ec_point_fprintf(out, *P, MONTGOMERY_xz, m);
        fputc('\n', out);
    }
}

/* Set P to zero (the neutral point): (0::0) */
template<typename layer>
void montgomery_point_set_zero(ec_point<layer> & P,
        typename layer::Modulus const & m)
{
    m.set0(P.x);
    m.set0(P.z);
}

/* Set Q to the same point as P but with z = 1.
 * Return 1 if it worked, 0 if the computation of the modular inverse of P.z
 * failed.
 * P and Q can be the same variables
 */
template<typename layer>
int montgomery_point_to_affine(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m)
{
    typename layer::Residue t(m);

    int const ret = m.inv(t, P.z);
    if (ret) {
        m.mul(Q.x, P.x, t);
        m.set1(Q.z);
    }
    return ret;
}

/* Convert the point P in affine before printing it (mostly used for debug) */
template<typename layer>
void montgomery_point_fprintf_affine(FILE * out, ec_point<layer> const & P,
        typename layer::Modulus const & m)
{
    ec_point<layer> Paff(m);
    if (montgomery_point_to_affine(Paff, P, m))
        ec_point_fprintf(out, Paff, MONTGOMERY_xz, m);
    else
        fmt::print(out, "(not a affine point)");
}

template<typename layer>
void montgomery_point_from_edwards_point(ec_point<layer> & PM,
        ec_point<layer> const & PE, int const out_aff,
        typename layer::Modulus const & m)
{
    m.add(PM.x, PE.z, PE.y);
    m.sub(PM.z, PE.z, PE.y);
    if (out_aff)
        montgomery_point_to_affine(PM, PM, m);
}

/* montgomery_dbl (Q, P)
 *     Q <- 2*P
 * It is permissible to let P and Q use the same memory.
 * Cost:
 *    3M + 2S + 4add/sub
 */
template<typename layer>
ATTRIBUTE_ALWAYS_INLINE
inline void montgomery_dbl_inl(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
#if ELLM_SAFE_ADD
    if (m.is0(P.z)) {
        ASSERT(m.is0(P.x));
    }
#endif

#ifdef ECM_COUNT_OPS
    _count_montgomery_dbl++;
#endif

    typename layer::Residue u(m), v(m), w(m);

    m.add(u, P.x, P.z);
    m.sqr(u, u);            /* u = (x + z)^2 */
    m.sub(v, P.x, P.z);
    m.sqr(v, v);            /* v = (x - z)^2 */
    m.mul(Q.x, u, v);       /* x2 = (x^2 - z^2)^2 */
    m.sub(w, u, v);         /* w = 4 * x * z */
    m.mul(u, w, b);         /* u = x * z * (A + 2) */
    m.add(u, u, v);         /* u = x^2 + x * z * A + z^2 */
    m.mul(Q.z, w, u);       /* Q_z = (4xz) * (x^2 + xzA + z^2) */

#if ELLM_SAFE_ADD
    if (m.is0(Q.z))
        m.set0(Q.x);
#endif
}

/* montgomery_dadd (R, P, Q, D)
 *     R <- P+Q with D=P-Q or D=Q-P
 * R may be identical to P, Q and/or D.
 * This function assumes that P !~= Q, i.e. that there is
 * no t!=0 so that P.x = t*Q.x and P.z = t*Q.z, for otherwise the result
 * is (0:0) although it shouldn't be (which actually is good for factoring!).
 * Cost:
 *    4M + 2S + 6add/sub
 */
template<typename layer>
ATTRIBUTE_ALWAYS_INLINE
inline void montgomery_dadd_inl(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, ec_point<layer> const & D,
        typename layer::Residue const & b MAYBE_UNUSED,
        typename layer::Modulus const & m)
{
#ifdef ECM_COUNT_OPS
    _count_montgomery_dadd++;
#endif

#if ELLM_SAFE_ADD
    /* Handle case where at least one input point is point at infinity */
    if (m.is0(P.z)) {
        ASSERT(m.is0(P.x));
        ec_point_set(R, Q, m, MONTGOMERY_xz);
        return;
    }
    if (m.is0(Q.z)) {
        ASSERT(m.is0(Q.x));
        ec_point_set(R, P, m, MONTGOMERY_xz);
        return;
    }
#endif

    typename layer::Residue u(m), v(m), w(m);

    m.sub(u, P.x, P.z);
    m.add(v, Q.x, Q.z);
    m.mul(u, u, v);         /* u = (Px-Pz)*(Qx+Qz) */
    m.add(w, P.x, P.z);
    m.sub(v, Q.x, Q.z);
    m.mul(v, w, v);         /* v = (Px+Pz)*(Qx-Qz) */
    m.add(w, u, v);         /* w = 2*(Qx*Px - Qz*Pz)*/
    m.sub(v, u, v);         /* v = 2*(Qz*Px - Qx*Pz) */
#if ELLM_SAFE_ADD
    /* Check if v == 0, which happens if P=Q or P=-Q.
       If P=-Q, set result to point at infinity.
       If P=Q, use montgomery_dbl() instead.
       This test only works if P=Q on the pseudo-curve modulo N, i.e.,
       if N has several prime factors p, q, ... and P=Q or P=-Q on E_p but
       not on E_q, this test won't notice it. */
    if (m.is0(v)) {
        /* Test if difference is point at infinity */
        if (m.is0(D.z)) {
            ASSERT(m.is0(D.x));
            montgomery_dbl_inl(R, P, m, b); /* Yes, points are identical, use doubling */
        } else {
            montgomery_point_set_zero(R, m); /* Set result to point at infinity */
        }
        return;
    }
#endif
    m.sqr(w, w);            /* w = 4*(Qx*Px - Qz*Pz)^2 */
    m.sqr(v, v);            /* v = 4*(Qz*Px - Qx*Pz)^2 */
    m.set(u, D.x);          /* save D.x */
    m.mul(R.x, w, D.z);     /* may overwrite D.x */
    m.mul(R.z, u, v);
}

/* In the bytecode interpreters, which chain these operations, the
 * fixed-size layers need them inlined: gcc does not do it by itself, and
 * the calls cost 10% on 64 bits. For mod_mpz_new, inlining them costs 6%
 * instead, so they stay out of line there.
 */
template<typename layer>
constexpr bool montgomery_inline_ops = !std::is_same_v<typename layer::Integer, cxx_mpz>;

template<typename layer>
NO_INLINE
void montgomery_dbl_ool(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    montgomery_dbl_inl(Q, P, m, b);
}

template<typename layer>
NO_INLINE
void montgomery_dadd_ool(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, ec_point<layer> const & D,
        typename layer::Residue const & b,
        typename layer::Modulus const & m)
{
    montgomery_dadd_inl(R, P, Q, D, b, m);
}

template<typename layer>
ATTRIBUTE_ALWAYS_INLINE
inline void montgomery_dbl(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    if constexpr (montgomery_inline_ops<layer>)
        montgomery_dbl_inl(Q, P, m, b);
    else
        montgomery_dbl_ool(Q, P, m, b);
}

template<typename layer>
ATTRIBUTE_ALWAYS_INLINE
inline void montgomery_dadd(ec_point<layer> & R, ec_point<layer> const & P,
        ec_point<layer> const & Q, ec_point<layer> const & D,
        typename layer::Residue const & b,
        typename layer::Modulus const & m)
{
    if constexpr (montgomery_inline_ops<layer>)
        montgomery_dadd_inl(R, P, Q, D, b, m);
    else
        montgomery_dadd_ool(R, P, Q, D, b, m);
}

/* montgomery_smul_ul (R, Rp1, P, k)
 *     R <- k*P
 *     Rp1 <- (k+1)*P       if Rp1 != nullptr
 * R or Rp1 can be the same variable as P.
 * Cost:
 *    log2(k) DBL and log2(k)-1 dADD to compute R and Rp1
 *    computing only R save 1 DBL or 1 dADD depending on the parity of k
 *    If the second most significant bit of k is 0, we save 1 DBL
 */
template<typename layer>
void montgomery_smul_ul(ec_point<layer> & R, ec_point<layer> * Rp1,
        ec_point<layer> const & P, unsigned long const k,
        typename layer::Modulus const & m,
        typename layer::Residue const & b)
{
    ec_point<layer> T0(m), T1(m);

    if (k == 0UL) {
        montgomery_point_set_zero(R, m);
        if (Rp1)
            ec_point_set(*Rp1, P, m, MONTGOMERY_xz);
    } else if (k == 1UL) {
        ec_point_set(R, P, m, MONTGOMERY_xz);
        if (Rp1)
            montgomery_dbl(*Rp1, P, m, b);
    } else if (k == 2UL) {
        if (Rp1) {
            montgomery_dbl(T1, P, m, b);
            montgomery_dadd(*Rp1, T1, P, P, b, m);
            ec_point_set(R, T1, m, MONTGOMERY_xz);
        } else {
            montgomery_dbl(R, P, m, b);
        }
    } else { /* k >= 3 */
        /* Montgomery Ladder */
        unsigned long mask = ~(0UL);
        mask -= mask / 2; /* Now the most significant bit of i is set */
        while ((mask & k) == 0)
            mask >>= 1;

        /* Most significant bit of k is 1, do it outside the loop */
        ec_point_set(T0, P, m, MONTGOMERY_xz); /* starting value T0 = P */
        montgomery_dbl(T1, P, m, b);           /* starting value T1 = 2P */
        mask >>= 1;

        /* If the second most significant bit of k is 0, then we do the
         * iteration manually (to avoid to compute again 2P)
         * As k >= 3, we know that in this case k has at least 3 bits.
         */
        if (!(k & mask)) { /* (T0,T1) <- (2*P, 3*P) */
            ec_point_set(T0, T1, m, MONTGOMERY_xz);
            montgomery_dadd(T1, T1, P, P, b, m);
            mask >>= 1;
        }

        for (; mask > 1; mask >>= 1) { /* loop invariant: T1-T0 = P */
            if (k & mask) { /* (T0,T1) <- (T1+T0, 2*T1) */
                montgomery_dadd(T0, T1, T0, P, b, m);
                montgomery_dbl(T1, T1, m, b);
            } else { /* (T0,T1) <- (2*T0, T1+T0) */
                montgomery_dadd(T1, T1, T0, P, b, m);
                montgomery_dbl(T0, T0, m, b);
            }
        }

        /* Deal with least significant bit outside the loop */
        if (k & mask) {
            montgomery_dadd(R, T1, T0, P, b, m);
            if (Rp1)
                montgomery_dbl(*Rp1, T1, m, b);
        } else {
            montgomery_dbl(R, T0, m, b);
            if (Rp1)
                montgomery_dadd(*Rp1, T1, T0, P, b, m);
        }
    }
}

/* Return the order of a Montgomery curve, using the Jacobi symbol.
 * The curve is described by the curve coefficient A and a point P (which
 * gives the Jacobi symbol of the curve coefficient B).
 * P must not be a point of order 2 and P.z must be invertible modulo m.
 * This has complexity O(m): it is only meant for small moduli, which
 * must fit in 64 bits.
 */
template<typename layer>
uint64_t montgomery_curve_order(typename layer::Residue const & A,
        ec_point<layer> const & P,
        typename layer::Modulus const & m)
{
    typename layer::Residue x(m), t(m), one(m);

    /* Compute x = P.x/P.z mod m */
    if (!m.inv(x, P.z))
        return 0;
    m.mul(x, x, P.x);

    m.set1(one);

    /* Compute x^3 + A*x^2 + x and see if it is a square */
    m.set(t, x);
    m.add(t, t, A);
    m.mul(t, t, x);
    m.add(t, t, one);
    m.mul(t, t, x);
    int const bchar = m.jacobi(t);
    ASSERT(bchar != 0);

    uint64_t order = 2; /* One for (0, 0, 1), one for the point at infinity */
    uint64_t const p = uint64_t(m.getmod());
    for (uint64_t i = 1; i < p; i++) {
        m.set(x, i);
        m.set(t, x);
        m.add(t, t, A);
        m.mul(t, t, x);
        m.add(t, t, one);
        m.mul(t, t, x);
        if (bchar == 1)
            order = order + uint64_t(1L + long(m.jacobi(t)));
        else
            order = order + uint64_t(1L - long(m.jacobi(t)));
    }

    return order;
}

#endif /* CADO_SIEVE_ECM_EC_ARITH_MONTGOMERY_HPP */
