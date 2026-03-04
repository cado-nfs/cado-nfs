#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "galois_action.hpp"
#include "special-q.hpp"
#include "misc.h"
#include "arith/mod_ul.h"
#include "renumber.hpp"
#include "typedefs.h"
#include "runtime_numeric_cast.hpp"
#include "macros.h"

/* XXX several operations in this file are written with the "usual
 * conversion rules" of the C standard in mind (§6.3.1.8). In particular,
 * expressions such as a-b or CA*a where a is int64_t and b and CA are
 * both uint64_t will always entail an implicit cast of a to uint64_t.
 *
 * We also make sure that we only shift non-negative integers.
 *
 * The apply_ab functions here have been written with integer overflow
 * conditions in mind, so we hope we got them right. They're currently
 * written as ASSERT_ALWAYS, it's reasonably trivial to throw specific
 * exceptions, of course. Note that the runtime_numeric_cast<> calls can
 * also trigger runtime_numeric_cast_failure (either
 * runtime_numeric_cast_underflow or runtime_numeric_cast_overflow). We
 * _hope_ that all overflows are detected. There a a few cases where an
 * overflow is reported even though it could be resolved with some
 * extra effort.
 */

/* action: none
 * x -> x
 * matrix=[[1, 0], [0, 1]]
 */
class galois_action_none final : public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 1;
    }

    unsigned long apply(unsigned long r, unsigned long) const final
    {
        return r;
    }

    special_q apply(special_q const & q) const final { return q; }

    int apply_ab(int64_t &, uint64_t &) const final { return 1; }
    unsigned long apply_ab_cofactor(unsigned long, unsigned long) const final
    {
        return 1;
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        return CA * a + CB * b;
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        mpz_mul_uint64(r, a, CA);
        mpz_addmul_uint64(r, b, CB);
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> x";
    }

    void print_name(std::ostream& os) const final
    {
        os << "Id";
    }
};

/* action: 1/y or autom2.1 or autom2.1g
 * x -> 1/x
 * matrix=[[0, 1], [1, 0]]
 */
class galois_action_inv final: public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 2;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        if (r == 0UL) {
            return p;
        } else if (r == p) {
            return 0UL;
        } else {
            modulusul_t pp;
            residueul_t rr;
            modul_initmod_ul(pp, p);
            modul_init_noset0(rr, pp);
            modul_set_ul_reduced(rr, r, pp);

            modul_inv(rr, rr, pp);      /* 1/r */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        if (q.r == 0)
            return { q.p, q.p, q.side };
        else if (q.r == q.p)
            return { q.p, 0, q.side };
        else {
            return { q.p, q.r.invmod(q.p), q.side };
        }
    }

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long) const final
    {
        return r;
    }

    int apply_ab(int64_t & a, uint64_t & b) const final {
        if (a < 0) {
            auto ua = static_cast<uint64_t>(-a);
            /* since b is nonnegative, as soon as we know that it casts
             * to int64_t, then its opposite also does */
            auto sb = runtime_numeric_cast<int64_t>(b);
            a = -sb;
            b = ua;
            return 1;
        } else {
            auto const ua = static_cast<uint64_t>(a);
            a = runtime_numeric_cast<int64_t>(b);
            b = ua;
            return -1;
        }
    }


    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (-b, -a) ~ (b, a)
         *
         * Rule:
         *  the representative is the pair with the smallest first coefficient
         *  in absolute value.
         *
         * Algo:
         *  |a| < b, we compute the hash of (a,b)
         *  |a| > b, we compute the hash of (b, a)   = (b, |a|)     if a > 0
         *                                  (-b, -a) = (-b, |a|)    if a < 0
         *  |a| = b, it implies (a,b) = (1,1) whose image is itself
         *                   or (a,b) = (-1,1) whose image is itself
         */
        uint64_t const a_abs = safe_abs64(a);
        if (a_abs <= b) {
            return CA * a + CB * b;
        } else {
            return CA * ((a >= 0) ? b : -b) + CB * a_abs;
        }
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        if (mpz_cmpabs(a, b) <= 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CB);
        } else {
            mpz_mul_uint64(r, b, CA);
            mpz_addmul_uint64(r, a, CB);
            if (mpz_sgn(a) < 0)
                mpz_neg(r, r);
        }
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> 1/x";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom2.1g";
    }
};

/* action: _y or autom2.2 or autom2.2g
 * x -> -x
 * matrix=[[-1, 0], [0, 1]]
 */
class galois_action_neg final: public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 2;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        return (r == p || r == 0) ? r : p-r;
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        if (q.r == 0 || q.p == q.r)
            return q;
        else
            return { q.p, q.p-q.r, q.side };
    }

    unsigned long apply_ab_cofactor(unsigned long, unsigned long) const final
    {
        return 1;
    }

    int apply_ab(int64_t & a, uint64_t &) const final {
        ASSERT_ALWAYS(a > 0 || (-a) > 0);
        a = -a;
        return -1;
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (a, -b) ~ (-a, b)
         *
         * Rule:
         *  the representative is the pair with the smallest first coeff > 0.
         */
        return CA * safe_abs64(a) + CB * b;
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        mpz_mul_uint64(r, a, CA);
        mpz_abs(r, r);
        mpz_addmul_uint64(r, b, CB);
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> -x";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom2.2g";
    }
};

/* action: autom3.1 or autom3.1g
 * x -> 1-1/x = (x-1)/x
 * matrix=[[1, -1], [1, 0]]
 */
class galois_action_autom3_1 final: public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 3;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        if (r == 0UL) {
            return p;
        } else if (r == p) {
            return 1UL;
        } else {
            modulusul_t pp;
            residueul_t rr;
            modul_initmod_ul(pp, p);
            modul_init_noset0(rr, pp);
            modul_set_ul_reduced(rr, r, pp);

            modul_inv(rr, rr, pp);      /* 1/r */
            modul_neg(rr, rr, pp);      /* -1/r */
            modul_add1(rr, rr, pp);     /* 1-1/r */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        /* 0 -> infinity -> 1 -> 0 */
        if (q.r == 0)
            return { q.p, q.p, q.side };
        else if (q.r == 1)
            return { q.p, 0, q.side };
        else if (q.p == q.r)
            return { q.p, 1, q.side };
        else {
            auto ri = q.r.invmod(q.p);
            ri -= 1;
            return { q.p, q.p-ri, q.side };
        }
    }

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long) const final
    {
        return r;
    }

    int apply_ab(int64_t & a, uint64_t & b) const final {
        if (a <= 0) {
            /* INT64_MIN <= a <= 0 and 0 <= b <= UINT64_MAX */
            auto na = runtime_numeric_cast<int64_t>(b);
            /* compute b-a, but we don't want to overflow! */
            auto nb = b - a;
            ASSERT_ALWAYS(nb >= b && nb >= static_cast<uint64_t>(-a));
            a = na;
            b = nb;
            return 1;
        } else if (static_cast<uint64_t>(a) < b) {
            /* 0 < a <= INT64_MAX and 0 <= b <= UINT64_MAX and a < b*/
            auto na = runtime_numeric_cast<int64_t>(b);
            auto nb = b - a;
            /* 0 < b - a <= UINT64_MAX */
            a = na;
            b = nb;
            return 1;
        } else {
            /* 0 < a <= INT64_MAX and 0 <= b <= UINT64_MAX and a >= b */
            auto na = -runtime_numeric_cast<int64_t>(b);
            auto nb = a - b;
            /* 0 <= a - b <= -INT64_MIN */
            a = na;
            b = nb;
            return -1;
        }
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (b, b-a) ~ (-b, a-b)
         *  (b-a, -a) ~ (a-b, a)
         *
         * Rule:
         *  the representative is the pair with the smallest first coefficient
         *  (<-> the only pair with the first coeff < 0)
         *
         * Algo:
         *
         *  If a < 0 < b: (a, b) -> (b, b-a) -> (b-a, -a)
         *      the only negative first coeff is the first one: a
         *      we compute the hash of (a, b)
         *  If 0 < a < b: (a, b) -> (b, b-a) -> (b-a, -a) ~ (a-b, a)
         *      the only negative first coeff is the last one: a-b
         *      we compute the hash of (a-b, a)
         *  If 0 < b < a: (a, b) -> (-b, a-b) -> (a-b, a)
         *      the only negative first coeff is the second one: -b
         *      we compute the hash of (-b, a-b)
         */
        if (a <= 0) {
            return CA * a + CB * b;
        } else if (static_cast<uint64_t>(a) < b) {
            return CA * (a - b) + CB * a;
        } else {
            return CA * -b + CB * (a - b);
        }
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        if (mpz_sgn(a) <= 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CB);
        } else {
            mpz_mul_uint64(r, b, CA);
            mpz_neg(r, r);
            mpz_addmul_uint64(r, a, CB);
            if (mpz_cmp(a, b) < 0) {
                mpz_addmul_uint64(r, a, CA);
            } else {
                mpz_submul_uint64(r, b, CB);
            }
        }
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> (x-1)/x";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom3.1g";
    }
};

/* action: autom3.2 or autom3.2g
 * x -> -1-1/x = (-x-1)/x
 * matrix=[[-1, -1], [1, 0]]
 */
class galois_action_autom3_2 final : public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 3;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        if (r == 0UL) {
            return p;
        } else if (r == p) {
            return p-1UL;
        } else {
            modulusul_t pp;
            residueul_t rr;
            modul_initmod_ul(pp, p);
            modul_init_noset0(rr, pp);
            modul_set_ul_reduced(rr, r, pp);

            modul_inv(rr, rr, pp);      /* 1/r */
            modul_add1(rr, rr, pp);     /* 1+1/r */
            modul_neg(rr, rr, pp);      /* -1-1/r */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        /* 0 -> infinity -> -1 -> 0 */
        if (q.r == 0)
            return { q.p, q.p, q.side };
        else if (q.p == q.r)
            return { q.p, q.p-1, q.side };
        else {
            auto ri = q.r.invmod(q.p);
            ri += 1;
            /* at this point, ri is in [1..p], so p-ri is the reduced
             * representative */
            return { q.p, q.p-ri, q.side };
        }
    }

    int apply_ab(int64_t & a, uint64_t & b) const final {
        if (a >= 0) {
            /* 0 <= a <= INT64_MAX and 0 <= b <= UINT64_MAX */
            /* if b <= INT64_MAX, then -b >= INT64_MIN */
            auto na = -runtime_numeric_cast<int64_t>(b);
            /* compute a + b, but we don't want to overflow! */
            auto nb = a + b;
            ASSERT_ALWAYS(nb >= b && nb >= static_cast<uint64_t>(a));
            a = na;
            b = nb;
            return -1;
        } else if (b > static_cast<uint64_t>(-a)) {
            /* INT_MIN <= a < 0 and 0 < -a < b <= UINT64_MAX */
            /* So in particular, 0 < a + b <= UINT64_MAX */
            auto na = -runtime_numeric_cast<int64_t>(b);
            auto nb = a + b;
            a = na;
            b = nb;
            return -1;
        } else {
            /* INT64_MIN <= a < 0 and 0 <= b <= UINT64_MAX
             * and 0 < -a-b <= -INT64_MIN <= UINT64_MAX */
            auto na = runtime_numeric_cast<int64_t>(b);
            auto nb = -a - b;
            a = na;
            b = nb;
            return 1;
        }
    }

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long) const final
    {
        return r;
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (b, -a-b) ~ (-b, a+b)
         *  (-a-b, a) ~ (a+b, -a)
         *
         * Rule:
         *  the representative is the pair with the largest first coeff
         *  (<-> the only pair with the first coeff > 0)
         *
         * Algo:
         *
         *  If 0 < a: (a, b) -> (-b, a+b) -> (-a-b, a)
         *      the only positive coeff is the first one: a
         *      we compute the hash of (a, b)
         *  If -b < a < 0: (a, b) -> (-b, a+b) -> (a+b, -a)
         *      the only positive coeff is the last one: a+b
         *      we compute the hash of (a+b, -a)
         *  If a < -b < 0: (a,b) -> (b, -a-b) -> (a+b, -a)
         *      the only positive coeff is the second one: b
         *      we compute the hash of (b, -a-b)
         */
        if (a > 0) {
            return CA * a + CB * b;
        } else if (b > static_cast<uint64_t>(-a)) {
            /* a+b is converted to unsigned, so it makes little sense to
             * check its sign */
            return CA * (a + b) - CB * a;
        } else {
            return CA * b - CB * (a + b);
        }
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        if (mpz_sgn(a) > 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CB);
        } else {
            mpz_mul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
            if (mpz_cmpabs(a, b) < 0) { /* known: a <= 0 and b > 0*/
                mpz_addmul_uint64(r, a, CA);
            } else {
                mpz_submul_uint64(r, b, CB);
            }
        }
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> (-x-1)/x";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom3.2g";
    }
};

/* action: autom4.1 or autom4.1g
 * x -> -(x+1)/(x-1) = -1-2/(x-1)
 * matrix=[[-1, -1], [1, -1]]
 */
class galois_action_autom4_1 final : public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 4;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        if (r == 1UL) {
            return p;
        } else if (r == p) {
            return p-1UL;
        } else {
            modulusul_t pp;
            residueul_t rr;
            modul_initmod_ul(pp, p);
            modul_init_noset0(rr, pp);
            modul_set_ul_reduced(rr, r, pp);

            modul_sub_ul(rr, rr, 1, pp);    /* r-1 */
            modul_inv(rr, rr, pp);          /* 1/(r-1) */
            modul_add(rr, rr, rr, pp);      /* 2/(r-1) */
            modul_add1(rr, rr, pp);         /* 1+2/(r-1) */
            modul_neg(rr, rr, pp);          /* -1-2/(r-1) */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        /* 0 -> 1 -> infinity -> -1 -> 0 */
        if (q.r == 1)
            return { q.p, q.p, q.side };
        else if (q.p == q.r)
            return { q.p, q.p-1, q.side };
        else {
            auto ri = q.r - 1;
            mpz_invert(ri, ri, q.p);
            ri += ri;
            if (ri >= q.p)
                ri -= q.p;
            ri += 1;
            /* as above, ri is now in [1..p], so -ri is represented by
             * p-ri */
            return { q.p, q.p-ri, q.side };
        }
    }

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long p) const final
    {
        return r ? (r-1) : (p-1);
    }

    int apply_ab(int64_t & a, uint64_t & b) const final {
        int m;
        if (a >= 0) {
            /* the image is (a-b, a+b) with m = -1 */
            /* 0 <= a <= INT64_MAX and 0 <= b <= UINT64_MAX */
            /* compute a + b, but we don't want to overflow! */
            auto nb = a + b;
            ASSERT_ALWAYS(nb >= b && nb >= static_cast<uint64_t>(a));
            /* We have INT64_MIN <= -INT64_MAX <= b-a <= UINT64_MAX
             * So no underflow can occur, but overflow is possible
             */
            auto mna = b - a;
            if (b >= static_cast<uint64_t>(a)) {
                /* 0 <= b-a <= UINT64_MAX ; overflow is possible */
                a = -runtime_numeric_cast<int64_t>(mna);
            } else {
                /* -INT64_MAX <= b-a < 0 ; no underflow */
                a = -static_cast<int64_t>(mna);
            }
            b = nb;
            m = -1;
        } else if (b >= static_cast<uint64_t>(-a)) {
            /* the image is (a-b, a+b) with m = -1 */
            /* INT_MIN <= a < 0 and 0 < -a <= b <= UINT64_MAX */
            /* so that 0 <= a + b <= UINT64_MAX */
            auto nb = a + b;
            /* first compute b - a, make sure it doesn't wrap around with
             * uint64_t's, and then check against -INT64_MIN (well,
             * technically we're checking against INT64_MAX and we're
             * missing the possibility that the final a could be
             * INT64_MIN) */
            auto mna = b - a;
            ASSERT_ALWAYS(mna >= b);
            ASSERT_ALWAYS(mna >= static_cast<uint64_t>(-a));
            a = -runtime_numeric_cast<int64_t>(mna);
            b = nb;
            m = -1;
        } else {
            /* the image is (b-a, -a-b) with m = 1 */
            /* INT_MIN <= a < 0 and 0 <= b < -a <= -INT64_MIN <= UINT64_MAX */
            /* INT_MIN <= a + b < 0, so 0 < -a-b < -INT64_MIN <= UINT64_MAX */
            /* b-a might overflow the int64_t range:
             * 0 <= b < b - a <= UINT64_MAX - INT64_MIN */
            auto nb = -a - b;
            auto na = b - a;
            /* first make sure we don't wrap around in uint64_t, then
             * check the int64_t range */
            ASSERT_ALWAYS(na >= b);
            ASSERT_ALWAYS(na >= static_cast<uint64_t>(-a));
            a = runtime_numeric_cast<int64_t>(na);
            b = nb;
            m = 1;
        }
        /* if we have the same parity, we need to take out a factor
         * of 2. Note that this is simply observed by the parity of
         * any of the resulting values.
         *
         * XXX doing this check late causes a few extraneous overflows
         * that we could possibly avoid. However it's not entirely
         * trivial, and probably not worth the trouble. One way to go
         * could be to check nb&1 early on, and if even, compute (for
         * example) (a+b)/2 as (a/2)+(b/2)+(a&1), which doesn't overflow.
         */
        if (b % 2 == 0) {
            b /= 2;
            a /= 2;
            m *= 2;
        }
        return m;
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (-a+b, -a-b) ~ (a-b, a+b)
         *  (-2b, 2a) ~ (-b, a) ~ (b, -a)
         *  (b+a, b-a) ~ (-b-a, a-b)
         *
         * Rule:
         *  the representative is the pair with the largest first coeff
         *
         * Algo:
         *  If a < -b < 0: (a, b) -> (b-a, -(a+b)) -> (b, -a) -> (a+b, b-a)
         *      the largest coeff is the second one: b-a
         *      we compute the hash of (b-a, -(a+b))
         *  If -b < a < 0: (a, b) -> (a-b, a+b) -> (b, -a) -> (a+b, b-a)
         *      the largest coeff is the third one: b
         *      we compute the hash of (b, -a)
         *  If 0 < a < b: (a, b) -> (a-b, a+b) -> (-b, a) -> (a+b, b-a)
         *      the largest coeff is the last one: a+b
         *      we compute the hash of (a+b, b-a)
         *  If 0 < b < a: (a, b) -> (a-b, a+b) -> (-b, a) -> (-(a+b), a-b)
         *      the largest coeff is the first one: a
         *      we compute the hash of (a, b)
         *  Note: Do not forget to divide by 2 when necessary!
         */
        uint64_t const a_abs = safe_abs64(a);
        if (a > 0 && a_abs > b) {
            return CA * a + CB * b;
        } else if (a < 0 && a_abs < b) {
            return CA * b - CB * a;
        } else if (a > 0) { /* implies 0 < a < b */
            if (a_abs % 2 == b % 2) {
                return CA * ((a + b) >> 1) + CB * ((b - a) >> 1);
            } else {
                return CA * (a + b) + CB * (b - a);
            }
        } else { /* a < -b < 0 */
            if (a_abs % 2 == b % 2) {
                return CA * ((b - a) >> 1) + CB * ((-a - b) >> 1);
            } else {
                return CA * (b - a) - CB * (a + b);
            }
        }
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        if (mpz_sgn(a) > 0 && mpz_cmpabs(a, b) > 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CB);
        } else if (mpz_sgn(a) < 0 && mpz_cmpabs(a, b) < 0) {
            mpz_mul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
        } else if (mpz_sgn(a) > 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_addmul_uint64(r, b, CB);
            if (mpz_congruent_2exp_p(a, b, 1u)) {
                mpz_fdiv_q_2exp(r, r, 1u);
            }
        } else {
            mpz_mul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_submul_uint64(r, b, CB);
            if (mpz_congruent_2exp_p(a, b, 1u)) {
                mpz_fdiv_q_2exp(r, r, 1u);
            }
        }
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> -(x+1)/(x-1)";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom4.1g";
    }
};

/* action: autom6.1 or autom6.1g
 * x -> -(2*x+1)/(x-1)
 * matrix=[[-2, -1], [1, -1]]
 */
class galois_action_autom6_1 final : public galois_action_impl_base
{
public:
    unsigned int get_order() const final
    {
        return 6;
    }

    unsigned long apply(unsigned long r, unsigned long p) const final
    {
        if (r == 1UL) {
            return p;
        } else if (r == p) {
            return p-2UL;
        } else {
            modulusul_t pp;
            residueul_t rr, three, two;
            modul_initmod_ul(pp, p);
            modul_init_noset0(rr, pp);
            modul_init_noset0(three, pp);
            modul_init_noset0(two, pp);
            modul_set_ul_reduced(rr, r, pp);
            modul_set_ul(three, 3UL, pp);
            modul_set_ul(two, 2UL, pp);

            modul_sub_ul(rr, rr, 1, pp);    /* r-1 */
            modul_inv(rr, rr, pp);          /* 1/(r-1) */
            modul_mul(rr, three, rr, pp);   /* 3/(r-1) */
            modul_add(rr, two, rr, pp);     /* 2+3/(r-1) */
            modul_neg(rr, rr, pp);          /* -2-3/(r-1) */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(two, pp);
            modul_clear(three, pp);
            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
    }

    special_q apply(special_q const & q) const final {
        ASSERT_ALWAYS(q.is_prime());
        /* 0 -> 1 -> infinity -> -2 -> -1 -> -1/2 -> 0 */
        if (q.r == 1)
            return { q.p, q.p, q.side };
        else if (q.p == q.r)
            return { q.p, q.p-2, q.side };
        else {
            auto ri = q.r - 1;
            mpz_invert(ri, ri, q.p);
            ri *= 3;
            if (ri >= q.p) ri -= q.p;
            if (ri >= q.p) ri -= q.p;
            ri += 2;    /* ri == 2 + 3/(r-1) */
            /* at this point, ri is in [2..p+1]. We have two cases */
            if (ri > q.p) {
                /* ri is p+1, which occurs only if 2+3/(r-1)==1, so r=-2.
                 * In this case, we want to return -1.
                 */
                return { q.p, q.p-1, q.side };
            } else {
                /* ri is in [1..p], so -ri is represented by p-ri */
                return { q.p, q.p-ri, q.side };
            }
        }
    }

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long p) const final
    {
        return r ? (r-1) : (p-1);
    }

    int apply_ab(int64_t & a, uint64_t & b) const final {
        /* It's very tricky here to avoid all overflow conditions. We're
         * ever so slightly approximating them by assuming that two_b
         * doesn't overflow. It's probably possible to refine this a
         * little bit. Beyond that, it's quite similar to the previous
         * code, with a few key differences.
         */
        auto two_b = b + b;
        ASSERT_ALWAYS(two_b >= b);
        int m;
        if (a >= 0) {
            /* the image is (a-b, a+2b) with m = -1 */
            /* 0 <= a <= INT64_MAX and 0 <= 2b <= UINT64_MAX */
            /* compute a + 2b, but we don't want to overflow! */
            auto nb = a + two_b;
            ASSERT_ALWAYS(nb >= two_b && nb >= static_cast<uint64_t>(a));
            /* We have INT64_MIN <= -INT64_MAX <= 2b-a <= UINT64_MAX
             * So no underflow can occur, but overflow is possible
             */
            auto mna = b - a;
            if (b >= static_cast<uint64_t>(a)) {
                /* 0 <= b-a <= UINT64_MAX ; overflow is possible */
                a = -runtime_numeric_cast<int64_t>(mna);
            } else {
                /* -INT64_MAX <= b-a < 0 ; no underflow */
                a = -static_cast<int64_t>(mna);
            }
            b = nb;
            m = -1;
        } else if (two_b >= static_cast<uint64_t>(-a)) {
            /* the image is (a-b, a+2b) with m = -1 */
            /* INT_MIN <= a < 0 and 0 < -a <= 2b <= UINT64_MAX */
            /* so that 0 <= a + 2b <= UINT64_MAX */
            /* Note that there is a small possibility here, if we disable
             * the two_b overflow test, that 2b+a still fits. We don't
             * cover this case.
             */
            auto nb = a + two_b;
            /* first compute b - a, make sure it doesn't wrap around with
             * uint64_t's, and then check against -INT64_MIN (well,
             * technically we're checking against INT64_MAX and we're
             * missing the possibility that the final a could be
             * INT64_MIN) */
            auto mna = b - a;
            ASSERT_ALWAYS(mna >= b);
            ASSERT_ALWAYS(mna >= static_cast<uint64_t>(-a));
            a = -runtime_numeric_cast<int64_t>(mna);
            b = nb;
            m = -1;
        } else {
            /* the image is (b-a, -a-2b) with m = 1 */
            /* INT_MIN <= a < 0 and 0 <= 2b < -a <= -INT64_MIN <= UINT64_MAX */
            /* INT_MIN <= a + 2b < 0, so 0 < -a-2b < -INT64_MIN <= UINT64_MAX */
            /* b-a might overflow the int64_t range:
             * 0 <= b < b - a <= UINT64_MAX - INT64_MIN */
            auto nb = -a - two_b;
            auto na = b - a;
            /* first make sure we don't wrap around in uint64_t, then
             * check the int64_t range */
            ASSERT_ALWAYS(na >= b);
            ASSERT_ALWAYS(na >= static_cast<uint64_t>(-a));
            a = runtime_numeric_cast<int64_t>(na);
            b = nb;
            m = 1;
        }
        /* basically the same discussion as above remains valid here. */
        if (b % 3 == 0) {
            b /= 3;
            a /= 3;
            m *= 3;
        }
        return m;
    }

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        /* Orbit of (a,b):
         *  (a, b)
         *  (-a+b, -a-2*b) ~ (a-b, a+2*b)
         *  (-3*b, 3*a+3*b) ~ (-b, a+b) ~ (b, -a-b)
         *  (a+2*b, -2*a-b) ~ (-a-2*b, 2*a+b)
         *  (-3*a-3*b, 3*a) ~ (-a-b, a) ~ (a+b, -a)
         *  (2*a+b, -a+b) ~ (-2*a-b, a-b)
         *
         * Rule:
         *  the representative is the pair with the largest first coeff
         *
         * Algo:
         *  If a < -2b < 0: (a, b) -> (-a+b, -a-2*b) -> (b, -a-b)
         *                      -> (a+2*b, -2*a-b) -> (a+b, -a) -> (2*a+b, -a+b)
         *      the largest coeff is the second one: -a+b
         *      we compute the hash of (-a+b, -a-2*b)
         *  If -2b < a < -b < 0: (a, b) -> (a-b, a+2*b) -> (b, -a-b)
         *                      -> (a+2*b, -2*a-b) -> (a+b, -a) -> (2*a+b, -a+b)
         *      the largest coeff is the third one: (b, -a-b)
         *      we compute the hash of (b, -a-b)
         *  If -b < a < -b/2 < 0: (a, b) -> (a-b, a+2*b) -> (-b, a+b)
         *                      -> (a+2*b, -2*a-b) -> (a+b, -a) -> (2*a+b, -a+b)
         *      the largest coeff is the fourth one: (a+2*b)
         *      we compute the hash of (a+2*b, -2*a-b)
         *  If -b/2 < a < 0: (a, b) -> (a-b, a+2*b) -> (-b, a+b)
         *                      -> (-a-2*b, 2*a+b) -> (a+b, -a) -> (2*a+b, -a+b)
         *      the largest coeff is the fifth one: a+b
         *      we compute the hash of (a+b, -a)
         *  If 0 < a < b: (a, b) -> (a-b, a+2*b) -> (-b, a+b) -> (-a-2*b, 2*a+b)
         *                              -> (-a-b, a) -> (2*a+b, -a+b)
         *      the largest coeff is the last one: 2*a+b
         *      we compute the hash of (2*a+b, -a+b)
         *  If 0 < b < a: (a, b) -> (a-b, a+2*b) -> (-b, a+b) -> (-a-2*b, 2*a+b)
         *                              -> (-a-b, a) -> (-2*a-b, a-b)
         *      the largest coeff is the first one: a
         *      we compute the hash of (a, b)
         *  Note: Do not forget to divide by 3 when necessary!
         */
        uint64_t const a_abs = safe_abs64(a);
        if (a > 0 && a_abs > b) { /* 0 < b < a */
            return CA * a + CB * b;
        } else if (a > 0) { /* 0 < a < b */
            if (a_abs % 3 == b % 3) {
              return CA * ((2*a+b)/3) + CB * ((b - a)/3);
            } else {
              return CA * (2*a+b) + CB * (b - a);
            }
        } else if (2*a_abs < b) { /* -b/2 < a < 0 */
            return CA * (a + b) - CB * a;
        } else if (a_abs < b) { /* -b < a < -b/2 < 0 */
            if ((3 - (a_abs % 3)) % 3 == b % 3) {
                return CA * (((b << 1) + a)/3) + CB * (((a_abs << 1)-b)/3);
            } else {
                return CA * ((b << 1) + a) + CB * ((a_abs << 1)-b);
            }
        } else if (a_abs < 2*b) { /* -2b < a < -b < 0 */
            return CA * b - CB * (a + b);
        } else { /* a < -2b < 0 */
            if ((3 - (a_abs % 3)) % 3 == b % 3) {
                /* note that negation does not commute with division by 3.
                 * For example if y = a+(b<<1) == (uint64_t) -51 then we
                 * have -(y/3) != (y/3)
                 */
                return CA * ((b - a)/3) + CB * ((-a - (b << 1))/3);
            } else {
                return CA * (b - a) - CB * (a + (b << 1));
            }
        }
    }

    uint64_t
    hash_ab(mpz_srcptr a, mpz_srcptr b, uint64_t CA, uint64_t CB) const final
    {
        ASSERT(mpz_sgn(b) >= 0);
        cxx_mpz r;
        cxx_mpz s;
        mpz_add(s, a, b);
        if (mpz_sgn(a) > 0 && mpz_cmpabs(a, b) > 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CB);
        } else if (mpz_sgn(a) > 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_mul_2exp(r, r, 1u);
            mpz_addmul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_addmul_uint64(r, b, CB);
            if (mpz_fdiv_ui(a, 3u) == mpz_fdiv_ui(b, 3u)) {
                mpz_divexact_ui(r, r, 3u);
            }
        } else if (mpz_cmpabs(a, s) < 0) {
            mpz_mul_uint64(r, a, CA);
            mpz_addmul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
        } else if (mpz_cmpabs(a, b) < 0) {
            mpz_mul_uint64(r, b, CA);
            mpz_mul_2exp(r, r, 1u);
            mpz_addmul_uint64(r, a, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_submul_uint64(r, a, CB);
            mpz_submul_uint64(r, b, CB);
            if (mpz_fdiv_ui(a, 3u) == mpz_fdiv_ui(b, 3u)) {
                mpz_divexact_ui(r, r, 3u);
            }
        } else if (mpz_sgn(s) < 0 && -mpz_cmpabs(b, s) < 0) {
            mpz_mul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_submul_uint64(r, b, CB);
        } else {
            mpz_mul_uint64(r, a, CA);
            mpz_neg(r, r);
            mpz_addmul_uint64(r, b, CA);
            mpz_submul_uint64(r, a, CB);
            mpz_submul_uint64(r, b, CB);
            mpz_submul_uint64(r, b, CB);
            if (mpz_fdiv_ui(a, 3u) == mpz_fdiv_ui(b, 3u)) {
                mpz_divexact_ui(r, r, 3u);
            }
        }
        return mpz_sgn(r) * mpz_tdiv_uint64(r, 0xffffffffffffffff);
    }

    void print_action(std::ostream& os) const final
    {
        os << "x -> -(2*x+1)/(x-1)";
    }

    void print_name(std::ostream& os) const final
    {
        os << "autom6.1g";
    }
};

galois_action::galois_action(const std::string &action)
{
    if (action == "none" || action.empty() || action == "id" || action == "Id")
        impl = std::make_unique<galois_action_none>();
    else if (action == "autom2.1" || action == "autom2.1g" || action == "1/y")
        impl = std::make_unique<galois_action_inv>();
    else if (action == "autom2.2" || action == "autom2.2g" || action == "_y")
        impl = std::make_unique<galois_action_neg>();
    else if (action == "autom3.1" || action == "autom3.1g")
        impl = std::make_unique<galois_action_autom3_1>();
    else if (action == "autom3.2" || action == "autom3.2g")
        impl = std::make_unique<galois_action_autom3_2>();
    else if (action == "autom4.1" || action == "autom4.1g")
        impl = std::make_unique<galois_action_autom4_1>();
    else if (action == "autom6.1" || action == "autom6.1g")
        impl = std::make_unique<galois_action_autom6_1>();
    else
        throw std::runtime_error("invalid action string");
}

unsigned int galois_action::get_order () const
{
    return impl->get_order();
}

unsigned long galois_action::apply (unsigned long r, unsigned long p) const
{
    return impl->apply(r, p);
}

unsigned long galois_action::apply_ab_cofactor(unsigned long r, unsigned long p) const
{
    return impl->apply_ab_cofactor(r, p);
}
int galois_action::apply_ab(int64_t & a, uint64_t & b) const
{
    return impl->apply_ab(a, b);
}

special_q galois_action::apply(special_q const & q) const
{
    return impl->apply(q);
}

uint64_t galois_action::hash_ab(int64_t a, uint64_t b,
                                uint64_t CA, uint64_t CB) const
{
    return impl->hash_ab(a, b, CA, CB);
}

uint64_t galois_action::hash_ab(mpz_srcptr a, mpz_srcptr b,
                                uint64_t CA, uint64_t CB) const
{
    return impl->hash_ab(a, b, CA, CB);
}

size_t galois_action::compute_action_on_index(std::vector<index_t> &sigma,
                                              renumber_t const & tab) const
{
    size_t norbits = 0;
    unsigned long const order = get_order();
    /* Map of all r corresponding to the current (p, side) to its index. */
    std::unordered_map<p_r_values_t, index_t> index_of_r;

    sigma.resize(tab.size());

    index_t i = 0;
    for (auto it = tab.begin() ; it != tab.end() ; ) {
        auto x = *it;
        if (x.p == 0) {
            sigma[i] = i; /* extra columns are unchanged by galois action */
            ++i, ++it;
        } else if (tab.is_bad(x)) {
            /* Bad ideals are left unchanged.
             * It should be okay in most cases.
             *  - for factorization, filter_galois is not used
             *  - for dl, it only means that columns corresponding to badideals
             *    are left unchanged instead of keeping one for each orbit. It
             *    should not change the correctness of the matrix. It only
             *    means that we could have reduced the number of columns a
             *    little bit more. But the number of badideals is very small so
             *    it should be negligeable (a few dozens columns that we could
             *    have deleted in the absolute worst case)
             *  - for class group computation using quadratic sieve, it is the
             *    same as for dl. And the number of badideals should be even
             *    smaller.
             */
            /* XXX Note on how to (maybe) handle bad ideals:
             *  Let (p, r) be a tuple corresponding to nbad bad ideals and let
             *  sigma be the galois action of order order.
             *  - If sigma(r) = r mod p. Does it imply order == nbad ?
             *    If it it the case, all bad ideals above (p, r) are in the same
             *    orbit and one can choose any one as the representative.
             *  - If sigma(r) = r' mod p, r' != r. Does it imply that (p, r')
             *    also corresponds to nbad bad ideals ? If it is the case, one
             *    should also check that the branches of the two are compatible
             *    (i.e., a branch (p^k, rb) with rb = r mod p should have the
             *    same rule for handling exponents as the corresponding branch
             *    (p^k, rb') with sigma(rb) = rb' mod p^k and rb' = r' mod p).
             *
             *    Unanswered question: how does one rewrite exponents ?
             */
            sigma[i] = i;
            ++i, ++it;
        } else {
            index_of_r.clear();
            renumber_t::p_r_side idc = x;
            renumber_t::p_r_side const id0 = idc;
            do
            {
                index_of_r.emplace(idc.r, i);
                ++it;
                ++i;
                if (it == tab.end())
                    break;
                idc = *it;
            } while (idc.same_p(id0));

            /* The map index_of_r now contains values that we need
             * to group into orbits of size order or orbits of size 1.
             */
            for (auto const &elt: index_of_r) {
                /* compute the orbit of r0=elt.first, using the index
                 * i0=elt.second as representative for the whole orbit
                 */
                unsigned int n;
                p_r_values_t const r0 = elt.first;
                index_t const i0 = elt.second;

                p_r_values_t sigma_r = apply(r0, id0.p);

                sigma[i0] = i0;
                for (n = 1; sigma_r != r0; n++) {
                    if (index_of_r.count(sigma_r) == 0) {
                        std::cerr << "Error, while computing orbit of "
                                  << "ideal (" << id0.p << ", " << r0
                                  << ") on side " << id0.side << ", got r="
                                  << sigma_r << " which does not correspond"
                                  << " to an ideal" << "\n";
                        throw std::runtime_error("wrong orbit computation");
                    }
                    sigma[index_of_r[sigma_r]] = i0;
                    index_of_r.erase(sigma_r);
                    sigma_r = apply(sigma_r, id0.p);
                }
                if (n != 1 && n != order) {
                    std::cerr << "Error, orbit of ideal (" << id0.p << ", "
                              << r0 << ") on side " << id0.side
                              << " has an orbit of length " << n
                              << "instead of 1 or " << order << "\n";
                    throw std::runtime_error("wrong orbit size");
                }
                norbits += n == order;
            }
        }
    }
    return norbits;
}

std::ostream& operator<<(std::ostream &os, const galois_action &g)
{
    os << "Galois action ";
    g.impl->print_name(os);
    os << " (";
    g.impl->print_action(os);
    os << ", order=" << g.impl->get_order() << ")";
    return os;
}
