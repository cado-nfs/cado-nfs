#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <iostream>         // for std::cerr
#include <stdexcept>        // for std::runtime_error
#include <string>
#include <unordered_map>    // for unordered_map
#include <vector>           // for vector

#include "galois_action.hpp"
#include "misc.h"           // for safe_abs64
#include "arith/mod_ul.h"         // for modul_clear, modul_clearmod, modul_get_ul
#include "renumber.hpp"
#include "typedefs.h"

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

    uint64_t hash_ab(int64_t a, uint64_t b,
                     uint64_t CA, uint64_t CB) const final
    {
        return CA * (uint64_t) a + CB * b;
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
        if (a_abs < b) {
            return CA * (uint64_t) a + CB * b;
        } else {
            return CA * ((a >= 0) ? b : -b) + CB * a_abs;
        }
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
        auto const ua = (uint64_t) a;
        if (a <= 0) {
            return CA * ua + CB * b;
        } else if (ua < b) {
            return CA * (ua-b) + CB * ua;
        } else {
            return CA * -b + CB * (ua-b);
        }
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
        uint64_t const a_abs = safe_abs64(a);
        if (a > 0) {
            return CA * (uint64_t) a + CB * b;
        } else if (a_abs < b) {
            return CA * (b-a_abs) + CB * a_abs;
        } else {
            return CA * b + CB * (a_abs-b);
        }
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

#if 0
    uint64_t next_ab (int64_t *ar, uint64_t *br, int64_t a, uint64_t b) const
    {
        uint64_t a_abs = safe_abs64(a);
        if (a < 0 && a_abs > b) /* a+b < 0 */
        {
            *ar = (int64_t) (b + a_abs);    /* b-a */
            *br = a_abs - b;                /* -a-b */
        } else {
            *ar = a - (int64_t) b;          /* a-b */
            *br = a+b;
        }
    }
#endif

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
            return CA * (uint64_t) a + CB * b;
        } else if (a < 0 && a_abs < b) {
            return CA * b + CB * a_abs;
        } else if (a > 0) { /* implies 0 < a < b */
            if (a_abs % 2 == b % 2) {
                return CA * ((a_abs+b) >> 1) + CB * ((b-a_abs) >> 1);
            } else {
                return CA * (a_abs+b) + CB * (b-a_abs);
            }
        } else { /* a < -b < 0 */
            if (a_abs % 2 == b % 2) {
                return CA * ((b+a_abs) >> 1) + CB * ((a_abs-b) >> 1);
            } else {
                return CA * (b+a_abs) + CB * (a_abs-b);
            }
        }
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
            modul_add(rr, two, rr, pp);     /* 1+3/(r-1) */
            modul_neg(rr, rr, pp);          /* -2-3/(r-1) */
            unsigned long const sigma_r = modul_get_ul(rr, pp);

            modul_clear(two, pp);
            modul_clear(three, pp);
            modul_clear(rr, pp);
            modul_clearmod(pp);
            return sigma_r;
        }
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

         *  (a, b)
         *  -> (-a+b, -a-2*b) ~ (a-b, a+2*b)
         *  -> (-3*b, 3*a+3*b) ~ (-b, a+b) ~ (b, -a-b)
         *  -> (a+2*b, -2*a-b) ~ (-a-2*b, 2*a+b)
         *  -> (-3*a-3*b, 3*a) ~ (-a-b, a) ~ (a+b, -a)
         *  -> (2*a+b, -a+b) ~ (-2*a-b, a-b)
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
            return CA * (uint64_t) a + CB * b;
        } else if (a > 0) { /* 0 < a < b */
            if (a_abs % 3 == b % 3) {
              return CA * ((2*a_abs+b)/3) + CB * ((b - a_abs)/3);
            } else {
              return CA * (2*a_abs+b) + CB * (b - a_abs);
            }
        } else if (2*a_abs < b) { /* -b/2 < a < 0 */
            return CA * (b-a_abs) + CB * a_abs;
        } else if (a_abs < b) { /* -b < a < -b/2 < 0 */
            if ((UINT64_C(3) - (a_abs % 3)) % 3 == b % 3) {
                return CA * (((b << 1)-a_abs)/3) + CB * (((a_abs << 1)-b)/3);
            } else {
                return CA * ((b << 1)-a_abs) + CB * ((a_abs << 1)-b);
            }
        } else if (a_abs < 2*b) { /* -2b < a < -b < 0 */
            return CA * b + CB * (a_abs-b);
        } else { /* a < -2b < 0 */
            if ((UINT64_C(3) - (a_abs % 3)) % 3 == b % 3) {
                return CA * ((b+a_abs)/3) + CB * ((a_abs-(b << 1))/3);
            } else {
                return CA * (b+a_abs) + CB * (a_abs-(b << 1));
            }
        }
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
    : impl(nullptr)
{
    if (action == "none" || action.empty() || action == "id" || action == "Id")
        impl = static_cast<impl_ptr>(new galois_action_none());
    else if (action == "autom2.1" || action == "autom2.1g" || action == "1/y")
        impl = static_cast<impl_ptr>(new galois_action_inv());
    else if (action == "autom2.2" || action == "autom2.2g" || action == "_y")
        impl = static_cast<impl_ptr>(new galois_action_neg());
    else if (action == "autom3.1" || action == "autom3.1g")
        impl = static_cast<impl_ptr>(new galois_action_autom3_1());
    else if (action == "autom3.2" || action == "autom3.2g")
        impl = static_cast<impl_ptr>(new galois_action_autom3_2());
    else if (action == "autom4.1" || action == "autom4.1g")
        impl = static_cast<impl_ptr>(new galois_action_autom4_1());
    else if (action == "autom6.1" || action == "autom6.1g")
        impl = static_cast<impl_ptr>(new galois_action_autom6_1());
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

uint64_t galois_action::hash_ab(int64_t a, uint64_t b,
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

    sigma.resize(tab.get_size());

    for (index_t i = 0; i < tab.get_size(); ) {
        if (tab.is_additional_column(i)) {
            sigma[i] = i; /* extra columns are unchanged by galois action */
            i++;
        } else if (tab.is_bad(i)) {
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
            i++;
        } else {
            index_of_r.clear();
            renumber_t::p_r_side idc = tab.p_r_from_index(i);
            renumber_t::p_r_side const id0 = idc;
            do
            {
                index_of_r.emplace(idc.r, i);
                i++;
                if (i == tab.get_size())
                    break;
                idc = tab.p_r_from_index(i);
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

galois_action::~galois_action()
{
    delete impl;
}
