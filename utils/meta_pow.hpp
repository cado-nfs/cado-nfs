#ifndef UTILS_META_POW_HPP_
#define UTILS_META_POW_HPP_

/* This is a generic adapter that provides a pow function for a type T,
 * given a function object that performs the reduction that is to be
 * applied after each multiply.
 *
 * We assume the following.
 *  - type T is constructible from the constant 1
 *  - type T is copy-constructible
 *  - comparison of T to zero is possible.
 *  - for a reducer r of type R, r(T&, T&&) is defined.
 *  - T::mul(T &, T const &, T const &) is defined as a static function.
 *  This last requirement should probably be changed to having a global
 *  mul() overload in the global namespace.
 */

#include <gmp.h>

#include "macros.h"

namespace cado {

template <typename T, typename R> struct meta_pow { /*{{{*/
    R const & reducer;
    explicit meta_pow(R const & reducer)
        : reducer(reducer)
    {
    }

    private:
    void wrapped(T & B, T const & A, unsigned long n) const
    {
        ASSERT_ALWAYS(&B != &A);
        ASSERT_ALWAYS(n > 0);
        unsigned long k = ((~0UL) >> 1) + 1;
        for (; k > n; k >>= 1)
            ;
        B = A;
        for (; k >>= 1;) {
            T::mul(B, B, B);
            reducer(B, B);
            if (n & k) {
                T::mul(B, B, A);
                reducer(B, B);
            }
        }
    }
    void wrapped(T & B, T const & A, mpz_srcptr n) const
    {
        ASSERT_ALWAYS(&B != &A);
        ASSERT_ALWAYS(mpz_sgn(n) > 0);
        unsigned long k = mpz_sizeinbase(n, 2) - 1;
        B = A;
        for (; k--;) {
            T::mul(B, B, B);
            reducer(B, B);
            if (mpz_tstbit(n, k)) {
                T::mul(B, B, A);
                reducer(B, B);
            }
        }
    }

    public:
    void operator()(T & B, T const & A, unsigned long n) const
    {
        if (n == 0) {
            B = 1;
            return;
        }
        if (A == 0) {
            B = 0;
            return;
        }
        if (&B == &A) {
            wrapped(B, T(A), n);
        } else {
            wrapped(B, A, n);
        }
    }
    void operator()(T & B, T const & A, mpz_srcptr n) const
    {
        if (mpz_sgn(n) == 0) {
            B = 1;
            return;
        }
        if (A.degree() < 0) {
            B = 0;
            return;
        }
        if (&B == &A) {
            wrapped(B, T(A), n);
        } else {
            wrapped(B, A, n);
        }
    }
}; /*}}}*/

} /* namespace cado */

#endif	/* UTILS_META_POW_HPP_ */
