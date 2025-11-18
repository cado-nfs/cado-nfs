#ifndef CADO_UTILS_ARITHXX_API_HPP
#define CADO_UTILS_ARITHXX_API_HPP

#include <cstddef>
#include <cstdint>

#include <array>
#include <memory>
#include <new>

#include <utility>

#include "utils_cxx.hpp"
#include "arithxx_residue_std_op.hpp"
#include "macros.h"
#include "misc.h"

namespace arithxx_details {
    using std::index_sequence;
    using std::make_index_sequence;

    /* The functions here are defined for every instantiation of the api.
     * This type is always a base class of the Modulus class, so that
     * downcasting via static_cast is OK (this is what CRTP is all
     * about).
     */
    template <typename layer>
        struct api {
            using Modulus = typename layer::Modulus;
            using Residue = typename layer::Residue;
            using Integer = typename layer::Integer;
            
            /* Data members */
            Integer m;

            Integer const & getmod() const { return m; }

            explicit api(Integer const & m)
                : m(m)
            {}

            // we can't do this static assert yet (because Modulus is
            // incomplete), even though we'd like to.
            // static_assert(std::is_base_of<api<layer>, Modulus>::value, "CRTP condition violated");

            Modulus & downcast() { return static_cast<Modulus&>(*this); }
            Modulus const & downcast() const { return static_cast<Modulus const &>(*this); }

            using ResidueOp = arithxx_details::ResidueStdOp<layer>;

            std::unique_ptr<Residue[]>
                newArray(Modulus const * mm, size_t const len)
                {
                    void * t = operator new[](len * sizeof(Residue));
                    if (!t)
                        throw std::bad_alloc();
                    auto * ptr = static_cast<Residue *>(t);
                    for (size_t i = 0; i < len; i++) {
                        new (&ptr[i]) Residue(*mm);
                    }
                    return std::unique_ptr<Residue[]>(ptr);
                }

            private:
            template <size_t... Is>
                static std::array<Residue, sizeof...(Is)>
                make_array_impl(Modulus const & m, index_sequence<Is...>)
                {
                    return {((void)Is, Residue(m))...};
                }

            public:
            template <size_t N>
                std::array<Residue, N> make_array() const
                {
                    return make_array_impl(downcast(), make_index_sequence<N>());
                }

            Residue operator()(Integer const& a) const
            {
                Residue r(downcast());
                downcast().set(r, a);
                return r;
            }
            template <typename T>
                Residue operator()(T const a) const
                requires cado::converts_via<T, uint64_t>
                {
                    Residue r(downcast());
                    downcast().set(r, uint64_t(a));
                    return r;
                }
            template <typename T>
                Residue operator()(T const a) const
                requires cado::converts_via<T, int64_t>
                {
                    Residue r(downcast());
                    downcast().set(r, int64_t(a));
                    return r;
                }

            bool is_strong_pseudoprime_base2() const;
            bool is_strong_lucas_pseudoprime() const;
            bool is_prime() const;
            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            constexpr bool sprp2_is_enough() const { return false; }

            bool batchinv(Residue * r, Residue const * a, size_t n, Residue const * c) const;

            void gcd (Integer &, const Residue &) const;
            int jacobi(Residue const &) const;


            bool inv(Residue &, Residue const &) const;
            bool intinv(Integer & r, Integer const & a) const;
            bool inv_odd(Residue & r, Residue const & a) const;
            bool inv_powerof2(Residue & r, Residue const & a) const;
#if 0
            /* 
             * This is an interesting alternative. Currently, a function
             * of the api that isn't implemented produces a linking
             * failure, which is hard to debug. By adding a static assert
             * like here, we force a compiler error earlier on. The
             * consequence is that we _must_ override at the terminal
             * class level, and we can't do special instantiations (since
             * that would lead to multiply defined methods). It is
             * generally not a problem, and to this point both approaches
             * have been an option. However, losing this flexibility
             * becomes annoying when we deal with overloads.
             *
             * One thing that might not work, though: forcing
             * instantiation of api_bysize<blah> does not seem to trigger
             * instantiation of api<blah>: if I remove the explicit
             * instantiation at the end of mod64.cpp I get linking
             * errors.
             * Very weird.
             */

            {
                static_assert(false, "must be overridden");
            }
#endif

            void pow(Residue &, Residue const &, uint64_t const *, size_t) const;
            void pow2(Residue &r, const uint64_t *e, size_t e_nrwords) const;

            void pow(Residue & r, Residue const & b, uint64_t e) const;
            void pow2(Residue &r, uint64_t e) const;

            void pow(Residue & r, Residue const & b, Integer const & e) const;
            void pow2(Residue &r, const Integer &e) const;

            /* {{{ V_dadd and V_dbl for Lucas sequences.
             *
             * No implementation has specialized Lucas sequences code, so
             * it's fine to define implementations right here. We compute
             * V_n = V_n(b, 2), as defined by V_0=2, V_1=b, V_2=b^2-2
             * (or, if b=x+y with xy=1, V_n=x^n+y^n).
             */

            /* Given a = V_n (x), b = V_m (x) and d = V_{n-m} (x), the
             * differential addition V_dadd computes V_{m+n} (x).  r can
             * be the same variable as a or b but must not be the same
             * variable as d.
             */
            void V_dadd(Residue & r, Residue const & a, Residue const & b,
                    Residue const & d) const
            {
                auto const & me = downcast();
                ASSERT(&r != &d);
                me.mul(r, a, b);
                me.sub(r, r, d);
            }

            /* Given a = V_n (x) and two = 2, compute V_{2n} (x).  r can
             * be the same variable as a but must not be the same
             * variable as two.
             */
            void V_dbl(Residue & r, Residue const & a, Residue const & two) const
            {
                auto const & me = downcast();
                ASSERT(&r != &two);
                me.sqr(r, a);
                me.sub(r, r, two);
            }
            /* }}} */

            void V(Residue & r, Residue * rp1, Residue const & b, Integer const & k) const;

            private:
            template<int n>
                bool divn(Residue &, Residue const &) const;

            public:
            bool div3(Residue &, Residue const &) const;
            bool div5(Residue &, Residue const &) const;
            bool div7(Residue &, Residue const &) const;
            bool div11(Residue &, Residue const &) const;
            bool div13(Residue &, Residue const &) const;

            /* {{{ set(*2), set_reduced(*1), set0 */
            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            void set(Residue & r, Residue const & s) const { r = s; }
            void set(Residue & r, int64_t const s) const
            {
                auto const & me = downcast();
                me.set(r, safe_abs64(s));
                if (s < 0)
                    me.neg(r, r);
            }

            /* Sets the residueredc2ul2_t to the class represented by the integer s.
               Assumes that s is reduced (mod m), i.e. 0 <= s < m */
            void set_reduced(Residue & r, uint64_t const s) const {
                auto const & me = downcast();
                me.set(r, s);
            }

            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            void set0(Residue & r) const { r.r = 0; }

            /* }}} */

            /* {{{ equal is0 */

            /* do we really want to keep these two, or should we use operator== ?
             * comparison to 1 in montgomery form is tricky, though.
             */
            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            bool equal(Residue const & a, Residue const & b) const
            {
                return a.r == b.r;
            }

            // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
            bool is0(Residue const & a) const { return a.r == 0; }

            /* }}} */

            /* {{{ set(*2), set_reduced(*1), set1, get, is1: possibly overloaded by redc */
            void set(Residue & r, uint64_t const s) const { r.r = s; }
            void set(Residue & r, Integer const & s) const {
                auto const & me = downcast();
                if (s < me.m)
                    me.set_reduced(r, s);
                else
                    me.set_reduced(r, s % me.m);
            }
            void set_reduced(Residue & r, Integer const & s) const { r.r = s; }
            void set1(Residue & r) const { r.r = (m != 1); }
            Integer get(Residue const & r) const { return r.r; }
            bool is1(Residue const& a) const { return a.r == 1; }
            /* }}} */
            void neg(Residue & r, Residue const & a) const
            {
                auto const & me = downcast();
                me.assertValid(a);
                if (me.is0(a))
                    me.set0(r);
                else
                    r.r = m - a.r;
            }
        };


    template <typename layer, typename Integer = typename layer::Integer>
        struct api_bysize;
} /* namespace arithxx_details */


#endif	/* UTILS_ARITHXX_API_HPP_ */
