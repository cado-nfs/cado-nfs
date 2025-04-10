#ifndef UTILS_ARITHXX_API_HPP_
#define UTILS_ARITHXX_API_HPP_

#include <cstddef>
#include <cstdint>

#include <array>
#include <memory>
#include <new>

#if __cplusplus >= 201402L
#include <utility>
#endif

#include "arithxx_residue_std_op.hpp"

namespace arithxx_details {

#if __cplusplus < 201402L
    /* we need to import equivalents of c++14's std::index_sequence and
     * std::make_index_sequence with only c++11. It's just clutter,
     * really.
     */
    template<std::size_t... Is> struct index_sequence{};
    template<std::size_t N, std::size_t... Is>
        struct make_index_sequence : make_index_sequence<N-1, N-1, Is...>{};
    template<std::size_t... Is>
        struct make_index_sequence<0, Is...> : index_sequence<Is...>{};
#else
    using std::index_sequence;
    using std::make_index_sequence;
#endif

    /* The functions here are defined for every instantiation of the api.
     * This type is always a base class of the Modulus class, so that
     * downcasting via static_cast is OK (this is what CRTP is all
     * about).
     */
    template <typename layer>
        struct api {
            typedef typename layer::Modulus Modulus;
            typedef typename layer::Residue Residue;
            typedef typename layer::Integer Integer;
            
            // we can't do this static assert yet (because Modulus is
            // incomplete), even though we'd like to.
            // static_assert(std::is_base_of<api<layer>, Modulus>::value, "CRTP condition violated");

            Modulus & downcast() { return static_cast<Modulus&>(*this); }
            Modulus const & downcast() const { return static_cast<Modulus const &>(*this); }

            typedef arithxx_details::ResidueStdOp<layer> ResidueOp;

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

            void pow(Residue &, Residue const &, uint64_t) const;
            void pow(Residue &, Residue const & b, Integer const & e) const;
            void pow(Residue &, Residue const &, uint64_t const *, size_t) const;
            bool sprp(Residue const &) const;
            bool sprp2() const;
            bool isprime() const;

            bool div3(Residue &, Residue const &) const;
            bool div5(Residue &, Residue const &) const;
            bool div7(Residue &, Residue const &) const;
            bool div11(Residue &, Residue const &) const;
            bool div13(Residue &, Residue const &) const;

            void V(Residue & r, Residue * rp1, Residue const & b, uint64_t k) const;
            bool batchinv(Residue * r, Residue const * a, size_t n, Residue const * c) const;

            void gcd (Integer &, const Residue &) const;
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
             * instantiation of api64<blah> does not seem to trigger
             * instantiation of api<blah>: if I remove the explicit
             * instantiation at the end of mod64.cpp I get linking
             * errors.
             * Very weird.
             */

            {
                static_assert(false, "must be overridden");
            }
#endif

            void pow2 (Residue &r, uint64_t e) const;
            void pow2 (Residue &r, const uint64_t *e, size_t e_nrwords) const;
            void pow2 (Residue &r, const Integer &e) const;
            void pow3 (Residue &r, uint64_t e) const;

            protected:
            bool find_minus1 (Residue &r1, const Residue &minusone, int po2) const;
        };
}

#endif	/* UTILS_ARITHXX_API_HPP_ */
