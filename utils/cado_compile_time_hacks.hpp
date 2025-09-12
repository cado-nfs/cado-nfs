#ifndef UTILS_CADO_COMPILE_TIME_HACKS_HPP_
#define UTILS_CADO_COMPILE_TIME_HACKS_HPP_

#include <type_traits>
#include <limits>
#include <array>

namespace cado_math_aux {
    /* {{{ compile-time left shift */
    namespace details {
    template<typename T, int n>
        struct multiply_by_poweroftwo_impl {
            constexpr T operator()(T const & x) const {
                return multiply_by_poweroftwo_impl<T, n-1>()(2*x);
            }
        };

    template<typename T>
        struct multiply_by_poweroftwo_impl<T, 0> {
            constexpr T operator()(T const & x) const {
                return x;
            }
        };
    } /* namespace details */
    template<int n, typename T>
    static T constexpr multiply_by_poweroftwo(T const & x) {
        return details::multiply_by_poweroftwo_impl<T, n>()(x);
    }
    /* }}} */

    template<typename T>
        class constant_time_square_root {
            static constexpr T mid(T a, T b) { return (a+b)/2; }
            static constexpr T above(T m, T x) { return m*m > x; }
            static constexpr T recurse(T a, T b, T m, T x)
            { return above(m, x) ? sqrt(a, m, x) : sqrt(m, b, x); }
            static constexpr T sqrt(T a, T b, T m, T x)
            { return m == a ? a : recurse(a, b, m, x); }
            static constexpr T sqrt(T a, T b, T x)
            { return sqrt(a, b, mid(a, b), x); }
            public:
            static constexpr T
                sqrt(T x)
                requires std::is_integral_v<T>
            { return sqrt(T(0), (T(1) << (std::numeric_limits<T>::digits/2)), x); }
        };

    template<typename T>
        static constexpr T
        constant_sqrt(T x)
        requires std::is_integral_v<T>
        {
            return constant_time_square_root<T>::sqrt(x);
        }

    /* some fun with compile time code */
    /* {{{ compile time evaluation of 2^k mod n */
    template<int k, int n> struct pow2_mod {
        template<int x> struct sq_mod {
            static constexpr int value = (x * x) % n;
        };
        template<int x> struct dbl_mod {
            static constexpr int value = (x + x) % n;
        };
        struct even {
            static constexpr int value = sq_mod<pow2_mod<k/2, n>::value>::value;
        };
        struct odd {
            static constexpr int value = dbl_mod<even::value>::value;
        };
        static constexpr int value = (k % 2) ? odd::value : even::value;
    };
    template<int n> struct pow2_mod<0, n>
    {
        static constexpr int value = 1;
    };
    /* }}} */
    /* {{{ compile time bezout relation */
    template<int a, int b,
        int u0 = 1, int v0 = 0, int g0 = a,
        int u1 = 0, int v1 = 1, int g1 = b>
            struct bezout_relation;
    template<int a, int b,
        int u0, int v0, int g0,
        int u1, int v1, int g1>
            struct bezout_relation : public bezout_relation<
                                     a, b,
                                     u1, v1, g1, 
                                     u0 - (g0/g1) * u1, v0 - (g0 / g1) * v1, g0 % g1>
    {};
    template<int a, int b,
        int u0, int v0, int g0,
        int u1, int v1>
            struct bezout_relation<a, b, u0, v0, g0, u1, v1, 0> {
                static constexpr int gcd() { return g0; }
                static constexpr int u() { return u0; }
                static constexpr int v() { return v0; }
            };
    /* }}} */
    /* {{{ compile-time map on an fixed-size array */
    template<int N, template <int, int> class F, int left=N, int... Rest>
        struct array_map : public array_map<N, F, left - 1, F<N, left - 1>::value(), Rest...> {
        };

    template<int N, template <int, int> class F, int... Rest>
        struct array_map<N, F, 0, Rest...> {
            static constexpr
                std::array<int, N>
                value { Rest... };
        };
    /* }}} */

    /* {{{ compute -1/i mod n at compile time */
    template<int n, int k> struct minus_inverse_mod_n {
        using B = bezout_relation<n, n-k>;
        static constexpr int value() {
            return B::gcd() == 1 ? (B::v() < 0 ? n + B::v() : B::v()) : 0;
        }
    };
    /* }}} */
    /* {{{ inverse modulo powers of two */
    template<int k, typename T, T current = (3*k)^2, int c = 5, bool done = c >= std::numeric_limits<T>::digits>
        struct invmod;

    template<int k, typename T, T current, int c>
        struct invmod<k, T, current, c, false>
        : public invmod<k, T, current * (2 - current * k), c*2>
        {
            static_assert(k & 1, "k must be odd");
            static_assert(!std::is_signed_v<T>, "T must be an unsigned type");
        };
    template<int k, typename T, T current, int c>
        struct invmod<k, T, current, c, true> {
            static constexpr T value = current;
        };
    /* }}} */

    /* for any type that defines the <<= and += operators, provide code
     * that computes a constant multiple of a given value (of any
     * acceptable input type for +=).
     * Optionally, a chooser template (such as at_most<8>) can be passed
     * in order to control when to fall back on multiplication at
     * runtime. If this evaluates to false, then a *= operator on the
     * destination type is emitted.
     */

    namespace details {
        template<unsigned int n, bool is_odd = n & 1>
            struct mul_c_impl;

        template<unsigned int n>
            struct mul_c_impl<n, false> {
                using super = mul_c_impl<n/2>;
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        super::mul(t, a);
                        t <<= 1;
                    }
                static constexpr int number_of_shifts() {
                    return 1 + super::number_of_shifts();
                }
                static constexpr int number_of_additions() {
                    return super::number_of_shifts();
                }
            };
        template<unsigned int n>
            struct mul_c_impl<n, true> {
                using super = mul_c_impl<n/2>;
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        super::mul(t, a);
                        t <<= 1;
                        t += a;
                    }
                static constexpr int number_of_shifts() {
                    return 1 + super::number_of_shifts();
                }
                static constexpr int number_of_additions() {
                    return 1 + super::number_of_shifts();
                }
            };

        template<>
            struct mul_c_impl<1> {
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        t = a;
                    }
                static constexpr int number_of_shifts() {
                    return 0;
                }
                static constexpr int number_of_additions() {
                    return 0;
                }
            };

        template<>
            struct mul_c_impl<0> {
                template<typename T, typename U>
                    static void mul(T & t, U const &) {
                        t = 0;
                    }
                static constexpr int number_of_shifts() {
                    return 0;
                }
                static constexpr int number_of_additions() {
                    return 0;
                }
            };

        template<unsigned int n, bool is_small>
            struct mul_c_impl_choice
            : public mul_c_impl<n> {
            };

        template<unsigned int n>
            struct mul_c_impl_choice<n, false> {
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        t = a;
                        t *= n;
                    }
            };


        /* addmul is different because it doesn't use double-and-add.
         * It's not very useful, so I'm commenting it out.
         */
#if 0
        template<unsigned int n>
            struct addmul_c_impl {
                template<typename T, typename U>
                    static void addmul(T & t, U const & a) {
                        addmul_c_impl<n-1>::addmul(t, a);
                        t += a;
                    }
            };
        template<>
            struct addmul_c_impl<0> {
                template<typename T, typename U>
                    static void addmul(T &, U const &) {
                    }
            };

        template<unsigned int n, bool is_very_small, typename chooser_mul>
            struct addmul_c_impl_choice
            : public addmul_c_impl<n> {};

        template<unsigned int n, typename chooser_mul>
            struct addmul_c_impl_choice<n, false, chooser_mul> {
                template<typename T, typename U>
                    static void addmul(T & t, U const & a) {
                        T t2;
                        mul_c<n, chooser_mul>(t2, a);
                        t += t2;
                    }
            };
#endif
    }   /* namespace details */

    /* Three chooser types. Evaluation to true means "use an addition
     * chain".
     */
    template<unsigned int max>
        struct at_most {
            template<unsigned int n> struct choose : public
                                                     std::integral_constant<bool, (n <= max)> {
                                                     };
        };

    struct always {
        template<unsigned int n> struct choose : public
                                                 std::true_type {
                                                 };
    };

    struct never {
        template<unsigned int n> struct choose : public
                                                 std::false_type {
                                                 };
    };

    template<unsigned int n, typename chooser = at_most<8>, typename T, typename U>
        void mul_c(T & t, U const & a, chooser const & = chooser{})
        {
            details::mul_c_impl_choice<n, chooser::template choose<n>::value>::mul(t, a);
        }
#if 0
    /* see above. the code is fine, but I'm yet to be convinced that it's
     * useful
     */
    template<unsigned int n, typename chooser = at_most<2>, typename chooser_mul, typename T, typename U>
        void addmul_c(T & t, U const & a, chooser const & = {}, chooser_mul const & = {})
        {
            details::addmul_c_impl_choice<n, chooser::template choose<n>::value, chooser_mul>::addmul(t, a);
        }
#endif

} /* namespace cado_math_aux */

#endif	/* UTILS_CADO_COMPILE_TIME_HACKS_HPP_ */
