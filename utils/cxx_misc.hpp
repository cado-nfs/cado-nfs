#ifndef CXX_MISC_HPP_
#define CXX_MISC_HPP_

#include <limits>
#include <type_traits>

/* miscellaneous stuff for C++ */

/* This is handy for range-based for loops that also need a counter.
 */
template<typename T>
struct increment_counter_on_dtor {
    T & a;
    increment_counter_on_dtor(T & a) : a(a) {}
    ~increment_counter_on_dtor() { ++a; }
};

/* usage example:
 *
    size_t counter = 0;
    for(auto x : foo) {
        increment_counter_on_dtor<size_t> _dummy(counter);
        if (x % 2 == 0) continue;
        std::cout << "[" << counter << "]: " << x << "\n";
    }

 * of course, &x-begin(foo) is also acceptable, but that works only for
 * random-access iterators.
 */


#if __cplusplus < 201402L
namespace std {
template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
}
#endif

/* A type trait that checks whether an integral type T can be cast losslessly
   to an integral type U.

   In particular, both T and U must be integral types, must have the same
   signedness, and the maximal permissible value of type T must be no greater
   than that of U. (We assume value ranges of signed types to be essentially
   symmetric around 0).

   Example use:

   template <typename T, integral_fits_t<T, long> = 0 >
   void print(T v) {printf("%ld\n", (long) v);}
   template <typename T, integral_fits_t<T, unsigned long> = 0 >
   void print(T v) {printf("%lu\n", (unsigned long) v);}
*/

/* Note: with gcc 9.2.1, a debug build can't instantiate the full check
 * "both integral + same sign + compatible maxval" lazily, and therefore
 * we get a warning. So we have to resort to an ugly workaround.
 */

template <typename T, typename U>
struct integral_fits_pre_ {
    static constexpr bool value = std::is_integral<T>::value && std::is_integral<U>::value &&
                                  std::is_signed<T>::value == std::is_signed<U>::value;
};

template<bool pre_flag, typename T, typename U>
struct integral_fits_post;

template<typename T, typename U>
struct integral_fits_post<true, T, U> {
    static constexpr bool value = std::numeric_limits<T>::max() <= std::numeric_limits<U>::max();
};
template<typename T, typename U>
struct integral_fits_post<false, T, U> {
    static constexpr bool value = false;
};

template <typename T, typename U>
struct integral_fits_ {
    static constexpr bool value_pre = integral_fits_pre_<T, U>::value;
    static constexpr bool value = integral_fits_post<value_pre, T, U>::value;
};


template<bool> struct integral_fits_final : std::false_type {};
template<> struct integral_fits_final<true> : std::true_type { typedef bool type; };

template <typename T, typename U>
struct integral_fits : integral_fits_final<integral_fits_<T, U>::value> {};

template <typename T, typename U >
using integral_fits_t = typename integral_fits<T, U>::type;

#endif	/* CXX_MISC_HPP_ */
