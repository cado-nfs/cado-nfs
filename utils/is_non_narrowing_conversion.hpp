#ifndef CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP
#define CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP

#include <type_traits>

namespace cado_math_aux {
    /* https://stackoverflow.com/questions/36270158/avoiding-narrowing-conversions-with-c-type-traits
     */
    namespace is_narrowing_conversion_detail
    {
        template<typename...> struct void_type { typedef void type; };

        template<typename From, typename To, typename = void>
            struct is_narrowing_conversion_impl : std::true_type {};
        template<typename From, typename To, typename = void>
            struct is_non_narrowing_conversion_impl : std::false_type {};

        template<typename From, typename To>
            struct is_narrowing_conversion_impl<From, To, typename void_type<decltype(To{std::declval<From>()})>::type> : std::false_type {};
        template<typename From, typename To>
            struct is_non_narrowing_conversion_impl<From, To, typename void_type<decltype(To{std::declval<From>()})>::type> : std::true_type {};
    }  // namespace detail

    template<typename From, typename To>
        struct is_narrowing_conversion : is_narrowing_conversion_detail::is_narrowing_conversion_impl<From, To> {};
    template<typename From, typename To>
        struct is_non_narrowing_conversion : is_narrowing_conversion_detail::is_non_narrowing_conversion_impl<From, To> {};
}

#endif	/* CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP_ */
