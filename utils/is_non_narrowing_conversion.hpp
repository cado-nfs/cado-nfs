#ifndef CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP
#define CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP

#include <type_traits>

/* XXX how does it compare to integral_fits defined in utils_cxx.hpp ? */

namespace cado_math_aux {
    /* https://stackoverflow.com/questions/36270158/avoiding-narrowing-conversions-with-c-type-traits
     */
    namespace is_narrowing_conversion_detail
    {
        template<typename From, typename To, typename = void>
            struct is_narrowing_conversion_impl : std::true_type {};
        template<typename From, typename To, typename = void>
            struct is_non_narrowing_conversion_impl : std::false_type {};

        template<typename From, typename To>
            struct is_narrowing_conversion_impl<From, To, std::void_t<decltype(To{std::declval<From>()})>> : std::false_type {};
        template<typename From, typename To>
            struct is_non_narrowing_conversion_impl<From, To, std::void_t<decltype(To{std::declval<From>()})>> : std::true_type {};
    }  /* namespace is_narrowing_conversion_detail */

    template<typename From, typename To>
        struct is_narrowing_conversion : is_narrowing_conversion_detail::is_narrowing_conversion_impl<From, To> {};
    template<typename From, typename To>
        struct is_non_narrowing_conversion : is_narrowing_conversion_detail::is_non_narrowing_conversion_impl<From, To> {};

    template<typename From, typename To>
    inline constexpr bool is_narrowing_conversion_v = is_narrowing_conversion<From, To>::value;
    template<typename From, typename To>
    inline constexpr bool is_non_narrowing_conversion_v = is_non_narrowing_conversion<From, To>::value;
} /* namespace cado_math_aux */

#endif	/* CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP_ */
