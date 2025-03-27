#ifndef CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP_
#define CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP_

#include <type_traits>
#include <limits>

namespace cado_math_aux {
    namespace is_non_narrowing_integral_conversion_details {
        template<bool b, typename From, typename To>
        struct checks;

        template<typename T, typename U,
            bool b = std::is_integral<T>::value && std::is_integral<U>::value>
        struct impl {
            static constexpr bool value = checks<b, T, U>::value;
        };

        /* We don't want this to be instantiated at all if From or To are not
         * integral types
         */
        template<bool b, typename From, typename To>
        struct checks
        {
            static constexpr bool value =
            std::is_signed<From>::value == std::is_signed<To>::value &&
            std::is_convertible<From, To>::value &&
            std::numeric_limits<From>::min() >= std::numeric_limits<To>::min() &&
            std::numeric_limits<From>::max() <= std::numeric_limits<To>::max();
        };

        template<typename From, typename To>
        struct checks<false, From, To>
        {
            static constexpr bool value = false;
        };
    }

#if 1
    template<typename From, typename To>
    struct is_non_narrowing_integral_conversion
    : public is_non_narrowing_integral_conversion_details::impl<From, To>
    {}
    ;
#else
    template<typename From, typename To>
    struct is_non_narrowing_integral_conversion
    {
        static constexpr bool value =
            std::is_integral<From>::value &&
            std::is_integral<To>::value &&
            std::is_signed<From>::value == std::is_signed<To>::value;

    }
    ;
#endif
}

#endif	/* CADO_UTILS_IS_NON_NARROWING_INTEGRAL_CONVERSION_HPP_ */
