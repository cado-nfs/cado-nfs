#ifndef CADO_RUNTIME_NUMERIC_CAST_HPP
#define CADO_RUNTIME_NUMERIC_CAST_HPP

#include <type_traits>
#include <limits>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include "fmt/format.h"

namespace runtime_numeric_cast_details {
struct runtime_numeric_cast_failure : public std::runtime_error {
    explicit runtime_numeric_cast_failure(std::string const & s)
        : std::runtime_error(s)
    {}
};
struct runtime_numeric_cast_underflow : public runtime_numeric_cast_failure {
    explicit runtime_numeric_cast_underflow(std::string const & s)
    : runtime_numeric_cast_failure(s) {}
};
struct runtime_numeric_cast_overflow : public runtime_numeric_cast_failure {
    explicit runtime_numeric_cast_overflow(std::string const & s)
    : runtime_numeric_cast_failure(s) {}
};

template<typename TO_TYPE>
struct runtime_numeric_cast_impl {
    struct cast_underflow : public runtime_numeric_cast_underflow {
        template<typename FROM_TYPE>
        explicit cast_underflow(FROM_TYPE x)
        : runtime_numeric_cast_underflow(fmt::format("underflow while casting {} of type {} to type {}", x, typeid(FROM_TYPE).name(), typeid(TO_TYPE).name())) {}
    };
    struct cast_overflow : public runtime_numeric_cast_overflow {
        template<typename FROM_TYPE>
        explicit cast_overflow(FROM_TYPE x)
        : runtime_numeric_cast_overflow(fmt::format("overflow while casting {} of type {} to type {}", x, typeid(FROM_TYPE).name(), typeid(TO_TYPE).name())) {}
    };
    static_assert(std::is_integral<TO_TYPE>::value
            /* we can probably get along with floating point destination
             * types, too. However it's not totally trivial, given that
             * we rely on the max() valued of a narrowed-to type target
             * to be cast-able to the the source type. Which isn't the
             * case for long->double
             */
            // || std::is_floating_point<TO_TYPE>::value
            , "destination type TO_TYPE must be integral");
    template<typename FROM_TYPE> struct traits {
    static_assert(std::is_integral<FROM_TYPE>::value,
            "source type FROM_TYPE must be integral");
    static constexpr bool T_signed = std::is_signed<FROM_TYPE>::value;
    static constexpr bool U_signed = std::is_signed<TO_TYPE>::value;
    static constexpr int T_digits = std::numeric_limits<FROM_TYPE>::digits;
    static constexpr int U_digits = std::numeric_limits<TO_TYPE>::digits;
    static constexpr bool same_sign = T_signed == U_signed;
    static constexpr bool both_signed = T_signed && U_signed;
    static constexpr bool both_unsigned = !T_signed && !U_signed;
    static constexpr bool widening = U_digits >= T_digits;
    static constexpr bool same_size = U_digits == T_digits;
    static constexpr bool strict_widening = widening && !same_size;
    static constexpr bool narrowing = U_digits < T_digits;
    static constexpr bool to_signed = U_signed && ! T_signed;
    static constexpr bool to_unsigned = T_signed && ! U_signed;
    static_assert(!(same_sign && widening) || (std::numeric_limits<FROM_TYPE>::max() <= std::numeric_limits<TO_TYPE>::max()), "same-sign, widening cast expects compatible max values");
    static_assert(!(to_unsigned && widening) || (std::numeric_limits<FROM_TYPE>::max() <= std::numeric_limits<TO_TYPE>::max()), "to-unsigned, widening cast expects compatible max values");
    static_assert(!(to_signed && strict_widening) || (std::numeric_limits<FROM_TYPE>::max() <= std::numeric_limits<TO_TYPE>::max()), "to-signed, strict widening cast expects compatible max values");
    static_assert(!(same_sign && widening) || (std::numeric_limits<FROM_TYPE>::min() >= std::numeric_limits<TO_TYPE>::min()), "same-sign, widening cast expects compatible min values");
    };
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::same_sign && traits<FROM_TYPE>::widening>::type * = nullptr) { return FROM_TYPE(x); }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::to_signed && traits<FROM_TYPE>::strict_widening>::type * = nullptr) { return FROM_TYPE(x); }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::both_signed && traits<FROM_TYPE>::narrowing>::type * = nullptr) {
        if (x < FROM_TYPE(std::numeric_limits<TO_TYPE>::min()))
            throw cast_underflow(x);
        if (x > FROM_TYPE(std::numeric_limits<TO_TYPE>::max()))
            throw cast_overflow(x);
        return FROM_TYPE(x);
    }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::both_unsigned && traits<FROM_TYPE>::narrowing>::type * = nullptr) {
        if (x > FROM_TYPE { std::numeric_limits<TO_TYPE>::max() })
            throw cast_overflow(x);
        return FROM_TYPE(x);
    }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::to_signed && traits<FROM_TYPE>::narrowing>::type * = nullptr) {
        if (x > FROM_TYPE { std::numeric_limits<TO_TYPE>::max() })
            throw cast_overflow(x);
        return FROM_TYPE(x);
    }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::to_signed && traits<FROM_TYPE>::same_size>::type * = nullptr) {
        if (x > FROM_TYPE { std::numeric_limits<TO_TYPE>::max() })
            throw cast_overflow(x);
        return FROM_TYPE(x);
    }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::to_unsigned && traits<FROM_TYPE>::narrowing>::type * = nullptr) {
        if (x < FROM_TYPE { 0 })
            throw cast_underflow(x);
        if (x > FROM_TYPE { std::numeric_limits<TO_TYPE>::max() })
            throw cast_overflow(x);
        return FROM_TYPE(x);
    }
    template<typename FROM_TYPE>
    TO_TYPE operator()(FROM_TYPE x, typename std::enable_if<traits<FROM_TYPE>::to_unsigned && traits<FROM_TYPE>::widening>::type * = nullptr) {
        if (x < FROM_TYPE { 0 })
            throw cast_underflow(x);
        return FROM_TYPE(x);
    }
};
}

template<typename TO_TYPE>
struct runtime_numeric_cast {
    // NOLINTBEGIN(hicpp-explicit-conversions)
    TO_TYPE x;
    template<typename FROM_TYPE>
        runtime_numeric_cast(FROM_TYPE x)
            : x(runtime_numeric_cast_details::runtime_numeric_cast_impl<TO_TYPE>()(x))
        {}
    operator TO_TYPE() { return x; }
    // NOLINTEND(hicpp-explicit-conversions)
};



#endif	/* CADO_RUNTIME_NUMERIC_CAST_HPP */
