#include "cado.h" // IWYU pragma: keep

#include <memory>
#include <stdexcept>

#include <gmp.h>

#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "utils_cxx.hpp"

auto fmt::formatter<cxx_mpz>::format(cxx_mpz const & z, format_context& ctx) const
-> format_context::iterator
{
    auto out = ctx.out();
    using astring = std::unique_ptr<char[], free_delete<char>>;

    // see https://fmt.dev/dev/syntax/#:~:text=integer%20presentation%20types
    switch(specs_.type()) {
        case presentation_type::none:
        case presentation_type::dec:
            {
                const astring str(mpz_get_str(nullptr, 10, z.x));
                return format_to(out, "{}", str.get());
            }
        case presentation_type::hex:
            {
                const int b = specs_.upper() ? -16 : 16;
                const astring str(mpz_get_str(nullptr, b, z.x));
                if (specs_.alt()) {
                    const char * prefix = specs_.upper() ? "0X" : "0x";
                    if (str[0] == '-') {
                        return format_to(out, "-{}{}", prefix, str.get()+1);
                    } else {
                        return format_to(out, "{}{}", prefix, str.get());
                    }
                } else {
                    return format_to(out, "{}", str.get());
                }
            }
        case presentation_type::oct:
            {
                const int b = 8;
                const astring str(mpz_get_str(nullptr, b, z.x));
                const char * prefix = "0";
                if (z == 0)
                    return format_to(out, "0");

                if (str[0] == '-') {
                    return format_to(out, "-{}{}", prefix, str.get()+1);
                } else {
                    return format_to(out, "{}{}", prefix, str.get());
                }
            }
        case presentation_type::bin:
            {
                const int b = 2;
                const astring str(mpz_get_str(nullptr, b, z.x));
                if (specs_.alt()) {
                    const char * prefix = specs_.upper() ? "0B" : "0b";
                    if (str[0] == '-') {
                        return format_to(out, "-{}{}", prefix, str.get()+1);
                    } else {
                        return format_to(out, "{}{}", prefix, str.get());
                    }
                } else {
                    return format_to(out, "{}", str.get());
                }
            }
        case presentation_type::chr:
        default:
            throw std::runtime_error("not implemented");
            break;
    }
}
