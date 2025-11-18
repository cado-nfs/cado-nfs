#include "cado.h"       // IWYU pragma: keep

#include <ostream>

#include "fmt/base.h"
#include <mpfr.h>
#include <mpc.h>

#include "cxx_mpfr.hpp"
#include "cxx_mpc.hpp"

std::ostream & operator<<(std::ostream & os, cxx_mpc const & x)
{
    return os << "(" << mpc_realref(x) << "," << mpc_imagref(x) << ")";
}

auto fmt::formatter<cxx_mpc>::format(cxx_mpc const & x, format_context& ctx) const -> format_context::iterator
{
    /* Is this optimization desirable or not? */
    if (mpfr_zero_p(mpc_imagref(x.x)))
        return formatter<cxx_mpfr>::format(mpc_realref(x.x), ctx);

    /* From here on, we follow the same convention as fmt::format for
     * std::complex.  */
    if (!mpfr_zero_p(mpc_realref(x.x))) {
        format_to(ctx.out(), "(");
        formatter<cxx_mpfr>::format(mpc_realref(x.x), ctx);
        /* an mpfr with the sign bit set will (or should) print with a
         * starting minus */
        if (!mpfr_signbit(mpc_imagref(x.x)))
            format_to(ctx.out(), "+");
    }

    formatter<cxx_mpfr>::format(mpc_imagref(x.x), ctx);
    if (mpfr_number_p(mpc_imagref(x.x))) {
        format_to(ctx.out(), "i");
    } else {
        format_to(ctx.out(), " i");
    }

    if (!mpfr_zero_p(mpc_realref(x.x)))
        format_to(ctx.out(), ")");

    return ctx.out();
}

