#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cmath>
#include <cstring>

#include <ios>
#include <istream>
#include <ostream>
#include <string>
#include <memory>

#include "fmt/base.h"
#include "fmt/format.h"
#include <mpfr.h>

#include "cado_expression_parser.hpp"
#include "cxx_mpfr.hpp"
#include "runtime_numeric_cast.hpp"
#include "macros.h"

/* code here _tries_ to be consistent with libc++ / libfmt in the
 * following way, to the extent that is possible.
 *  - if a double (d) prints as (td) in some given context (a certain
 *  choice of width, precision, presentation type, whether it's libc++ or
 *  libfmt, etc), then if we create a cxx_mpfr from (d), printing should be
 *  the same if the printing context is unchanged.
 *  - in the same context, if we create a cxx_mpfr from the printed
 *  string (td), then printing this cxx_mpfr should look the same as (d).
 *  In some occasions, this can only be achieved by printing with limited
 *  precision (corresponding to the precision of double).
 *
 * Hex printing with mpfr is known to be inconsistent, there isn't much
 * we plan to do about it.
 *
 * XXX A caveat is that we don't have a convenient way to print a
 * cxx_mpfr with the maximum precision without specifying the precision
 * exactly, which is a bit sad.
 */


/* It's surprisingly difficult */
std::ostream & operator<<(std::ostream & os, cxx_mpfr const & x)
{
    char * s = nullptr;

    const std::string fill(1, os.fill());
    const auto align_left = (os.flags() & std::ios_base::left) != 0;
    const auto w = os.width();
    const auto p = os.precision();
    const auto pt = os.flags() & std::ios_base::floatfield;
    constexpr decltype(pt) feag[] = {
        std::ios_base::fixed,
        std::ios_base::scientific,
        std::ios_base::fixed | std::ios_base::scientific,
        {},
    };

    std::string fs = "%";
    // #0-+[space]
    if (os.flags() & std::ios_base::showpoint)
        fs += "#";
    if (os.flags() & std::ios_base::showpos)
        fs += "+";

    if (fill == "0")
        fs += "0";
    if (align_left)
        fs += "-";

    if (w > 0)
        fs += fmt::format("{}", w);
    if (p >= 0)
        fs += fmt::format(".{}", p);

    if (pt == feag[0]) {
        /* as with std::fixed */
        fs += "Rf";
    } else if (pt == feag[1]) {
        /* as with std::scientific */
        fs += "Re";
    } else if (pt == feag[2]) {
        /* as with std::hexfloat */
        fs += "Ra";
    } else {
        /* as with std::defaultfloat */
        fs += "Rg";
    }
    mpfr_asprintf(&s, fs.c_str(), x.x);
    os << s;
    free(s);
    return os;
}


auto fmt::formatter<cxx_mpfr>::format(cxx_mpfr const & x, format_context& ctx) const
-> format_context::iterator {
    char * s = nullptr;
    std::string fs = "%";
    /*
    if (os.flags() & std::ios_base::showpoint)
        fs += "#";
    if (os.flags() & std::ios_base::showpos)
        fs += "+";
        */
    auto specs = specs_;

    const std::string fill(specs.fill<char>());
    const bool align_left = specs.align() == align::left;
    const auto w = specs.width;
    auto p = specs.precision;
    auto pt = specs.type();
    constexpr decltype(pt) feag[] = {
        presentation_type::fixed,
        presentation_type::exp,
        presentation_type::hexfloat,
        presentation_type::general };

    auto pt2 = pt;

    if (p < 0) {
        if (pt != presentation_type::none) {
            p = 6;
        }
    }
    if (pt == presentation_type::none) {
        if (mpfr_regular_p(x.x)) {
            // p = 6;
            auto str = std::make_unique<char[]>(mpfr_get_str_ndigits(10, x.prec()));
            mpfr_exp_t e;
            mpfr_get_str(str.get(), &e, 10, 0, x.x, MPFR_RNDN);
            /* the implicit decimal point is to the _left_ of the first
             * digit */
            e--;
            if (e < -4 || e >= (p > 0 ? p : 16)) {
                pt2 = presentation_type::exp;
            } else {
                pt2 = presentation_type::fixed;
                if (p < 0)
                    p = runtime_numeric_cast<int>(strlen(str.get())) - 1 - (str[0] == '-');
            }
        } else {
            /* doesn't realy matter which one we choose. It's only about
             * printing 0
             */
            pt2 = presentation_type::fixed;
        }
    }

    if (fill == "0")
        fs += "0";
    if (align_left)
        fs += "-";

    if (w > 0)
        fs += fmt::format("{}", w);
    if (p >= 0)
        fs += fmt::format(".{}", p);

    if (pt == feag[0]) {
        /* as with std::fixed */
        fs += "Rf";
    } else if (pt == feag[1]) {
        /* as with std::scientific */
        fs += "Re";
    } else if (pt == feag[2]) {
        /* as with std::hexfloat */
        fs += "Ra";
    } else if (pt == feag[3]) {
        /* presentation_type::general */
        fs += "Rg";
    } else {
        ASSERT_ALWAYS(pt == presentation_type::none);
        if (pt2 == presentation_type::exp)
            fs += "Re";
        else if (pt2 == presentation_type::fixed)
            fs += "Rf";
        else
            ASSERT_ALWAYS(0);
    }
    mpfr_asprintf(&s, fs.c_str(), x.x);
    std::string ss(s);
    free(s);

    if (pt == presentation_type::none && ss.find('.') != std::string::npos) {
        auto e = ss.find('e');
        if (e == std::string::npos)
            e = ss.size();
        auto t = e;
        for( ; ss[t-1] == '0' ; )
            t--;
        if (ss[t-1] == '.')
            t--;
        ss = ss.substr(0, t) + ss.substr(e);
    }
    format_to(ctx.out(), "{}", ss);
    return ctx.out();
}

std::istream & operator>>(std::istream & is, cxx_mpfr::input_with_precision xp)
{
    using cado_expression_parser_details::number_literal;
    using traits = cado_expression_parser_details::number_traits<cxx_mpfr>;
    number_literal N;
    if (is >> N) {
        try {
            xp.x = traits::from_number_literal(N, xp.p);
            if (mpfr_nan_p(xp.x.x) && !N.full.empty() && N.full[0] == '-')
                mpfr_setsign(xp.x, xp.x, 1, MPFR_RNDN);
        } catch (cado_expression_parser_details::parse_error const &) {
            is.setstate(std::ios::failbit);
        }
    }
    return is;
}
