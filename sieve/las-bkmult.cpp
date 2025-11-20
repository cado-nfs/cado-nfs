#include "cado.h" // IWYU pragma: keep

#include <ios>
#include <sstream>

#include "fmt/format.h"

#include "utils_cxx.hpp"
#include "las-bkmult.hpp"

#include "macros.h"

buckets_are_full::buckets_are_full(bkmult_specifier::key_type const & key,
                                   int b, int r, int t)
    : key(key)
    , bucket_number(b)
    , reached_size(r)
    , theoretical_max_size(t)
{
    std::ostringstream os;
    os << "Fullest level-" << bkmult_specifier::printkey(key) << " bucket #"
       << b << ", wrote " << reached_size << "/" << theoretical_max_size << "";
    message = os.str();
}

bkmult_specifier::bkmult_specifier(char const * specifier)
{
    char const * p = specifier;
    for (char const * q; *p; p = q + (*q != '\0')) {
        char const * colon = nullptr;
        for (q = p; *q && *q != ','; q++)
            if (*q == ':')
                colon = q;
        /* parse from p to q */
        if (colon) {
            std::string const vs(colon + 1, q);
            std::istringstream is(vs);
            double v;
            is >> v;
            ASSERT_ALWAYS(!(is.rdstate() & std::ios_base::failbit));
            ASSERT_ALWAYS(colon - p == 2);
            ASSERT_ALWAYS(p[0] >= '1' && p[0] <= '9');
            ASSERT_ALWAYS(p[1] == 's' || p[1] == 'l');
            dict_t::key_type const key(p[0] - '0', p[1]);
            dict[key] = v;
        } else {
            std::string const vs(p, q);
            std::istringstream is(vs);
            is >> base;
            ASSERT_ALWAYS(!(is.rdstate() & std::ios_base::failbit));
        }
    }
}

std::string bkmult_specifier::print_all() const
{
    return fmt::format("{},{}",
            base, join(dict.begin(), dict.end(), " ",
                [](auto const& item) {
                    auto const & [ key, mult ] = item;
                    auto const & [ lev, hint ] = key;
                    return fmt::format("{}{}:{:.3f}", lev, hint, mult);
                }));
}
