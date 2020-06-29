#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <sstream>

#include "las-bkmult.hpp"

#include "macros.h"

buckets_are_full::buckets_are_full(bkmult_specifier::key_type const& key, int b, int r, int t) : key(key), bucket_number(b), reached_size(r), theoretical_max_size(t) {
    std::ostringstream os;
    os << "Fullest level-"<<bkmult_specifier::printkey(key)<<" bucket #"<<b<<", wrote "<<reached_size<<"/"<<theoretical_max_size<<"";
    message = os.str();
}

/* provide these here in order to avoid emitting this code for all users.
 */
buckets_are_full::~buckets_are_full() = default;
buckets_are_full::buckets_are_full(buckets_are_full const &) = default;

bkmult_specifier::bkmult_specifier(const char * specifier)
{
    const char * p = specifier;
    for(const char * q ; *p ; p = q + (*q != '\0') ) {
        const char * colon = NULL;
        for(q = p ; *q && *q != ',' ; q++)
            if (*q == ':')
                colon = q;
        /* parse from p to q */
        if (colon) {
            std::string vs(colon+1, q);
            std::istringstream is(vs);
            double v;
            is >> v;
            ASSERT_ALWAYS(!(is.rdstate() & std::ios_base::failbit));
            ASSERT_ALWAYS(colon - p == 2);
            ASSERT_ALWAYS(p[0] >= '1' && p[0] <= '9');
            ASSERT_ALWAYS(p[1] == 's' || p[1] == 'l');
            dict_t::key_type key(p[0]-'0', p[1]);
            dict.insert(std::make_pair(key, v));
        } else {
            std::string vs(p, q);
            std::istringstream is(vs);
            is >> base;
            ASSERT_ALWAYS(!(is.rdstate() & std::ios_base::failbit));
        }
    }
}

std::string bkmult_specifier::print_all() const
{
    std::ostringstream os;
    os << base;
    for(auto const & x : dict) {
        key_type const& key(x.first);
        os << "," << key.first << key.second << ":" << x.second;
    }
    return os.str();
}
