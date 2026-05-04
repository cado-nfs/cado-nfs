#ifndef UTILS_PRIME_POWER_FACTORIZATION_HPP_
#define UTILS_PRIME_POWER_FACTORIZATION_HPP_

#include <map>
#include <utility>

#include "cxx_mpz.hpp"
#include "params.hpp"

namespace cado {

struct prime_power : public std::pair<cxx_mpz, int> {};
struct prime_power_factorization : public std::map<cxx_mpz, int> {};

} /* namespace cado */


namespace cado::params {
template<> struct parser<cado::prime_power> {
    bool operator()(std::string const & s, cado::prime_power & value) const {
        std::string pstr;
        const size_t pos = s.find('^');
        if (pos == std::string::npos) {
            pstr = s;
            value.second = 1;
        } else {
            pstr = s.substr(0, pos);
            if (!parse(s.substr(pos + 1), value.second))
                return false;
        }
        if (!parse(pstr, value.first))
            return false;
        return true;
    }
};

template<> struct parser<cado::prime_power_factorization> {
    bool operator()(std::string const & s, cado::prime_power_factorization & value) const {
        std::vector<cado::prime_power> tmp;
        /* also tolerate commas, they're easier to the shell */
        if (!parse(s, tmp, "*") && !parse(s, tmp, ","))
            return false;
        value.clear();
        for(auto const & [ p, e ] : tmp)
            value[p] += e;
        return true;
    }
};

} /* namespace cado::params */

#endif	/* UTILS_PRIME_POWER_FACTORIZATION_HPP_ */
