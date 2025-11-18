#ifndef CADO_POWERS_OF_P_HPP
#define CADO_POWERS_OF_P_HPP

#include <cstddef>

#include <map>
#include <vector>
#include <mutex>

#include <gmp.h>

#include "cxx_mpz.hpp"

struct power_lookup_table {
    mutable std::mutex mx;
    unsigned long p;
    std::map<int, size_t> m;
    std::vector<cxx_mpz> z;
    power_lookup_table(unsigned long p) : p(p) { }
    private:
    power_lookup_table(power_lookup_table&& o, std::lock_guard<std::mutex> const &)
        : p(o.p)
        , m(std::move(o.m))
        , z(std::move(o.z))
    {}
    public:
    power_lookup_table(power_lookup_table&& o) noexcept
        : power_lookup_table(std::move(o), std::lock_guard<std::mutex>(mx))
        {}

    cxx_mpz const & operator()(int i);
    cxx_mpz const & operator()(int i) const;
    private: cxx_mpz const & inside(int i);
};

#endif	/* CADO_POWERS_OF_P_HPP */
