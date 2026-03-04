#ifndef CADO_LAS_COFAC_STANDALONE_HPP
#define CADO_LAS_COFAC_STANDALONE_HPP

#include <cstdio>
#include <cstdint>

#include <iostream>
#include <list>
#include <vector>

#include "fmt/ostream.h"

#include "cxx_mpz.hpp"
#include "ecm/batch.hpp"
#include "las-divide-primes.hpp"
#include "relation.hpp"

class nfs_work_cofac; // IWYU pragma: keep
struct special_q; // IWYU pragma: keep
struct special_q_data_base; // IWYU pragma: keep
template <typename T> struct lock_guarded_container; // IWYU pragma: keep

struct cofac_standalone {
    std::vector<uint8_t> S;
    std::vector<cxx_mpz> norm;
    std::vector<factor_list_t> factors;
    std::vector<std::vector<cxx_mpz>> lps;
    int64_t a;
    uint64_t b;

#ifdef SUPPORT_LARGE_Q
    cxx_mpz az, bz;
#endif
    explicit operator relation_ab() const
    { 
#ifdef SUPPORT_LARGE_Q
        return { az, bz };
#else
        return { a, b };
#endif
    }
    cofac_standalone();
    cofac_standalone(int nsides, int N, size_t x, int logI,
                     special_q_data_base const & Q);
    bool trace_on_spot() const;
    /* TODO. Hmmm. How important is this ? We don't want to expose
     * dependence on a compile flag in a header */
    bool both_even() const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        return ((((a | b) & 1) == 0));
#else
        return ((mpz_even_p(az) && mpz_even_p(bz)));
#endif
    }/*}}}*/
    bool gcd_coprime_with_q(special_q const & E) const;
    bool ab_coprime() const;
    void print_as_survivor(FILE * f);
    relation get_relation(special_q const & doing) const;
    void transfer_to_cofac_list(lock_guarded_container<std::list<cofac_candidate>> & L);
    int factor_leftover_norms(nfs_work_cofac & wc);
    friend std::ostream& operator<<(
            std::ostream& os,
            cofac_standalone const & cur);
};

namespace fmt {
    template <> struct formatter<cofac_standalone>: ostream_formatter {};
}

#endif	/* CADO_LAS_COFAC_STANDALONE_HPP */
