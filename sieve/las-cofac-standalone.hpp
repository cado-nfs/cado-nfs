#ifndef CADO_LAS_COFAC_STANDALONE_HPP
#define CADO_LAS_COFAC_STANDALONE_HPP

#include <cstdio>                // for FILE, size_t
#include <cstdint>                // for uint8_t, int64_t, uint64_t
#include <vector>                 // for vector

#include "cxx_mpz.hpp"
#include "ecm/batch.hpp"              // for cofac_list
#include "las-divide-primes.hpp"  // for factor_list_t
#include "relation.hpp"           // for relation
class nfs_work_cofac; // IWYU pragma: keep
struct las_todo_entry; // IWYU pragma: keep
struct qlattice_basis; // IWYU pragma: keep
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
    explicit operator relation_ab() const { return { az, bz }; }
#else
    explicit operator relation_ab() const { return { a, b }; }
#endif
    cofac_standalone();
    cofac_standalone(int nsides, int N, size_t x, int logI, qlattice_basis const & Q);
    bool trace_on_spot() const;
    /* TODO. Hmmm. How important is this ? We don't want to expose
     * dependence on a compile flag in a header */
    inline bool both_even() const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        return ((((a | b) & 1) == 0));
#else
        return ((mpz_even_p(az) && mpz_even_p(bz)));
#endif
    }/*}}}*/
    bool gcd_coprime_with_q(las_todo_entry const & E) const;
    bool ab_coprime() const;
    void print_as_survivor(FILE * f);
    relation get_relation(las_todo_entry const & doing) const;
    void transfer_to_cofac_list(lock_guarded_container<cofac_list> & L, las_todo_entry const & doing);
    int factor_both_leftover_norms(nfs_work_cofac & wc);
};


#endif	/* CADO_LAS_COFAC_STANDALONE_HPP */
