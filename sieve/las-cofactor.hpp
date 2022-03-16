#ifndef LAS_COFACTOR_HPP_
#define LAS_COFACTOR_HPP_

// IWYU pragma: no_include <ext/alloc_traits.h>

#include <cstdint>               // for uint32_t
#include <cstdio>                // for FILE, NULL
#include <array>                  // for array, array<>::value_type
#include <mutex>                  // for mutex
#include <vector>                 // for vector
#include <gmp.h>                  // for mpz_sizeinbase
#include "cxx_mpz.hpp"  // for cxx_mpz
#include "ecm/facul.hpp"          // for facul_make_strategies, facul_strate...
#include "las-siever-config.hpp"  // for siever_config::side_config, siever_...
#include "params.h"     // param_list_ptr

class cofactorization_statistics {
    FILE * file;
    std::vector<std::vector<uint32_t>> cof_call;
    std::vector<std::vector<uint32_t>> cof_success;
    std::mutex lock;
public:
    cofactorization_statistics(param_list_ptr pl);
    bool active() { return file != NULL; }
    void call(int bits0, int bits1);
    void print();
    void call(std::array<cxx_mpz, 2> const & norm, std::array<int, 2> & cof_bitsize) {
        if (!active()) return;
        cof_bitsize[0] = mpz_sizeinbase(norm[0], 2);
        cof_bitsize[1] = mpz_sizeinbase(norm[1], 2);
        call(cof_bitsize[0], cof_bitsize[1]);
    }
    void success(std::array<int, 2> const & cof_bitsize) {
        if (!active()) return;
        cof_success[cof_bitsize[0]][cof_bitsize[1]]++;
    }
    void success(int bits0, int bits1)
    {
        if (!file) return;
        cof_success[bits0][bits1]++;
    }
    ~cofactorization_statistics();
    static void declare_usage(cxx_param_list & pl);
};

int check_leftover_norm (cxx_mpz const & n, siever_config::side_config const & sc);

int factor_both_leftover_norms(
        std::array<cxx_mpz, 2> & norms,
        std::array<std::vector<cxx_mpz>, 2> &,
        std::array<unsigned long, 2> const &,
        facul_strategies const &);

/* handy shortcut. Can't have it defined at the facul.hpp level because
 * facul does not know about las stuff. */
static inline facul_strategies * facul_make_strategies (siever_config const & conf, FILE* file, const int verbose);
static inline facul_strategies * facul_make_strategies (siever_config const & conf, FILE* file, const int verbose)
{
    std::array<unsigned long, 2> lim;
    std::array<unsigned int, 2> lpb;
    std::array<unsigned int, 2> mfb;
    std::array<int, 2> ncurves;
    auto plim = lim.begin();
    auto plpb = lpb.begin();
    auto pmfb = mfb.begin();
    auto pncurves = ncurves.begin();
    for(auto const & s : conf.sides) {
        *plim++ = s.lim;
        *plpb++ = s.lpb;
        *pmfb++ = s.mfb;
        *pncurves++ = s.ncurves;
    }

    if (file) {
        return new facul_strategies(
                lim,
                lpb,
                mfb,
                (conf.sublat_bound == 0), // with sublat, some primes are skipped.
                file, verbose);
    } else {
        return new facul_strategies(
                lim,
                lpb,
                mfb,
                ncurves,
                (conf.sublat_bound == 0), // with sublat, some primes are skipped.
                verbose);
    }
}

#endif
