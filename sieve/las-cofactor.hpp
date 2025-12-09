#ifndef CADO_LAS_COFACTOR_HPP
#define CADO_LAS_COFACTOR_HPP

// IWYU pragma: no_include <ext/alloc_traits.h>

#include <cstdio>
#include <cstdint>

#include <mutex>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "ecm/facul.hpp"
#include "las-siever-config.hpp"
#include "params.h"
#include "ecm/facul_strategies.hpp"

class cofactorization_statistics {
    /* We rarely use this, if ever. Do we want to generalize this at all?
     * How? For the moment we just bail out if more than two sides.
     */
    FILE * file;
    std::vector<std::vector<uint32_t>> cof_call;
    std::vector<std::vector<uint32_t>> cof_success;
    std::mutex lock;
public:
    cofactorization_statistics(param_list_ptr pl);
    bool active() { return file != nullptr; }
    void call(int bits0, int bits1);
    void print();
    void call(std::vector<cxx_mpz> const & norm, std::vector<int> & cof_bitsize) {
        if (!active()) return;
        ASSERT_ALWAYS(norm.size() == 1 || norm.size() == 2);
        cof_bitsize.assign(norm.size(), 0);
        cof_bitsize[0] = mpz_sizeinbase(norm[0], 2);
        if (norm.size() == 2) {
            cof_bitsize[1] = mpz_sizeinbase(norm[1], 2);
        } else {
            cof_bitsize.push_back(0);
        }
        call(cof_bitsize[0], cof_bitsize[1]);
    }
    void success(std::vector<int> const & cof_bitsize) {
        if (!active()) return;
        ASSERT_ALWAYS(cof_bitsize.size() == 2);
        cof_success[cof_bitsize[0]][cof_bitsize[1]]++;
    }
    void success(int bits0, int bits1)
    {
        if (!active()) return;
        cof_success[bits0][bits1]++;
    }
    ~cofactorization_statistics();
    static void declare_usage(cxx_param_list & pl);
};

int check_leftover_norm (cxx_mpz const & n, siever_side_config const & sc);

facul_status factor_leftover_norms(
        std::vector<cxx_mpz> const & norms,
        std::vector<std::vector<cxx_mpz>> &,
        std::vector<unsigned long> const &,
        facul_strategies const &);

/* handy shortcut. Can't have it defined at the facul.hpp level because
 * facul does not know about las stuff. */
static inline facul_strategies * facul_make_strategies (siever_config const & conf, FILE* file, int verbose);
static inline facul_strategies * facul_make_strategies (siever_config const & conf, FILE* file, int verbose)
{
    std::vector<unsigned long> lim(conf.sides.size());
    std::vector<unsigned int> lpb(conf.sides.size());
    std::vector<unsigned int> mfb(conf.sides.size());
    std::vector<int> ncurves(conf.sides.size());
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
