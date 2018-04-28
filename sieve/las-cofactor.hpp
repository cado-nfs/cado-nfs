#ifndef LAS_COFACTOR_HPP_
#define LAS_COFACTOR_HPP_

#include "gmp.h"
#include "las-siever-config.hpp"
#include "cxx_mpz.hpp"
#include "ecm/facul.hpp"
#include <array>
#include <vector>

class cofactorization_statistics {
    FILE * file;
    std::vector<std::vector<uint32_t>> cof_call;
    std::vector<std::vector<uint32_t>> cof_success;
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
};

int check_leftover_norm (cxx_mpz const & n, siever_config::side_config const & sc);

int factor_both_leftover_norms(
        std::array<cxx_mpz, 2> & norms,
        std::array<std::vector<cxx_mpz>, 2> &,
        std::array<unsigned long, 2> const &,
        facul_strategies_t const *);

#endif
