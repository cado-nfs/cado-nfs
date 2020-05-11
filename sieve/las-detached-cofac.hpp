#ifndef LAS_DETACHED_COFAC_HPP_
#define LAS_DETACHED_COFAC_HPP_

#include <cstdint>
#include <array>
#include <vector>
#include <memory>
#include "cxx_mpz.hpp"
#include "las-divide-primes.hpp"
#include "las-todo-entry.hpp"
#include "las-qlattice.hpp"
#include "relation.hpp"
#include "las-threads-work-data.hpp"

struct cofac_standalone {
    std::array<uint8_t, 2> S;
    std::array<cxx_mpz, 2> norm;
    std::array<factor_list_t, 2> factors;
    std::array<std::vector<cxx_mpz>, 2> lps;
    int64_t a;
    uint64_t b;
#ifdef SUPPORT_LARGE_Q
    cxx_mpz az, bz;
#endif
    cofac_standalone();
    cofac_standalone(int N, size_t x, int logI, qlattice_basis const & Q);
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
    bool gcd_coprime_with_q(las_todo_entry const & E);
    bool ab_coprime() const;
    void print_as_survivor(FILE * f);
    relation get_relation(las_todo_entry const & doing);
    void transfer_to_cofac_list(lock_guarded_container<cofac_list> & L, las_todo_entry const & doing);
    int factor_both_leftover_norms(nfs_work_cofac & wc);
};

struct detached_cofac_parameters : public cofac_standalone, public task_parameters {
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    detached_cofac_parameters(std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, cofac_standalone&& C) : cofac_standalone(C), wc_p(wc_p), aux_p(aux_p) {}
};

struct detached_cofac_result : public task_result {
    std::shared_ptr<relation> rel_p;
};

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int);

#endif	/* LAS_DETACHED_COFAC_HPP_ */
