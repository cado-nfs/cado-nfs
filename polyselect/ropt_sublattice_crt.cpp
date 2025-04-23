#include "cado.h" // IWYU pragma: keep

#include <cmath>

#include "ropt_single_sublattice_priority_queue_impl.hpp"
#include "ropt_sublattice_crt.h"
#include "cxx_mpz.hpp"


struct ropt_sublattice_crt_combiner {
    unsigned int nprimes;
    const unsigned int * primes;

    /* tops[0] ... tops[nprimes-1] is our main input */
    single_sublattice_priority_queue const * tops;

    ropt_bound_srcptr bound;

    /* pqueue is our only output */
    sublattice_priority_queue_ptr pqueue;

    /* pointers to within tops[0] ... to tops[nprimes-1] */
    std::vector<single_sublattice_info const *> infos;

    /* CRT multipliers */
    std::vector<cxx_mpz> multipliers;

    std::vector<float> logp;

    std::vector<unsigned int> pe;       // holds p^e
    cxx_mpz modulus;

    ropt_sublattice_crt_combiner(unsigned int nprimes, const unsigned int * primes, single_sublattice_priority_queue const * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue)
        : nprimes(nprimes)
        , primes(primes)
        , tops(tops)
        , bound(bound)
        , pqueue(pqueue)
        , infos(nprimes, nullptr)
        , multipliers(nprimes)
        , logp(nprimes, 0)
        , pe(nprimes, 0)
    {
        for(size_t i = 0 ; i < nprimes ; i++)
            logp[i] = log( (double) primes[i] );
    }

    void finalize_multipliers() {
        /* I wonder whether it is costly to recompute the modulus each
         * time around. I would guess not, but haven't checked.
         */
        mpz_set_ui(modulus, 1);

        for(size_t i = 0 ; i < nprimes ; i++) {
            pe[i] = primes[i];
            for(unsigned int j = 1 ; j < (unsigned int) infos[i]->e ; j++)
                pe[i] *= primes[i];
            mpz_mul_ui(modulus, modulus, pe[i]);
        }
        /* It's fine with very small sets of primes, but this approach is
         * sub-par when nprimes becomes large (which is not the case here
         * anyway) */
        cxx_mpz tmp;
        for (size_t i = 0; i < nprimes; i ++) {
            mpz_divexact_ui (multipliers[i], modulus, pe[i]);
            mpz_set_ui (tmp, pe[i]);
            mpz_invert (tmp, multipliers[i], tmp);
            mpz_mul (multipliers[i], multipliers[i], tmp);
        }
    }

    /**
     * Compute crt and add (u, v) to queue.
     *
     * infos[i]: lattice modulo tprimes[i], to power infos[i]->e
     */
    unsigned int combine () {
        cxx_mpz sum, tmp;
        cxx_mpz u_crt, v_crt;

        finalize_multipliers();

        /* compute u */
        mpz_set_ui (sum, 0);
        for (size_t i = 0; i < nprimes; i ++) {
            mpz_addmul_ui(sum, multipliers[i], infos[i]->u);
        }
        mpz_mod (u_crt, sum, modulus);
        mpz_sub (tmp, modulus, u_crt); // tmp > 0

        unsigned int count = 0;

        /* check if a signed, centered representative of u falls within
         * global_u_boundl .. global_u_boundr
         */
        /* if u is good, compute v */
        if ( mpz_cmp_si (u_crt, bound->global_u_boundr) <= 0 ||
                mpz_cmp_si (tmp, -bound->global_u_boundl) <= 0 ) {

            /* compute v */
            float val = 0.0;
            mpz_set_ui (sum, 0);
            for (size_t i = 0; i < nprimes; i ++) {
                mpz_addmul_ui(sum, multipliers[i], infos[i]->v);
                val += logp[i] * infos[i]->val;
            }
            mpz_mod (v_crt, sum, modulus);

            mpz_add(tmp, bound->global_v_boundl, modulus);

            if (mpz_cmp(v_crt, bound->global_v_boundr) <= 0 ||
                    mpz_cmp(v_crt, tmp) >= 0 ) {

                /* insert this node */
                sublattice_priority_queue_push(pqueue, u_crt, v_crt, modulus, val );
                count = 1;
            }
        }

        return count;
    }

    unsigned int operator()(unsigned int i = 0)
    {
        if (i == nprimes)
            return combine();
        /* iterate through all elements of tops[i], put them one by one in
         * infos[i], and recurse on that.
         */
        auto qi = single_sublattice_priority_queue_impl::cast(tops[i]);
        unsigned int count = 0;
        for(auto const & q : *qi) {
            infos[i] = &q;
            count += (*this)(i+1);
        }
        return count;
    }
};

unsigned int ropt_sublattice_combine_all_crt(unsigned int nprimes, const unsigned int * primes, single_sublattice_priority_queue const * tops, ropt_bound_srcptr bound, sublattice_priority_queue_ptr pqueue)
{
    return ropt_sublattice_crt_combiner(nprimes, primes, tops, bound, pqueue)();
}
