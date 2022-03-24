#ifndef LAS_SIEVE_SHARED_DATA_HPP_
#define LAS_SIEVE_SHARED_DATA_HPP_

#include <array>                       // for array
#include <map>                         // for map, swap
#include <memory>                      // for shared_ptr
#include <utility>                     // for pair

#include "cado_poly.h"
#include "mpz_poly.h"
#include "ecm/facul.hpp"                   // for facul_strategies_t
#include "fb-types.h"                  // for fbprime_t
#include "fb.hpp"                      // for fb_factorbase, fb_factorbase::...
#include "las-siever-config.hpp"       // for siever_config, siever_config::...
#include "lock_guarded_container.hpp"  // for lock_guarded_container
#include "trialdiv.hpp"                // for trialdiv_data
#include "gmp_aux.h"

struct j_divisibility_helper; // IWYU pragma: keep
struct unsieve_data; // IWYU pragma: keep
struct cxx_param_list;

/*
 * A sieve_shared_data struct is only a thin layer above the factor base. It
 * provides the factor base, and the possibility to access some data
 * structures that depend on the current configuration (e.g. on I).
 *
 * We have a singleton of this type, contained in the las_info structure.

 *
 * It must be regarded as a (mostly) const repository from where we fetch
 * data that is of interest to many special-q's.
 *
 * The data at the sieve_shared_data level should not be attached to a
 * particular config, nor to particular q.
 *
 * In the multi-las context, sieve_shared_data is *shared* between the different
 * sub-processes.
 *
 * TODO: at some point, the query mechanism for the cached entries will
 * depend on the memory binding.
 */
struct sieve_shared_data {
    cxx_cado_poly cpoly;

    /* {{{ side_data */
    struct side_data {
        cxx_mpz_poly f;     /* handy. Alias to cpoly.pols[side] */

        /* fb depends on
         *      cpoly
         *      conf.sides[side].lim
         *      conf.sides[side].powlim
         *
         */
        fb_factorbase fb;

        inline bool no_fb() const { return fb.empty(); }

        fb_factorbase::slicing const * get_factorbase_slicing(fb_factorbase::key_type fbK) {
            /* implemented in fb_factorbase::operator[] (fb.hpp), with
             * locking.
             */
            return &fb[fbK];
        }


        /* trialdiv_data depends on
         *      fbK.td_thresh
         * and marginally on
         *      fbK.thresholds[0]
         */

        private:
        struct equivalent_fbK_for_td {
            bool operator()(fb_factorbase::key_type const & a, fb_factorbase::key_type const & b) const {
                if (a.thresholds[0] < b.thresholds[0]) return true;
                if (a.thresholds[0] > b.thresholds[0]) return false;
                return a.td_thresh < b.td_thresh;
            }
        };
        lock_guarded_container<
            std::map<
                fb_factorbase::key_type,
                trialdiv_data,
                equivalent_fbK_for_td
            >
        > trialdiv_data_cache;
        cxx_gmp_randstate rstate;       /* use is protected by trialdiv_data_cache.mutex() */
        public:
        /* in las-trialdiv.cpp */
        trialdiv_data const * get_trialdiv_data(fb_factorbase::key_type fbK, fb_factorbase::slicing const * fbs);

        side_data(int side, cxx_cado_poly const & cpoly, cxx_param_list & pl, int nthreads = 1);
        side_data() = default;
        side_data(side_data &&) = default;
        side_data& operator=(side_data &&) = default;
    };

    /* }}} */

    side_data sides[2];

    /* Most of the member functions below which take a side argument
     * should be members of the side_data structure instead.
     *
     * Note that have have only init_* functions below, no clear_*. In
     * fact, it's mostly a misnomer, and I apologize for it. The
     * automatic dtor will work fine for all these fields, hence there is
     * no need for specific clear_* functions. The magic lies most of the
     * time in the use of shared_ptr's. Maybe the init_* ones should have
     * been named setup_* or something like this.
     */

    /* in las-unsieve.cpp */
    /* Data for unsieving locations where gcd(i,j) > 1 */
    /* This gets initialized only when I is finally decided */
    private:
    lock_guarded_container<
        std::map<
            std::pair<int, int>,
            unsieve_data
        >
    > us_cache;
    public:
    unsieve_data const * get_unsieve_data(siever_config const & conf);

    /* in las-unsieve.cpp */
    /* Data for divisibility tests p|i in lines where p|j */
    private:
    lock_guarded_container<
        std::map<
            unsigned int,
            j_divisibility_helper
        >
    > jdiv_cache;
    public:
    j_divisibility_helper const * get_j_divisibility_helper(int J);

    /* We use a shared_ptr only to provide the custom dtor.
     *
     * We have a cache because during the descent, several cofactoring
     * strategies are possible.
     */
    private:
    const char *cofactfilename;
    lock_guarded_container<
        std::map<
            siever_config,
            std::shared_ptr<facul_strategies>,
            siever_config::has_same_cofactoring::comparison
        >
    > facul_strategies_cache;
    public:
    facul_strategies const * get_strategies(siever_config const & conf);

    ~sieve_shared_data();
    sieve_shared_data(cxx_cado_poly const & cpoly, cxx_param_list & pl);
    void load_factor_base(cxx_param_list & pl, int nthreads = 1);
    // sieve_shared_data() = default;
    sieve_shared_data(sieve_shared_data&&) = default;
    sieve_shared_data& operator=(sieve_shared_data&&) = default;
    static void declare_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list &, int nsides = 2);
};

#endif	/* LAS_SIEVE_SHARED_DATA_HPP_ */
