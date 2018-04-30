#ifndef LAS_SIEVE_INFO_HPP_
#define LAS_SIEVE_INFO_HPP_


/* various includes that are needed to define the sieve_info type */
#include "cado_poly.h"
#include "las-siever-config.hpp"
#include "bucket.hpp"  /* for bkmult_specifier */
#include "las-todo-entry.hpp"
#include "las-qlattice.hpp"
#include "las-norms.hpp"
#include "trialdiv.h"
#include "fb.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-smallsieve-types.hpp"
#include "las-unsieve.hpp"
#include "ecm/facul.hpp"

#include <memory>
#ifdef HAVE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace std { using boost::shared_ptr; using boost::make_shared; }
#endif

/* General information about the siever, based on some input-dependent
 * configuration data with an impact on the output (as opposed to e.g.
 * file names, or verbosity flags, which do not affect the output).
 *
 * A sieve_info struct depends on the special-q.
 */
struct sieve_info {
    cxx_cado_poly const * cpoly_ptr; /* The polynomial pair */
    inline cxx_cado_poly const & cpoly() { return *cpoly_ptr; }

    /* This conditions the validity of the sieve_info_side members
     * sides[0,1], as well as some other members. If config fields match,
     * we can basically use the same sieve_info structure. */
    siever_config conf;

    /* This field gets initialized only via get_sieve_info_from_config */
    const bkmult_specifier * bk_multiplier = nullptr;

    las_todo_entry doing;

    uint32_t J;
    uint32_t I;

    // description of the q-lattice. The values here should remain
    // compatible with those in ->conf (this concerns notably the bit
    // size as well as the special-q side).
    //
    // Within las, this is set via sieve_info::recover_per_sq_values(),
    // which is quite normal since skew gauss is done within
    // sieve_range_adjust.
    qlattice_basis qbasis;

    // parameters for bucket sieving
    /* Actual number of buckets at toplevel used by current special-q */
    uint32_t nb_buckets[FB_MAX_PARTS];

    /* Largest level for which the corresponding fb_part is not empty.
     * This is set via sieve_info::update */
    int toplevel;

    /* {{{ sieve_side_info */
    struct side_info {
        /* lognorms is set by sieve_info::update_norm_data(), and depends on
         *      conf.logA
         *      conf.sides[side].lpb
         *      conf.sides[side].sublat
         *      cpoly
         *      qbasis
         *      conf.logI
         *      conf.J
         *
         * We recompute the lognorms for each q.
         */
        std::shared_ptr<lognorm_base> lognorms;

        /* fb depends on
         *      cpoly
         *      conf.sides[side].lim
         *      conf.sides[side].powlim
         *
         * Morally there should be one single factor base. Currently, in
         * las_descent, we have several. See comment in
         * sieve_info::sieve_info
         *
         * sieve_info::sieve_info is also the place from where the factor
         * base is created.
         */
        std::shared_ptr<fb_factorbase> fb;

        /* fbK depends on
         *      conf.logI
         *      conf.sides[side].lim
         *      nthreads
         * and is set by sieve_info::update */
        fb_factorbase::key_type fbK;

        /* *fbs is always (*fb)[fbK]. It's kept as a pointer for speedy
         * access and so that the structure remains
         * default-constructible. This is set by sieve_info::update
         */
        fb_factorbase::slicing * fbs = NULL;

        /* trialdiv_data depends on fbK.td_thresh, and marginally on
         * fbK.thresholds[0].
         *
         * TODO
         * Currently it is set by sieve_info::update, but this looks
         * completely wrong. There's absolutely no reason why we can't
         * cache it based on fbK.td_thresh.
         */
        std::shared_ptr<trialdiv_divisor_t> trialdiv_data;


        /* TODO
         * fb_smallsieved is a middle man which we should kill. Really,
         * fbs->small_sieve_entries has everything we need, and we should
         * use that instead.
         *
         * This is currently set in sieve_info::update
         */
        std::shared_ptr<std::vector<fb_entry_general> > fb_smallsieved;
        /* Offsets into fb_smallsieved; the entries between these
           two offsets are to be small-sieved, the others are not. */
        size_t resieve_start_offset, resieve_end_offset;

        /* precomp_plattice_dense: caching of the FK-basis in sublat mode.
         * (for the toplevel only)
         *
         * It's quite unholy. A vector of pointers to vectors of FK bases.
         *
         * Allocation is made in las-fill-in-buckets.cpp
         *
         * There's some cleanup that is done directly in las.
         * */
        precomp_plattice_dense_t precomp_plattice_dense;

        /* This is updated by applying the special-q lattice transform to
         * the factor base. This is a "current status" that gets updated
         * as we sieve. It's initialized by small_sieve_init
         */
        small_sieve_data_t ssd;
    };
    /* }}} */

    side_info sides[2];

    /* Most of the member functions below which take a side argument
     * should be members of the side_info structure instead.
     *
     * Note that have have only init_* functions below, no clear_*. In
     * fact, it's mostly a misnomer, and I apologize for it. The
     * automatic dtor will work fine for all these fields, hence there is
     * no need for specific clear_* functions. The magic lies most of the
     * time in the use of shared_ptr's. Maybe the init_* ones should have
     * been named setup_* or something like this.
     */

    /* in las-fb.cpp */
    // void init_factor_bases( param_list_ptr pl);
    // void share_factor_bases(sieve_info & other);

    /* in las-trialdiv.cpp */
    void init_trialdiv(int side);

    /* in las-unsieve.cpp */
    /* Data for unsieving locations where gcd(i,j) > 1 */
    /* This gets initialized only when I is finally decided */
    unsieve_data us;
    void init_unsieve_data() { us = unsieve_data(conf.logI, conf.logA); }

    /* in las-unsieve.cpp */
    /* Data for divisibility tests p|i in lines where p|j */
    /* This gets initialized only when J is finally decided. */
    j_divisibility_helper j_div;
    void init_j_div() { j_div = j_divisibility_helper(J); }

    /* in las-cofactor.cpp */
    std::shared_ptr<facul_strategies_t> strategies;
    void init_strategies(param_list_ptr pl);

    /* These functions must be called before actually sieving */
    void update_norm_data ();
    void update (unsigned int nr_workspaces);

    sieve_info(siever_config const & sc, cxx_cado_poly const & cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl, bool try_fbc = false);

    sieve_info() {
        cpoly_ptr = NULL;
        I = J = 0;
        memset(nb_buckets, 0, sizeof(nb_buckets));
        toplevel = 0;
    }

    /* This constructor is meant to cover the case where we want to store
     * all that to a common registry. It would be cleaner if the registry
     * were a hidden feature of the class, buried in las-types.cpp. The
     * danger is with the static initialization order fiasco which lurks
     * near, so I've been refraining so far.
     */
    static sieve_info & get_sieve_info_from_config(siever_config const & sc, cxx_cado_poly const & cpoly, std::list<sieve_info> & registry, cxx_param_list & pl, bool try_fbc = false);

    void set_for_new_q(las_todo_entry const & doing, qlattice_basis const & Q, int J) {
        this->doing = doing;
        qbasis = Q;
        qbasis.set_q(doing.p, doing.prime_sq);
        if (!qbasis.prime_sq) {
            qbasis.prime_factors = doing.prime_factors;
        }
        this->J = J;
    }
};

/*  */

#endif	/* LAS_SIEVE_INFO_HPP_ */
