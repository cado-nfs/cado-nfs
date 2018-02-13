#ifndef LAS_TYPES_HPP_
#define LAS_TYPES_HPP_

#include <stdint.h>
#include "las-config.h"
#include "las-base.hpp"
#include "las-todo-entry.hpp"
#include "las-siever-config.hpp"
#include "las-norms.hpp"
#ifdef DLP_DESCENT
#include "las-dlog-base.hpp"
#endif
#include "fb.hpp"
#include "trialdiv.h"
// #include "bucket.hpp"
#include "cado_poly.h"
#include "ecm/facul.hpp"
#include "fb-types.h"
#include "las-plattice.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-forwardtypes.hpp"
#include "las-report-stats.hpp"
#include "las-unsieve.hpp"
#include "las-qlattice.hpp"
#include "las-smallsieve-types.hpp"
#include "ecm/batch.h"
#include <list>
#include <vector>
#include <stack>

#include <memory>
#ifdef HAVE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace std { using boost::shared_ptr; using boost::make_shared; }
#endif

/* This one wants to have siever_config defined */
#include "las-descent-trees.hpp"

/* {{{ sieve_info
 *
 * General information about the siever, based on some input-dependent
 * configuration data with an impact on the output (as opposed to e.g.
 * file names, or verbosity flags, which do not affect the output).
 */
struct sieve_info {
    cado_poly_srcptr cpoly; /* The polynomial pair */

    /* This conditions the validity of the sieve_info_side members
     * sides[0,1], as well as some other members */
    siever_config conf;

    /* This field gets initialized only via get_sieve_info_from_config */
    const bkmult_specifier * bk_multiplier = nullptr;

    las_todo_entry doing;

    // sieving area. Note that in the (conf) member struct, we find
    // logI_adjusted (will be renamed eventually), as well as logA.
    uint32_t J;
    uint32_t I;

    // description of the q-lattice. The values here should remain
    // compatible with those in ->conf (this concerns notably the bit
    // size as well as the special-q side).
    qlattice_basis qbasis;

    // parameters for bucket sieving
    /* Actual number of buckets at toplevel used by current special-q */
    uint32_t nb_buckets[FB_MAX_PARTS];

    /* Largest level for which the corresponding fb_part is not empty */
    int toplevel;

    /* {{{ sieve_side_info */
    struct side_info {
        std::shared_ptr<lognorm_base> lognorms;

        std::shared_ptr<trialdiv_divisor_t> trialdiv_data;
        std::shared_ptr<std::vector<fb_general_entry> > fb_smallsieved;
        /* Iterators that point into fb_smallsieved; the entries between these
           two iterators are to be small-sieved, the others are not. */
        size_t resieve_start_offset, resieve_end_offset;

        /* The reading, mapping or generating the factor base all create the
         * factor base in several pieces: small primes, and large primes split
         * into one piece for each thread.
         */
        std::shared_ptr<fb_factorbase> fb;

        /* Caching of the FK-basis in sublat mode */
        precomp_plattice_dense_t precomp_plattice_dense;

        /* When threads pick up this sieve_info structure, they should check
         * their bucket allocation */
        double max_bucket_fill_ratio[FB_MAX_PARTS];

        /* This is updated by applying the special-q lattice transform to
         * the factor base. */
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
    void init_factor_bases( param_list_ptr pl);
    void share_factor_bases(sieve_info & other);

    void init_fb_smallsieved(int side);
    void print_fb_statistics(int side);

    /* in las-trialdiv.cpp */
    void init_trialdiv(int side);

    /* in las-unsieve.cpp */
    /* Data for unsieving locations where gcd(i,j) > 1 */
    /* This gets initialized only when I is finally decided */
    unsieve_data us;
    void init_unsieve_data() { us = unsieve_data(conf.logI_adjusted, conf.logA); }

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
    void update (size_t nr_workspaces);

    sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl);
    sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, cxx_param_list & pl);

    sieve_info() {
        cpoly = NULL;
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
    static sieve_info & get_sieve_info_from_config(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & registry, cxx_param_list & pl);

    void recover_per_sq_values(sieve_range_adjust const& Adj) {
        doing = Adj.doing;
        qbasis = Adj.Q;
        qbasis.set_q(doing.p, doing.prime_sq);
        if (!qbasis.prime_sq) {
            qbasis.prime_factors = doing.prime_factors;
        }
        ASSERT_ALWAYS(conf.logI_adjusted == Adj.logI);
        ASSERT_ALWAYS(I == (1UL << Adj.logI));
        J = Adj.J;
    }
};

/* }}} */

struct las_augmented_output_channel {
    int verbose;
    FILE *output;
    const char * outputname; /* keep track of whether it's gzipped or not */
    las_augmented_output_channel(cxx_param_list & pl);
    ~las_augmented_output_channel();
};


/* {{{ las_info
 *
 * las_info holds general data, mostly unrelated to what is actually
 * computed within a sieve. las_info also contains outer data, which
 * lives outside the choice of one particular way to configure the siever
 * versus another.
 */
struct las_info : private NonCopyable, public las_augmented_output_channel {
    // ----- general operational flags
    int nb_threads;
    const char * galois; /* a string to indicate which galois to use in las */
    int suppress_duplicates;

    /* It's not ``general operational'', but global enough to be here */
    cado_poly cpoly;
    gmp_randstate_t rstate;

    // ----- default config and adaptive configs
    siever_config_pool config_pool;

    bkmult_specifier bk_multiplier { 1.0 };

    void grow_bk_multiplier(bkmult_specifier::key_type const& key, double d) {
        bk_multiplier.grow(key, d);
    }


    /* There may be several configured sievers. This is used mostly for
     * the descent.
     * TODO: For now these different sievers share nothing of their
     * factor bases, which is a shame. We should examine ways to get
     * around this limitation, but it is hard. One could imagine some,
     * though. Among them, given the fact that only _one_ siever is
     * active at a time, it might be possible to sort the factor base
     * again each time a new siever is used (but maybe it's too
     * expensive). Another way could be to work only on sharing
     * bucket-sieved primes.
     */
    /* TODO: not sure std::list is the right container. The difficult
     * part is that we really don't want to copy the factor bases */
    std::list<sieve_info> sievers;

    // ----- todo list and various specification of what the siever will
    // be doing.
    std::stack<las_todo_entry> todo;
    unsigned int nq_pushed;
    unsigned int nq_max;
    int random_sampling;
    cxx_mpz todo_q0;
    cxx_mpz todo_q1;
    FILE * todo_list_fd;
 
    /* For composite special-q */
    int allow_composite_q;
    uint64_t qfac_min;
    uint64_t qfac_max;

    // ----- stuff roughly related to the descent
    unsigned int max_hint_bitsize[2];
    int * hint_lookups[2]; /* quick access indices into hint_table */
    /* This is an opaque pointer to C++ code. */
    void * descent_helper;
#ifdef  DLP_DESCENT
    las_dlog_base dlog_base;
#endif
    mutable descent_tree tree;
    void init_hint_table(param_list_ptr);
    void clear_hint_table();
    // ----- batch mode
    int batch; /* batch mode for cofactorization */
    int batch_print_survivors;
    cofac_list L; /* store (a,b) and corresponding cofactors in batch mode */

    /* ----- cofactorization statistics for the default config */
    FILE *cof_stats_file; // null means no stats.
    uint32_t **cof_call; /* cof_call[r][a] is the number of calls of the
                          * cofactorization routine with a cofactor of r
                          * bits on the rational side, and a bits on the
                          * algebraic side */
    uint32_t **cof_succ; /* cof_succ[r][a] is the corresponding number of
                          * successes, i.e., of call that lead to a
                          * relation */
    void init_cof_stats(param_list_ptr);
    void clear_cof_stats();
    void print_cof_stats();

    const char *dump_filename;
    dumpfile dumpfiles[2];

    las_info(cxx_param_list &);
    ~las_info();
};
/* }}} */

enum {
  OUTPUT_CHANNEL,
  ERROR_CHANNEL,
  STATS_CHANNEL,
  TRACE_CHANNEL,
  NR_CHANNELS /* This must be the last element of the enum */
};


#endif	/* LAS_TYPES_HPP_ */
