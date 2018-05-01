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
#include "las-sieve-info.hpp"
#include "las-cofactor.hpp"
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

// #define HILIGHT_START   "\e[01;31m"
// #define HILIGHT_END   "\e[00;30m"
#define HILIGHT_START   ""
#define HILIGHT_END   ""

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
    cxx_cado_poly cpoly;
    gmp_randstate_t rstate;

    // ----- default config and adaptive configs
    siever_config_pool config_pool;

    bkmult_specifier bk_multiplier { 1.0 };

    void grow_bk_multiplier(bkmult_specifier::key_type const& key, double d) {
        bk_multiplier.grow(key, d);
    }

    /* There may be several configured sievers. This is used mostly for
     * the descent.  */
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

    /* Would this rather go somewhere else ? In a global (not per-sq)
     * version of nfs_work_cofac perhaps ?
     * 
     * We're really over-using mutable modifiers with this struct. It's
     * annoying me, and more than a hint at the fact that we should think
     * the design a bit differently.
     */
    mutable cofac_list L; /* store (a,b) and corresponding cofactors in batch mode */

    /* ----- cofactorization statistics for the default config */
    mutable cofactorization_statistics cofac_stats;

    const char *dump_filename;
    mutable dumpfile dumpfiles[2];

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
