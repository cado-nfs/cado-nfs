#ifndef CADO_LAS_INFO_HPP
#define CADO_LAS_INFO_HPP

#include "cado_config.h"               // for HAVE_HWLOC

#include <condition_variable>
#include <cstdint>
#include <list>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cado_poly.h"
#include "ecm/batch.hpp"
#include "fb.hpp"
#include "gmp_aux.h"
#ifdef HAVE_HWLOC
#include "hwloc-aux.h"
#endif
#include "ecm/facul_strategies.hpp"
#include "las-bkmult.hpp"
#include "las-cofactor.hpp"
#include "las-dlog-base.hpp"
#include "las-memory.hpp"
#include "las-parallel.hpp"
#include "las-side-config.hpp"
#include "las-sieve-shared-data.hpp"
#include "las-siever-config.hpp"
#include "params.h"
#include "trialdiv.hpp"
#include "utils_cxx.hpp"

// scan-headers: stop here

/* forward decls of j_divisibility_helper and unsieve_data are not
 * sufficient.
 */
#include "las-unsieve.hpp"      // IWYU pragma: keep
#include "lock_guarded_container.hpp"  // for lock_guarded_container

/* This one wants to have siever_config defined */
#include "las-descent-trees.hpp"

// #define HILIGHT_START   "\e[01;31m"
// #define HILIGHT_END   "\e[00m"
#define HILIGHT_START   ""
#define HILIGHT_END   ""


/*  las_info
 *
 * las_info holds general data, mostly unrelated to what is actually
 * computed within a sieve. las_info also contains outer data, which
 * lives outside the choice of one particular way to configure the siever
 * versus another.
 */
struct las_info : public las_parallel_desc, private NonCopyable {
    // ----- general operational flags
    const char * galois; /* a string to indicate which galois to use in las */
    int suppress_duplicates;
    int adjust_strategy = 0;

    /* It's not ``general operational'', but global enough to be here */
    cxx_cado_poly cpoly;
    cxx_gmp_randstate rstate;

    // ----- default config and adaptive configs
    siever_config_pool config_pool;

    private:
    /* It's slightly unfortunate to have "mutable" here, obviously. The
     * root cause is the fetching of strategies for cofactoring in
     * duplicate suppression mode. As the call to relation_is_duplicate
     * happens really deep in the call chain, we have a const ref to las
     * at this point, and this is generally a good thing ! So the
     * get_strategies() method must be const.
     *
     * Other cache access calls, on the other hand, are all relatively
     * shallow and can access a non-const ref to las very close by. So
     * there is no compelling need to have the mutable keyword here, as
     * we can afford to call them with the non-const ref (see
     * nfs_work::prepare_for_new_q and nfs_work_cofac::nfs_work_cofac).
     */
#ifdef HAVE_HWLOC
    mutable std::map<cxx_hwloc_nodeset, sieve_shared_data> shared_structure_cache;
    sieve_shared_data & local_cache() const {
        return shared_structure_cache.at(current_memory_binding());
    }
    std::map<cxx_hwloc_nodeset, las_memory_accessor> las_memory_accessor_cache;
    public:
    las_memory_accessor & local_memory_accessor() {
        return las_memory_accessor_cache.at(current_memory_binding());
    }
    private:
#else
    mutable sieve_shared_data shared_structure_private;
    sieve_shared_data & local_cache() const {
        return shared_structure_private;
    }
    las_memory_accessor las_memory_accessor_private;
    public:
    las_memory_accessor & local_memory_accessor() {
        return las_memory_accessor_private;
    }
    private:
#endif

    public:
    /* These accessors are for everyone to use. */
    fb_factorbase::slicing const * get_factorbase_slicing(int side, fb_factorbase::key_type fbK) {
        sieve_shared_data::side_data & s(local_cache().sides[side]);
        return s.get_factorbase_slicing(fbK);
    }
    trialdiv_data const * get_trialdiv_data(int side, fb_factorbase::key_type fbK, fb_factorbase::slicing const * fbs) {
        return local_cache().sides[side].get_trialdiv_data(fbK, fbs);
    }
    unsieve_data const * get_unsieve_data(siever_config const & conf) {
        return local_cache().get_unsieve_data(conf);
    }
    j_divisibility_helper const * get_j_divisibility_helper(int J) {
        return local_cache().get_j_divisibility_helper(J);
    }
    facul_strategies const * get_strategies(siever_config const & conf) const {
        return local_cache().get_strategies(conf);
    }
    bool no_fb(int side) const {
        return local_cache().sides[side].no_fb();
    }

    private:
    bkmult_specifier bk_multiplier { 1.0 };
    mutable std::mutex mm;

    public:
    bool grow_bk_multiplier(bkmult_specifier::key_type const& key, double& ratio, double & new_value, double& old_value) {/*{{{*/
        const std::lock_guard<std::mutex> foo(mm);
        old_value = bk_multiplier.get(key);
        if (old_value > new_value) {
            return false;
        } else {
            /* no point in growing a lot more than new_value! */
            if (ratio * old_value > new_value) {
                ratio = new_value / old_value;
            }
            new_value = bk_multiplier.grow(key, ratio);
            return true;
        }
    }/*}}}*/
    bkmult_specifier get_bk_multiplier() const {/*{{{*/
        const std::lock_guard<std::mutex> foo(mm);
        return bk_multiplier;
    }/*}}}*/

    /* For composite special-q: note present both in las_info and
     * las_todo_list */
    bool allow_composite_q = false;
    uint64_t qfac_min = 1024;
    uint64_t qfac_max = UINT64_MAX;
    bool is_in_qfac_range(uint64_t p) const {
        return (p >= qfac_min) && (p <= qfac_max);
    }

    std::vector<unsigned long> dupqmin;   /* smallest q sieved, for dupsup */
    std::vector<unsigned long> dupqmax;   /* largest q sieved, for dupsup */
 
    // ----- stuff roughly related to the descent
    /* This is an opaque pointer to C++ code. */
    void * descent_helper;
    las_dlog_base dlog_base;
    mutable descent_tree tree;
    void init_hint_table(param_list_ptr);
    void clear_hint_table();

    // ----- special mode when we don't compute relations, but read them
    // from a relation cache instead.

    std::string relation_cache;
    void reproduce_relations_from_cache(las_todo_entry const & doing);
    
    // ----- batch mode
    int batch; /* batch mode for cofactorization */

    /* the batch_print_survivors holds several variables that are
     * attached to the process of printing the survivors to files. This
     * is activated by the batch_print_survivors option, and only if the
     * batch_print_survivors_filename option is set (otherwise printing
     * just goes to stdout, period). Printing
     * is typically handled by one or maybe several print-dedicated
     * threads. In most normal cases, one printing thread should be
     * enough. Most fields in here are accessed with mm acquired. cv is
     * here for printing threads waiting for new work notifications.
     */
    struct batch_print_survivors_t {
        mutable std::mutex mm;
        mutable std::condition_variable cv;
        bool done = false;
        std::list<cofac_list> todo;

        const char *filename; // basename for the files
        uint64_t    filesize; // number of survivors per file
        int         counter;  // current index of filename
        int number_of_printers;
        std::vector<std::thread>  printer_threads;     // id of the thread doing writing (if any)
        operator bool() const { return filename != nullptr; }

        void doit();
    } batch_print_survivors;

    std::vector<batch_side_config> bsides;

    /* Would this rather go somewhere else ? In a global (not per-sq)
     * version of nfs_work_cofac perhaps ?
     * 
     * We're really over-using mutable modifiers with this struct. It's
     * annoying me, and more than a hint at the fact that we should think
     * the design a bit differently.
     */
    mutable lock_guarded_container<cofac_list> L; /* store (a,b) and corresponding cofactors in batch mode */

    /* ----- cofactorization statistics for the default config */
    mutable cofactorization_statistics cofac_stats;

    const char *dump_filename;


    /* typicall call order is as follows */
    explicit las_info(cxx_param_list &);
    template<typename... Args> void set_parallel(cxx_param_list &pl, Args&& ...args) {
        (las_parallel_desc&)*this = las_parallel_desc(pl, std::forward<Args>(args)...);
        prepare_sieve_shared_data(pl);
    }
    void prepare_sieve_shared_data(cxx_param_list & pl);
    void load_factor_base(cxx_param_list & pl);

    static void declare_usage(cxx_param_list & pl);
    static void configure_switches(cxx_param_list & pl);
    static void configure_aliases(cxx_param_list & pl);
};

#endif	/* CADO_LAS_INFO_HPP */
