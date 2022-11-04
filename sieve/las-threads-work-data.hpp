#ifndef LAS_THREADS_WORK_DATA_HPP_
#define LAS_THREADS_WORK_DATA_HPP_

#include <cstddef>                     // for NULL
#include <cstdint>                     // for uint32_t
#include <array>                       // for array
#include <vector>                      // for vector
#include "ecm/batch.hpp"                   // for cofac_list
#include "ecm/facul.hpp"                   // for facul_strategies_t
#include "fb.hpp"                      // for fb_factorbase, fb_factorbase::...
#include "las-bkmult.hpp"              // for bkmult_specifier
#include "las-config.h"                // for FB_MAX_PARTS
#include "las-dumpfile.hpp"            // for dumpfile_t
#include "las-norms.hpp"               // for lognorm_smart
#include "las-plattice.hpp"            // for precomp_plattice_dense_t
#include "las-qlattice.hpp"            // for qlattice_basis
#include "las-siever-config.hpp"       // for siever_config
#include "las-smallsieve-types.hpp"    // for small_sieve_data_t
#include "las-threads.hpp"             // for reservation_group
#include "las-todo-entry.hpp"          // for las_todo_entry
#include "lock_guarded_container.hpp"  // for lock_guarded_container
#include "multityped_array.hpp"        // for multityped_array

class las_memory_accessor; // IWYU pragma: keep
class nfs_aux; // IWYU pragma: keep
class thread_pool; // IWYU pragma: keep
struct j_divisibility_helper; // IWYU pragma: keep
struct las_info; // IWYU pragma: keep
struct trialdiv_data; // IWYU pragma: keep
struct unsieve_data; // IWYU pragma: keep
template <int LEVEL, typename HINT> class bucket_array_t; // IWYU pragma: keep

#define NUMBER_OF_BAS_FOR_THREADS(n)    ((n) == 1 ? 1 : ((n) + 2))

/*
 * This structure holds the key algorithmic data that is used in las. It
 * is intentionally detached from the rest of the ``stats-like'' control
 * data (found in nfs_aux, defined in las-auxiliary-data.hpp)
 *
 * This is really sort of a "work space", with allocated buckets, bucket
 * regions, and so on.
 *
 * Two important aspects here:
 *
 *  - Everything here is attached to one special-q only. Concurrent
 *    access for various special-q's is not possible.
 *  - Allocated space for one structure may be reused for another
 *    special-q.
 *
 * We have here nb_threads threads that will work with nb_threads+1 (or 1
 * if nb_threads==1 anyway) reservation_arrays in each data member of the
 * two reservation_groups in the groups[] data member. This +1 is here to
 * allow work to spread somewhat more evenly.
 *
 * Thread-private memory areas such as bucket regions are allocated in
 * the thread_data fields.
 */
class nfs_work {
    public:
    las_info const & las;
    las_memory_accessor & local_memory;
    private:

    const int nr_workspaces;

    public:

    bkmult_specifier bk_multiplier;

    /* Largest level for which the corresponding fb_part is not empty (on
     * either side).
     * This is set via prepare_for_new_q() */
    int toplevel;

    /* This is set via prepare_for_new_q() */
    int nb_buckets[FB_MAX_PARTS];

    siever_config conf;

    qlattice_basis Q;

    /* These are fetched from the sieve_shared_data structure, which
     * caches them */
    j_divisibility_helper const * jd = NULL;
    unsieve_data const * us = NULL;
    uint32_t J = 0;

    /* This is used only in batch mode. The list of cofactorization
     * candidates will be transfered to the main list when we're done
     * with this special-q */
    lock_guarded_container<cofac_list> cofac_candidates;

    struct side_data {
        reservation_group group;

        /* lognorms is set by prepare_for_new_q(), and depends on
         *      conf.logA
         *      conf.sides[side].lpb
         *      conf.sides[side].sublat
         *      cpoly
         *      qbasis
         *      conf.logI
         *      conf.J
         *
         * We recompute the lognorms for each q.
         *
         */
        lognorm_smart lognorms;

        /* fbK depends on
         *      conf.logI
         *      conf.sides[side].lim
         *      nthreads
         * and is set by prepare_for_new_q() */
        fb_factorbase::key_type fbK;

        /* *fbs is always (*fb)[fbK]. It's kept as a pointer for speedy
         * access and so that the structure remains
         * default-constructible.
         *
         * It is set by prepare_for_new_q().
         */
        fb_factorbase::slicing const * fbs = NULL;

        inline bool no_fb() const { return fbs == NULL; }

        trialdiv_data const * td;

        /* precomp_plattice_dense: caching of the FK-basis in sublat mode.
         * (for the toplevel only). This is not the same as the
         * precomp_plattice that is done for lower levels.
         *
         * This is (obviously) recomputed for each special-q in sublat
         * mode. The only reason why it's here is because we would like
         * the storage to remain allocated, adnd avoid constant
         * malloc/free.
         */
        multityped_array<precomp_plattice_dense_t, 1, FB_MAX_PARTS> precomp_plattice_dense;
        void precomp_plattice_dense_clear();

        /* This is updated by applying the special-q lattice transform to
         * the factor base. This is a "current status" that gets updated
         * as we sieve. It's initialized by small_sieve_init
         *
         * Again, this is recomputed for each special-q, and is only put
         * here as an allocation optimization.
         */
        small_sieve_data_t ssd;

        /* the "group" member is not default-constructible,
         * unfortunately.
         */
        side_data(int nr_arrays)
            : group(nr_arrays)
        {
        }

        template <int LEVEL, typename HINT> void reset_all_pointers();

        template <int LEVEL, typename HINT>
            bucket_array_t<LEVEL, HINT> &
            reserve_BA(int wish) {
                return group.get<LEVEL, HINT>().reserve(wish);
            }

        template <int LEVEL, typename HINT>
            int rank_BA(bucket_array_t<LEVEL, HINT> const & BA) {
                return group.get<LEVEL, HINT>().rank(BA);
            }

        template <int LEVEL, typename HINT>
            void
            release_BA(bucket_array_t<LEVEL, HINT> &BA) {
                return group.get<LEVEL, HINT>().release(BA);
        }

        /*
         * not even needed. Better to expose only reserve() and release()
         template <int LEVEL, typename HINT>
         std::vector<bucket_array_t<LEVEL, HINT>> &
         bucket_arrays() {return group.get<LEVEL, HINT>().bucket_arrays();}
         */

        template <int LEVEL, typename HINT>
            std::vector<bucket_array_t<LEVEL, HINT>> const &
            bucket_arrays() const {return group.cget<LEVEL, HINT>().bucket_arrays();}

        dumpfile_t dumpfile;
    };

    std::array<side_data, 2> sides;

    /* All of this exists _for each thread_ */
    struct thread_data {
        struct side_data {
            /* The real array where we apply the sieve.
             * This has size BUCKET_REGION_0 and should be close to L1
             * cache size. */
            unsigned char *bucket_region = NULL;
        };

        nfs_work &ws;  /* a pointer to the parent structure, really */
        std::array<side_data, 2> sides;
        /* SS is used only in process_bucket region */
        unsigned char *SS = NULL;

        /* A note on SS versus sides[side].bucket_region.
         *
         * All of this in only used in process_bucket_region. However, as
         * these are pretty hammered areas, we prefer to have them
         * allocated once and for all.
         *
         * norm initialization goes to sides[side].bucket_region.  Some
         * tolerance is subtracted from these lognorms to account for
         * accepted cofactors.
         *
         * SS is the bucket region where buckets are applied, and where
         * the small sieve is done. The qualification test is SS[x] >=
         * sides[side].bucket_region[x].
         *
         * TODO: that makes three 64kb memory areas per thread, which
         * might be overkill ?
         */

        thread_data(nfs_work &);

        /* Currently we prevent thread_data to be copy-able, so that the
         * situation of what we do with the thread-private pointers on
         * the bucket region is more clear. But we could as well arrange
         * so that this is possible.
         */
        thread_data(thread_data const &) = delete;
        thread_data(thread_data &&);
        ~thread_data();
        friend class nfs_work;
        private:
        void allocate_bucket_regions();
    };

    std::vector<thread_data> th;

    nfs_work(las_info & _las);
    nfs_work(las_info & _las, int);
    private:
    void zeroinit_defaults();
    void compute_toplevel_and_buckets();        // utility
    public:
    /* This uses the same reference as this->las, except that we want it
     * non-const */
    void prepare_for_new_q(las_info &);

    void allocate_buckets(nfs_aux&, thread_pool&);
    private:
    void allocate_bucket_regions();
    public:
    void buckets_alloc();
    void buckets_free();


    private:
    template <int LEVEL, typename HINT> double buckets_max_full();

    public:
    double check_buckets_max_full();
    template <typename HINT> double check_buckets_max_full(int level);
};

/* Should it be made a shared pointer too ? Probably. */
class nfs_work_cofac {
    public:
    las_info const & las;
    siever_config sc;
    las_todo_entry doing;

    facul_strategies const * strategies;

    /* yes, the ctor takes a non-const reference, but this->las is const.
     * This is because the ctor wants to access the cache in the las
     * structure. */
    nfs_work_cofac(las_info & las, nfs_work const & ws);
};

#endif	/* LAS_THREADS_WORK_DATA_HPP_ */
