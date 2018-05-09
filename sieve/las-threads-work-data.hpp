#ifndef LAS_THREADS_WORK_DATA_HPP_
#define LAS_THREADS_WORK_DATA_HPP_

#include "las-forwardtypes.hpp"
#include "las-types.hpp"
#include "las-threads.hpp"
#include "ecm/facul.hpp"

/*
 * This structure holds the key algorithmic data that is used in las. It
 * is intentionally detached from the rest of the ``stats-like'' control
 * data (found in nfs_aux, defined in las-auxiliary-data.hpp)
 *
 * Two important aspects here:
 *
 *  - this structure remains the same when sieve_info changes
 *  - no concurrent access with two different sieve_info structures is
 *    possible.
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
    private:
    const int nr_workspaces;
    reservation_group groups[2]; /* one per side */

    public:
    /* All of this exists _for each thread_ */
    struct thread_data {
        struct side_data {
            /* The real array where we apply the sieve.
             * This has size BUCKET_REGION_0 and should be close to L1
             * cache size. */
            unsigned char *bucket_region = NULL;

            ~side_data();

            private:
            void allocate_bucket_region();
            friend struct thread_data;
            public:
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
        thread_data(thread_data const &);
        ~thread_data();
        void allocate_bucket_regions();
    };

    std::vector<thread_data> th;

    nfs_work(las_info const & _las);
    nfs_work(las_info const & _las, int);

    void allocate_buckets(sieve_info const & si);
    void allocate_bucket_regions();
    void buckets_alloc();
    void buckets_free();

    private:
    template <int LEVEL, typename HINT>
        double buckets_max_full();

    public:

    double check_buckets_max_full();

    template <typename HINT>
        double check_buckets_max_full(int level, HINT const & hint);

    template <int LEVEL, typename HINT> void reset_all_pointers(int side);

    template <int LEVEL, typename HINT>
        bucket_array_t<LEVEL, HINT> &
        reserve_BA(const int side, int wish) {
            return groups[side].get<LEVEL, HINT>().reserve(wish);
        }

    template <int LEVEL, typename HINT>
        int rank_BA(const int side, bucket_array_t<LEVEL, HINT> const & BA) {
            return groups[side].get<LEVEL, HINT>().rank(BA);
        }

    template <int LEVEL, typename HINT>
        void
        release_BA(const int side, bucket_array_t<LEVEL, HINT> &BA) {
            return groups[side].get<LEVEL, HINT>().release(BA);
        }

    /*
     * not even needed. Better to expose only reserve() and release()
     template <int LEVEL, typename HINT>
     std::vector<bucket_array_t<LEVEL, HINT>> &
     bucket_arrays(int side) {return groups[side].get<LEVEL, HINT>().bucket_arrays();}
     */

    template <int LEVEL, typename HINT>
        std::vector<bucket_array_t<LEVEL, HINT>> const &
        bucket_arrays(int side) const {return groups[side].cget<LEVEL, HINT>().bucket_arrays();}

    void diagnosis(int level, std::array<fb_factorbase::slicing const *, 2> fbs) const {
        groups[0].diagnosis(0, level, *fbs[0]);
        groups[1].diagnosis(1, level, *fbs[1]);
    }
};

/* Should it be made a shared pointer too ? Probably. */
class nfs_work_cofac {
    public:
    las_info const & las;
    siever_config const & sc;
    las_todo_entry doing;

    std::shared_ptr<facul_strategies_t> strategies;

    nfs_work_cofac(las_info const& las, sieve_info const & si);
};

#endif	/* LAS_THREADS_WORK_DATA_HPP_ */
