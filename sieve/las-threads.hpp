#ifndef CADO_LAS_THREADS_HPP
#define CADO_LAS_THREADS_HPP

#include <cstddef>

#include <algorithm>
#include <array>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

#include "bucket.hpp"
#include "fb.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "threadpool.hpp"
#include "macros.h"
#include "verbose.hpp"

class las_memory_accessor;
class nfs_aux;

/* A set of n bucket arrays, all of the same type, and methods to reserve one
   of them for exclusive use and to release it again. */
template <bucket_array_type T>
class reservation_array_base {
    public:
    static constexpr int level = T::level;
    using update_t = T::update_t;
    using hint_t = update_t::hint_t;
    static constexpr bool has_longhint_v = hint_t::is_long_v;

    static_assert(level <= MAX_TOPLEVEL);

    protected:
    /* typically, T is here bucket_array_t<LEVEL, HINT>. It's a
     * non-copy-able object. Yet, it's legit to use std::vectors's on
     * such objects in c++11, provided that we limit ourselves to the
     * right constructor, and compiled code never uses allocation
     * changing modifiers.
     */
    std::vector<T> BAs;

    public:

    explicit reservation_array_base(size_t n) : BAs(n) { }

    reservation_array_base(reservation_array_base const &) = delete;
    reservation_array_base& operator=(reservation_array_base const&) = delete;

    reservation_array_base(reservation_array_base &&) = default;
    reservation_array_base& operator=(reservation_array_base &&) = default;


    /* Allocate enough memory to be able to store at least n_bucket buckets,
       each of size at least fill_ratio * bucket region size. */
    void allocate_buckets(las_memory_accessor & memory, int n_bucket, double fill_ratio, int logI, nfs_aux&, thread_pool&);

    ATTRIBUTE_NODISCARD
    std::vector<T> const& bucket_arrays() const { return BAs; }

    ATTRIBUTE_NODISCARD
    size_t rank(T const & BA) const { return &BA - BAs.data(); }

    void reset_all_pointers(std::unique_lock<std::mutex> &) {
        for(auto & A : BAs) A.reset_pointers();
    }

    void slice_statistics(int side, fb_factorbase::slicing const & fbs) const {
        verbose_fmt_print(0, 2,
                "# diagnosis for {} buckets on side {} ({} arrays defined)\n",
                level, hint_t::rtti[0], side, BAs.size());
        for(auto const & A : BAs) {
            /* Tell which slices have been processed using this array
             * exactly */
            A.slice_statistics(side, &A - &BAs[0], fbs);
        }
    }
};

/* bucket arrays with shorthints are filled by competing threads, and we
 * want a priority queue so that threads pick the least full array.
 */
template <bucket_array_type T, bool has_longhint_v = T::update_t::hint_t::is_long_v>
class reservation_array;

template<bucket_array_type T>
class reservation_array<T, false> : public reservation_array_base<T> {
    static constexpr bool has_longhint_v = false;
    using super = reservation_array_base<T>;
    std::unique_ptr<std::mutex> my_lock = std::make_unique<std::mutex>();
    auto get_lock() const { return std::unique_lock(*my_lock); }
    std::unique_ptr<std::condition_variable> cv = std::make_unique<std::condition_variable>();
    using available_bucket_t = std::pair<double, size_t>;
    struct prioritize_least_full_bucket {
        /* a priority queue takes the "top" element, so the comparator
         * element C must be such that C(others, the_one_we_want) is
         * always true. Which means that it must behave as std::greater<>
         */
        bool operator()(available_bucket_t const & a, available_bucket_t const & b) const {
            return a.first > b.first;
        }
    };
    std::priority_queue<available_bucket_t ,std::vector<available_bucket_t>, prioritize_least_full_bucket> available_buckets;

    struct acquired_BA {
        reservation_array<T> & parent;
        T & BA;
        explicit acquired_BA(reservation_array<T> & parent)
            : parent(parent)
              , BA(parent.inner_reserve())
        {}
        T & access() { return BA; }
        ~acquired_BA() { parent.release(BA); }
        acquired_BA(acquired_BA const &) = delete;
        acquired_BA& operator=(acquired_BA const &) = delete;
        acquired_BA(acquired_BA &&) = delete;
        acquired_BA& operator=(acquired_BA &&) = delete;
    };

    T & inner_reserve();
    void release(T &BA);

    /* The occupancy that we record when an array is released goes stale
     * as soon as the array is emptied, so the queue must be rebuilt
     * whenever that happens.
     */
    void reset_queue() {
        available_buckets = decltype(available_buckets)();
        for(size_t i = 0 ; i < super::BAs.size() ; i++)
            available_buckets.emplace(0, i);
    }

    public:
    ~reservation_array() = default;
    reservation_array(reservation_array const &) = delete;
    reservation_array& operator=(reservation_array const&) = delete;

    reservation_array(reservation_array && o) = default;
    reservation_array& operator=(reservation_array &&) = default;

    explicit reservation_array(size_t n)
        : super(n)
    {
        for(size_t i = 0 ; i < super::BAs.size() ; i++)
            available_buckets.emplace(0, i);
    }

    reservation_array() = default;

    /* allocate_memory() resets the write pointers, so this empties all
     * arrays. It happens once per special-q.
     */
    void allocate_buckets(las_memory_accessor & memory, int n_bucket,
            double fill_ratio, int logI, nfs_aux & aux, thread_pool & pool)
    {
        super::allocate_buckets(memory, n_bucket, fill_ratio, logI, aux, pool);
        if (n_bucket <= 0) return;
        auto lock = get_lock();
        reset_queue();
    }

    void reset_all_pointers() {
        auto lock = get_lock();
        super::reset_all_pointers(lock);
        reset_queue();
    }

    acquired_BA reserve() { return acquired_BA(*this); }
};

/* buckets with longhints are only used in downsort, and we can use a
 * much simpler mechanism in that case.
 */
template <bucket_array_type T>
class reservation_array<T, true> : public reservation_array_base<T> {
    std::unique_ptr<std::mutex> my_lock = std::make_unique<std::mutex>();
    auto get_lock() const { return std::unique_lock(*my_lock); }
    static constexpr bool has_longhint_v = true;
    using super = reservation_array_base<T>;

    public:
    explicit reservation_array(size_t n) : super(n) { }

    ~reservation_array() = default;
    reservation_array(reservation_array const &) = delete;
    reservation_array& operator=(reservation_array const&) = delete;

    reservation_array(reservation_array && o) = default;
    reservation_array& operator=(reservation_array &&) = default;

    void reset_all_pointers() {
        auto lock = get_lock();
        super::reset_all_pointers(lock);
    }

    T & acquire(size_t rank) { return super::BAs[rank]; }
};

/* A group of reservation arrays, one for each possible update type.
   Also defines a getter function, templated by the desired type of
   update, that returns the corresponding reservation array, i.e.,
   it provides a type -> object mapping. */
class reservation_group {
    friend class nfs_work;
    private:
    template <int LEVEL, typename HINT>
    using target_array_t =
        std::vector<reservation_array<bucket_array_t<LEVEL, HINT>>>;

    static_assert(MAX_TOPLEVEL <= 3);

    using RAs_t = std::tuple<
          target_array_t<1, shorthint_t>
        , target_array_t<1, emptyhint_t>
#if MAX_TOPLEVEL >= 2
        , target_array_t<2, shorthint_t>
        , target_array_t<2, emptyhint_t>
        , target_array_t<1, longhint_t>
        , target_array_t<1, logphint_t>
#endif
#if MAX_TOPLEVEL >= 3
        , target_array_t<3, shorthint_t>
        , target_array_t<3, emptyhint_t>
        , target_array_t<2, longhint_t>
        , target_array_t<2, logphint_t>
#endif
    >;

    RAs_t RAs;

    public:

    /* the bucket_batch_size argument indicates how many level-1
     * reservation arrays are to be considered simultaneously during
     * downsorting.
     *
     * It is 1 by default, so that by
     * we simultaneously process as many level-1 regions as we can find inside a
     * level-2 region.
     * When bucket_batch_size is increased to 2 or more (it has to be a
     * power of two), several level-1 reservation_arrays are stored in the
     * reservation_group, and several level-2 buckets are processed
     * simultaneously during downsorting in order to fill them.
     *
     * This might even trickle to higher levels if
     * bucket_batch_size exceeds (1 << LOG_BUCKET_REGION_step). In that
     * case, we would downsort several level-3 buckets simultanously in
     * order to fill several level-2 reservation_arrays.
     *
     * nslots(LEVEL, bucket_batch_size) is the number of reservation
     * arrays at a given level.
     */
    static int nslots(int level, int bucket_batch_size) {
        ASSERT_ALWAYS(!(bucket_batch_size & (bucket_batch_size - 1)));
        ASSERT_ALWAYS(level > 0);
        int b = bucket_batch_size;
        for(int i = 1 ; i < level ; i++)
            b = iceildiv(b, 1 << LOG_BUCKET_REGION_step);
        return b;
    }
    template<std::size_t LEVEL>
    static int nslots(int bucket_batch_size) {
        return nslots(LEVEL, bucket_batch_size);
    }
    private:
    template <typename Vector_t>
    static auto make_vector(int bucket_batch_size, int nr_workspaces) {
        Vector_t vec;
        static constexpr int LEVEL = Vector_t::value_type::level;
        const int count = nslots<LEVEL>(bucket_batch_size);
        vec.reserve(count);
        for (int i = 0; i < count; ++i) {
            vec.emplace_back(nr_workspaces);
        }
        return vec;
    }

    template <std::size_t... Is>
    reservation_group(int bucket_batch_size, int nr_workspaces, std::index_sequence<Is...>)
        : RAs(make_vector<std::tuple_element_t<Is, RAs_t>>(bucket_batch_size, nr_workspaces)...) {}

public:
    template <int LEVEL, typename HINT>
    [[nodiscard]] auto& get_all_slots() {
        return std::get<target_array_t<LEVEL, HINT>>(RAs);
    }

    template <int LEVEL, typename HINT>
    [[nodiscard]] auto const& get_all_slots() const {
        return std::get<target_array_t<LEVEL, HINT>>(RAs);
    }

    template <int LEVEL, typename HINT>
        requires (LEVEL>1)
    [[nodiscard]] auto& get() {
        return std::get<target_array_t<LEVEL, HINT>>(RAs)[0];
    }

    template <int LEVEL, typename HINT>
        requires (LEVEL>1)
    [[nodiscard]] auto const& get() const {
        return std::get<target_array_t<LEVEL, HINT>>(RAs)[0];
    }

    /* we expect that only level-1 buckets will effectively be
     * multiplexed. So it only makes sense to accept an extra slot
     * argument for those.
     *
     * For the moment, we'll take a default parameter.
     */
    template <int LEVEL, typename HINT>
        requires (LEVEL==1)
    [[nodiscard]] auto& get(int slot) {
        return std::get<target_array_t<LEVEL, HINT>>(RAs)[slot];
    }

    template <int LEVEL, typename HINT>
        requires (LEVEL==1)
    [[nodiscard]] auto const& get(int slot) const {
        return std::get<target_array_t<LEVEL, HINT>>(RAs)[slot];
    }

public:
    /* Reserve the required number of bucket arrays. For shorthint BAs, we
     * need at least as many as there are threads filling them (or more, for
     * balancing). This is controlled by the nr_workspaces field in
     * nfs_work.  For longhint, the parallelization scheme is a bit
     * different, hence we specify directly here the number of threads that
     * will fill these bucket arrays by downsosrting. Older code had that
     * downsorting single-threaded.
     *
     * Note that a reservation group is technically a tuple of _vectors_
     * of reservation arrays, which in turn are vectors of bucket arrays
     * because of multiplexing.
     */

    /* call the private ctor to initialize all RA members */
    explicit reservation_group(int bucket_batch_size, int nr_workspaces)
        : reservation_group(
              bucket_batch_size,
              nr_workspaces,
              std::make_index_sequence<std::tuple_size_v<RAs_t>>{}
    ) {}

    void allocate_buckets(
            las_memory_accessor & memory,
            const int *n_bucket,
            bkmult_specifier const& mult,
            std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI,
            nfs_aux & aux,
            thread_pool & pool,
            bool with_hints);

    void slice_statistics(int side, int level, fb_factorbase::slicing const & fbs) const {
        switch(level) {
            case 1:
                get<1, shorthint_t>(0).slice_statistics(side, fbs);
                get<1, longhint_t>(0).slice_statistics(side, fbs);
                break;
            case 2:
                get<2, shorthint_t>().slice_statistics(side, fbs);
                get<2, longhint_t>().slice_statistics(side, fbs);
                break;
            case 3:
                get<3, shorthint_t>().slice_statistics(side, fbs);
                break;
            default:
                ASSERT_ALWAYS(0);
        }
    }

    private:
    template<bool> void allocate_buckets(
            las_memory_accessor & memory,
            const int *n_bucket,
            bkmult_specifier const& mult,
            std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI,
            nfs_aux & aux,
            thread_pool & pool);
};

extern template class reservation_array<bucket_array_t<1, shorthint_t> >;
extern template class reservation_array<bucket_array_t<1, emptyhint_t> >;

#if MAX_TOPLEVEL >= 2
extern template class reservation_array<bucket_array_t<2, shorthint_t> >;
extern template class reservation_array<bucket_array_t<2, emptyhint_t> >;
extern template class reservation_array<bucket_array_t<1, longhint_t> >;
extern template class reservation_array<bucket_array_t<1, logphint_t> >;
#endif

#if MAX_TOPLEVEL >= 3
extern template class reservation_array<bucket_array_t<3, shorthint_t> >;
extern template class reservation_array<bucket_array_t<3, emptyhint_t> >;
extern template class reservation_array<bucket_array_t<2, longhint_t> >;
extern template class reservation_array<bucket_array_t<2, logphint_t> >;
#endif

static_assert(MAX_TOPLEVEL == 3);

#endif
