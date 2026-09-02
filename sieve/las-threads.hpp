#ifndef CADO_LAS_THREADS_HPP
#define CADO_LAS_THREADS_HPP

#include <cstddef>

#include <condition_variable>
#include <array>
#include <vector>
#include <queue>
#include <utility>
#include <mutex>
#include <tuple>

#include "bucket.hpp"
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
    mutable std::mutex my_lock;
    auto get_lock() const { return std::unique_lock(my_lock); }
    std::condition_variable cv;
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

    public:
    ~reservation_array() = default;
    reservation_array(reservation_array const &) = delete;
    reservation_array& operator=(reservation_array const&) = delete;

    /* even moves are unholy, of course. We only ever need them in places
     * where no locking is needed, and fortunately so. Because there is
     * no such thing as _moving_ a mutex, obviously. So there is ample
     * potential to shoot yourself in the foot if you use these, really.
     */
    reservation_array(reservation_array && o, std::unique_lock<std::mutex> &&) noexcept
        : super(std::move(o))
        {}

    reservation_array(reservation_array && o) noexcept
        : reservation_array(std::move(o), get_lock())
        {}

    reservation_array& operator=(reservation_array && o) noexcept {
        auto me = get_lock();
        auto them = o.get_lock();
        static_cast<super&>(*this) = std::move(o);
        return *this;
    }

    explicit reservation_array(size_t n)
        : super(n)
    {
        for(size_t i = 0 ; i < super::BAs.size() ; i++)
            available_buckets.emplace(0, i);
    }

    reservation_array() = default;

    void reset_all_pointers() {
        auto lock = get_lock();
        super::reset_all_pointers(lock);
        available_buckets = decltype(available_buckets)();
        for(size_t i = 0 ; i < super::BAs.size() ; i++)
            available_buckets.emplace(0, i);
    }

    acquired_BA reserve() { return acquired_BA(*this); }
};

/* buckets with longhints are only used in downsort, and we can use a
 * much simpler mechanism in that case.
 */
template <bucket_array_type T>
class reservation_array<T, true> : public reservation_array_base<T> {
    mutable std::mutex my_lock;
    auto get_lock() const { return std::unique_lock(my_lock); }
    static constexpr bool has_longhint_v = true;
    using super = reservation_array_base<T>;

    public:
    explicit reservation_array(size_t n) : super(n) { }

    ~reservation_array() = default;
    reservation_array(reservation_array const &) = delete;
    reservation_array& operator=(reservation_array const&) = delete;
    reservation_array(reservation_array && o, std::unique_lock<std::mutex> &&) noexcept
        : super(std::move(o))
        {}

    reservation_array(reservation_array && o) noexcept
        : reservation_array(std::move(o), get_lock())
        {}

    reservation_array& operator=(reservation_array && o) noexcept {
        auto me = get_lock();
        auto them = o.get_lock();
        static_cast<super&>(*this) = std::move(o);
        return *this;
    }

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
    using target_array_t = reservation_array<bucket_array_t<LEVEL, HINT>>;

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

    template <std::size_t... Is>
        reservation_group(int nr_bucket_arrays, std::index_sequence<Is...>)
            : RAs(((void)Is, nr_bucket_arrays)...) {}

public:
    template <int LEVEL, typename HINT>
        [[nodiscard]] auto& get() {
            return std::get<target_array_t<LEVEL, HINT>>(RAs);
        }

    template <int LEVEL, typename HINT>
        [[nodiscard]] auto const& get() const {
            return std::get<target_array_t<LEVEL, HINT>>(RAs);
        }

public:
    /* Reserve the required number of bucket arrays. For shorthint BAs, we
     * need at least as many as there are threads filling them (or more, for
     * balancing). This is controlled by the nr_workspaces field in
     * nfs_work.  For longhint, the parallelization scheme is a bit
     * different, hence we specify directly here the number of threads that
     * will fill these bucket arrays by downsosrting. Older code had that
     * downsorting single-threaded.
     */

    /* call the private ctor to initialize all RA members */
    explicit reservation_group(int nr_bucket_arrays)
        : reservation_group(
              nr_bucket_arrays,
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
                get<1, shorthint_t>().slice_statistics(side, fbs);
                get<1, longhint_t>().slice_statistics(side, fbs);
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
