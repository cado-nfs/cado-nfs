#ifndef CADO_LAS_THREADS_HPP
#define CADO_LAS_THREADS_HPP

#include <cstddef>

#include <condition_variable>
#include <array>
#include <vector>
#include <queue>
#include <utility>
#include <mutex>
#include <type_traits>

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
    mutable std::mutex my_lock;
    protected:
    auto get_lock() const { return std::unique_lock(my_lock); }
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
    reservation_array(reservation_array const &) = delete;
    reservation_array& operator=(reservation_array const&) = delete;

    /* I think that moves are ok */
    reservation_array(reservation_array &&) noexcept = default;
    reservation_array& operator=(reservation_array &&) noexcept = default;

    explicit reservation_array(size_t n)
        : super(n)
    {
        for(size_t i = 0 ; i < super::BAs.size() ; i++)
            available_buckets.emplace(0, i);
    }

    void reset_all_pointers() {
        auto lock = super::get_lock();
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
    static constexpr bool has_longhint_v = true;
    using super = reservation_array_base<T>;

    public:
    explicit reservation_array(size_t n) : super(n) { }

    void reset_all_pointers() {
        auto lock = super::get_lock();
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
  reservation_array<bucket_array_t<1, shorthint_t> > RA1_short;
  reservation_array<bucket_array_t<1, emptyhint_t> > RA1_empty;
#if MAX_TOPLEVEL >= 2
  reservation_array<bucket_array_t<2, shorthint_t> > RA2_short;
  reservation_array<bucket_array_t<2, emptyhint_t> > RA2_empty;
  reservation_array<bucket_array_t<1, longhint_t> > RA1_long;
  reservation_array<bucket_array_t<1, logphint_t> > RA1_logp;
#endif
#if MAX_TOPLEVEL >= 3
  reservation_array<bucket_array_t<3, shorthint_t> > RA3_short;
  reservation_array<bucket_array_t<3, emptyhint_t> > RA3_empty;
  reservation_array<bucket_array_t<2, longhint_t> > RA2_long;
  reservation_array<bucket_array_t<2, logphint_t> > RA2_logp;
#endif
  static_assert(MAX_TOPLEVEL == 3);
protected:
#define RA_NAME(LEVEL, HINTTYPE) RA ## LEVEL ## _ ## HINTTYPE
#define RA_ACCESSOR(LEVEL, HINTTYPE)					\
  template<int LEV, typename HINT>					\
  reservation_array<bucket_array_t<LEV, HINT> > &			\
  get()									\
  requires (LEV == (LEVEL) && std::is_same_v<HINT, CPP_PAD(HINTTYPE, hint_t)>)\
  { return RA_NAME(LEVEL,HINTTYPE); }                                   \
  template<int LEV, typename HINT>					\
  reservation_array<bucket_array_t<LEV, HINT> > const &			\
  cget() const								\
  requires (LEV == (LEVEL) && std::is_same_v<HINT, CPP_PAD(HINTTYPE, hint_t)>)\
  { return RA_NAME(LEVEL,HINTTYPE); }

  RA_ACCESSOR(1, short)
  RA_ACCESSOR(1, empty)
#if MAX_TOPLEVEL >= 2
  RA_ACCESSOR(1, long)
  RA_ACCESSOR(1, logp)
  RA_ACCESSOR(2, short)
  RA_ACCESSOR(2, empty)
#endif
#if MAX_TOPLEVEL >= 3
  RA_ACCESSOR(2, long)
  RA_ACCESSOR(2, logp)
  RA_ACCESSOR(3, short)
  RA_ACCESSOR(3, empty)
#endif
  static_assert(MAX_TOPLEVEL == 3);
public:
  reservation_group(int nr_bucket_arrays);
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
              RA1_short.slice_statistics(side, fbs);
              RA1_long.slice_statistics(side, fbs);
              break;
          case 2:
              RA2_short.slice_statistics(side, fbs);
              RA2_long.slice_statistics(side, fbs);
              break;
          case 3:
              RA3_short.slice_statistics(side, fbs);
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
