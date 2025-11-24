#ifndef CADO_LAS_THREADS_HPP
#define CADO_LAS_THREADS_HPP

#include <cstddef>

#include <condition_variable>
#include <array>
#include <vector>

#include "bucket.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "threadpool.hpp"
#include "macros.h"

class las_memory_accessor;
class nfs_aux;

/* A set of n bucket arrays, all of the same type, and methods to reserve one
   of them for exclusive use and to release it again. */
template <typename T>
class reservation_array : public monitor {
  static_assert(T::level <= MAX_TOPLEVEL);
    /* typically, T is here bucket_array_t<LEVEL, HINT>. It's a
     * non-copy-able object. Yet, it's legit to use std::vectors's on
     * such objects in c++11, provided that we limit ourselves to the
     * right constructor, and compiled code never uses allocation
     * changing modifiers.
     */
    std::vector<T> BAs;
    std::vector<bool> in_use;
    std::condition_variable cv;
  /* Return the index of the first entry that's currently not in use, or the
     first index out of array bounds if all are in use */
  ATTRIBUTE_NODISCARD
  size_t find_free() const {
    return std::find(in_use.begin(), in_use.end(), false) - in_use.begin();
  }
  T& use_(int i) {
      ASSERT_ALWAYS(!in_use[i]);
      in_use[i]=true; 
      return BAs[i];
  }
public:
  reservation_array(reservation_array const &) = delete;
  reservation_array& operator=(reservation_array const&) = delete;
  
  /* I think that moves are ok */
  reservation_array(reservation_array &&) noexcept = default;
  reservation_array& operator=(reservation_array &&) noexcept = default;

  typedef typename T::update_t update_t;

  explicit reservation_array(size_t n)
      : BAs(n)
      , in_use(n, false)
  { }

  /* Allocate enough memory to be able to store at least n_bucket buckets,
     each of size at least fill_ratio * bucket region size. */
  void allocate_buckets(las_memory_accessor & memory, int n_bucket, double fill_ratio, int logI, nfs_aux&, thread_pool&);
  // typename std::vector<T>::const_iterator cbegin() const {return BAs.cbegin();}
  // typename std::vector<T>::const_iterator cend() const {return BAs.cend();}
  // std::vector<T>& arrays() { return BAs; }

  ATTRIBUTE_NODISCARD
  std::vector<T> const& bucket_arrays() const { return BAs; }

  ATTRIBUTE_NODISCARD
  int rank(T const & BA) const { return &BA - BAs.data(); }

  void reset_all_pointers() { for(auto & A : BAs) A.reset_pointers(); }

  T &reserve(int);
  void release(T &BA);
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
