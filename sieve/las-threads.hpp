#ifndef LAS_THREADS_HPP_
#define LAS_THREADS_HPP_

#include <pthread.h>
#include <algorithm>
#include <vector>
#include "threadpool.hpp"
#include "las-forwardtypes.hpp"
#include "bucket.hpp"
#include "fb.hpp"
#include "las-report-stats.hpp"
#include "las-base.hpp"
#include "tdict.hpp"

/* A set of n bucket arrays, all of the same type, and methods to reserve one
   of them for exclusive use and to release it again. */

template <typename T>
class reservation_array : public buckets_are_full::callback_base, private monitor {
    /* typically, T is here bucket_array<LEVEL, HINT>. It's a
     * non-copy-able object. Yet, it's legit to use std::vectors's on
     * such objects in c++11, provided that we limit ourselves to the
     * right constructor, and compiled code never uses allocation
     * changing modifiers.
     */
    std::vector<T> BAs;
    std::vector<bool> in_use;
  condition_variable cv;
  /* Return the index of the first entry that's currently not in use, or the
     first index out of array bounds if all are in use */
  size_t find_free() const {
    return std::find(in_use.begin(), in_use.end(), false) - in_use.begin();
  }
  reservation_array(reservation_array const &) = delete;
  reservation_array& operator=(reservation_array const&) = delete;
public:
  typedef typename T::update_t update_t;
  reservation_array(reservation_array &&) = default;
  reservation_array(size_t n) : BAs(n), in_use(n, false) { }

  /* Allocate enough memory to be able to store at least n_bucket buckets,
     each of size at least fill_ratio * bucket region size. */
  void allocate_buckets(int n_bucket, double fill_ratio, int logI, nfs_aux&, thread_pool&);
  // typename std::vector<T>::const_iterator cbegin() const {return BAs.cbegin();}
  // typename std::vector<T>::const_iterator cend() const {return BAs.cend();}
  // std::vector<T>& arrays() { return BAs; }
  std::vector<T> const& bucket_arrays() const { return BAs; }
  inline int rank(T const & BA) const { return &BA - &BAs.front(); }

  void reset_all_pointers() { for(auto & A : BAs) A.reset_pointers(); }

  T &reserve(int);
  void release(T &BA);
  void diagnosis(int side, fb_factorbase::slicing const & fbs) const override {
      int LEVEL = T::level;
      typedef typename T::hint_type HINT;
      verbose_output_print(0, 2, "# diagnosis for %d%c buckets on side %d (%zu arrays defined)\n",
              LEVEL, HINT::rtti[0], side, BAs.size());
      for(auto const & A : BAs) {
          /* Tell which slices have been processed using this array
           * exactly */
          A.diagnosis(side, &A - &BAs[0], fbs);
      }
    }
};

class nfs_work;

/* A group of reservation arrays, one for each possible update type.
   Also defines a getter function, templated by the desired type of
   update, that returns the corresponding reservation array, i.e.,
   it provides a type -> object mapping. */
class reservation_group {
  friend class nfs_work;
  reservation_array<bucket_array_t<1, shorthint_t> > RA1_short;
  reservation_array<bucket_array_t<2, shorthint_t> > RA2_short;
  reservation_array<bucket_array_t<3, shorthint_t> > RA3_short;
  reservation_array<bucket_array_t<1, longhint_t> > RA1_long;
  reservation_array<bucket_array_t<2, longhint_t> > RA2_long;
protected:
  template<int LEVEL, typename HINT>
  reservation_array<bucket_array_t<LEVEL, HINT> > &
  get();

  template <int LEVEL, typename HINT>
  const reservation_array<bucket_array_t<LEVEL, HINT> > &
  cget() const;
public:
  reservation_group(int nr_bucket_arrays);
  void allocate_buckets(const int *n_bucket,
          bkmult_specifier const& multiplier,
          std::array<double, FB_MAX_PARTS> const &
          fill_ratio, int logI, nfs_aux&, thread_pool&);
  void diagnosis(int side, int level, fb_factorbase::slicing const & fbs) const {
      switch(level) {
          case 1:
              RA1_short.diagnosis(side, fbs);
              RA1_long.diagnosis(side, fbs);
              break;
          case 2:
              RA2_short.diagnosis(side, fbs);
              RA2_long.diagnosis(side, fbs);
              break;
          case 3:
              RA3_short.diagnosis(side, fbs);
              break;
          default:
              ASSERT_ALWAYS(0);
      }
  }
};

#endif
