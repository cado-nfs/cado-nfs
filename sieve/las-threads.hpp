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

class thread_workspaces;

struct thread_data;

/* All of this exists _for each thread_ */
struct thread_side_data : private NonCopyable {

  /* For small sieve */
  std::vector<spos_t> ssdpos;
  std::vector<spos_t> rsdpos;

  /* The real array where we apply the sieve.
     This has size BUCKET_REGION_0 and should be close to L1 cache size. */
  unsigned char *bucket_region = NULL;
  sieve_checksum checksum_post_sieve;

  thread_side_data();
  ~thread_side_data();

  private:
  void allocate_bucket_region();
  friend struct thread_data;
  public:
  void update_checksum(){checksum_post_sieve.update(bucket_region, BUCKET_REGION);}
};

struct thread_data : private NonCopyable {
  const thread_workspaces *ws;  /* a pointer to the parent structure, really */
  int id;
  thread_side_data sides[2];
  las_info const * plas;
  sieve_info * psi;
  las_report rep;       /* XXX obsolete, will be removed once the
                           timetree_t things absorbs everything. We still
                           need to decide what we do with the non-timer
                           (a.k.a. counter) data. */
  /* SS is used only in process_bucket region */
  unsigned char *SS;
  bool is_initialized;
  uint32_t first_region0_index;
  thread_data();
  ~thread_data();
  void init(const thread_workspaces &_ws, int id, las_info const & las);
  void pickup_si(sieve_info& si);
  void update_checksums();
};

struct thread_data_task_wrapper : public task_parameters {
    thread_data * th;
    void * (*f)(timetree_t&, thread_data *);
};

/* A set of n bucket arrays, all of the same type, and methods to reserve one
   of them for exclusive use and to release it again. */

template <typename T>
class reservation_array : public buckets_are_full::callback_base, private monitor {
  T * const BAs;
  bool * const in_use;
  const size_t n;
  condition_variable cv;
  /* Return the index of the first entry that's currently not in use, or the
     first index out of array bounds if all are in use */
  size_t find_free() const {
    return std::find(&in_use[0], &in_use[n], false) - &in_use[0];
  }
  reservation_array(reservation_array const &) = delete;
  reservation_array& operator=(reservation_array const&) = delete;
public:
  typedef typename T::update_t update_t;
  reservation_array(reservation_array &&) = default;
  reservation_array(size_t n)
    : BAs(new T[n]), in_use(new bool[n]), n(n)
  {
    for (size_t i = 0; i < n; i++) {
      in_use[i] = false;
    }
  }
  ~reservation_array() {
    delete[] BAs;
    delete[] in_use;
  }

  /* Allocate enough memory to be able to store at least n_bucket buckets,
     each of size at least fill_ratio * bucket region size. */
  void allocate_buckets(const uint32_t n_bucket, double fill_ratio, int logI);
  const T* cbegin() const {return &BAs[0];}
  const T* cend() const {return &BAs[n];}

  void reset_all_pointers() {
      for (T * it = &BAs[0]; it != &BAs[n]; it++) {
          it->reset_pointers();
      }
  }

  T &reserve();
  void release(T &BA);
  void diagnosis(int side, fb_factorbase::slicing const & fbs) const override {
      int LEVEL = T::level;
      typedef typename T::hint_type HINT;
      verbose_output_print(2, 3, "# overfull diagnosis for %d%c buckets on side %d (%zu arrays defined)\n",
              LEVEL, HINT::rtti[0], side, n);
      for(const T * t = cbegin() ; t != cend() ; ++t) {
          /* Tell which slices have been processed using this array
           * exactly */
          t->diagnosis(side, t - cbegin(), fbs);
      }
    }
};

/* A group of reservation arrays, one for each possible update type.
   Also defines a getter function, templated by the desired type of
   update, that returns the corresponding reservation array, i.e.,
   it provides a type -> object mapping. */
class reservation_group {
  friend class thread_workspaces;
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
  reservation_group(size_t nr_bucket_arrays);
  void allocate_buckets(const uint32_t *n_bucket,
          bkmult_specifier const& multiplier,
          std::array<double, FB_MAX_PARTS> const &
          fill_ratio, int logI);
  void diagnosis(int side, fb_factorbase::slicing const & fbs) const {
       RA1_short.diagnosis(side, fbs);
       RA2_short.diagnosis(side, fbs);
       RA3_short.diagnosis(side, fbs);
       RA1_long.diagnosis(side, fbs);
       RA2_long.diagnosis(side, fbs);
  }
};

class thread_workspaces : private NonCopyable {
  const size_t nr_workspaces;
  sieve_info * psi;
  const unsigned int nr_sides; /* Usually 2 */
  reservation_group groups[2]; /* one per side */

public:
  // FIXME: thrs should be private!
  thread_data *thrs;

  thread_workspaces(size_t nr_workspaces, unsigned int nr_sides, las_info& _las);
  ~thread_workspaces();
  void pickup_si(sieve_info& si);
  void thread_do_using_pool(thread_pool&, void * (*) (timetree_t&, thread_data *));
  // void thread_do(void * (*) (thread_data *));
  void buckets_alloc();
  void buckets_free();
  template <int LEVEL, typename HINT>
  double buckets_max_full();

  void accumulate_and_clear(las_report_ptr, sieve_checksum *);

  template <int LEVEL, typename HINT>
  void reset_all_pointers(int side);

  template <int LEVEL, typename HINT>
  bucket_array_t<LEVEL, HINT> &
  reserve_BA(const int side) {return groups[side].get<LEVEL, HINT>().reserve();}

  template <int LEVEL, typename HINT>
  void
  release_BA(const int side, bucket_array_t<LEVEL, HINT> &BA) {
    return groups[side].get<LEVEL, HINT>().release(BA);
  }

  /* Iterator over all the bucket arrays of a given type on a given side */
  template <int LEVEL, typename HINT>
  const bucket_array_t<LEVEL, HINT> *
  cbegin_BA(const int side) const {return groups[side].cget<LEVEL, HINT>().cbegin();}

  template <int LEVEL, typename HINT>
  const bucket_array_t<LEVEL, HINT> *
  cend_BA(const int side) const {return groups[side].cget<LEVEL, HINT>().cend();}

  /* Realistically, this can't be useful. Or maybe only with no
   * middle-primes being sieved. The reason for that is that the info we
   * want to print has to be printed after fill_in_buckets, but before
   * apply_buckets. And with middle primes, these get interleaved */
  void diagnosis(std::array<fb_factorbase::slicing const *, 2> fbs) const {
      groups[0].diagnosis(0, *fbs[0]);
      groups[1].diagnosis(1, *fbs[1]);
  }
};

#endif
