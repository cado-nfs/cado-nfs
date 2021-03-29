#include "cado.h" // IWYU pragma: keep

#include <memory>                  // for allocator_traits<>::value_type

#include "las-threads.hpp"

#include "macros.h"                // for ASSERT_ALWAYS

#include "las-report-stats.hpp"    // for TIMER_CATEGORY
#include "las-auxiliary-data.hpp"         // for nfs_aux
#include "tdict.hpp"               // for timetree_t, slot
#include "verbose.h"    // verbose_output_print
class las_memory_accessor; // IWYU pragma: keep


template <typename T>
void
reservation_array<T>::allocate_buckets(las_memory_accessor & memory, int n_bucket, double fill_ratio, int logI, nfs_aux & aux, thread_pool & pool)
{
    if (n_bucket <= 0) return;

  /* We estimate that the updates will be evenly distributed among the n
     different bucket arrays, so each gets fill_ratio / n.
     However, for a large number of threads, we need a bit of margin.
     In principle, one should check that the number of threads asked by the user
     is not too large compared to the number of slices (i.e. the size of the
     factor bases).
     */
  double ratio = fill_ratio;

  size_t n = BAs.size();
  for (size_t i = 0; i < n; i++) {
      auto & B(BAs[i]);
      /* Arrange so that the largest allocations are done first ! */
      double cost = ratio/n * BUCKET_REGIONS[T::level] * n_bucket * sizeof(typename T::update_t);
      pool.add_task_lambda([=,&B,&aux,&memory](worker_thread * worker,int){
            timetree_t & timer(aux.th[worker->rank()].timer);
            ENTER_THREAD_TIMER(timer);
#ifndef DISABLE_TIMINGS
            timetree_t::accounting_sibling dummy(timer, tdict_slot_for_alloc_buckets);
#endif
            TIMER_CATEGORY(timer, bookkeeping());
            B.allocate_memory(memory, n_bucket, ratio / n, logI);
              }, i, 2, cost);
      /* queue 2. Joined in nfs_work::allocate_buckets */
  }
}

template <typename T>
T &reservation_array<T>::reserve(int wish)
{
  enter();
  const bool verbose = false;
  const bool choose_least_full = true;
  size_t i;

  size_t n = BAs.size();

  if (wish >= 0) {
      ASSERT_ALWAYS(!in_use[wish]);
      i=wish;
      goto happy;
  }

  while ((i = find_free()) == n)
    wait(cv);

  if (choose_least_full) {
    /* Find the least-full bucket array. A bucket array that has one, or
     * maybe several full buckets, but isn't full on average may still be
     * used. We'll prefer the least full bucket arrays anyway.
     */
    if (verbose)
      verbose_output_print(0, 3, "# Looking for least full bucket array\n");
    double least_full = 1000; /* any large value */
    int least_full_index = -1;
    for (i = 0; i < n; i++) {
      if (in_use[i])
        continue;
      double full = BAs[i].average_full();
      if (full < least_full) {
          least_full = full;
          least_full_index = i;
      }
    }
    if (least_full_index >= 0) {
        if (verbose)
            verbose_output_print(0, 3, "# Bucket %d is %.0f%% full\n",
                    least_full_index, least_full * 100.);
        i = least_full_index;
        goto happy;
    }
    /*
     * Now all bucket arrays are full on average. We're going to scream
     * and throw an exception. Now where we ought to go from here is not
     * really decided upon based on our analysis, but rather on the check
     * that is done in check_buckets_max_full. Here we'll just throw a
     * mostly phony exception that will maybe be caught and acted upon,
     * or maybe not.
     */
    leave(); /* important ! */

    auto k = bkmult_specifier::getkey<typename T::update_t>();
    verbose_output_print(0, 1, "# Error: %s buckets are full (least avg %f), throwing exception\n",
            bkmult_specifier::printkey(k).c_str(),
            least_full);
    throw buckets_are_full(k, -1, least_full * 1e6, 1 * 1e6); 
  }
happy:
  in_use[i] = true;
  leave();
  return BAs[i];
}

template <typename T>
void reservation_array<T>::release(T &BA) {
  enter();
  ASSERT_ALWAYS(&BA >= BAs.data());
  ASSERT_ALWAYS(&BA < BAs.data() + BAs.size());
  size_t i = &BA - BAs.data();
  in_use[i] = false;
  signal(cv);
  leave();
}

/* Reserve the required number of bucket arrays. For shorthint BAs, we
 * need at least as many as there are threads filling them (or more, for
 * balancing). This is controlled by the nr_workspaces field in
 * nfs_work.  For longhint, the parallelization scheme is a bit
 * different, hence we specify directly here the number of threads that
 * will fill these bucket arrays by downsosrting. Older code had that
 * downsorting single-threaded.
 */
reservation_group::reservation_group(int nr_bucket_arrays)
  : RA1_short(nr_bucket_arrays),
    RA2_short(nr_bucket_arrays),
    RA3_short(nr_bucket_arrays),
    /* currently the parallel downsort imposes restrictions on the number
     * of bucket arrays we must have here and there. In particular #2s ==
     * #1l.
     */
    RA1_long(nr_bucket_arrays),
    RA2_long(nr_bucket_arrays),
    RA1_empty(nr_bucket_arrays),
    RA2_empty(nr_bucket_arrays),
    RA3_empty(nr_bucket_arrays),
    RA1_logp(nr_bucket_arrays),
    RA2_logp(nr_bucket_arrays)
{
}


/* TODO: we may expose two distinct runtime functions that trigger
 * allocation either on the bare or non-bare bucket arrays.
 */
void
reservation_group::allocate_buckets(
        las_memory_accessor & memory,
        const int *n_bucket,
        bkmult_specifier const& mult,
        std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI,
        nfs_aux & aux,
        thread_pool & pool,
        bool with_hints)
{
  /* Short hint updates are generated only by fill_in_buckets(), so each BA
     gets filled only by its respective FB part */
  typedef typename decltype(RA1_short)::update_t T1s;
  typedef typename decltype(RA2_short)::update_t T2s;
  typedef typename decltype(RA3_short)::update_t T3s;
  typedef typename decltype(RA1_long)::update_t T1l;
  typedef typename decltype(RA2_long)::update_t T2l;

  /* We use the same multiplier definitions for both "with" and "without
   * hints".
   */
  if (with_hints) {
      RA1_short.allocate_buckets(memory, n_bucket[1], mult.get<T1s>()*fill_ratio[1], logI, aux, pool);
      RA2_short.allocate_buckets(memory, n_bucket[2], mult.get<T2s>()*fill_ratio[2], logI, aux, pool);
      RA3_short.allocate_buckets(memory, n_bucket[3], mult.get<T3s>()*fill_ratio[3], logI, aux, pool);

      /* Long hint bucket arrays get filled by downsorting. The level-2 longhint
         array gets the shorthint updates from level 3 sieving, and the level-1
         longhint array gets the shorthint updates from level 2 sieving as well
         as the previously downsorted longhint updates from level 3 sieving. */
      RA1_long.allocate_buckets(memory, n_bucket[1], mult.get<T1l>()*(fill_ratio[2] + fill_ratio[3]), logI, aux, pool);
      RA2_long.allocate_buckets(memory, n_bucket[2], mult.get<T2l>()*fill_ratio[3], logI, aux, pool);
  } else {
      RA1_empty.allocate_buckets(memory, n_bucket[1], mult.get<T1s>()*fill_ratio[1], logI, aux, pool);
      RA2_empty.allocate_buckets(memory, n_bucket[2], mult.get<T2s>()*fill_ratio[2], logI, aux, pool);
      RA3_empty.allocate_buckets(memory, n_bucket[3], mult.get<T3s>()*fill_ratio[3], logI, aux, pool);
      RA1_logp.allocate_buckets(memory, n_bucket[1], mult.get<T1l>()*(fill_ratio[2] + fill_ratio[3]), logI, aux, pool);
      RA2_logp.allocate_buckets(memory, n_bucket[2], mult.get<T2l>()*fill_ratio[3], logI, aux, pool);
  }
}


/* 
   We want to map the desired bucket_array type to the appropriate
   reservation_array in reservation_group, which we do by explicit
   specialization. Endless copy-paste here...  maybe use a virtual
   base class and an array of base-class-pointers, which then get
   dynamic_cast to the desired return type?
*/
template<>
reservation_array<bucket_array_t<1, shorthint_t> > &
reservation_group::get()
{
  return RA1_short;
}
template<>
reservation_array<bucket_array_t<2, shorthint_t> > &
reservation_group::get()
{
  return RA2_short;
}
template<>
reservation_array<bucket_array_t<3, shorthint_t> > &
reservation_group::get() {
  return RA3_short;
}
template<>
reservation_array<bucket_array_t<1, longhint_t> > &
reservation_group::get()
{
  return RA1_long;
}
template<>
reservation_array<bucket_array_t<2, longhint_t> > &
reservation_group::get()
{
  return RA2_long;
}

/* mapping types to objects is a tricky business. The code below looks
 * like red tape. But sophisticated means to avoid it would be even
 * longer (the naming difference "get" versus "cget" is not the annoying
 * part here -- it's just here as a decoration. The real issue is that we
 * want the const getter to return a const reference) */
template<>
const reservation_array<bucket_array_t<1, shorthint_t> > &
reservation_group::cget() const
{
  return RA1_short;
}
template<>
const reservation_array<bucket_array_t<2, shorthint_t> > &
reservation_group::cget() const
{
  return RA2_short;
}
template<>
const reservation_array<bucket_array_t<3, shorthint_t> > &
reservation_group::cget() const
{
  return RA3_short;
}
template<>
const reservation_array<bucket_array_t<1, longhint_t> > &
reservation_group::cget() const
{
  return RA1_long;
}
template<>
const reservation_array<bucket_array_t<2, longhint_t> > &
reservation_group::cget() const
{
  return RA2_long;
}
template <>
const reservation_array<bucket_array_t<3, longhint_t> > &
reservation_group::cget<3, longhint_t>() const
{
    ASSERT_ALWAYS(0);
}


/* And ditto for empty or near-empty hints */
template<>
reservation_array<bucket_array_t<1, emptyhint_t> > &
reservation_group::get()
{
  return RA1_empty;
}
template<>
reservation_array<bucket_array_t<2, emptyhint_t> > &
reservation_group::get()
{
  return RA2_empty;
}
template<>
reservation_array<bucket_array_t<3, emptyhint_t> > &
reservation_group::get() {
  return RA3_empty;
}
template<>
reservation_array<bucket_array_t<1, logphint_t> > &
reservation_group::get()
{
  return RA1_logp;
}
template<>
reservation_array<bucket_array_t<2, logphint_t> > &
reservation_group::get()
{
  return RA2_logp;
}
template<>
const reservation_array<bucket_array_t<1, emptyhint_t> > &
reservation_group::cget() const
{
  return RA1_empty;
}
template<>
const reservation_array<bucket_array_t<2, emptyhint_t> > &
reservation_group::cget() const
{
  return RA2_empty;
}
template<>
const reservation_array<bucket_array_t<3, emptyhint_t> > &
reservation_group::cget() const
{
  return RA3_empty;
}
template<>
const reservation_array<bucket_array_t<1, logphint_t> > &
reservation_group::cget() const
{
  return RA1_logp;
}
template<>
const reservation_array<bucket_array_t<2, logphint_t> > &
reservation_group::cget() const
{
  return RA2_logp;
}
template <>
const reservation_array<bucket_array_t<3, logphint_t> > &
reservation_group::cget<3, logphint_t>() const
{
    ASSERT_ALWAYS(0);
}

template class reservation_array<bucket_array_t<1, shorthint_t> >;
template class reservation_array<bucket_array_t<2, shorthint_t> >;
template class reservation_array<bucket_array_t<3, shorthint_t> >;
template class reservation_array<bucket_array_t<1, longhint_t> >;
template class reservation_array<bucket_array_t<2, longhint_t> >;
template class reservation_array<bucket_array_t<1, emptyhint_t> >;
template class reservation_array<bucket_array_t<2, emptyhint_t> >;
template class reservation_array<bucket_array_t<3, emptyhint_t> >;
template class reservation_array<bucket_array_t<1, logphint_t> >;
template class reservation_array<bucket_array_t<2, logphint_t> >;
