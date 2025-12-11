#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <climits>

#include <array>

#include "bucket.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-report-stats.hpp"
#include "las-threads.hpp"
#include "macros.h"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "verbose.h"

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
  const double ratio = fill_ratio;

  const size_t n = BAs.size();
  for (size_t i = 0; i < n; i++) {
      auto & B(BAs[i]);
      /* Arrange so that the largest allocations are done first ! */
      const auto cost = (double) (ratio/n * BUCKET_REGIONS[T::level] * n_bucket * sizeof(typename T::update_t));
      pool.add_task_lambda([=,&B,&aux,&memory](worker_thread * worker,int){
            timetree_t & timer(aux.th[worker->rank()].timer);
            ENTER_THREAD_TIMER(timer);
#ifndef DISABLE_TIMINGS
            const timetree_t::accounting_sibling dummy(timer, tdict_slot_for_alloc_buckets);
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
  my_unique_lock u(*this);
  const bool verbose = false;
  const bool choose_least_full = true;
  size_t i;

  const size_t n = BAs.size();

  if (wish >= 0)
      return use_(wish);

  while ((i = find_free()) == n)
      wait(cv, u);

  if (choose_least_full) {
    /* Find the least-full bucket array. A bucket array that has one, or
     * maybe several full buckets, but isn't full on average may still be
     * used. We'll prefer the least full bucket arrays anyway.
     */
    if (verbose)
      verbose_fmt_print(0, 3, "# Looking for least full bucket array\n");
    double least_full = 1000; /* any large value */
    size_t least_full_index = SIZE_MAX;
    for (i = 0; i < n; i++) {
      if (in_use[i])
        continue;
      const double full = BAs[i].average_full();
      if (full < least_full) {
          least_full = full;
          least_full_index = i;
      }
    }
    if (least_full_index != SIZE_MAX) {
        if (verbose)
            verbose_fmt_print(0, 3, "# Bucket {} is {:.0f}% full\n",
                    least_full_index, least_full * 100.);
        i = least_full_index;
        return use_(i);
    }
    /*
     * Now all bucket arrays are full on average. We're going to scream
     * and throw an exception. Now where we ought to go from here is not
     * really decided upon based on our analysis, but rather on the check
     * that is done in check_buckets_max_full. Here we'll just throw a
     * mostly phony exception that will maybe be caught and acted upon,
     * or maybe not.
     */

    auto k = bkmult_specifier::getkey<typename T::update_t>();
    verbose_fmt_print(0, 1, "# Error: {} buckets are full (least avg {}), throwing exception\n",
            bkmult_specifier::printkey(k),
            least_full);
    throw buckets_are_full(k, -1, least_full * 1e6, 1 * 1e6); 
  }
  return use_(i);
}

template <typename T>
void reservation_array<T>::release(T &BA) {
    const my_unique_lock u(*this);
    ASSERT_ALWAYS(&BA >= BAs.data());
    ASSERT_ALWAYS(&BA < BAs.data() + BAs.size());
    const size_t i = &BA - BAs.data();
    in_use[i] = false;
    signal(cv);
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
  : RA1_short(nr_bucket_arrays)
  , RA1_empty(nr_bucket_arrays)
#if MAX_TOPLEVEL >= 2
    /* currently the parallel downsort imposes restrictions on the number
     * of bucket arrays we must have here and there. In particular #2s ==
     * #1l.
     */
  , RA2_short(nr_bucket_arrays)
  , RA2_empty(nr_bucket_arrays)
  , RA1_long(nr_bucket_arrays)
  , RA1_logp(nr_bucket_arrays)
#endif
#if MAX_TOPLEVEL >= 3
  , RA3_short(nr_bucket_arrays)
  , RA3_empty(nr_bucket_arrays)
  , RA2_long(nr_bucket_arrays)
  , RA2_logp(nr_bucket_arrays)
#endif
{
    static_assert(MAX_TOPLEVEL == 3);
}


template<bool with_hints>
void
reservation_group::allocate_buckets(
        las_memory_accessor & memory,
        const int *n_bucket,
        bkmult_specifier const& mult,
        std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI,
        nfs_aux & aux,
        thread_pool & pool)
{
  /* We use the same multiplier definitions for both "with" and "without
   * hints".
   */

  /* Short hint updates are generated only by fill_in_buckets(), so each BA
     gets filled only by its respective FB part */
  auto & r1s = get<1, typename hints_proxy<with_hints>::s>();
  using T1s = bucket_update_t<1, typename hints_proxy<with_hints>::s>;
  r1s.allocate_buckets(memory, n_bucket[1], mult.get<T1s>()*fill_ratio[1], logI, aux, pool);

  /* Long hint bucket arrays get filled by downsorting. The level-2
   * longhint array gets the shorthint updates from level 3 sieving,
   * and the level-1 longhint array gets the shorthint updates from
   * level 2 sieving as well as the previously downsorted longhint
   * updates from level 3 sieving. */

#if MAX_TOPLEVEL >= 2
  auto & r2s = get<2, typename hints_proxy<with_hints>::s>();
  auto & r1l = get<1, typename hints_proxy<with_hints>::l>(); 
  using T2s = bucket_update_t<2, typename hints_proxy<with_hints>::s>;
  using T1l = bucket_update_t<1, typename hints_proxy<with_hints>::l>;
  r2s.allocate_buckets(memory, n_bucket[2], mult.get<T2s>()*fill_ratio[2], logI, aux, pool);
  {
      double s = 0;
      for(int level = 2 ; level <= MAX_TOPLEVEL ; level++)
          s += fill_ratio[level];
      r1l.allocate_buckets(memory, n_bucket[1], mult.get<T1l>() * s, logI, aux, pool);
  }
#endif

#if MAX_TOPLEVEL >= 3
  auto & r3s = get<3, typename hints_proxy<with_hints>::s>();
  auto & r2l = get<2, typename hints_proxy<with_hints>::l>(); 
  using T3s = bucket_update_t<3, typename hints_proxy<with_hints>::s>;
  using T2l = bucket_update_t<2, typename hints_proxy<with_hints>::l>;
  r3s.allocate_buckets(memory, n_bucket[3], mult.get<T3s>()*fill_ratio[3], logI, aux, pool);
  {
      double s = 0;
      for(int level = 3 ; level <= MAX_TOPLEVEL ; level++)
          s += fill_ratio[level];
      r2l.allocate_buckets(memory, n_bucket[2], mult.get<T2l>() * s, logI, aux, pool);
  }
#endif
  static_assert(MAX_TOPLEVEL == 3);
}

void reservation_group::allocate_buckets(
        las_memory_accessor & memory,
        const int *n_bucket,
        bkmult_specifier const& mult,
        std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI,
        nfs_aux & aux,
        thread_pool & pool,
        bool with_hints)
{
    if (with_hints)
        allocate_buckets<true>(memory, n_bucket, mult, fill_ratio, logI, aux, pool);
    else
        allocate_buckets<false>(memory, n_bucket, mult, fill_ratio, logI, aux, pool);
}

template class reservation_array<bucket_array_t<1, shorthint_t> >;
template class reservation_array<bucket_array_t<1, emptyhint_t> >;
#if MAX_TOPLEVEL >= 2
template class reservation_array<bucket_array_t<2, shorthint_t> >;
template class reservation_array<bucket_array_t<2, emptyhint_t> >;
template class reservation_array<bucket_array_t<1, longhint_t> >;
template class reservation_array<bucket_array_t<1, logphint_t> >;
#endif

#if MAX_TOPLEVEL >= 3
template class reservation_array<bucket_array_t<3, shorthint_t> >;
template class reservation_array<bucket_array_t<3, emptyhint_t> >;
template class reservation_array<bucket_array_t<2, longhint_t> >;
template class reservation_array<bucket_array_t<2, logphint_t> >;
#endif
static_assert(MAX_TOPLEVEL == 3);
