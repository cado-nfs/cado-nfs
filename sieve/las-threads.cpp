#include "cado.h" // IWYU pragma: keep

#include <cstddef>

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
#include "verbose.hpp"

class las_memory_accessor; // IWYU pragma: keep


template <typename T>
void
reservation_array_base<T>::allocate_buckets(las_memory_accessor & memory, int n_bucket, double fill_ratio, int logI, nfs_aux & aux, thread_pool & pool)
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
      pool.add_task_lambda([=,&B,&aux,&memory](worker_thread * worker){
            timetree_t & timer(aux.th[worker->rank()].timer);
            ENTER_THREAD_TIMER(timer);
#ifndef DISABLE_TIMINGS
            const timetree_t::accounting_sibling dummy(timer, tdict_slot_for_alloc_buckets);
#endif
            TIMER_CATEGORY(timer, bookkeeping());
            auto tt = timer.trace(worker->rank(), chronograms::ALLOC {});

            B.allocate_memory(memory, n_bucket, ratio / n, logI);
              }, i, thread_pool::QUEUE_MISC, cost);
      /* queue 2. Joined in nfs_work::allocate_buckets */
  }
}

template <typename T>
T & reservation_array<T, false>::inner_reserve()
{
    auto lock = super::get_lock();

    while (available_buckets.empty())
        cv.wait(lock);

    auto [ ratio, i ] = available_buckets.top();
    available_buckets.pop();

    verbose_fmt_print(0, 3, "# Bucket {} is {:.0f}% full\n",
            i, ratio * 100.);

    /* We used to have a mechanism that detected the situation where no
     * bucket array was claiming any room available. This turned out to
     * be a no-op since average_full() alwayrs returns something anyway.
     */
    return super::BAs[i];
}

template <typename T>
void reservation_array<T, false>::release(T &BA) {
    auto lock = super::get_lock();
    const double ratio = BA.average_full();
    available_buckets.emplace(ratio, super::rank(BA));
    cv.notify_one();
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
