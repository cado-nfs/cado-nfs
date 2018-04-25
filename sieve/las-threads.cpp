#include "cado.h"
#include "memory.h"
#include "las-threads.hpp"
#include "las-types.hpp"
#include "las-config.h"

void thread_side_data::allocate_bucket_region()
{
  /* Allocate memory for each side's bucket region. Our intention is to
   * avoid doing it in case we have no factor base. Not that much because
   * of the spared memory, but rather because it's a useful way to trim
   * the irrelevant parts of the code in that case.
   */
  if (!bucket_region)
  bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
}

thread_side_data::thread_side_data() { }

thread_side_data::~thread_side_data()
{
  if (bucket_region) contiguous_free(bucket_region);
  bucket_region = NULL;
}

thread_data::thread_data() : is_initialized(false)
{
  /* Allocate memory for the intermediate sum (only one for both sides) */
  SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
  las_report_init(rep);
}

thread_data::~thread_data()
{
  ASSERT_ALWAYS(is_initialized);
  ASSERT_ALWAYS(SS != NULL);
  contiguous_free(SS);
  SS = NULL;
  las_report_clear(rep);
}

void thread_data::init(const thread_workspaces &_ws, const int _id, las_info const & _las)
{
  ASSERT_ALWAYS(!is_initialized);
  ws = &_ws;
  id = _id;
  plas = & _las;
  is_initialized = true;
}

void thread_data::pickup_si(sieve_info & _si)
{
  psi = & _si;
  for (int side = 0 ; side < 2 ; side++)
      sides[side].allocate_bucket_region();
  /*
  sieve_info & si(*psi);
  for (int side = 0 ; side < 2 ; side++) {
      if (si.sides[side].fb)
          sides[side].set_fbs(si.sides[side].fb.get());
  }
  */
}

void thread_data::update_checksums()
{
  for(int s = 0 ; s < 2 ; s++)
    sides[s].update_checksum();
}

template <typename T>
void
reservation_array<T>::allocate_buckets(const uint32_t n_bucket, const double fill_ratio, int logI)
{
  /* We estimate that the updates will be evenly distributed among the n
     different bucket arrays, so each gets fill_ratio / n.
     However, for a large number of threads, we need a bit of margin.
     In principle, one should check that the number of threads asked by the user
     is not too large compared to the number of slices (i.e. the size of the
     factor bases).
     */
  double ratio = fill_ratio;
  size_t n = BAs.size();
  for (size_t i = 0; i < n; i++)
    BAs[i].allocate_memory(n_bucket, ratio / n, logI);
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
      in_use[wish]=true;
      leave();
      return BAs[wish];
  }

  while ((i = find_free()) == n)
    wait(cv);

  if (choose_least_full) {
    /* Find the least-full bucket array */
    double least_full = 1.;
    double most_full = 0;
    /* We want indices too. For most full buckets, we even want to report
     * the exact index of the overflowing buckets */
    size_t least_full_index = n;
    std::pair<size_t, unsigned int> most_full_index;
    if (verbose)
      verbose_output_print(0, 3, "# Looking for least full bucket\n");
    for (size_t j = 0; j < n; j++) {
      if (in_use[j])
        continue;
      double full = BAs[j].average_full();
      if (verbose)
        verbose_output_print(0, 3, "# Bucket %zu is %.0f%% full\n",
                             j, full * 100.);
      if (full > most_full) {
          most_full = full;
          most_full_index = std::make_pair(j, 0);
      }
      if (full < least_full) {
        least_full = full;
        least_full_index = j;
      }
    }
    if (least_full_index == n || least_full >= 1.) {
        verbose_output_print(0, 1, "# Error: buckets are full, throwing exception\n");
        ASSERT_ALWAYS(most_full > 1);
        size_t j = most_full_index.first;
        unsigned int i = most_full_index.second;
        /* important ! */
        leave();
        /* We don't know which side we are, yet, but our caller knows */
        throw buckets_are_full(this, -1, bkmult_specifier::getkey<typename T::update_t>(), i,
                BAs[j].nb_of_updates(i), BAs[j].room_allocated_for_updates(i));
    }
    i = least_full_index;
  }
  in_use[i] = true;
  leave();
  return BAs[i];
}

template <typename T>
void reservation_array<T>::release(T &BA) {
  enter();
  ASSERT_ALWAYS(&BA >= &BAs.front());
  ASSERT_ALWAYS(&BA < &BAs[BAs.size()]);
  size_t i = &BA - &BAs[0];
  in_use[i] = false;
  signal(cv);
  leave();
}

template class reservation_array<bucket_array_t<1, shorthint_t> >;
template class reservation_array<bucket_array_t<2, shorthint_t> >;
template class reservation_array<bucket_array_t<3, shorthint_t> >;
template class reservation_array<bucket_array_t<1, longhint_t> >;
template class reservation_array<bucket_array_t<2, longhint_t> >;

/* Reserve the required number of bucket arrays. For shorthint BAs, we
 * need at least as many as there are threads filling them (or more, for
 * balancing). This is controlled by the nr_workspaces field in
 * thread_workspaces.  For longhint, the parallelization scheme is a bit
 * different, hence we specify directly here the number of threads that
 * will fill these bucket arrays by downsosrting. Older code had that
 * downsorting single-threaded.
 */
reservation_group::reservation_group(int nr_bucket_arrays, int)
  : RA1_short(nr_bucket_arrays),
    RA2_short(nr_bucket_arrays),
    RA3_short(nr_bucket_arrays),
    /* currently the parallel downsort imposes restrictions on the number
     * of bucket arrays we must have here and there. In particular #2s ==
     * #1l.
     */
    RA1_long(nr_bucket_arrays),
    RA2_long(nr_bucket_arrays)
{
}

void
reservation_group::allocate_buckets(const uint32_t *n_bucket,
        bkmult_specifier const& mult,
        std::array<double, FB_MAX_PARTS> const & fill_ratio, int logI)
{
  /* Short hint updates are generated only by fill_in_buckets(), so each BA
     gets filled only by its respective FB part */
  typedef typename decltype(RA1_short)::update_t T1s;
  typedef typename decltype(RA2_short)::update_t T2s;
  typedef typename decltype(RA3_short)::update_t T3s;
  typedef typename decltype(RA1_long)::update_t T1l;
  typedef typename decltype(RA2_long)::update_t T2l;
  RA1_short.allocate_buckets(n_bucket[1], mult.get<T1s>()*fill_ratio[1], logI);
  RA2_short.allocate_buckets(n_bucket[2], mult.get<T2s>()*fill_ratio[2], logI);
  RA3_short.allocate_buckets(n_bucket[3], mult.get<T3s>()*fill_ratio[3], logI);

  /* Long hint bucket arrays get filled by downsorting. The level-2 longhint
     array gets the shorthint updates from level 3 sieving, and the level-1
     longhint array gets the shorthint updates from level 2 sieving as well
     as the previously downsorted longhint updates from level 3 sieving. */
  RA1_long.allocate_buckets(n_bucket[1], mult.get<T1l>()*(fill_ratio[2] + fill_ratio[3]), logI);
  RA2_long.allocate_buckets(n_bucket[2], mult.get<T2l>()*fill_ratio[3], logI);
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



thread_workspaces::thread_workspaces(int n, int _nr_sides, las_info & las)
  : nb_threads(n),
    nr_workspaces(n + 2),
    groups { {nr_workspaces, nb_threads}, {nr_workspaces, nb_threads} }
{
    /* Well, groups is an array of side 2 anyway... */
    ASSERT_ALWAYS(_nr_sides == 2);

    thrs = new thread_data[nb_threads];
    ASSERT_ALWAYS(thrs != NULL);

    for(int i = 0 ; i < nb_threads; i++) {
        thrs[i].init(*this, i, las);
    }
}

thread_workspaces::~thread_workspaces()
{
    delete[] thrs;
}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void
thread_workspaces::pickup_si(sieve_info & _si)
{
    psi = & _si;
    sieve_info & si(*psi);
    for (int i = 0; i < nb_threads; ++i) {
        thrs[i].pickup_si(_si);
    }
    /* Always allocate the max number of buckets (i.e., as if we were using the
       max value for J), even if we use a smaller J due to a poor q-lattice
       basis */
    /* Take some margin depending on parameters */
    /* Multithreading perturbates the fill-in ratio */
    bkmult_specifier const & multiplier = * si.bk_multiplier;
    verbose_output_print(0, 2, "# Reserving buckets with a multiplier of %s\n",
            multiplier.print_all().c_str());
    for (unsigned int i_side = 0; i_side < nr_sides; i_side++) {
        if (!_si.sides[i_side].fb) continue;
        groups[i_side].allocate_buckets(si.nb_buckets,
                multiplier,
                si.sides[i_side].fbs->stats.weight,
                si.conf.logI);
    }
}

task_result * thread_do_task_wrapper(const worker_thread * worker, const task_parameters *_param)
{
    const thread_data_task_wrapper *th = static_cast<const thread_data_task_wrapper *>(_param);
    (*th->f)(worker->timer, th->th);
    return new task_result;
}

void
thread_workspaces::thread_do_using_pool(thread_pool &pool, void * (*f)(timetree_t&, thread_data *))
{
    thread_data_task_wrapper * ths = new thread_data_task_wrapper[nb_threads];
    for (int i = 0; i < nb_threads; ++i) {
        ths[i].th = &thrs[i];
        ths[i].f = f;
    }

    for (int i = 0; i < nb_threads; ++i) {
        pool.add_task(thread_do_task_wrapper, &ths[i], i);
    }
    for (int i = 0; i < nb_threads; ++i) {
        task_result *result = pool.get_result();
        delete result;
    }
    delete[] ths;
}

template <int LEVEL, typename HINT>
double
thread_workspaces::buckets_max_full()
{
    double mf0 = 0;
    typedef bucket_array_t<LEVEL, HINT> BA_t;
    for(unsigned int side = 0 ; side < nr_sides ; side++) {
      for (auto const & BA : bucket_arrays<LEVEL, HINT>(side)) {
        unsigned int fullest;
        const double mf = BA.max_full(&fullest);
        if (mf > mf0) mf0 = mf;
        if (mf > 1)
            throw buckets_are_full(
                    &(groups[side].cget<LEVEL,HINT>()),
                    side,
                    bkmult_specifier::getkey<typename BA_t::update_t>(),
                    fullest,
                    BA.nb_of_updates(fullest),
                    BA.room_allocated_for_updates(fullest));
      }
    }
    return mf0;
}
template double thread_workspaces::buckets_max_full<1, shorthint_t>();
template double thread_workspaces::buckets_max_full<2, shorthint_t>();
template double thread_workspaces::buckets_max_full<3, shorthint_t>();
template double thread_workspaces::buckets_max_full<1, longhint_t>();
template double thread_workspaces::buckets_max_full<2, longhint_t>();


void
thread_workspaces::accumulate_and_clear(las_report_ptr rep, sieve_checksum *checksum)
{
    for (int i = 0; i < nb_threads; ++i) {
        las_report_accumulate_and_clear(rep, thrs[i].rep);
        for (unsigned int side = 0; side < nr_sides; side++)
            checksum[side].update(thrs[i].sides[side].checksum_post_sieve);
    }
}

template <int LEVEL, typename HINT>
void
thread_workspaces::reset_all_pointers(int side) {
    groups[side].get<LEVEL, HINT>().reset_all_pointers();
}

template void thread_workspaces::reset_all_pointers<1, shorthint_t>(int);
template void thread_workspaces::reset_all_pointers<2, shorthint_t>(int);
template void thread_workspaces::reset_all_pointers<3, shorthint_t>(int);
