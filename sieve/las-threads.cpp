#include "cado.h"
#include "memory.h"
#include "las-threads.h"
#include "las-types.h"
#include "las-config.h"

thread_side_data::thread_side_data()
{
  /* Allocate memory for each side's bucket region */
  bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
}

thread_side_data::~thread_side_data()
{
  contiguous_free(bucket_region);
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

void thread_data::init(const thread_workspaces &_ws, const int _id, las_info_srcptr _las)
{
  ASSERT_ALWAYS(!is_initialized);
  ws = &_ws;
  id = _id;
  las = _las;
  is_initialized = true;
}

void thread_data::pickup_si(sieve_info_ptr _si)
{
  si = _si;
  for (int side = 0 ; side < si->cpoly->nb_polys ; side++) {
    sides[side].set_fb(si->sides[side]->fb);
  }
}

void thread_data::update_checksums()
{
  for(int s = 0 ; s < si->cpoly->nb_polys ; s++)
    sides[s].update_checksum();
}

template <typename T>
void
reservation_array<T>::allocate_buckets(const uint32_t n_bucket, const double fill_ratio)
{
  /* We estimate that the updates will be evenly distrubuted among the n
     different bucket arrays, so each gets fill_ratio / n.
     However, for a large number of threads, we need a bit of margin.
     In principle, one should check that the number of threads asked by the user
     is not too large compared to the number of slices (i.e. the size of the
     factor bases).
     */
  double ratio = fill_ratio;
  if (n > 1) {
      ratio *= 1.2;
  } else if (n > 8) {
      ratio *= 1.3;
  }
  for (size_t i = 0; i < n; i++)
    BAs[i].allocate_memory(n_bucket, ratio / n);
}

template <typename T>
T &reservation_array<T>::reserve()
{
  enter();
  const bool verbose = false;
  const bool choose_least_full = true;
  size_t i;

  while ((i = find_free()) == n)
    wait(cv);

  if (choose_least_full) {
    /* Find the least-full bucket array */
    size_t best_i = n;
    double best_full = 1.;
    if (verbose)
      verbose_output_print(0, 3, "# Looking for least full bucket\n");
    for (size_t j = 0; j < n; j++) {
      if (in_use[j])
        continue;
      double full = BAs[j].max_full();
      if (verbose)
        verbose_output_print(0, 3, "# Bucket %zu is %.0f%% full\n",
                             j, full * 100.);
      if (full < best_full) {
        best_full = full;
        best_i = j;
      }
    }
    if (!(best_i != n && best_full < 1.)) {
        fprintf(stderr, "Error: buckets are full!\n"
                "This may occurr if you have set too many threads compared to the sizes of the factor bases.\n"
                "Please try again with less threads, and report bug if the problem is still there.\n");
        ASSERT_ALWAYS(0);
    }
    i = best_i;
  }
  in_use[i] = true;
  leave();
  return BAs[i];
}

template <typename T>
void reservation_array<T>::release(T &BA) {
  enter();
  size_t i;
  for (i = 0; i < n; i++) {
    /* Compare object identity (=address), not contents */
    if (&BAs[i] == &BA)
      break;
  }
  ASSERT_ALWAYS(i < n);
  in_use[i] = false;
  signal(cv);
  leave();
}

template class reservation_array<bucket_array_t<1, shorthint_t> >;
template class reservation_array<bucket_array_t<2, shorthint_t> >;
template class reservation_array<bucket_array_t<3, shorthint_t> >;
template class reservation_array<bucket_array_t<1, longhint_t> >;
template class reservation_array<bucket_array_t<2, longhint_t> >;

/* Reserve the required number of bucket arrays. For shorthint BAs, we need 
   as many as there are threads filling them (or more, for balancing).
   For longhint, we need only one each, as those are filled only by
   downsorting, and all the downsorting of one level n bucket into a level
   n-1 bucket array is done as a single task, by a single thread. */
reservation_group::reservation_group(const size_t nr_bucket_arrays)
  : RA1_short(nr_bucket_arrays),
    RA2_short(nr_bucket_arrays),
    RA3_short(nr_bucket_arrays),
    RA1_long(1),
    RA2_long(1)
{
}

void
reservation_group::allocate_buckets(const uint32_t *n_bucket, const double *fill_ratio)
{
  /* Short hint updates are generated only by fill_in_buckets(), so each BA
     gets filled only by its respective FB part */
  RA1_short.allocate_buckets(n_bucket[1], fill_ratio[1]);
  RA2_short.allocate_buckets(n_bucket[2], fill_ratio[2]);
  RA3_short.allocate_buckets(n_bucket[3], fill_ratio[3]);

  /* Long hint bucket arrays get filled by downsorting. The level-2 longhint
     array gets the shorthint updates from level 3 sieving, and the level-1
     longhint array gets the shorthint updates from level 2 sieving as well
     as the previously downsorted longhint updates from level 3 sieving. */
  RA1_long.allocate_buckets(n_bucket[1], fill_ratio[2] + fill_ratio[3]);
  RA2_long.allocate_buckets(n_bucket[2], fill_ratio[3]);
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

/* Fricken C++. There has to be a better way of handling const getter methods
   than duplicating all the code. */
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



thread_workspaces::thread_workspaces(const size_t _nr_workspaces,
  const unsigned int _nr_sides, las_info_ptr las)
  : nr_workspaces(_nr_workspaces), nr_sides(_nr_sides)
{
    thrs = new thread_data[_nr_workspaces];
    ASSERT_ALWAYS(thrs != NULL);
    groups = new reservation_group_ptr[_nr_sides];
    ASSERT_ALWAYS(groups != NULL);

    for(size_t i = 0 ; i < nr_workspaces; i++) {
        thrs[i].init(*this, i, las);
    }
    for (unsigned int i = 0; i < nr_sides; i++)
      groups[i] = new reservation_group(_nr_workspaces);
}

thread_workspaces::~thread_workspaces()
{
    delete[] thrs;
    for (unsigned int i = 0; i < nr_sides; i++)
      delete groups[i];
    delete[] groups;
}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void
thread_workspaces::pickup_si(sieve_info_ptr _si)
{
    si = _si;
    for (size_t i = 0; i < nr_workspaces; ++i) {
        thrs[i].pickup_si(_si);
    }
    /* Always allocate the max number of buckets (i.e., as if we were using the
       max value for J), even if we use a smaller J due to a poor q-lattice
       basis */
    for (unsigned int i_side = 0; i_side < nr_sides; i_side++)
      groups[i_side]->allocate_buckets(si->nb_buckets_max, si->sides[i_side]->max_bucket_fill_ratio);
}

void
thread_workspaces::thread_do(void * (*f) (thread_data *))
{
    if (nr_workspaces == 1) {
        /* Then don't bother with pthread calls */
        (*f)(&thrs[0]);
        return;
    }
    pthread_t * th = (pthread_t *) malloc(nr_workspaces * sizeof(pthread_t)); 
    ASSERT_ALWAYS(th);

#if 0
    /* As a debug measure, it's possible to activate this branch instead
     * of the latter. In effect, this causes las to run in a
     * non-multithreaded way, albeit strictly following the code path of
     * the multithreaded case.
     */
    for (size_t i = 0; i < nr_workspaces; ++i) {
        (*f)(&thrs[i]);
    }
#else
    for (size_t i = 0; i < nr_workspaces; ++i) {
        int ret = pthread_create(&(th[i]), NULL, 
		(void * (*)(void *)) f,
                (void *)(&thrs[i]));
        ASSERT_ALWAYS(ret == 0);
    }
    for (size_t i = 0; i < nr_workspaces; ++i) {
        int ret = pthread_join(th[i], NULL);
        ASSERT_ALWAYS(ret == 0);
    }
#endif

    free(th);
}

template <int LEVEL, typename HINT>
double
thread_workspaces::buckets_max_full()
{
    double mf0 = 0;
    
    for(unsigned int side = 0 ; side < nr_sides ; side++) {
      for (const bucket_array_t<LEVEL, HINT> *it_BA = cbegin_BA<LEVEL, HINT>(side);
           it_BA != cend_BA<LEVEL, HINT>(side); it_BA++) {
        const double mf = it_BA->max_full();
        if (mf > mf0) mf0 = mf;
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
thread_workspaces::accumulate(las_report_ptr rep, sieve_checksum *checksum)
{
    for (size_t i = 0; i < nr_workspaces; ++i) {
        las_report_accumulate(rep, thrs[i].rep);
        for (unsigned int side = 0; side < nr_sides; side++)
            checksum[side].update(thrs[i].sides[side].checksum_post_sieve);
    }
}

template <int LEVEL, typename HINT>
void
thread_workspaces::reset_all_pointers(int side) {
    groups[side]->get<LEVEL, HINT>().reset_all_pointers();
}

template void thread_workspaces::reset_all_pointers<1, shorthint_t>(int);
template void thread_workspaces::reset_all_pointers<2, shorthint_t>(int);
template void thread_workspaces::reset_all_pointers<3, shorthint_t>(int);
