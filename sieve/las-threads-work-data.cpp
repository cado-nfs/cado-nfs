#include "cado.h"
#include "las-types.hpp"
#include "las-threads-work-data.hpp"

void nfs_work::thread_data::side_data::allocate_bucket_region()
{
  /* Allocate memory for each side's bucket region. Our intention is to
   * avoid doing it in case we have no factor base. Not that much because
   * of the spared memory, but rather because it's a useful way to trim
   * the irrelevant parts of the code in that case.
   */
  if (!bucket_region)
  bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
}

nfs_work::thread_data::thread_data(thread_data const & o)
    : ws(o.ws), sides {o.sides}
{
    SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
    memcpy(SS, o.SS, BUCKET_REGION);
}

nfs_work::thread_data::side_data::~side_data()
{
  if (bucket_region) contiguous_free(bucket_region);
  bucket_region = NULL;
}

nfs_work::thread_data::thread_data(nfs_work & ws)
    : ws(ws)
{
  /* Allocate memory for the intermediate sum (only one for both sides) */
  SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
}

nfs_work::thread_data::~thread_data()
{
    contiguous_free(SS);
}

void nfs_work::thread_data::allocate_bucket_regions()
{
    for (auto & S : sides)
        S.allocate_bucket_region();
}

nfs_work::nfs_work(las_info const & _las)
    : nfs_work(_las, _las.nb_threads + 2)
{}
nfs_work::nfs_work(las_info const & _las, int nr_workspaces)
    : las(_las),
    nr_workspaces(nr_workspaces),
    groups { {nr_workspaces}, {nr_workspaces} },
    th(_las.nb_threads, thread_data(*this))
{ }

nfs_work_cofac::nfs_work_cofac(las_info const& las, sieve_info const & si) :
    las(las),
    sc(si.conf),
    doing(si.doing),
    strategies(si.strategies)
{}
