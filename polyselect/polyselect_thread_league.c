#include "cado.h"
#include "polyselect_thread_league.h"
#include "polyselect_primes_table.h"

void polyselect_thread_league_init(polyselect_thread_league_ptr league, polyselect_main_data_srcptr main, unsigned int league_index MAYBE_UNUSED)
{
    league->main = main;
    /* We'll do our best to never use it. But it some occasions, it might
     * be quite hard. */
    // league->main_nonconst = main;
    dllist_init_head(&league->async_jobs);
#ifdef HAVE_HWLOC
    if (league->main->bind) {
        hwloc_topology_t topo = main->topology;
        int mdepth = hwloc_get_type_depth(topo, HWLOC_OBJ_NUMANODE);
        /* hwloc 2 will answer HWLOC_TYPE_DEPTH_NUMANODE but hwloc 1.x
         * will return something in the topology, or possibly UNKNOWN */
        if (mdepth == HWLOC_TYPE_DEPTH_UNKNOWN) {
            league->membind_set = NULL;
        } else {
            hwloc_obj_t node = hwloc_get_obj_by_depth(topo, mdepth, league_index);
            league->membind_set = hwloc_bitmap_alloc();
            hwloc_bitmap_copy(league->membind_set, node->nodeset);
        }
    }
#endif
    pthread_mutex_init(&league->lock, NULL);

    polyselect_primes_table_init(league->pt, main->P);
}

void polyselect_thread_league_clear(polyselect_thread_league_ptr league)
{
    pthread_mutex_destroy(&league->lock);
#ifdef HAVE_HWLOC
    if (league->main->bind && league->membind_set)
        hwloc_bitmap_free(league->membind_set);
#endif
    polyselect_primes_table_clear(league->pt);
}
