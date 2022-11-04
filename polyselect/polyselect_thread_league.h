#ifndef POLYSELECT_THREAD_LEAGUE_H_
#define POLYSELECT_THREAD_LEAGUE_H_

#include "polyselect_main_data.h"
#include "polyselect_primes_table.h"
#include "dllist.h"
#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

/* polyselect "leagues" are collection of "teams" of threads, all of
 * which work together on the same NUMA node. We typically have exactly
 * one league per NUMA node. Among a league, all thread teams, and thus
 * all threads, share the same list of primes, in particular. They also
 * share the same set of asynchronous jobs (if these are enabled).
 *
 */
#ifdef __cplusplus
extern "C" {
#endif


struct polyselect_thread_league_s {
    /* We do keep a link to the parent main data, but only for the
     * purpose of fetching some performance-irrelevant things in
     * read-only mode. Anything that involves non-const access has to be
     * done at the higher level, so that we know what we're doing, really
     * */
    polyselect_main_data_srcptr main;

    /* well, okay. sometimes we have to make compromises */
    // polyselect_main_data_ptr main_nonconst;

    polyselect_primes_table pt;

    /* list of async jobs. Not necessarily used. */
    struct dllist_head async_jobs;

    pthread_mutex_t lock;

    unsigned int league_index;

#ifdef HAVE_HWLOC
    hwloc_nodeset_t membind_set;
    // it's not really necessary to enforce the cpu binding here, as it's
    // rather done at the lower thread level (as opposed to thread
    // league or thread team).
    // hwloc_cpuset_t  cpubind_set;
#endif
};

typedef struct polyselect_thread_league_s polyselect_thread_league[1];
typedef struct polyselect_thread_league_s * polyselect_thread_league_ptr;
typedef const struct polyselect_thread_league_s * polyselect_thread_league_srcptr;

extern void polyselect_thread_league_init(polyselect_thread_league_ptr grp, polyselect_main_data_srcptr main, unsigned int league_index);
extern void polyselect_thread_league_clear(polyselect_thread_league_ptr league);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_THREAD_LEAGUE_H_ */
