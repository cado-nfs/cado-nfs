#ifndef POLYSELECT_THREAD_H_
#define POLYSELECT_THREAD_H_

#include <pthread.h>
#include <gmp.h>

#include "polyselect_thread_league.h"
#include "polyselect_thread_team.h"
#include "polyselect_poly_header.h"
#include "polyselect_stats.h"
#include "dllist.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_thread_s {
    unsigned int thread_index;

#ifdef HAVE_HWLOC
    hwloc_cpuset_t cpubind_set;
#endif

    /* the ctor will initialize it to pthread_self so that it's not an
     * undefined value, but beyond that, the moment from which this
     * starts to make sense is after pthread_create, of course.
     */
    pthread_t tid;

    polyselect_thread_team_ptr team;

    /* yes, it's ugly. We definitely want that for some rare occasions */
    /* XXX Is it restricted to just grabbing the main lock every now and
     * then ? */
    // polyselect_main_data_ptr main_nonconst_use_rarely;
    pthread_mutex_t * main_lock;

    enum wait_cause why_wait;
    int is_in_sync_group;
    int is_unlocked;

    unsigned int index_in_sync_zone;

    /*******************************************
     * THE FIELDS BELOW ARE INITIALIZED LATE ! *
     *                                         *
     * (see polyselect_thread_late_init)       *
     *******************************************/

    /* Are there any asynchronous tasks that this local thread in
     * particular has set aside, and that should be picked up for further
     * processing ?
     */
    struct dllist_head async_jobs;

    /* This points to a pool of available polyselect_match_info_t ' s that
     * can be reused if needed, and put into the async_jobs list.
     *
     * Whenever a new match must be processed,
     * polyselect_shash2_find_collision_multi can
     * take the opportunity to reuse data from this list, or allocate a
     * new one if it is empty.
     *
     * Conversely, the matching function has the choice to dispose the
     * structure right away (with free), or to put it aside in this list
     */
    struct dllist_head empty_job_slots;

    /* These stats receive things that are processed asynchronously.
     */
    polyselect_stats stats;

    /* reuse our private randstate from the stats structure. */

    gmp_randstate_ptr rstate;
};
typedef struct polyselect_thread_s polyselect_thread[1];
typedef struct polyselect_thread_s * polyselect_thread_ptr;
typedef const struct polyselect_thread_s * polyselect_thread_srcptr;

extern void polyselect_thread_init(polyselect_thread_ptr thread, polyselect_thread_team_ptr team, polyselect_main_data_ptr main, unsigned int thread_index);
extern void polyselect_thread_late_init(polyselect_thread_ptr);
extern void polyselect_thread_clear(polyselect_thread_ptr);
extern void polyselect_thread_bind(polyselect_thread_ptr);

/* These functions are not attached to a polyselect_thread object in
 * particular. */
extern void polyselect_thread_chronogram_init(const char * filename);
extern void polyselect_thread_chronogram_clear();


extern void polyselect_thread_chronogram_chat(polyselect_thread_srcptr arg, const char * fmt, ...);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_THREAD_H_ */
