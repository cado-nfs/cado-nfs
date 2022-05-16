#ifndef POLYSELECT_THREAD_TEAM_H_
#define POLYSELECT_THREAD_TEAM_H_

#include <pthread.h>

#include "polyselect_main_data.h"
#include "polyselect_thread_league.h"
#include "polyselect_proots.h"
#include "polyselect_shash.h"
#include "barrier.h"

#define xxxDEBUG_POLYSELECT_THREADS

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_thread_s;

enum signal_cause {
    S_NONE,
    S_NEW_JOB,
    S_LEAVING_BARRIER,
    S_PLEASE_ENTER_BARRIER,
    S_LEFT_BARRIER,
    S_REACHED_BARRIER,
    S_NO_ROAMING,
};


enum wait_cause {
    W_NONE,
    W_NEW_JOB,
    W_NO_ROAMING_AND_NO_BARRIER_TAIL,
    W_LEAVING_BARRIER,
    W_REACHING_BARRIER,
};



/* threads in a given team engage **collectively** in operations,
 * notably in the collision_on_p and collision_on_sq passes.
 *
 * notice that much of what used to be defined at the thread level is now
 * defined at the team level instead.
 */
struct polyselect_thread_team_s {
    unsigned int team_index;    /* this is a _global_ index: 0 <=
                                   team_index < league->main->number_of_teams */

    int idx;
    mpz_t ad;
    polyselect_stats stats;

    /* used by several routines. Actually a replica of data which exists
     * elsewhere, mostly */
    polyselect_poly_header_t header;

    /* This is just a copy (const pointer) of stats->rstate. Since we
     * have a random state there, let's just use it...
     */
    gmp_randstate_ptr rstate;

    /* not entirely clear that this deserves here. I think that it could
     * probably live on the stack somehere in the collisions_ functions
     * now */
    polyselect_proots_t R;

    unsigned int size;          /* number of threads below */
    // unsigned int sync_busy;     /* 0 <= sync_busy <= size */

    int done;   /* when there's no further sync task to do */

    polyselect_thread_league_ptr league;

    /* We take this NON-CONST pointer for two reasons.
     *  1 - idx (which governs which ad is looked at next) is handled
     *      globally.
     *  2 - stats are committed to the main stats structure.
     */
    polyselect_main_data_ptr main_nonconst;

    pthread_mutex_t lock;

    struct {
        /* At any point inside the main loop, we have the invariant:
         *
         * count.ready + count.sync + count.async == size
         *
         * However, this invariant does not hold at the beginning of the
         * thread team setup.
         *
         * we always has roaming <= sync2 <= sync
         *
         * The only transitions that are allowed (and always with the
         * lock held) are:
         *
         * (not in group) -> ready++
         * ready-- -> async++ (when picking an async task)
         * ready-- -> sync++ (when entering a group doing sync tasks)
         * sync -> sync2++   (when entering a sync task)
         * sync2 -> roaming++ (begin_roaming)
         * roaming-- -> sync2 (end_roaming)
         * sync2-- -> sync   (end_sync_task)
         * sync-- -> ready   (leave_sync_group)
         *
         * transitions are only allowed with the team lock held and while
         * waiting on the condition that is exactly matching the current
         * state.
         *
         */
        unsigned int ready;
        unsigned int sync;
        unsigned int async;
        unsigned int sync2;
        unsigned int roaming;
        pthread_cond_t w_ready_full;
        pthread_cond_t w_ready_empty;
        pthread_cond_t w_sync_empty;
        pthread_cond_t w_job;
        pthread_cond_t w_async_empty;
        pthread_cond_t w_sync2;
        pthread_cond_t w_roaming;
        /* This reuses the team lock ! */
        barrier_t roaming_barrier;
    } count[1];

    struct {
        /* All access to this struct must be protected by the team lock */

        /* This holds the info on the precise task that is to be done by the
         * workers. Since we pass a pointer to the polyselect_thread
         * structure, the called function f has access to all relevant
         * information, and its prototype can remain quite lean.
         *
         * IMPORTANT: function f will always have the team lock held upon
         * entry. If it chooses to release it (which it should), it must
         * re-take it before exiting.
         */
        unsigned int expected;

        void (*f)(struct polyselect_thread_s *);
        /* This is passed on a case-by-case basis depending on the function.
         * Type-punning is deliberately the way to go here.
         */
        void * arg;

        /* This barrier is only used for synchronizations _within the
         * synchronous task_ (and outside lockstep)
         */
        barrier_t barrier;
    } task[1];

    /* This is used for rapid collision checking.
     *
     * We have one table that is used by each thread for filling (so it
     * could conceivably live with the thread itself), but the detection
     * works transversally, which is why we prefer to have this structure
     * at the team level
     */
    polyselect_shash_t * SH;
};

typedef struct polyselect_thread_team_s polyselect_thread_team[1];
typedef struct polyselect_thread_team_s * polyselect_thread_team_ptr;
typedef const struct polyselect_thread_team_s * polyselect_thread_team_srcptr;

extern void polyselect_thread_team_init(polyselect_thread_team_ptr team, polyselect_thread_league_ptr league, polyselect_main_data_ptr, unsigned int team_index);
extern void polyselect_thread_team_late_init(polyselect_thread_team_ptr team);
extern void polyselect_thread_team_clear(polyselect_thread_team_ptr team);
extern void polyselect_thread_team_set_idx(polyselect_thread_team_ptr team, unsigned int idx);

/* Yes, this takes the thread pointer as well, because this is posted
 * only from the leader */
extern void polyselect_thread_team_post_work(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread, void (*f)(struct polyselect_thread_s *), void * arg);
extern void polyselect_thread_team_post_work_stop(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

/* All these transitions must be called with team->lock held.
 */
//

// ready++
extern void polyselect_thread_team_i_am_ready(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

// sync -> sync2++   (when entering a sync task)
extern void polyselect_thread_team_enter_sync_task(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);
// sync2-- -> sync   (end_sync_task)
extern void polyselect_thread_team_leave_sync_task(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

// ready-- -> async++ (when picking an async task)
extern void polyselect_thread_team_enter_async(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);
extern void polyselect_thread_team_leave_async(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

// ready-- -> sync++ (when entering a group doing sync tasks)
extern void polyselect_thread_team_enter_sync_zone(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);
// sync-- -> ready   (leave_sync_zone)
extern void polyselect_thread_team_leave_sync_zone(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

// sync2 -> roaming++ (begin_roaming)
extern void polyselect_thread_team_enter_roaming(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);
// roaming-- -> sync2 (end_roaming)
extern void polyselect_thread_team_leave_roaming(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

extern void polyselect_thread_team_roaming_barrier(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread);

/* temporary helpers */
extern void cond_helper_broadcast(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread, enum signal_cause s, int ignore_async, ...);
extern void cond_helper_wait(polyselect_thread_team_ptr team, struct polyselect_thread_s * thread, enum wait_cause w, ...);




#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_THREAD_TEAM_H_ */
