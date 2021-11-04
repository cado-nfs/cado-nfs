#include "cado.h"
#include <stdarg.h>
#include "polyselect_thread_team.h"
#include "polyselect_thread.h"
#include "barrier.h"

void polyselect_thread_team_init(polyselect_thread_team_ptr team, polyselect_thread_league_ptr league, polyselect_main_data_ptr main, unsigned int team_index)
{
    memset(team, 0, sizeof(polyselect_thread_team));

    team->team_index = team_index;
    team->size = league->main->finer_grain_threads;
    team->league = league;
    team->main_nonconst = main;
    team->done = 0;
    pthread_mutex_init(&team->lock, NULL);
    pthread_cond_init(&team->count->w_ready_full, NULL);
    pthread_cond_init(&team->count->w_ready_empty, NULL);
    pthread_cond_init(&team->count->w_sync_empty, NULL);
    pthread_cond_init(&team->count->w_job, NULL);
    pthread_cond_init(&team->count->w_async_empty, NULL);
    pthread_cond_init(&team->count->w_sync2, NULL);
    pthread_cond_init(&team->count->w_roaming, NULL);
    team->count->async = 0;
    team->count->sync = 0;
    team->count->sync2 = 0;
    team->count->roaming = 0;
    team->count->ready = 0;
    team->task->expected = 0;
    team->task->f = NULL;
    barrier_init(&team->count->roaming_barrier, &team->lock, 0);
}

/* polyselect_thread_team_late_init is called _after_ thread
 * creation and binding, so that the memory that we touch here is
 * allocated in the proper place.
 */
void polyselect_thread_team_late_init(polyselect_thread_team_ptr team)
{
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main_data = league->main;
    mpz_init(team->ad);
    polyselect_poly_header_init(team->header);
    polyselect_proots_init(team->R, main_data->d, league->pt->lenPrimes);
    polyselect_stats_init(team->stats, main_data->keep);
    team->rstate = team->stats->rstate;

    team->SH = malloc(main_data->finer_grain_threads * sizeof(polyselect_shash_t));
    size_t size_hint = POLYSELECT_SHASH_ALLOC_RATIO * league->pt->lenPrimes;
    polyselect_shash_init_multi(team->SH, size_hint, main_data->finer_grain_threads);

    //// pthread_cond_init(&team->sync_task->wait_begintask, NULL);
    //// pthread_cond_init(&team->sync_task->wait_endtask, NULL);
    //// barrier_init(&team->sync_task->barrier, &team->lock, 0);
}

void polyselect_thread_team_set_idx(polyselect_thread_team_ptr team, unsigned int idx)
{
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main_data = league->main;
    mpz_set_ui(team->ad, main_data->incr);
    mpz_mul_si(team->ad, team->ad, idx);
    mpz_add(team->ad, team->ad, main_data->admin);
    polyselect_poly_header_set_Nd(team->header, main_data->N, main_data->d);
    polyselect_poly_header_set_ad(team->header, team->ad);
}

void polyselect_thread_team_clear(polyselect_thread_team_ptr team)
{
    //// barrier_destroy(&team->sync_task->barrier, &team->lock);
    //// pthread_cond_destroy(&team->sync_task->wait_begintask);
    //// pthread_cond_destroy(&team->sync_task->wait_endtask);
    polyselect_shash_clear_multi(team->SH, team->league->main->finer_grain_threads);
    free(team->SH);
   
    polyselect_poly_header_clear(team->header);
    polyselect_proots_clear(team->R);
    polyselect_stats_clear(team->stats);
    mpz_clear(team->ad);

    pthread_cond_destroy(&team->count->w_ready_full);
    pthread_cond_destroy(&team->count->w_ready_empty);
    pthread_cond_destroy(&team->count->w_sync_empty);
    pthread_cond_destroy(&team->count->w_job);
    pthread_cond_destroy(&team->count->w_async_empty);
    pthread_cond_destroy(&team->count->w_sync2);
    pthread_cond_destroy(&team->count->w_roaming);
    barrier_destroy(&team->count->roaming_barrier, &team->lock);
    pthread_mutex_destroy(&team->lock);

}

#if 0
void accept_wait_causes(polyselect_thread_ptr thread, int ignore_async, ...)
{
    va_list ap;
    va_start(ap, ignore_async);
    polyselect_thread_team_ptr team = thread->team;
    /* make sure that all are waiting as we expect they are.  */
    unsigned int ti = thread->thread_index;
    unsigned int team_base = ti - (ti % team->size);
    for(unsigned int i = 0 ; i < team->size ; i++) {
        if (thread->thread_index == (team_base + i))
            continue;
        polyselect_thread_ptr T = &team->main_nonconst->threads[team_base + i];
        if (ignore_async && !T->is_in_sync_group)
            continue;
        enum wait_cause w;
        for( ; (w = va_arg(ap, enum wait_cause)) != W_NONE ; ) {
            if (w == W_NONE)
                break;
            if (T->why_wait == w)
                break;
        }
        if (w == W_NONE)
            fprintf(stderr, "ERROR, thread %d sees that thread %d is waiting for something unexpected\n", thread->thread_index, team_base + i);
    }
    va_end(ap);
}

void accept_signal_causes(polyselect_thread_ptr thread, ...)
{
    va_list ap;
    va_start(ap, thread);
    enum signal_cause s;
    for( ; (s = va_arg(ap, enum signal_cause)) != S_NONE ; ) {
        /* Check that we're awaken by a reason we expect */
        if (thread->team->why_signal == va_arg(ap, enum signal_cause))
            return;
    }
    fprintf(stderr, "ERROR, thread %d waits for roaming or barrier to end, but gets awaken by something unexpected\n", thread->thread_index);
}
#endif

const char * signal_string[] = {
    [S_NEW_JOB] = "new job",
    [S_NO_ROAMING] = "no roaming",
    [S_LEFT_BARRIER] = "left barrier",
    [S_REACHED_BARRIER] = "reached barrier",
    [S_PLEASE_ENTER_BARRIER] = "please enter barrier",
    [S_NONE] = "**NO SIGNAL, BUG**",
};

const char * wait_string[] = {
    [W_NONE] = "**NO WAIT CAUSE, BUG**",
    [W_NEW_JOB] = "new job",
    [W_LEAVING_BARRIER] = "leaving barrier",
    [W_REACHING_BARRIER] = "reaching barrier",
    [W_NO_ROAMING_AND_NO_BARRIER_TAIL] = "end of roaming and barrier tail",
};

#if 0
void cond_helper_broadcast(polyselect_thread_team_ptr team, polyselect_thread_ptr thread, enum signal_cause s, int ignore_async, ...)
{
    fprintf(stderr, "thread %d signals %s\n", thread->thread_index,
            signal_string[s]);
    team->why_signal = s;
    // accept_wait_causes(thread, 1, W_NEW_JOB, W_NONE);
    va_list ap, aq;
    va_start(ap, ignore_async);
    va_copy(aq, ap);
    // polyselect_thread_team_ptr team = thread->team;
    /* make sure that all are waiting as we expect they are.  */
    unsigned int ti = thread->thread_index;
    unsigned int team_base = ti - (ti % team->size);
    for(unsigned int i = 0 ; i < team->size ; i++) {
        if (thread->thread_index == (team_base + i))
            continue;
        polyselect_thread_ptr T = &team->main_nonconst->threads[team_base + i];
        if (ignore_async && !T->is_in_sync_group)
            continue;
        enum wait_cause w;
        for( ; (w = va_arg(ap, enum wait_cause)) != W_NONE ; ) {
            if (T->why_wait == w)
                break;
        }
        if (w == W_NONE)
            fprintf(stderr, "ERROR, thread %d signals %s and sees that thread %d is waiting for unexpected: %s\n", thread->thread_index, signal_string[s], team_base + i, wait_string[T->why_wait]);
        for( ; (w = va_arg(aq, enum wait_cause)) != W_NONE ; ) {
            fprintf(stderr, "  (expected: %s)\n", wait_string[w]);
        }
    }
    va_end(ap);
    va_end(aq);
    pthread_cond_broadcast(&team->wait);
}

void cond_helper_wait(polyselect_thread_team_ptr team, polyselect_thread_ptr thread, enum wait_cause w, ...)
{
    fprintf(stderr, "thread %d waits for %s\n", thread->thread_index,
            wait_string[w]);
    thread->why_wait = w;
    pthread_cond_wait(&team->wait, &team->lock);
    thread->why_wait = W_NONE;
    va_list ap;
    va_list aq;
    va_start(ap, w);
    va_copy(aq, ap);
    enum signal_cause s;
    for( ; (s = va_arg(ap, enum signal_cause)) != S_NONE ; ) {
        /* Check that we're awaken by a reason we expect */
        if (thread->team->why_signal == s)
            break;
    }
    if (s == S_NONE) {
        fprintf(stderr, "ERROR, thread %d waits for %s and is awaken by unexpected %s\n", thread->thread_index, wait_string[w], signal_string[thread->team->why_signal]);
        for( ; (s = va_arg(aq, enum signal_cause)) != S_NONE ; ) {
            fprintf(stderr, "  (expected: %s)\n", signal_string[s]);
        }
    }
    va_end(ap);
    va_end(aq);
}
#endif

void polyselect_thread_team_post_work(polyselect_thread_team_ptr team, polyselect_thread_ptr thread, void (*f)(polyselect_thread_ptr), void * arg)
{
    /* This is called with the team lock held ! */
    ASSERT_ALWAYS(team == thread->team);

    /*
    fprintf(stderr, "thread %d posts work for %d sync thread in team %d\n",
            thread->thread_index, team->count->sync, team->team_index);
            */

    team->task->expected = team->count->sync;
    team->task->f = f;
    team->task->arg = arg;
    // team->task->done = 0;
    
    barrier_resize_unlocked(&team->count->roaming_barrier, team->task->expected);
 
    // for performance reasons, we prefer if the caller takes care of
    // making this call _when it is needed_, that is when the SH
    // structure is used underneath.
    //
    // polyselect_shash_reset_multi(team->SH, tk->expected);

    pthread_cond_broadcast(&team->count->w_job);

    /* do my share ! */
    polyselect_thread_team_enter_sync_task(team, thread);
    (*f)(thread);
    /* same as polyselect_thread_team_leave_sync_task, but we're the one
     * that is interested in the result, so we want to make sure that
     * everyone's done.
     */
    for(team->count->sync2-- ; team->count->sync2 ; )
        pthread_cond_wait(&team->count->w_sync2, &team->lock);
}

void polyselect_thread_team_post_work_stop(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    ASSERT_ALWAYS(team == thread->team);

    /* We want to make sure that there's no thread in sync2 state at
     * the moment.
     */
    for( ; team->count->sync2 ; ) {
        pthread_cond_wait(&team->count->w_sync2, &team->lock);
    }

    /*
    fprintf(stderr, "thread %d posts STOP for %d sync thread in team %d\n",
            thread->thread_index, team->count->sync, team->team_index);
            */

    /* At this point no thread is reading the barrier state, we may
     * modify its characteristics.
     *
     * There's no roaming thread either, so that we're sure that all are
     * waiting on the condition wait. (note that the counters are checked
     * after the cond_wait).
     */

    team->task->expected = 0;
    // team->task->f = NULL;
    // team->task->arg = NULL;
    // team->task->done = 0;

    pthread_cond_broadcast(&team->count->w_job);
}

void polyselect_thread_team_i_am_ready(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED) {
    if (++team->count->ready == team->size) {
        pthread_cond_broadcast(&team->count->w_ready_full);
    } else {
        pthread_cond_wait(&team->count->w_ready_full, &team->lock);
    }
}

/* called with team lock held.
 *
 * as all transition functions, this function must be called with the
 * team lock held
 */
void polyselect_thread_team_enter_async(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    team->count->async++;
    --team->count->ready;
    if (team->count->ready == 0)
        pthread_cond_broadcast(&team->count->w_ready_empty);
#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d has entered async group (current: sync=%d async=%d ready=%d\n",
            thread->thread_index, team->count->sync, team->count->async, team->count->ready);
#endif
}

/* called with team lock held.
 *
 * as all transition functions, this function must be called with the
 * team lock held
 */
void polyselect_thread_team_leave_async(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    team->count->ready++;
    if (--team->count->async == 0)
        pthread_cond_broadcast(&team->count->w_async_empty);
}

/* called with team lock held.
 */
void polyselect_thread_team_enter_sync_zone(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d wants to enter sync group (current: sync=%d ready=%d\n",
            thread->thread_index, team->count->sync, team->count->ready);
#endif
    --team->count->ready;
    /*
            */
    /* If all threads have taken a decision as to what they're going to
     * do, allow the first thread in the sync zone to orchestrate the
     * work.
     */

    /* If we don't wait, this can happen:
     *
# thread 5 completed ad=9005760 at time=53.70s ; ad: 21.28s
thread 5 leaves sync group (current: sync=4 ready=0
thread 6 leaves sync group (current: sync=3 ready=1
thread 7 leaves sync group (current: sync=2 ready=2
thread 4 leaves sync group (current: sync=1 ready=3
thread 4 wants to enter sync group (current: sync=0 ready=4
thread 4 has entered sync group (current: sync=1 ready=3
thread 4 is 0-th sync thread in team 1
thread 5 wants to enter sync group (current: sync=1 ready=3
thread 5 has entered sync group (current: sync=2 ready=2
thread 5 is 1-th sync thread in team 1
thread 6 has entered async group (current: sync=2 async=1 ready=1
1636062335.338436 6 enter match
thread 7 has entered async group (current: sync=2 async=2 ready=0
1636062335.338457 7 enter match

    if (team->count->ready == 0)
        pthread_cond_broadcast(&team->count->w_ready_empty);
    else
        pthread_cond_wait(&team->count->w_ready_empty, &team->lock);
     */


    /* If we do sync++ before the wait, we'll violate our promise that a
     * thread in sync is necessarily waiting on w_job */
    thread->index_in_sync_zone = team->count->sync++;
#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d has entered sync group (current: sync=%d ready=%d\n",
            thread->thread_index, team->count->sync, team->count->ready);
#endif

}

/* called with team lock held.
 *
 * This is the terminal call for the sync workers. The workers goes from
 * state "sync" to state "ready".
 */
void polyselect_thread_team_leave_sync_zone(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d leaves sync group (current: sync=%d ready=%d\n",
            thread->thread_index, team->count->sync, team->count->ready);
#endif
    team->count->ready++;
    /*
     * There's a subtle catch here, since we want the full
     * wind-down sequence to terminate before we can start entering
     * threads in a sync zone again for a new ad.
     */
    if (--team->count->sync == 0)
        pthread_cond_broadcast(&team->count->w_sync_empty);
    else
        pthread_cond_wait(&team->count->w_sync_empty, &team->lock);

}

/* called with team lock held.
 *
 * This is called after each sync task, and in particular before the
 * leave_sync_group call.
 */
void polyselect_thread_team_leave_sync_task(polyselect_thread_team_ptr team MAYBE_UNUSED, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    if (--team->count->sync2 == 0)
        pthread_cond_broadcast(&team->count->w_sync2);
}

void polyselect_thread_team_enter_sync_task(polyselect_thread_team_ptr team MAYBE_UNUSED, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    ++team->count->sync2;
}

/* This must be called with the lock acquired ! */
void polyselect_thread_team_enter_roaming(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    team->count->roaming++;
}

/* This must be called with the lock released ! */
void polyselect_thread_team_leave_roaming(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    if (--team->count->roaming == 0)
        pthread_cond_broadcast(&team->count->w_roaming);
}

void polyselect_thread_team_roaming_barrier(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    /* This acquires team->lock */
    barrier_wait(&team->count->roaming_barrier, NULL, NULL, NULL);
}

