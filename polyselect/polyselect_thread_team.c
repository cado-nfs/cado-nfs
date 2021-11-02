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
    pthread_cond_init(&team->wait, NULL);
    team->sync_ready = 0;
    team->sync_roaming = 0;
    team->leaving_barrier = 0;
    team->sync_task->expected = 0;
    team->sync_task->in_barrier = 0;
    team->sync_task->f = NULL;
    team->why_signal = S_NONE;

    /* This is not the whole story. We defer some of the initialization
     * to polyselect_thread_team_late_init
     */
    pthread_barrier_init(&team->barrier, NULL, team->size);
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

    pthread_barrier_destroy(&team->barrier);

    pthread_cond_destroy(&team->wait);
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

void polyselect_thread_team_post_work(polyselect_thread_team_ptr team, polyselect_thread_ptr thread, void (*f)(polyselect_thread_ptr), void * arg)
{
    /* This is called with the team lock held ! */
    ASSERT_ALWAYS(thread->index_in_sync_team == 0);
    ASSERT_ALWAYS(team == thread->team);

    struct polyselect_thread_team_sync_task * tk = team->sync_task;

    fprintf(stderr, "thread %d posts work for %d sync thread in team %d\n",
            thread->thread_index, team->sync_ready, team->team_index);

    tk->expected = team->sync_ready;
    tk->f = f;
    tk->arg = arg;
    tk->done = 0;
    
    // barrier_resize_unlocked(&tk->barrier, tk->expected);
 
    // for performance reasons, we prefer if the caller takes care of
    // making this call _when it is needed_, that is when the SH
    // structure is used underneath.
    //
    // polyselect_shash_reset_multi(team->SH, tk->expected);


    cond_helper_broadcast(team, thread, S_NEW_JOB, 1, W_NEW_JOB, W_NONE);

    /* do my share ! */
    (*f)(thread);
}

void polyselect_thread_team_post_work_stop(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    /* This is called with the team lock held !
     *
     * We're **NOT** going to release the team lock, so that any team
     * member that enters it anew will have a fresh situation to decide
     * upon
     */
    ASSERT_ALWAYS(thread->index_in_sync_team == 0);
    ASSERT_ALWAYS(team == thread->team);

    struct polyselect_thread_team_sync_task * tk = team->sync_task;

    /* Wait until the previous barrier is cleared by everyone */
    /* XXX or could we / should we check on team->sync_task->in_barrier <
     * team->sync_task->expected  instead?
     */
    for( ; team->sync_roaming || team->leaving_barrier ; ) {
        cond_helper_wait(team, thread, W_NO_ROAMING_AND_NO_BARRIER_TAIL, S_LEFT_BARRIER,S_NO_ROAMING, S_NONE);
    }

    fprintf(stderr, "thread %d posts STOP for %d sync thread in team %d\n",
            thread->thread_index, team->sync_ready, team->team_index);

    /* At this point no thread is reading the barrier state, we may
     * modify its characteristics.
     *
     * There's no roaming thread either, so that we're sure that all are
     * waiting on the condition wait. (note that the counters are checked
     * after the cond_wait).
     */

    tk->expected = 0;
    // tk->f = NULL;
    // tk->arg = NULL;
    tk->done = 0;

    polyselect_thread_team_sync_group_leave(team, thread);
}

/* This must be called with the lock acquired ! */
void polyselect_thread_team_sync_group_leave(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    /* This is really just a barrier wait */
    struct polyselect_thread_team_sync_task * tk = team->sync_task;
    /* wait on the previous barrier */
    for( ; team->sync_roaming || team->leaving_barrier ; ) {
        cond_helper_wait(team, thread, W_NO_ROAMING_AND_NO_BARRIER_TAIL, S_LEFT_BARRIER,S_NO_ROAMING, S_NONE);
    }
    if (++tk->in_barrier == team->sync_ready) {
        team->leaving_barrier = 1;
        cond_helper_broadcast(team, thread, S_REACHED_BARRIER, 1, W_REACHING_BARRIER, W_NONE);
    } else {
        /* This call is slightly special. Not all threads are poised to
         * enter the barrier from this call, since they might be simply
         * waiting on a new task. It seems really clumsy, I'm pretty sure
         * it's possible to do better, but presently my impression is
         * that I have to try and signal them once, so that there's
         * necessarily going to enter the barrier.
         */
        cond_helper_broadcast(team, thread, S_PLEASE_ENTER_BARRIER, 1, W_NEW_JOB, W_NONE);
        do {
            cond_helper_wait(team, thread, W_REACHING_BARRIER, S_REACHED_BARRIER,S_NO_ROAMING, S_NONE);
        } while (!team->leaving_barrier);
    }
    team->sync_ready--;
    if (--tk->in_barrier == 0) {
        team->leaving_barrier = 0;
        cond_helper_broadcast(team, thread, S_LEFT_BARRIER, 1, W_LEAVING_BARRIER, W_NONE);
    }
}

static void common_barrier_wait_lockstep(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    /* This is really just a barrier wait */
    struct polyselect_thread_team_sync_task * tk = team->sync_task;
    /* wait on the previous barrier */
    for( ; team->sync_roaming || team->leaving_barrier ; ) {
        cond_helper_wait(team, thread, W_NO_ROAMING_AND_NO_BARRIER_TAIL, S_LEFT_BARRIER,S_NO_ROAMING, S_NONE);
    }
    if (++tk->in_barrier == tk->expected) {
        team->leaving_barrier = 1;
        cond_helper_broadcast(team, thread, S_REACHED_BARRIER, 1, W_REACHING_BARRIER, W_NONE);
    } else {
        do {
            cond_helper_wait(team, thread, W_REACHING_BARRIER, S_REACHED_BARRIER,S_NO_ROAMING, S_NONE);
        } while (!team->leaving_barrier);
    }
    if (--tk->in_barrier == 0) {
        team->leaving_barrier = 0;
        cond_helper_broadcast(team, thread, S_LEFT_BARRIER, 1, W_LEAVING_BARRIER, W_NONE);
    }
}

/* This must be called with the lock acquired ! */
void polyselect_thread_team_end_subtask(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    fprintf(stderr, "team %d reaches end_subtask with roaming=%d leaving_barrier=%d in_barrier=%d sync_ready=%d\n",
            team->team_index,
            team->sync_roaming,
            team->leaving_barrier,
            team->sync_task->in_barrier,
            team->sync_ready);
    common_barrier_wait_lockstep(team, thread);
    fprintf(stderr, "team %d completes end_subtask with roaming=%d leaving_barrier=%d in_barrier=%d sync_ready=%d\n",
            team->team_index,
            team->sync_roaming,
            team->leaving_barrier,
            team->sync_task->in_barrier,
            team->sync_ready);
}

/* This must be called with the lock acquired ! */
void polyselect_thread_team_sync_group_enter(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    team->sync_ready++;
}

/* This must be called with the lock acquired ! */
void polyselect_thread_team_sync_group_begin_roaming(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    team->sync_roaming++;
    pthread_mutex_unlock(&team->lock);
}

/* This must be called with the lock released ! */
void polyselect_thread_team_sync_group_end_roaming(polyselect_thread_team_ptr team, polyselect_thread_ptr thread MAYBE_UNUSED)
{
    pthread_mutex_lock(&team->lock);
    if (--team->sync_roaming == 0) {
        cond_helper_broadcast(team, thread, S_NO_ROAMING, 1, W_NO_ROAMING_AND_NO_BARRIER_TAIL, W_NONE);
    }
}

void polyselect_thread_team_sync_group_roaming_barrier(polyselect_thread_team_ptr team, polyselect_thread_ptr thread)
{
    pthread_mutex_lock(&team->lock);
    if (--team->sync_roaming == 0) {
        cond_helper_broadcast(team, thread, S_NO_ROAMING, 1, W_NO_ROAMING_AND_NO_BARRIER_TAIL, W_NONE);
    }
    common_barrier_wait_lockstep(team, thread);
    team->sync_roaming++;
    pthread_mutex_unlock(&team->lock);
}

