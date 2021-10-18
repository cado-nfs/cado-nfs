#include "cado.h"
#include "polyselect_thread_team.h"
#include "polyselect_thread.h"
#include "barrier.h"

void polyselect_thread_team_init(polyselect_thread_team_ptr team, polyselect_thread_league_ptr league, polyselect_main_data_ptr main, unsigned int team_index)
{
    memset(team, 0, sizeof(polyselect_thread_team));

    team->team_index = team_index;
    team->size = league->main->finer_grain_threads;
    team->sync_busy = 0;
    team->league = league;
    team->main_nonconst = main;
    team->done = 0;
    pthread_mutex_init(&team->lock, NULL);
    team->sync_task->expected_participants = 0;
    team->sync_task->f = NULL;

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
    polyselect_main_data_srcptr main = league->main;
    mpz_init(team->ad);
    polyselect_poly_header_init(team->header);
    polyselect_proots_init(team->R, main->d, league->pt->lenPrimes);
    polyselect_stats_init(team->stats, main->keep);
    team->rstate = team->stats->rstate;

    team->SH = malloc(main->finer_grain_threads * sizeof(polyselect_shash_t));
    size_t size_hint = POLYSELECT_SHASH_ALLOC_RATIO * league->pt->lenPrimes;
    polyselect_shash_init_multi(team->SH, size_hint, main->finer_grain_threads);

    pthread_cond_init(&team->sync_task->wait_begintask, NULL);
    pthread_cond_init(&team->sync_task->wait_endtask, NULL);
    barrier_init(&team->sync_task->barrier, &team->lock, 0);
}

void polyselect_thread_team_set_idx(polyselect_thread_team_ptr team, unsigned int idx)
{
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main = league->main;
    mpz_set_ui(team->ad, main->incr);
    mpz_mul_si(team->ad, team->ad, idx);
    mpz_add(team->ad, team->ad, main->admin);
    polyselect_poly_header_set_Nd(team->header, main->N, main->d);
    polyselect_poly_header_set_ad(team->header, team->ad);
}

void polyselect_thread_team_clear(polyselect_thread_team_ptr team)
{
    barrier_destroy(&team->sync_task->barrier, &team->lock);
    pthread_cond_destroy(&team->sync_task->wait_begintask);
    pthread_cond_destroy(&team->sync_task->wait_endtask);
    polyselect_shash_clear_multi(team->SH, team->league->main->finer_grain_threads);
    free(team->SH);
   
    polyselect_poly_header_clear(team->header);
    polyselect_proots_clear(team->R);
    polyselect_stats_clear(team->stats);
    mpz_clear(team->ad);

    pthread_barrier_destroy(&team->barrier);

    pthread_mutex_destroy(&team->lock);

}

void polyselect_thread_team_post_work(polyselect_thread_team_ptr team, polyselect_thread_ptr thread, void (*f)(polyselect_thread_ptr), void * arg)
{
    /* This is called with the team lock held ! */
    ASSERT_ALWAYS(thread->index_in_sync_team == 0);
    ASSERT_ALWAYS(team == thread->team);

    struct polyselect_thread_team_sync_task * tk = team->sync_task;

    tk->expected_participants = team->sync_busy;
    tk->f = f;
    tk->arg = arg;
    tk->done = 0;
    barrier_resize_unlocked(&tk->barrier, tk->expected_participants);
 
    // for performance reasons, we prefer if the caller takes care of
    // making this call _when it is needed_, that is when the SH
    // structure is used underneath.
    //
    // polyselect_shash_reset_multi(team->SH, tk->expected_participants);


    /* Important: the only moment when a thread of this team releases
     * the team lock after having increased sync_busy is when waiting on
     * the sync_wait condition! */

    pthread_cond_broadcast(&tk->wait_begintask);

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

    barrier_finish_unlocked(&tk->barrier);

    tk->expected_participants = 0;
    tk->f = NULL;
    tk->arg = NULL;
    tk->done = 0;
    team->sync_busy = 0;
    
    /* I think it's enough. */
    pthread_cond_broadcast(&tk->wait_begintask);
}


void polyselect_thread_team_end_subtask(polyselect_thread_team_ptr team)
{
    struct polyselect_thread_team_sync_task * tk = team->sync_task;

    /* called with team_lock held */
    barrier_wait_unlocked(&tk->barrier, NULL, NULL, NULL);
}
