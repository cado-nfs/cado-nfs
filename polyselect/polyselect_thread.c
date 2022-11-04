#include "cado.h"
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include "polyselect_thread.h"
#ifdef HAVE_HWLOC
#include "hwloc.h"
#endif
#include "polyselect_match.h"
#include "dllist.h"

void polyselect_thread_bind(polyselect_thread_ptr thread MAYBE_UNUSED)
{
#ifdef HAVE_HWLOC
    polyselect_thread_team_ptr team = thread->team;
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main_data = league->main;

    if (!main_data->bind)
        return;

    pthread_mutex_t * mlock = thread->main_lock;
    /* Do the binding. Temporarily acquire the main_data lock in order to do
     * so */
    char * cstr = NULL;
    char * mstr = NULL;
    if (thread->cpubind_set) {
        int rc = hwloc_bitmap_asprintf(&cstr, thread->cpubind_set);
        ASSERT_ALWAYS(rc >= 0);
    }
    if (league->membind_set) {
        int rc = hwloc_bitmap_asprintf(&mstr, league->membind_set);
        ASSERT_ALWAYS(rc >= 0);
    }

    pthread_mutex_lock(mlock);
    if (thread->cpubind_set) {
        hwloc_set_cpubind(main_data->topology, thread->cpubind_set, HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
        printf("# cpu-binding thread %d to [%s]\n", thread->thread_index, cstr);
    }
    if (league->membind_set) {
        hwloc_set_membind(main_data->topology, league->membind_set, HWLOC_MEMBIND_BIND,
                HWLOC_MEMBIND_THREAD |
                HWLOC_MEMBIND_STRICT |
                HWLOC_MEMBIND_BYNODESET);
        printf("# mem-binding thread %d to [%s]\n", thread->thread_index, mstr);
    }

    pthread_mutex_unlock(mlock);
    if (cstr)
        free(cstr);
    if (mstr)
        free(mstr);
#endif
}

void polyselect_thread_init(polyselect_thread_ptr thread, polyselect_thread_team_ptr team, polyselect_main_data_ptr main, unsigned int thread_index)
{
    thread->thread_index = thread_index;
    thread->team = team;
    thread->tid = pthread_self();
#ifdef HAVE_HWLOC
    thread->cpubind_set = NULL;
#endif
    thread->main_lock = &main->lock;
    thread->index_in_sync_zone = UINT_MAX;
    thread->why_wait = W_NONE;
    thread->is_in_sync_group = 0;
    thread->is_unlocked = 1;
}

FILE * chronogram = NULL;
pthread_mutex_t chronogram_lock = PTHREAD_MUTEX_INITIALIZER;

void polyselect_thread_chronogram_init(const char * filename)
{
    if (!filename) {
        chronogram = NULL;
        return;
    }
    chronogram = fopen(filename, "w");
    setbuf(chronogram, NULL);
}

void polyselect_thread_chronogram_clear()
{
    if (!chronogram)
        return;
    fclose(chronogram);
}

void polyselect_thread_chronogram_chat(polyselect_thread_srcptr thread, const char * fmt, ...)
{
    if (!chronogram) return;

    // TODO: do we want to restrict printing within groups? It's not
    // clear, in fact.
    // if (!polyselect_thread_is_leader(thread)) return;

    va_list ap;
    va_start(ap, fmt);
    char * msg;
    int rc = vasprintf(&msg, fmt, ap);
    ASSERT_ALWAYS(rc >= 0);
    struct timeval tv[1];
    gettimeofday(tv, NULL);
    pthread_mutex_lock(&chronogram_lock);
    // OSX has its special type for tv_usec...
    fprintf(chronogram, "%lu.%06lu %d %s\n", (unsigned long) tv->tv_sec, (unsigned long) tv->tv_usec, thread->thread_index, msg);
    pthread_mutex_unlock(&chronogram_lock);
    free(msg);
    va_end(ap);
}

static void polyselect_thread_provision_job_slots(polyselect_thread_ptr thread, size_t n)
{
    for( ; n-- ; ) {
        polyselect_match_info_ptr job = malloc(sizeof(polyselect_match_info_t));
        polyselect_match_info_init_trivial(job);
        dllist_push_back(&thread->empty_job_slots, &job->queue);
    }
}


static void polyselect_thread_unprovision_job_slots(polyselect_thread_ptr thread)
{
    dllist_for_each_safe(&thread->empty_job_slots, ptr) {
        dllist_pop(ptr);
        polyselect_match_info_ptr job = dllist_entry(ptr, struct polyselect_match_info_s, queue);
        polyselect_match_info_clear(job);
        free(job);
    }
}

void polyselect_thread_late_init(polyselect_thread_ptr thread)
{
    polyselect_thread_team_ptr team = thread->team;
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main_data = league->main;

    polyselect_stats_init(thread->stats, main_data->keep);
    dllist_init_head(&thread->async_jobs);
    dllist_init_head(&thread->empty_job_slots);

    polyselect_thread_provision_job_slots(thread, 32);

    thread->rstate = thread->stats->rstate;
}


void polyselect_thread_clear(polyselect_thread_ptr thread)
{
    /* thread->rstate is just a pointer */

    ASSERT_ALWAYS(dllist_is_empty(&thread->async_jobs));
    polyselect_stats_clear(thread->stats);
    polyselect_thread_unprovision_job_slots(thread);
}


