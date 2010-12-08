#include <stdlib.h>
#include "worker-threads.h"
#include "macros.h"

static void * thread_job(void * ti)
{
    struct worker_threads_group * tg = ((struct worker_thread_info *) ti) -> tg;
    int i = ((struct worker_thread_info *) ti) -> i;

    pthread_mutex_lock(&tg->mu);
    tg->working++;
    pthread_cond_signal(&tg->cond_out);
    /* The mutex is released later by pthread_cond_wait */

    unsigned int res = 0;

    for( ; ; ) {
        int rc = pthread_cond_wait(&tg->cond_in,&tg->mu);
	ASSERT_ALWAYS(rc == 0);
        tg->working++;
        pthread_mutex_unlock(&tg->mu);

        /* This is called with the mutex unlocked. */
        if (tg->f == NULL)
            break;

        (*tg->f)(tg, i, tg->arg);

        /* resume serialized behaviour */
        pthread_mutex_lock(&tg->mu);
        res = ++(tg->done);

        if (res == tg->n) {
            tg->done = 0;
            tg->working = 0;
            rc = pthread_cond_signal(&tg->cond_out);
	    ASSERT_ALWAYS(rc == 0);
        }
    }

    pthread_mutex_unlock(&tg->mu);

    return NULL;
}

void worker_threads_do(struct worker_threads_group * tg, worker_func_t f, void *arg)
{
    tg->f = f;
    tg->arg = arg;
    /* This mutex is associated to the cond_in and cond_out conditions.  */
    pthread_mutex_lock(&tg->mu);
    int rc = pthread_cond_broadcast(&tg->cond_in);
    ASSERT_ALWAYS(rc == 0);
    /* All threads are now working */
    /* Now switch to the outgoing condition, wait for the last thread to
     * release us */
    rc = pthread_cond_wait(&tg->cond_out,&tg->mu);
    ASSERT_ALWAYS(rc == 0);
    pthread_mutex_unlock(&tg->mu);
}


struct worker_threads_group * worker_threads_init(unsigned int n)
{
    struct worker_threads_group * tg;
    tg = malloc(sizeof(struct worker_threads_group));

    tg->n = n;
    tg->threads = malloc(n * sizeof(pthread_t));
    tg->thread_table = malloc(n * sizeof(struct worker_thread_info));
    pthread_cond_init(&tg->cond_in, NULL);
    pthread_cond_init(&tg->cond_out, NULL);
    pthread_mutex_init(&tg->mu, NULL);
    tg->working = 0;
    tg->done = 0;
    for(unsigned int i = 0 ; i < n ; i++) {
        int rc;
        tg->thread_table[i].i = i;
        tg->thread_table[i].tg = tg;
        void * ti = &(tg->thread_table[i]);
        rc = pthread_create(&(tg->threads[i]), NULL, &thread_job, ti);
        ASSERT_ALWAYS(rc == 0);
    }

    /* We must make sure that all threads have started properly. Since
     * we've got everything available for writing a barrier, we spare the
     * hassle of defining an extra data member for it. */
    pthread_mutex_lock(&tg->mu);
    for( ; tg->working != n ; ) {
        int rc = pthread_cond_wait(&tg->cond_out, &tg->mu);
	ASSERT_ALWAYS(rc == 0);
    }
    pthread_mutex_unlock(&tg->mu);

    return tg;
}

void worker_threads_clear(struct worker_threads_group * tg)
{
    /* tg->f == NULL is used as a termination condition */
    tg->f = NULL;
    pthread_mutex_lock(&tg->mu);
    pthread_cond_broadcast(&tg->cond_in);
    pthread_mutex_unlock(&tg->mu);
    for(unsigned int i = 0 ; i < tg->n ; i++) {
        void * res;
        int rc = pthread_join(tg->threads[i], &res);
        ASSERT_ALWAYS(rc == 0);
    }
    free(tg->threads);
    free(tg->thread_table);
    pthread_mutex_destroy(&tg->mu);
    pthread_cond_destroy(&tg->cond_in);
    pthread_cond_destroy(&tg->cond_out);
    free(tg);
}

