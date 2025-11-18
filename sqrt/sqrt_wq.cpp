#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <pthread.h>

#include "sqrt_wq.hpp"
#include "macros.h"

// puts a task on the pending list.
// returns a handle which might be waited for by a join.
struct wq_task * wq_push(struct work_queue * wq, wq_func_t f, void * arg)
{

    struct wq_task * t = (struct wq_task *) malloc(sizeof(struct wq_task));
    t->f = f;
    t->arg = arg;
    t->done = 0;
    if (t->f) {
        // t->f == NULL is used as a special marker. In this case the
        // mutexes are not initialized
        pthread_mutex_init(t->m_, NULL);
        pthread_cond_init(t->c_, NULL);
    }

    pthread_mutex_lock(wq->m);
    t->next = wq->head;
    wq->head = t;
    pthread_cond_signal(wq->c);
    pthread_mutex_unlock(wq->m);
    return t;
}

// grabs a task as soon as one is available.
struct wq_task * wq_pop_wait(struct work_queue * wq)
{
    struct wq_task * t;
    pthread_mutex_lock(wq->m);
    for( ; wq->head == NULL ; ) {
        pthread_cond_wait(wq->c, wq->m);
    }
    t = wq->head;
    wq->head = t->next;
    t->next = NULL;
    pthread_mutex_unlock(wq->m);
    return t;
}

void * wq_waiter(void * wq)
{
    for( ; ; ) {
        struct wq_task * t = wq_pop_wait((struct work_queue *) wq);
        if (t->f == NULL) {
            /* t == NULL is an indication that we're reaching the
             * end-of-work signal. For this special case, t->m_ and t->c_
             * have not been initialized.  */
            free(t);
            break;
        }
        void * res = (*t->f)(t->arg);
        pthread_mutex_lock(t->m_);
        t->done = 1;
        t->arg = res;   /* We reuse the arg field */
        pthread_cond_signal(t->c_);     // signal waiters.
        pthread_mutex_unlock(t->m_);
    }
    /* we have to obey the pthread_create proto. */
    return NULL;
}

void * wq_join(struct wq_task * t)
{
    pthread_mutex_lock(t->m_);
    for( ; t->done == 0 ; ) {
        pthread_cond_wait(t->c_, t->m_);
    }
    void * res = t->arg;
    pthread_mutex_unlock(t->m_);
    pthread_cond_destroy(t->c_);
    pthread_mutex_destroy(t->m_);
    free(t);
    return res;
}

void wq_init(struct work_queue * wq, unsigned int n)
{
    pthread_mutex_init(wq->m, NULL);
    pthread_cond_init(wq->c, NULL);
    wq->head = NULL;
    wq->clients = (pthread_t *) malloc(n * sizeof(pthread_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_create(wq->clients + i, NULL, &wq_waiter, wq);
    }
}

void wq_clear(struct work_queue * wq, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        // schedule end-of-work for everybody.
        wq_push(wq, NULL, NULL);
    }
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_join(wq->clients[i], NULL);
    }

    ASSERT_ALWAYS(wq->head == NULL);
    pthread_mutex_destroy(wq->m);
    pthread_cond_destroy(wq->c);
}

/*  */

