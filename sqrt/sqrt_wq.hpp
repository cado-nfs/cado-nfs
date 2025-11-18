#ifndef SQRT_SQRT_WQ_HPP_
#define SQRT_SQRT_WQ_HPP_

#include <pthread.h>

/* POSIX threads stuff : work queues */

typedef void * (*wq_func_t)(void *);

struct wq_task {
    wq_func_t f;
    void * arg;

    // the remaining fields are reserved. In particular, access to done
    // must be mutex protected.
    int done;

    pthread_mutex_t m_[1];
    pthread_cond_t c_[1];

    struct wq_task * next;
    // struct wq_task * prev;
};

struct work_queue {
    pthread_mutex_t m[1];
    pthread_cond_t c[1];   // used by waiters.
    pthread_t * clients;
    // TODO: Use a single tasklist item for a doubly linked list.
    struct wq_task * head;
};

extern struct wq_task * wq_push(struct work_queue * wq, wq_func_t f, void * arg);
extern struct wq_task * wq_pop_wait(struct work_queue * wq);
extern void * wq_waiter(void * wq);
extern void * wq_join(struct wq_task * t);
extern void wq_init(struct work_queue * wq, unsigned int n);
extern void wq_clear(struct work_queue * wq, unsigned int n);




#endif	/* SQRT_SQRT_WQ_HPP_ */
