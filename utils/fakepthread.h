#ifndef FAKEPTHREAD_H_
#define FAKEPTHREAD_H_

#error "Sorry, cado-nfs now _requires_ posix threads. some parts of the code cannot be made to easily cope with the situation where threads are not there"

/* conditions are the trouble maker here */

#include <signal.h>
#include "macros.h"

// See select_mpi.h for info on my_* stuff.
#ifdef __cplusplus
extern "C" {
#endif


/* {{{ create/join */
typedef void * my_pthread_t;
typedef int my_pthread_attr_t;

static inline int my_pthread_create(my_pthread_t * /* restrict */ t, const my_pthread_attr_t * /* restrict */ a MAYBE_UNUSED, void *(*s)(void*), void * arg)
{
    void * u=(*s)(arg);
    *t = u;
    return 0;
}

static inline int my_pthread_join(my_pthread_t t, void ** vp)
{
    if (vp)
        *vp = t;
    return 0;
}
/* }}} */

/* {{{ posix pthread barriers */
typedef int my_pthread_barrier_t;
typedef int my_pthread_barrierattr_t;

static inline int my_pthread_barrier_init(my_pthread_barrier_t * /* restrict */ b MAYBE_UNUSED, const my_pthread_barrierattr_t * /* restrict */ attr MAYBE_UNUSED, unsigned int c)
{
    ASSERT_ALWAYS(c == 1);
    return 0;
}

#define MY_PTHREAD_BARRIER_SERIAL_THREAD   -1
static inline int my_pthread_barrier_wait(my_pthread_barrier_t * b MAYBE_UNUSED)
{
    return MY_PTHREAD_BARRIER_SERIAL_THREAD;
}
static inline int my_pthread_barrier_destroy(my_pthread_barrier_t * b MAYBE_UNUSED)
{
    return 0;
}
/* }}} */

/* {{{ advanced barriers */
typedef int barrier_t;

static inline int barrier_wait(barrier_t * b MAYBE_UNUSED,
        void (*in)(int, void *),
        void (*out)(int, void *), void * arg)
{
    if (in) (*in)(0, arg);
    if (out) (*out)(0, arg);
    return 1;   /* waker returns 1 */
}

static inline int
barrier_init (barrier_t *b MAYBE_UNUSED, int c)
{ ASSERT_ALWAYS(c==1); return 0; }

static inline int
barrier_destroy (barrier_t *b MAYBE_UNUSED)
{ return 0; }
/* }}} */

/* {{{ mutexes */
typedef int my_pthread_mutex_t;
typedef int my_pthread_mutexattr_t;

static inline int
my_pthread_mutex_init(my_pthread_mutex_t *m MAYBE_UNUSED,
        const my_pthread_mutexattr_t * attr MAYBE_UNUSED)
{ return 0; }
static inline int
my_pthread_mutex_lock(my_pthread_mutex_t * m MAYBE_UNUSED)
{ return 0; }
static inline int
my_pthread_mutex_unlock(my_pthread_mutex_t * m MAYBE_UNUSED)
{ return 0; }
static inline int
my_pthread_mutex_destroy(my_pthread_mutex_t * m MAYBE_UNUSED)
{ return 0; }
/* }}} */



/* {{{ rwlocks */
typedef int my_pthread_rwlock_t;
typedef int my_pthread_rwlockattr_t;

static inline int
my_pthread_rwlock_init(my_pthread_rwlock_t * m MAYBE_UNUSED,
        const my_pthread_rwlockattr_t *attr MAYBE_UNUSED)
{ return 0; }

static inline int
my_pthread_rwlock_wrlock(my_pthread_rwlock_t * m MAYBE_UNUSED) { return 0; }
static inline int
my_pthread_rwlock_rdlock(my_pthread_rwlock_t * m MAYBE_UNUSED) { return 0; }
static inline int
my_pthread_rwlock_unlock(my_pthread_rwlock_t * m MAYBE_UNUSED) { return 0; }
static inline int
my_pthread_rwlock_destroy(my_pthread_rwlock_t * m MAYBE_UNUSED) { return 0; }
/* }}} */

static inline int
my_pthread_sigmask(int how, const sigset_t * set, sigset_t * oset)
{ return sigprocmask(how, set, oset); }


#ifdef __cplusplus
}
#endif

#endif	/* FAKEPTHREAD_H_ */
