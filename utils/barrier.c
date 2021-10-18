#include "cado.h" // IWYU pragma: keep
#include <stddef.h> // NULL
#include <errno.h>
#include <pthread.h>
#include "barrier.h"

int barrier_init (barrier_t *barrier, pthread_mutex_t * lock_reuse, int count)
{
    int rc;

    /* count == 0 might be used to pre-initialize a barrier, meant to be
     * re-tweaked at little cost later on with barrier_resize */

    if (lock_reuse) {
        barrier->lock = lock_reuse;
    } else {
        rc = pthread_mutex_init (&barrier->lock_private, NULL);
        if (rc != 0) return rc;
        barrier->lock = &barrier->lock_private;
    }

    barrier->left = barrier->count = count;
    barrier->event = 0;

    rc = pthread_cond_init (&barrier->cv, NULL);
    if (rc != 0 && !lock_reuse) {
        pthread_mutex_destroy (&barrier->lock_private);
        return rc;
    }

    return 0;
}

int barrier_destroy (barrier_t *barrier, pthread_mutex_t * lock_reuse)
{
    int rc = EBUSY;

    rc = pthread_mutex_lock (barrier->lock);
    if (rc != 0) return rc;

    int ok = barrier->left == barrier->count;
    rc = pthread_mutex_unlock (barrier->lock);

    if (!ok) return -EBUSY;
    if (rc != 0) return rc;

    int r = 0;

    if (!lock_reuse)
        r = pthread_mutex_destroy (&barrier->lock_private);
    rc = pthread_cond_destroy (&barrier->cv);
    if (rc && !r) r = rc;

    return r;
}

int barrier_resize_unlocked(barrier_t * barrier, int count)
{
    int rc = 0;

    for ( ; rc == 0 && (barrier->event & 1) ; ) {
        rc = pthread_cond_wait (&barrier->cv, barrier->lock);
    }
    if (rc != 0) return rc;

    barrier->left = barrier->count = count;

    return 0;
}

int barrier_resize(barrier_t * barrier, int count)
{
    int rc;

    rc = pthread_mutex_lock (barrier->lock);
    if (rc != 0) return rc;

    barrier_resize_unlocked(barrier, count);

    rc = pthread_mutex_unlock (barrier->lock);
    return rc;
}

int barrier_finish_unlocked(barrier_t * barrier)
{
    int rc = 0;

    /* It could be that not all threads have exited the previous barrier. As
     * usual, only the contended case matters here. The uncontended case won't
     * even see this loop.
     */
    for ( ; rc == 0 && (barrier->event & 1) ; ) {
        rc = pthread_cond_wait (&barrier->cv, barrier->lock);
    }

    return rc;
}

int barrier_wait_unlocked(barrier_t * barrier, 
        void (*in)(int, void *),
        void (*out)(int, void *), void * arg)
{
    int rc = 0;

    rc = barrier_finish_unlocked(barrier);

    --barrier->left;

    /* Call the (*in) function with mutex locked, and sequential value,
     * in order. */
    if (in) (*in)(barrier->left, arg);

    if (barrier->left) {
        int event = barrier->event;

        /* protect against possible spurious wakeups */
        do {
            rc = pthread_cond_wait (&barrier->cv, barrier->lock);
            /* Error codes are returned as negative numbers */
            if (rc != 0) break;
        } while (event == barrier->event);
    } else {
        ++barrier->event;

        /* Wake up everybody. */
        rc = pthread_cond_broadcast (&barrier->cv);

        if (rc == 0)
            rc = BARRIER_SERIAL_THREAD;
    }

    /* This has the mutex locked */
    if (out) (*out)(barrier->left, arg);
    ++barrier->left;
    
    if (barrier->left == barrier->count) {
        /* We're leaving last. Increase barrier->event, so that its low bit
         * can be used as an indicator of previous barrier completion. */
        barrier->event++;
        pthread_cond_broadcast (&barrier->cv);
    }

    /* error: negative number, -(error code)
     * waker: return BARRIER_SERIAL_THREAD
     * other: return 0 */
    return rc;
}

int barrier_wait(barrier_t * barrier, 
        void (*in)(int, void *),
        void (*out)(int, void *), void * arg)
{
    int rc = 0;

    rc = pthread_mutex_lock (barrier->lock);
    if (rc != 0) return rc;

    rc = barrier_wait_unlocked(barrier, in, out, arg);

    pthread_mutex_unlock (barrier->lock);

    return rc;
}



