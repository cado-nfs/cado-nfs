#include "cado.h" // IWYU pragma: keep
#include <stddef.h> // NULL
#include <errno.h>
#include <pthread.h>
#include "barrier.h"

int barrier_init (barrier_t *barrier, int count)
{
    int rc;

    if (count == 0)
        return -EINVAL;

    rc = pthread_mutex_init (&barrier->lock, NULL);
    if (rc != 0) return rc;

    barrier->left = barrier->count = count;
    barrier->event = 0;

    rc = pthread_cond_init (&barrier->cv, NULL);
    if (rc != 0) {
        pthread_mutex_destroy (&barrier->lock);
        return rc;
    }

    return 0;
}

int barrier_destroy (barrier_t *barrier)
{
    int rc = EBUSY;

    rc = pthread_mutex_lock (&barrier->lock);
    if (rc != 0) return rc;

    int ok = barrier->left == barrier->count;
    rc = pthread_mutex_unlock (&barrier->lock);

    if (!ok) return -EBUSY;
    if (rc != 0) return rc;

    int r = 0;

    r = pthread_mutex_destroy (&barrier->lock);
    rc = pthread_cond_destroy (&barrier->cv);
    if (rc && !r) r = rc;

    return r;
}

int barrier_wait(barrier_t * barrier, 
        void (*in)(int, void *),
        void (*out)(int, void *), void * arg)
{
    int rc;

    rc = pthread_mutex_lock (&barrier->lock);
    if (rc != 0) return rc;

    /* It could be that not all threads have exited the previous barrier. As
     * usual, only the contended case matters here. The uncontended case won't
     * even see this loop.
     */
    for ( ; rc == 0 && (barrier->event & 1) ; ) {
        rc = pthread_cond_wait (&barrier->cv, &barrier->lock);
    }

    --barrier->left;

    /* Call the (*in) function with mutex locked, and sequential value,
     * in order. */
    if (in) (*in)(barrier->left, arg);

    if (barrier->left) {
        int event = barrier->event;

        /* protect against possible spurious wakeups */
        do {
            rc = pthread_cond_wait (&barrier->cv, &barrier->lock);
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
    pthread_mutex_unlock (&barrier->lock);

    /* error: negative number, -(error code)
     * waker: return BARRIER_SERIAL_THREAD
     * other: return 0 */
    return rc;
}
