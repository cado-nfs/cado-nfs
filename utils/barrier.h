#ifndef BARRIER_H_
#define BARRIER_H_

/* This is the interface to barrier synchronization waits.
 *
 * XXX Note that pthread_barrier_wait does only part of the stuff we
 * allow here.
 *
 * A barrier object must first be initialized with barrier_init. The
 * second argument specifies the number of threads expected to join.
 *
 * barrier_wait implements the wait. The function only returns when all
 * the threads have reached it. The second and third arguments provide
 * additional capabilities:
 *
 * If the second argument is non null, the function specified is called
 * for each thread with a mutex held, prior to reaching the barrier
 * point. The int argument passed to this function in an ordinal
 * representing the order in which threads reach the barrier. First
 * thread to reach calls the function with this argument set to 0, and
 * last to reach with the argument set to (n-1).
 *
 * If the third argument is not null, the function specified is called
 * with a mutex held _after_ reaching the barrier point.
 *
 * notes:
 *
 * It is safe to call barrier_wait several times in a row ; that is, in a
 * loop, you don't have to intermix barrier waits on two different
 * barriers. Even if a fast thread reaches the next barrier while the
 * slow one hardly escaped the previous one, the barrier data will not be
 * messed up.
 *
 * Positive return values indicate errors.
 *
 * BARRIER_SERIAL_THREAD (which is a negative number) is returned for
 * exactly one thread. 0 is returned for others.
 */

#include <pthread.h>

typedef struct barrier_tag {
    pthread_mutex_t     * lock;
    pthread_mutex_t     lock_private;
    pthread_cond_t      cv;
    int                 left;
    int                 count;
    int                 event;
} barrier_t;

#ifdef	__cplusplus
extern "C" {
#endif
#define BARRIER_SERIAL_THREAD   -1
extern int barrier_init (barrier_t *, pthread_mutex_t * lock, int);
extern int barrier_destroy (barrier_t *, pthread_mutex_t * lock);
extern int barrier_resize(barrier_t * barrier, int count);
extern int barrier_wait (barrier_t *,
        void (*)(int, void*), void (*)(int, void*), void *);

/* These calls are meant to allow more consistent behaviour when the
 * caller has already acquired the barrier's defined lock. It's only
 * useful when barrier_init was called with a custom lock, of course */
extern int barrier_resize_unlocked(barrier_t * barrier, int count);
extern int barrier_wait_unlocked (barrier_t *,
        void (*)(int, void*), void (*)(int, void*), void *);
extern int barrier_finish_unlocked(barrier_t * barrier);

#ifdef	__cplusplus
}
#endif

#endif	/* BARRIER_H_ */
