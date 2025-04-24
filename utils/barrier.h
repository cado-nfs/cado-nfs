#ifndef CADO_BARRIER_H
#define CADO_BARRIER_H

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

#ifdef	__cplusplus
#include <cstddef>
#endif

#include <pthread.h>

#include "macros.h"

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
#define BARRIER_SERIAL_THREAD   (-1)
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

#ifdef	__cplusplus
namespace cado_nfs {
/* interface is a subset of the c++20 std::barrier */
class barrier {
    barrier_t b[1];
    public:
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
    explicit inline barrier(std::ptrdiff_t expected) {
        barrier_init(b, nullptr, (int) expected);
    }
    inline ~barrier() { barrier_destroy(b, nullptr); }
    barrier(barrier const &) = delete;
    barrier(barrier &&) = delete;
    barrier& operator=(barrier const &) = delete;
    barrier& operator=(barrier &&) = delete;
    inline void arrive_and_wait() { barrier_wait(b, nullptr, nullptr, nullptr); }
};
}
#endif

#ifndef  HAVE_PTHREAD_BARRIER_WAIT
/* Also enable proxies for emulating the standard stuff with our more
 * sophisticated barriers */
#define my_pthread_barrier_t               barrier_t
#define my_pthread_barrierattr_t           int
#define MY_PTHREAD_BARRIER_SERIAL_THREAD   BARRIER_SERIAL_THREAD

#ifdef __cplusplus
extern "C" {
#endif

static inline int my_pthread_barrier_init(barrier_t * /* restrict */ barrier,
        const int * /* restrict */ attr MAYBE_UNUSED, unsigned count)
{
    return barrier_init(barrier, NULL, count);
}
static inline int my_pthread_barrier_wait(barrier_t * b)
{
    return barrier_wait(b, NULL, NULL, NULL);
}
static inline int my_pthread_barrier_destroy(barrier_t * b)
{
    return barrier_destroy(b, NULL);
}

#ifdef __cplusplus
}
#endif
#else
#define my_pthread_barrier_t               pthread_barrier_t
#define my_pthread_barrierattr_t           pthread_barrierattr_t
#define my_pthread_barrier_init            pthread_barrier_init
#define MY_PTHREAD_BARRIER_SERIAL_THREAD   PTHREAD_BARRIER_SERIAL_THREAD
#define my_pthread_barrier_wait            pthread_barrier_wait
#define my_pthread_barrier_destroy         pthread_barrier_destroy
#endif

#endif	/* CADO_BARRIER_H */
