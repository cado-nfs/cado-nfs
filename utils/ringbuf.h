#ifndef ROLLBUF_H_
#define ROLLBUF_H_

#include <sys/types.h>  // ssize_t
#include <stdio.h>      // for size_t, FILE

/* Example of a rotating, and reallocating buffer. A separate thread has
 * to fetch the data from the real source and fill the buffer. This way,
 * reads block as little as possible -- sort of a userland pipe, only we
 * get more control. Not clear we gain anything, but it's less fragile
 * than forking an external program.
 */
#include <pthread.h>

/* Getting data from a rolling buffer can be done via two interfaces. One
 * interface (ringbuf_get) provides a pointer, while the other
 * (ringbuf_get2) does not, and asks the ringbuf interface to provide a
 * pointer.
 *
 * The area pointed to by the pointer returned from ringbuf_get2 is free
 * to be read from by the calling thread until the next call to
 * ringbuf_get. This implies in particular that we assume that there is
 * exactly one threading getting data, no more.
 *
 * Note though that this is not implemented by means of a zero-copy
 * mechanism (doing so would make rollback_put break the active area of
 * data being read, in case a realloc() occurs). The data is copied to a
 * buffer exclusively dedicated to reading.
 *
 * (the reading buffer is allocated only if needed, and freed by
 * ringbuf_clear eventually)
 */
struct ringbuf_s {
    char * p;
    size_t alloc;
    size_t avail_to_read;
    size_t avail_to_write;
    const char * rhead;
    char * whead;
    char * rbuf;        /* Only for ringbuf_get2 */
    pthread_mutex_t mx[1];
    pthread_cond_t bored[1];
    int empty_count;
    int full_count;
    int done:1;
};

typedef struct ringbuf_s ringbuf[1];
typedef struct ringbuf_s * ringbuf_ptr;
typedef const struct ringbuf_s * ringbuf_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void ringbuf_init(ringbuf_ptr r, size_t initial_size);
extern void ringbuf_clear(ringbuf_ptr r);
extern int ringbuf_put(ringbuf_ptr r, char * p, size_t s);
extern void ringbuf_mark_done(ringbuf_ptr r);
extern int ringbuf_is_done(ringbuf_ptr r);

/* see above for the distinction between these two get() calls */
extern int ringbuf_get(ringbuf_ptr r, char * p, size_t s);
extern int ringbuf_get2(ringbuf_ptr r, void ** p, size_t s);

extern int ringbuf_strchr(ringbuf_ptr r, int c, size_t offset);

/* Equivalent of doing ringbuf_put for all bytes from the stdio stream f.
 *
 * The only difference is that this call does not automatically enlarge
 * the ring buffer, so an appropriate initial_size must have been
 * provided on initialization
 */
extern ssize_t ringbuf_feed_stream(ringbuf_ptr r, FILE * f);

extern int ringbuf_skip_get(ringbuf_ptr r, size_t s);

/* A quick accessor macro which does a 1-byte fetch from the ring buffer.
 * We must be sure that the get will succeed, and we must provide an
 * auxiliary pointer (here s_) which will be updated by the macro. s_ has
 * to be set to r_->rhead originally.
 *
 * n successive calls to RINGBUF_GET_ONE_BYTE must be followed by a call to
 * ringbug_skip_get(r_, n)
 */
#define RINGBUF_GET_ONE_BYTE(c_, r_, s_) do {				\
    c_ = *(s_)++;							\
    if (s_ >= r_->p + r_->alloc) {					\
        s_ = r_->p;							\
    }									\
} while (0)

extern int ringbug_skip_get(ringbuf_ptr r, size_t s);


#ifdef __cplusplus
}
#endif

#endif	/* ROLLBUF_H_ */
