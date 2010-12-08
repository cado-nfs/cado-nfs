#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>
#include <unistd.h>
#include <stdint.h>
#include "balancing_rollbuf.h"

/* This is a hack. Define to 1 to disable */
#define ROLLBUF_ALIGNED_RETURNS sizeof(uint32_t)

/* This has been tested with small buffer sizes, so as to stress the
 * pause/resume mechanism (based on pthread conditions). It seems to
 * work. Note though that this code heavily relies on the fact that there
 * is only _one_ thread reading data.
 */

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

/* refuse to fill memory with incoming data beyond this size. Incoming
 * data which can only be stored by expanding the rollbuf capacity beyond
 * this size are paused */
#define ROLLBUF_MAX_SIZE        (1 << 26)

/* Getting data from a rolling buffer can be done via two interfaces. One
 * interface (rollbuf_get) provides a pointer, while the other
 * (rollbuf_get2) does not, and asks the rollbuf interface to provide a
 * pointer.
 *
 * The area pointed to by the pointer returned from rollbuf_get2 is free
 * to be read from by the calling thread until the next call to
 * rollbuf_get. This implies in particular that we assume that there is
 * exactly one threading getting data, no more.
 *
 * Note though that this is not implemented by means of a zero-copy
 * mechanism (doing so would make rollback_put break the active area of
 * data being read, in case a realloc() occurs). The data is copied to a
 * buffer exclusively dedicated to reading.
 *
 * (the reading buffer is allocated only if needed)
 */
#define ROLLBUF_READING_BUFFER_SIZE     (1 << 16)

void rollbuf_init(rollbuf_ptr r)
{
    memset(r, 0, sizeof(rollbuf));
    pthread_mutex_init(r->mx, NULL);
    pthread_cond_init(r->bored, NULL);
}

void rollbuf_clear(rollbuf_ptr r)
{
    /*
    // fprintf(stderr, "rollbuf: %d times full, %d times empty\n",
            r->full_count, r->empty_count);
            */
    pthread_mutex_destroy(r->mx);
    pthread_cond_destroy(r->bored);
    free(r->rbuf);
}


/* must be called with mutex locked !!! */
static void rollbuf_grow__(rollbuf_ptr r)
{
    size_t newalloc = r->alloc;
    newalloc += newalloc / 2;
    newalloc |= sysconf(_SC_PAGESIZE) - 1;
    newalloc++;
    if (r->avail_to_read) {
        char * newp = malloc(newalloc);
        size_t tail = r->alloc - (r->rhead - r->p);
        ptrdiff_t head = r->avail_to_read - tail;
        // tail + head == r->avail_to_read
        if (head > 0) {
            memcpy(newp, r->rhead, tail);
            memcpy(newp + tail, r->p, head);
        } else {
            memcpy(newp, r->rhead, r->avail_to_read);
        }
        free(r->p);
        r->p = newp;
        r->rhead = r->p;
        r->whead = r->p + r->avail_to_read;
    } else {
        r->p = realloc(r->p, newalloc);
        r->rhead = r->p;
        r->whead = r->p;
    }
    r->avail_to_write += newalloc - r->alloc;
    r->alloc = newalloc;
}

int rollbuf_put(rollbuf_ptr r, char * p, size_t s)
{
    pthread_mutex_lock(r->mx);
    // fprintf(stderr, "put(%zu): (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    for( ; s > r->avail_to_write ; ) {
        if (r->alloc >= ROLLBUF_MAX_SIZE) {
            if (s < r->alloc) {
                /* Then we want to drain our pipe first. */
                // fprintf(stderr, "put(%zu): on hold (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
                r->full_count++;
                pthread_cond_wait(r->bored, r->mx);
                // fprintf(stderr, "put(%zu): resuming (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
                continue;
            } else {
                /* Here there is no hope to drain our pipe, so we must
                 * exceed our desired window size. This is not much of a
                 * problem though, since the curl reading buffer is not
                 * expected to grow too large (or so we hope...)
                 */
                // fprintf(stderr, "Warning: buffer growing beyond max size ! (previous=%zu)\n", r->alloc);
            }
        }

        rollbuf_grow__(r);
    }

    size_t tail = r->alloc - (r->whead - r->p);
    if (r->avail_to_write <= tail || s <= tail) {
        assert(s <= r->avail_to_write);
        // s = MIN(r->avail_to_write, s);
        assert(s <= tail);
        memcpy(r->whead, p, s);
        r->whead += s;
    } else {
        assert(tail > 0);
        assert(s > tail);
        memcpy(r->whead, p, tail);
#ifndef NDEBUG
        ptrdiff_t head = r->avail_to_write - tail;
        assert(head > 0);
        assert(s-tail <= (size_t) head);
#endif
        // s = tail + MIN((size_t) head, s - tail);
        memcpy(r->p, p + tail, s - tail);
        r->whead = r->p + (s-tail);
    }
    if ((size_t) (r->whead - r->p) == r->alloc) {
        r->whead = r->p;
    }
    r->avail_to_write -= s;
    r->avail_to_read += s;
    /* Could be that someone is waiting for data to be read */
    pthread_cond_signal(r->bored);
    pthread_mutex_unlock(r->mx);
    return s;
}

void rollbuf_mark_done(rollbuf_ptr r)
{
    pthread_mutex_lock(r->mx);
    assert(!r->done);
    r->done = 1;
    pthread_cond_signal(r->bored);
    pthread_mutex_unlock(r->mx);
}

int rollbuf_is_done(rollbuf_ptr r)
{
    pthread_mutex_lock(r->mx);
    int x = r->done;
    pthread_mutex_unlock(r->mx);
    return x;
}

int rollbuf_get(rollbuf_ptr r, char * p, size_t s)
{
    pthread_mutex_lock(r->mx);
    // fprintf(stderr, "get(%zu): (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    if (!r->done && r->avail_to_read < ROLLBUF_ALIGNED_RETURNS) {
        // fprintf(stderr, "get(%zu): on hold (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
        r->empty_count++;
        pthread_cond_wait(r->bored, r->mx);
        // fprintf(stderr, "get(%zu): resumed (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    }
    assert(r->done || r->avail_to_read >= ROLLBUF_ALIGNED_RETURNS);
    if (r->done && !r->avail_to_read) {
        pthread_mutex_unlock(r->mx);
        return 0;
    }
    if (r->avail_to_read < ROLLBUF_ALIGNED_RETURNS)
        assert(r->done);
    size_t tail = r->alloc - (r->rhead - r->p);
    assert(s >= ROLLBUF_ALIGNED_RETURNS);
    s = MIN(r->avail_to_read, s);
    if (s >= ROLLBUF_ALIGNED_RETURNS) s -= s % ROLLBUF_ALIGNED_RETURNS;
    if (r->avail_to_read <= tail || s <= tail) {
        assert(s <= tail);
        memcpy(p, r->rhead, s);
        r->rhead += s;
    } else {
        assert(tail > 0);
        assert(s > tail);
        memcpy(p, r->rhead, tail);
#ifndef NDEBUG
        ptrdiff_t head = r->avail_to_read - tail;
        assert(head > 0);
        assert((size_t) (s - tail) <= head);
        // s = tail + MIN((size_t) head, s - tail);
#endif
        assert(s >= tail);
        memcpy(p + tail, r->p, s - tail);
        r->rhead = r->p + (s-tail);
    }
    if ((size_t) (r->rhead - r->p) == r->alloc) {
        r->rhead = r->p;
    }
    r->avail_to_read -= s;
    r->avail_to_write += s;
    /* Could be that someone is waiting for room to write data */
    pthread_cond_signal(r->bored);
    pthread_mutex_unlock(r->mx);
    return s;
}

int rollbuf_get2(rollbuf_ptr r, void ** p, size_t s)
{
    if (*p) {
        return rollbuf_get(r, *p, s);
    }
    assert(s == 0);     // does not really make sense otherwise
    pthread_mutex_lock(r->mx);
    if (!r->rbuf) {
        r->rbuf = malloc(ROLLBUF_READING_BUFFER_SIZE);
    }
    pthread_mutex_unlock(r->mx);
    /*
       if (s > ROLLBUF_READING_BUFFER_SIZE || s == 0)
       s = ROLLBUF_READING_BUFFER_SIZE;
       */
    *p = r->rbuf;
    return rollbuf_get(r, r->rbuf, ROLLBUF_READING_BUFFER_SIZE);
}

