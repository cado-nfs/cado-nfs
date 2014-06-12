#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>
#include <unistd.h>
#include <stdint.h>
#include "macros.h"
#include "ringbuf.h"
#include "memory.h"
#include "portability.h"

/* This is a hack. Define to 1 to disable */
#define RINGBUF_ALIGNED_RETURNS sizeof(uint32_t)

/* Length of one write in preempt buffer. Between 64 and 1024 Ko
   seems the best. */
#define PREEMPT_ONE_READ        (1UL<<20)

/* This has been tested with small buffer sizes, so as to stress the
 * pause/resume mechanism (based on pthread conditions). It seems to
 * work. Note though that this code heavily relies on the fact that there
 * is only _one_ thread reading data.
 */

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

/* refuse to fill memory with incoming data beyond this size. Incoming
 * data which can only be stored by expanding the ringbuf capacity beyond
 * this size are paused */
#define RINGBUF_MAX_SIZE        (1 << 26)

/* This is only for ringbuf_get2, for the auxiliary malloc()'ed area
 * which it returns. Note that this is not guided by efficiency needs. */
#define RINGBUF_READING_BUFFER_SIZE     (1 << 16)

/* must be called with mutex locked !!! */
static void ringbuf_grow__(ringbuf_ptr r, size_t claim)
{
    size_t newalloc = r->alloc;
    if (!claim) {
        newalloc += newalloc / 2;
        newalloc |= pagesize(); /* home-made wrapper */
        newalloc++;
    } else if (claim <= r->alloc) {
        return;
    } else {
        newalloc = claim;
        /* round to page size */
        newalloc--;
        newalloc |= pagesize();
        newalloc++;
    }

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

void ringbuf_init(ringbuf_ptr r, size_t claim)
{
    memset(r, 0, sizeof(ringbuf));
    pthread_mutex_init(r->mx, NULL);
    pthread_cond_init(r->bored, NULL);
    if (claim) {
        ringbuf_grow__(r, claim);       /* mutex unneeded within init */
    }
}

void ringbuf_clear(ringbuf_ptr r)
{
    /*
    // fprintf(stderr, "ringbuf: %d times full, %d times empty\n",
            r->full_count, r->empty_count);
            */
    pthread_mutex_destroy(r->mx);
    pthread_cond_destroy(r->bored);
    free(r->rbuf);
    free(r->p);
}



int ringbuf_put(ringbuf_ptr r, char * p, size_t s)
{
    pthread_mutex_lock(r->mx);
    // fprintf(stderr, "put(%zu): (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    for( ; s > r->avail_to_write ; ) {
        if (r->alloc >= RINGBUF_MAX_SIZE) {
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

        ringbuf_grow__(r, 0);
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

void ringbuf_mark_done(ringbuf_ptr r)
{
    pthread_mutex_lock(r->mx);
    assert(!r->done);
    r->done = 1;
    pthread_cond_signal(r->bored);
    pthread_mutex_unlock(r->mx);
}

int ringbuf_is_done(ringbuf_ptr r)
{
    pthread_mutex_lock(r->mx);
    int x = r->done;
    pthread_mutex_unlock(r->mx);
    return x;
}

int ringbuf_get(ringbuf_ptr r, char * p, size_t s)
{
    pthread_mutex_lock(r->mx);
    // fprintf(stderr, "get(%zu): (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    if (!r->done && r->avail_to_read < RINGBUF_ALIGNED_RETURNS) {
        // fprintf(stderr, "get(%zu): on hold (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
        r->empty_count++;
        pthread_cond_wait(r->bored, r->mx);
        // fprintf(stderr, "get(%zu): resumed (ravail: %zu, wavail: %zu)\n", s, r->avail_to_read, r->avail_to_write);
    }
    assert(r->done || r->avail_to_read >= RINGBUF_ALIGNED_RETURNS);
    if (r->done && !r->avail_to_read) {
        pthread_mutex_unlock(r->mx);
        return 0;
    }
    if (r->avail_to_read < RINGBUF_ALIGNED_RETURNS)
        assert(r->done);
    size_t tail = r->alloc - (r->rhead - r->p);
    assert(s >= RINGBUF_ALIGNED_RETURNS);
    s = MIN(r->avail_to_read, s);
    if (s >= RINGBUF_ALIGNED_RETURNS) s -= s % RINGBUF_ALIGNED_RETURNS;
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
        assert((size_t) (s - tail) <= (size_t) head);
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

int ringbuf_get2(ringbuf_ptr r, void ** p, size_t s)
{
    if (*p) {
        return ringbuf_get(r, *p, s);
    }
    assert(s == 0);     // does not really make sense otherwise
    pthread_mutex_lock(r->mx);
    if (!r->rbuf) {
        r->rbuf = malloc(RINGBUF_READING_BUFFER_SIZE);
    }
    pthread_mutex_unlock(r->mx);
    /*
       if (s > RINGBUF_READING_BUFFER_SIZE || s == 0)
       s = RINGBUF_READING_BUFFER_SIZE;
       */
    *p = r->rbuf;
    return ringbuf_get(r, r->rbuf, RINGBUF_READING_BUFFER_SIZE);
}

ssize_t ringbuf_feed_stream(ringbuf_ptr r, FILE * f)
{
    ssize_t nread = 0;

    /* We are the only thread decreasing the avail_to_write counter in
     * rb. So we may keep a copy of its value, which will always be a
     * lower bound, provided that we accurately report our decreases both
     * to our local value and  to the global counter.  */
    pthread_mutex_lock(r->mx);
    size_t local_w_avail = r->avail_to_write;
    pthread_mutex_unlock(r->mx);

    for( ; ; ) {
        /* Make sure our writing space in the buffer is not empty */
        if (local_w_avail == 0) {
            pthread_mutex_lock(r->mx);
            for( ; ! r->avail_to_write ; ) {
                pthread_cond_wait(r->bored, r->mx);
            }
            local_w_avail = r->avail_to_write;
            pthread_mutex_unlock(r->mx);
        }
        /* We may now fread() from f, but only up to the _contiguous_
         * amount which is available in the buffer. This entails some
         * intimate dialogue with the ringbuf internals, which
         * obviously isn't cool (well, in fact, this whole thing
         * could probably be considered within the ringbuf API, after
         * all ?) */
        /* We are the only thread likely to call ringbuf_grow__,
         * which is the only (internal) call tinkering with r->p (and
         * hence the validity of r->whead */
        size_t tail = r->alloc - (r->whead - r->p);

        tail = MIN(tail, local_w_avail);

        /* restrict to reads of some maximum size, or we'll be too
         * long delivering data to our customers */
        tail = MIN(tail, PREEMPT_ONE_READ);

        size_t s = fread(r->whead, 1, tail, f);
        nread += s;

        if (s) {
            pthread_mutex_lock(r->mx);
            r->whead += s;
            r->avail_to_read += s;
            local_w_avail = r->avail_to_write -= s;
            if ((size_t) (r->whead - r->p) == r->alloc) {
                r->whead = r->p;
            }
            /* Could be that someone is waiting for data to be read */
            pthread_cond_signal(r->bored);
            pthread_mutex_unlock(r->mx);
        } else if (feof(f)) {
            return nread;
        } else {
            return -1;
        }
    }
}

/* Search for character c, from position (offset) bytes into the readable
 * part of the input buffer.  Return the offset in bytes from the initial
 * segment to the matching byte (thus an integer >= offset), or -1 if the
 * byte could not be found.
 *
 * This might be called by the reader thread with the mutex unlocked,
 * under the condition that the reading thread is unique.
 */
int ringbuf_strchr(ringbuf_ptr r, int c, size_t offset)
{
    pthread_mutex_lock(r->mx);
    ASSERT_ALWAYS(offset <= r->avail_to_read);
    pthread_mutex_unlock(r->mx);
    size_t tail = r->alloc - (r->rhead - r->p);
    int s = offset;
    for(; (size_t) s < tail ; s++) {
        if (r->rhead[s] == c)
            return s;
    }
    tail = r->avail_to_read - s;
    for(int t = 0 ; (size_t) t < tail ; s++,t++) {
        if (r->p[t] == c)
            return s;
    }
    return -1;
}

int ringbuf_skip_get(ringbuf_ptr r, size_t s)
{
    pthread_mutex_lock(r->mx);
    size_t d = (r->rhead-r->p) + s;
    if (d >= r->alloc) {
        d -= r->alloc;
    }
    r->rhead = r->p + d;
    r->avail_to_read -= s;
    r->avail_to_write += s;
    /* Could be that someone is waiting for room to write data */
    pthread_cond_signal(r->bored);
    pthread_mutex_unlock(r->mx);
    return 0;
}
