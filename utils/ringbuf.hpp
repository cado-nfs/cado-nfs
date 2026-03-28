#ifndef CADO_UTILS_RINGBUF_HPP
#define CADO_UTILS_RINGBUF_HPP

#include "cado_config.h"        // IWYU pragma: keep

/* Example of a rotating, and reallocating buffer. A separate thread has
 * to fetch the data from the real source and fill the buffer. This way,
 * reads block as little as possible -- sort of a userland pipe, only we
 * get more control. Not clear we gain anything, but it's less fragile
 * than forking an external program.
 *
 * This has been tested with small buffer sizes, so as to stress the
 * pause/resume mechanism (based on pthread conditions). It seems to
 * work.
 */

#include <cstdio>
#include <cstdint>

#include <mutex>
#include <condition_variable>
#include <memory>
#include <algorithm>

#include <sys/types.h>

#include "macros.h"
#include "portability.h"

struct ringbuf {
    std::unique_ptr<char[]> p;
    size_t alloc = 0;
    size_t avail_to_read = 0;
    size_t avail_to_write = 0;
    const char * rhead = nullptr;
    char * whead = nullptr;
    mutable std::mutex mx;
    std::condition_variable bored;
    int empty_count = 0;
    int full_count = 0;
    int done:1 = 0;

    explicit ringbuf(size_t claim)
    {
        if (claim == 0)
            claim = PREEMPT_BUF;

        if (claim)
            grow(claim);       /* mutex unneeded within init */
    }

    private:
    /* must be called with mutex locked !!! */
    void grow(size_t claim)
    {
        size_t newalloc = alloc;
        if (!claim) {
            newalloc += newalloc / 2;
            newalloc |= pagesize(); /* home-made wrapper */
            newalloc++;
        } else if (claim <= alloc) {
            return;
        } else {
            newalloc = claim;
            /* round to page size */
            newalloc--;
            newalloc |= pagesize();
            newalloc++;
        }

        auto newp = std::make_unique<char[]>(newalloc);
        if (avail_to_read) {
            const size_t tail = alloc - (rhead - p.get());
            if (tail < avail_to_read) {
                std::copy_n(rhead, tail, newp.get());
                std::copy_n(p.get(), avail_to_read - tail, newp.get() + tail);
            } else {
                std::copy_n(rhead, avail_to_read, newp.get());
            }
        }
        std::swap(p, newp);
        rhead = p.get();
        whead = p.get() + avail_to_read;
        avail_to_write += newalloc - alloc;
        alloc = newalloc;
    }

    /* refuse to fill memory with incoming data beyond this size. Incoming
     * data which can only be stored by expanding the ringbuf capacity beyond
     * this size are paused */
    static constexpr size_t RINGBUF_MAX_SIZE = 1 << 26;

    /* This is a hack. Define to 1 to disable */
    static constexpr size_t RINGBUF_ALIGNED_RETURNS = sizeof(uint32_t);

    /* Length of one write in preempt buffer. Between 64 and 1024 Ko
       seems the best. */
    static constexpr size_t PREEMPT_ONE_READ = 1UL<<20;

    /* Length of preempt buffer. Must be a power of 2. */
    static constexpr size_t PREEMPT_BUF = 1<<22;

    public:

    size_t put(const char * data, size_t s)
    {
        ASSERT_ALWAYS(data);
        std::unique_lock ux(mx);
        // fprintf(stderr, "put(%zu): (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
        for( ; s > avail_to_write ; ) {
            if (alloc >= RINGBUF_MAX_SIZE) {
                if (s < alloc) {
                    /* Then we want to drain our pipe first. */
                    // fprintf(stderr, "put(%zu): on hold (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
                    full_count++;
                    bored.wait(ux);
                    // fprintf(stderr, "put(%zu): resuming (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
                    continue;
                } else {
                    /* Here there is no hope to drain our pipe, so we must
                     * exceed our desired window size. This is not much of a
                     * problem though, since the curl reading buffer is not
                     * expected to grow too large (or so we hope...)
                     */
                    // fprintf(stderr, "Warning: buffer growing beyond max size ! (previous=%zu)\n", alloc);
                }
            }
            grow(0);
        }

        const size_t tail = alloc - (whead - p.get());
        if (avail_to_write <= tail || s <= tail) {
            ASSERT(s <= avail_to_write);
            // s = MIN(avail_to_write, s);
            ASSERT(s <= tail);
            whead = std::copy_n(data, s, whead);
        } else {
            ASSERT(tail > 0);
            ASSERT(s > tail);
            ASSERT_ALWAYS(whead);
            std::copy_n(data, tail, whead);
#ifndef NDEBUG
            ptrdiff_t head = avail_to_write - tail;
            ASSERT(head > 0);
            ASSERT(s-tail <= (size_t) head);
#endif
            // s = tail + MIN((size_t) head, s - tail);
            whead = std::copy_n(data + tail, s - tail, p.get());
        }
        if (whead == p.get() + alloc)
            whead = p.get();
        avail_to_write -= s;
        avail_to_read += s;
        /* Could be that someone is waiting for data to be read */
        bored.notify_one();
        return s;
    }

    void mark_done()
    {
        const std::scoped_lock dummy(mx);
        ASSERT(!done);
        done = -1;
        bored.notify_all();
    }

    int is_done() const
    {
        const std::scoped_lock dummy(mx);
        return done;
    }

    size_t get(char * data, size_t s)
    {
        std::unique_lock ux(mx);
        // fprintf(stderr, "get(%zu): (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
        while (!done && avail_to_read < RINGBUF_ALIGNED_RETURNS) {
            // fprintf(stderr, "get(%zu): on hold (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
            empty_count++;
            bored.wait(ux);
            // fprintf(stderr, "get(%zu): resumed (ravail: %zu, wavail: %zu)\n", s, avail_to_read, avail_to_write);
        }
        ASSERT(done || avail_to_read >= RINGBUF_ALIGNED_RETURNS);
        if (done && !avail_to_read)
            return 0;
        if (avail_to_read < RINGBUF_ALIGNED_RETURNS)
            ASSERT(done);
        const size_t tail = alloc - (rhead - p.get());
        ASSERT(s >= RINGBUF_ALIGNED_RETURNS);
        s = std::min(avail_to_read, s);
        if (s >= RINGBUF_ALIGNED_RETURNS) s -= s % RINGBUF_ALIGNED_RETURNS;
        if (avail_to_read <= tail || s <= tail) {
            ASSERT(s <= tail);
            std::copy_n(rhead, s, data);
            rhead += s;
        } else {
            ASSERT(tail > 0);
            ASSERT(s > tail);
            std::copy_n(rhead, tail, data);
#ifndef NDEBUG
            ptrdiff_t head = avail_to_read - tail;
            ASSERT(head > 0);
            ASSERT((size_t) (s - tail) <= (size_t) head);
            // s = tail + MIN((size_t) head, s - tail);
#endif
            ASSERT(s >= tail);
            std::copy_n(p.get(), s - tail, data + tail);
            rhead = p.get() + (s-tail);
        }
        if (rhead == p.get() + alloc)
            rhead = p.get();
        avail_to_read -= s;
        avail_to_write += s;
        /* Could be that someone is waiting for room to write data */
        bored.notify_one();
        return s;
    }

    size_t get_avail_to_write_safe() const {
        const std::scoped_lock dummy(mx);
        return avail_to_write;
    }

    size_t get_avail_to_read_safe() const {
        const std::scoped_lock dummy(mx);
        return avail_to_read;
    }

    /* Equivalent of doing put for all bytes from the stdio stream f.
     *
     * The only difference is that this call does not automatically enlarge
     * the ring buffer, so an appropriate initial_size must have been
     * provided on initialization
     */
    ssize_t feed_stream(FILE * f)
    {
        size_t nread = 0;

        /* We are the only thread decreasing the avail_to_write counter in
         * rb. So we may keep a copy of its value, which will always be a
         * lower bound, provided that we accurately report our decreases both
         * to our local value and  to the global counter.  */
        size_t local_w_avail = get_avail_to_write_safe();

        for( ; ; ) {
            /* Make sure our writing space in the buffer is not empty */
            if (local_w_avail == 0) {
                std::unique_lock ux(mx);
                for( ; ! avail_to_write ; ) {
                    full_count++;
                    bored.wait(ux);
                }
                local_w_avail = avail_to_write;
            }
            /* We may now fread() from f, but only up to the _contiguous_
             * amount which is available in the buffer. This entails some
             * intimate dialogue with the ringbuf internals, which
             * obviously isn't cool (well, in fact, this whole thing
             * could probably be considered within the ringbuf API, after
             * all ?) */
            /* We are the only thread likely to call grow(),
             * which is the only (internal) call tinkering with p (and
             * hence the validity of whead */
            size_t tail = alloc - (whead - p.get());

            tail = std::min(tail, local_w_avail);

            /* restrict to reads of some maximum size, or we'll be too
             * long delivering data to our customers */
            tail = std::min(tail, PREEMPT_ONE_READ);

            const size_t s = fread(whead, 1, tail, f);
            nread += s;

            if (s) {
                const std::scoped_lock dummy(mx);
                whead += s;
                avail_to_read += s;
                local_w_avail = avail_to_write -= s;
                if (whead == p.get() + alloc)
                    whead = p.get();
                /* Could be that someone is waiting for data to be read */
                bored.notify_one();
            } else if (feof(f)) {
                /* this interface returns ssize_t, which is a pity since it
                 * is not the same signedness as nread...
                 */
                return static_cast<ssize_t>(nread);
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
     *
     * FIXME: actually calling this function with the mutex lock can't
     * work because of get_avail_to_read_safe(). This all seems
     * contradictory, though. If we assume (mutex || single-thread), then
     * why do we even bother locking???
     */
    int strchr(int c, size_t offset) const
    {
        ASSERT_ALWAYS(offset <= get_avail_to_read_safe());
        size_t tail = alloc - (rhead - p.get());
        size_t s = offset;
        for(; s < tail ; s++) {
            if (rhead[s] == c)
                return static_cast<int>(s);
        }
        tail = avail_to_read - s;
        for(int t = 0 ; static_cast<size_t>(t) < tail ; s++,t++) {
            if (p[t] == c)
                return static_cast<int>(s);
        }
        return -1;
    }

    int skip_get(size_t s)
    {
        const std::scoped_lock dummy(mx);
        ASSERT_ALWAYS(s <= alloc);
        size_t d = (rhead - p.get()) + s;
        if (d >= alloc)
            d -= alloc;
        rhead = p.get() + d;
        avail_to_read -= s;
        avail_to_write += s;
        /* Could be that someone is waiting for room to write data */
        bored.notify_one();
        return 0;
    }

    /* The iterator does not *consume* the contents of the ring buffer,
     * which is in contrast with the get() method. n consecutive calls to
     * *it++ must be followed by a skip_get(n) call.
     */
    struct const_iterator {
        using value_type = char;
        ringbuf const & r;
        const char * p;
        const_iterator(ringbuf const & r, const char * p)
            : r(r)
            , p(p)
        {}
        explicit const_iterator(ringbuf const & r)
            : const_iterator(r, r.rhead)
        {}
        const_iterator(const_iterator const &) = default;
        const_iterator(const_iterator &&) = default;
        const_iterator& operator=(const_iterator const & a)
        {
            if (&a != this) {
                ASSERT_ALWAYS(&r == &a.r);
                p = a.p;
            }
            return *this;
        }
        /* the old RINGBUF_GET_ONE_BYTE macro returned an unsigned char */
        unsigned char operator*() const { return *p; }
        const_iterator& operator++() {
            p++;
            if (p >= r.p.get() + r.alloc)
                p = r.p.get();
            return *this;
        }
        const_iterator operator++(int) { auto c = *this; ++*this; return c; }
    };

    const_iterator begin() const {
        return { *this, rhead };
    }
    const_iterator end() const {
        return { *this, p.get() + ((rhead - p.get()) + avail_to_read) % alloc };
    }

};

typedef ringbuf * ringbuf_ptr;
typedef const ringbuf * ringbuf_srcptr;

/* A quick accessor macro which does a 1-byte fetch from the ring buffer.
 * We must be sure that the get will succeed, and we must provide an
 * auxiliary pointer (here s_) which will be updated by the macro. s_ has
 * to be set to r_->rhead originally.
 *
 * n successive calls to RINGBUF_GET_ONE_BYTE must be followed by a call to
 * ringbug_skip_get(r_, n)
 */
#define RINGBUF_GET_ONE_BYTE(c_, r_, s_) do {				\
    (c_) = (unsigned char) *(s_)++;					\
    if ((s_) >= (r_)->p.get() + (r_)->alloc) {				\
        (s_) = (r_)->p.get();						\
    }									\
} while (0)



#endif	/* CADO_UTILS_RINGBUF_HPP */
