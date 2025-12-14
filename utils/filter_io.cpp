#include "cado.h" // IWYU pragma: keep

#include <cerrno>                     // for errno
#include <climits>                    // for INT_MAX
#include <cstdio>                     // for fprintf, stderr, stdout, FILE
#include <cstdlib>                    // for abort, malloc, realloc, free
#include <cstring>                    // for memset, memcpy, strcmp, strerror
#include <cstdint>
#include <ctime>

#include <algorithm>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <condition_variable>
#if !defined(__x86_64) && !defined(__i386)
#include <atomic>
#endif

#include <pthread.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h>              // for rusage // IWYU pragma: keep
#endif
#include <gmp.h>

#include "barrier.h"                   // for barrier_destroy, barrier_init
#include "bit_vector.h"
#include "cado_popen.h"                // for cado_pclose2, cado_popen
#include "filter_io.h"
#include "gzip.h"                      // prepare_grouped_command_lines
#include "macros.h"                    // for ASSERT_ALWAYS, ASSERT, UNLIKELY
#include "ringbuf.h"                   // for ringbuf_s, ringbuf_ptr, RINGBU...
#include "stats.h"                     // stats_data_t
#include "portability.h" // sleep // IWYU pragma: keep
#include "runtime_numeric_cast.hpp"
#include "timing.h"
#include "typedefs.h"
#include "utils_cxx.hpp"

/* This is a configuration variable which may be set by the caller (it's
 * possible to bind it to a command-line argument)
 */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
int filter_rels_force_posix_threads = 0;

/*{{{ inflight buffer. See filter_io.tex for documentation. */

/* macros for the two locking models */

struct ifb_locking_posix {/*{{{*/
    static const int max_supported_concurrent = INT_MAX;
    template<typename T> struct critical_datatype {
        class t {
            T x;
            public:
            T load() const { return x; }
            void store(T a) { x = a; }
            explicit t(T const& a) : x(a) {}
            explicit t() : x(0) {}
            t& operator=(T const& a) { x = a; return *this; }
            T increment() { return x++; }
        };
    };
    using lock_t = std::mutex ;
    using cond_t = std::condition_variable  ;
    static void lock(lock_t * m) { m->lock(); }
    static void unlock(lock_t * m) { m->unlock(); }
    static void wait(cond_t * c, lock_t * m) {
        std::unique_lock<std::mutex> foo(*m, std::adopt_lock);
        c->wait(foo);
        foo.release();
    }
    static void signal(cond_t * c) { c->notify_one(); }
    static void signal_broadcast(cond_t * c) { c->notify_all(); }
    static int isposix() { return 1; }
};
/*}}}*/

/*{{{ define NANOSLEEP */
/* The realistic minimal non-CPU waiting with nanosleep is about 10 to 40
 * microseconds (1<<13 for nanosleep).  But all the I/O between the
 * threads have been buffered, and a thread does a nanosleep only if its
 * buffer is empty.  So I use here ~2ms (1<<21) to optimize CPU
 * scheduler.  Max pause is about 4 to 8ms (1<<22, 1<<23); above that,
 * the program is slowed down.
 */
#ifndef HAVE_NANOSLEEP
#ifdef HAVE_USLEEP
#define NANOSLEEP() usleep((unsigned long) (1<<21 / 1000UL))
#else
#define NANOSLEEP() sleep(0)
#endif
#else
static const struct timespec wait_classical = { 0, 1<<21 };
#define NANOSLEEP() nanosleep(&wait_classical, NULL)
#endif
/*}}}*/

struct ifb_locking_lightweight {/*{{{*/
    /* we don't support several threads wanting to write to the same
     * location (we could, if we were relying on atomic compare and swap,
     * for instance) */
    static const int max_supported_concurrent = 1;
    template<typename T> struct critical_datatype {
        /* See bug #30068
         *
         * The total store ordering on x86 implies that we can play very
         * dirty games with the completed[] and scheduled[] arrays. The
         * underlying assumptions need not be true in general, and
         * definitely do not hold on arm64.
         *
         * Ideally, there would be a way to qualify our operations on the
         * atomic type that resolve to no emitted code at all if the
         * hardware memory model is x86. But I can't find a way to do
         * that.
         */
#if !defined(__x86_64) && !defined(__i386)
        class t : private std::atomic<T> {
            using super = std::atomic<T>;
            public:
            T load() const { return super::load(std::memory_order_acquire); }
            explicit t(T const& a) : super(a) {}
            explicit t() : super(0) {}
            void store(T a) { super::store(a, std::memory_order_release); }
            t& operator=(T const& a) { store(a); return *this; }
            T increment() { return super::fetch_add(1, std::memory_order_acq_rel); }
        };
#else
        class t {
            volatile T x;
            public:
            T load() const { return x; }
            explicit t(T const& a) : x(a) {}
            explicit t() : x(0) {}
            void store(T a) { x = a; }
            t& operator=(T const& a) { store(a); return *this; }
            /* c++20 frowns upon volatile. Well, it kinda forces us to
             * use atomics, in fact. Let's silence the issue for the
             * moment. It may well be that the correct way to go is the
             * code branch above. But I still long for the optimal way to
             * write things so that there's a 1-to-1 correspondence with
             * the simple and easy code here.
             */
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-volatile"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#endif
            T increment() { return x++; }
#ifdef __clang__
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
        };
#endif
    };
    using lock_t = int;
    using cond_t = int;
    template<typename T> static T next(T a, int) { return a; }
    static void lock(lock_t *) {}
    static void unlock(lock_t *) {}
    static void wait(cond_t *, lock_t *) { NANOSLEEP(); }
    static void signal(cond_t *) {}
    static void signal_broadcast(cond_t *) {}
    static int isposix() { return 0; }
};
/*}}}*/

/* {{{ status table (utility for inflight_rels_buffer).
 *
 * In fact, when we use simple busy waits, we are restricted to
 * scheduled[k]==completed[k]+(0 or 1), and keeping track of the
 * processing level is useless. So we provide a trimmed-down
 * specialization for this case.
 *
 * the status table depends on the maximum number of threads per step. We
 * make it depend on the locking backend instead, for simplicity. */
template<typename locking>
struct status_table {
    using csize_t = typename locking::template critical_datatype<size_t>::t;
    typename locking::template critical_datatype<int8_t>::t x[SIZE_BUF_REL];
    /* {{{ ::catchup() (for ::schedule() termination) */
    void catchup(csize_t & last_completed, size_t last_scheduled, int level) {
        size_t c = last_completed.load();
        for( ; c < last_scheduled ; c++) {
            if (x[c & (SIZE_BUF_REL-1)].load() < level)
                break;
        }
        last_completed.store(c);
    }
    /*}}}*/
    /*{{{ ::catchup_until_mine_completed() (for ::complete()) */
    /* (me) is the absolute relation index of the relation I'm currently
     * processing (the value of schedule[k] when it was called prior to
     * giving me this relation to process).
     */
    void catchup_until_mine_completed(csize_t & last_completed, size_t me, int level) {
        const size_t slot = me & (SIZE_BUF_REL-1);
        size_t c = last_completed.load();
        ASSERT(x[slot].load() == (int8_t) (level-1));
        /* The big question is how far we should go. By not exactly answering
         * this question, we avoid the reading of scheduled[k], which is good
         * because it is protected by m[k-1]. And even if we could consider
         * doing a rwlock for reading this, it's too much burden. So we leave
         * open the possibility that many relation slots ahead of us already
         * have x[slot] set to k, yet we do not increment
         * completed[k] that far. This will be caught later on by further
         * processing at this level.
         *
         * This logic is problematic regarding termination, though. See the
         * termination code in ::complete()
         */
        for( ; c < me ; c++) {
            if (x[c & (SIZE_BUF_REL-1)].load() < level)
                break;
        }
        last_completed.store(c + (c == me));
        x[slot].increment();
        ASSERT(x[slot].load() == (int8_t) (level));
    }
    /*}}}*/
    void update_shouldbealreadyok(size_t slot, int level) {
        if (level < 0) {
            x[slot & (SIZE_BUF_REL-1)].store(level);
        } else {
            ASSERT(x[slot & (SIZE_BUF_REL-1)].load() == level);
        }
    }
};

template<>
struct status_table<ifb_locking_lightweight> {
    using csize_t = ifb_locking_lightweight::critical_datatype<size_t>::t;
    static void catchup(csize_t & last_completed, size_t last_scheduled, int) {
        ASSERT_ALWAYS(last_completed.load() == last_scheduled);
    }
    static void catchup_until_mine_completed(csize_t & last_completed, size_t, int) {
        last_completed.increment();
    }
    void update_shouldbealreadyok(size_t, int) {}
};
/* }}} */

/* {{{ inflight_rels_buffer: n-level buffer, with underyling locking
 * mechanism specified by the template class.  */
template<typename locking, int n>
struct inflight_rels_buffer {
    cado_nfs::barrier sync_point;
    std::unique_ptr<earlyparsed_relation[]> rels;        /* always malloc()-ed to SIZE_BUF_REL,
                                           which is a power of two */
    /* invariant:
     * scheduled_0 >= ... >= completed_{n-1} >= scheduled_0 - SIZE_BUF_REL
     */
    typename locking::template critical_datatype<size_t>::t completed[n];
    typename locking::template critical_datatype<size_t>::t scheduled[n];
    status_table<locking> status;
    typename locking::lock_t m[n];
    typename locking::cond_t bored[n];
    int active[n];     /* number of active threads */

    explicit inflight_rels_buffer(int nthreads_total);
    ~inflight_rels_buffer();

    inflight_rels_buffer(inflight_rels_buffer const&) = delete;
    inflight_rels_buffer(inflight_rels_buffer &&) = delete;
    inflight_rels_buffer& operator=(inflight_rels_buffer const&) = delete;
    inflight_rels_buffer& operator=(inflight_rels_buffer &&) = delete;

    void drain();
    earlyparsed_relation_ptr schedule(int);
    void complete(int, earlyparsed_relation_srcptr);
    
    /* computation threads joining the computation are calling these */
    void enter(int k) {
        locking::lock(m+k); active[k]++; locking::unlock(m+k);
        sync_point.arrive_and_wait();
    }
    /* leave() is a no-op, since active-- is performed as part of the
     * normal drain() call */
    void leave(int) { }

    /* The calling scenario is as follows.
     *
     * For the owner thread.
     *  - constructor
     *  - start workers.
     *  - enter(0)
     *  - some schedule(0) / complete(0) for relations which get fed in.
     *  - drain() once all are produced
     *  - leave(0)
     *
     * For the workers (there may be more at each level if
     * ifb_locking_posix is used):
     *  - enter(k)
     *  - a loop on with schedule(k) / complete(k), exiting when
     *    schedule() returns NULL.
     *  - leave(k)
     *
     * Currently the owner thread is weakly assumed to be the only
     * level-0 thread, but that does not seem to be an absolute necessity
     * from the design. Additional level-0 threads would induce a loop
     * similar to other worker threads, but the fine points haven't been
     * considered yet.
     *
     * The current implementation has leave() a no-op, and uses drain()
     * at the owner thread to to a shutdown. This could change.
     */
};
/* }}} */

/* {{{ instantiations for the locking buffer methods */

/*{{{ ::inflight_rels_buffer() */
template<typename locking, int n>
inflight_rels_buffer<locking, n>::inflight_rels_buffer(int nthreads_total)
    : sync_point(nthreads_total)
    , rels(new earlyparsed_relation[SIZE_BUF_REL])
{
    std::fill(completed, completed + n, 0UL);
    std::fill(scheduled, scheduled + n, 0UL);
    std::fill(active, active + n, 0);
    // std::fill(rels.get(), rels.get() + SIZE_BUF_REL, 0);
    memset(rels.get(), 0, SIZE_BUF_REL * sizeof(earlyparsed_relation));
}/*}}}*/
/*{{{ ::drain() */
/* This belongs to the buffer closing process.  The out condition of this
 * call is that all X(k) for k>0 terminate.  This call (as well as
 * init/clear) must be called on the producer side (step 0) (in a
 * multi-producer context, only one thread is entitled to call this) */
template<typename locking, int n>
void inflight_rels_buffer<locking, n>::drain()
{
    // size_t c = completed[0];
    active[0]--;

    for(int k = 0 ; k < n ; k++) {
        locking::lock(m + k);
        while(active[k]) {
            locking::wait(bored + k, m + k);
        }
        completed[k].store(SIZE_MAX);
        locking::signal_broadcast(bored + k);
        locking::unlock(m + k);
    }
}
/*}}}*/
/*{{{ ::~inflight_rels_buffer() */
/* must be called on producer side */
template<typename locking, int n>
inflight_rels_buffer<locking, n>::~inflight_rels_buffer()
{
    for(int i = 0 ; i < n ; i++) {
        ASSERT_ALWAYS_NOTHROW(active[i] == 0);
    }
    for(size_t i = 0 ; i < SIZE_BUF_REL ; i++) {
        if (rels[i]->primes != rels[i]->primes_data) {
            free(rels[i]->primes);
        }
        if (rels[i]->line) free(rels[i]->line);
        if (rels[i]->sm_alloc) {
            for(int j = 0 ; j < rels[i]->sm_alloc ; j++) {
                mpz_clear(rels[i]->sm[j]);
            }
            free(rels[i]->sm);
        }
        memset(rels[i], 0, sizeof(rels[i]));
    }
}
/*}}}*/

/*{{{ ::schedule() (generic) */
/* Schedule a new relation slot for processing at level k.
 *
 * This call may block until a relation is processed by level k-1 (or, if
 * k==0, until a slot is made available in the relation buffer).
 *
 * The relation is free for use by the current (consumer) thread until it
 * calls inflight_rels_buffer_completed. When the producing stream ends,
 * this function returns NULL. */
template<typename locking, int n>
earlyparsed_relation_ptr
inflight_rels_buffer<locking, n>::schedule(int k)
{
    int const prev = k ? (k-1) : (n-1);
    // coverity[result_independent_of_operands]
    ASSERT(active[k] <= locking::max_supported_concurrent);
    size_t s;
    size_t const a = k ? 0 : SIZE_BUF_REL;
    /* in 1-thread scenario, scheduled[k] == completed[k] */
    locking::lock(m + prev);
    if (locking::max_supported_concurrent == 1) {       /* static check */
        /* can't change */
        s = scheduled[k].load();
        while(s == a + completed[prev].load()) {
            locking::wait(bored + prev, m + prev);
        }
    } else {
        while((s=scheduled[k].load()) == a + completed[prev].load()) {
            locking::wait(bored + prev, m + prev);
        }
    }
    /* when completed[prev] == SIZE_MAX, the previous-level workers
     * are creating spuriouss relation created to trigger termination.
     * In this case, scheduled[prev] is safe to read now. we use it
     * as a marker to tell whether there's still work ahead of us, or
     * not.  */
    if (UNLIKELY(completed[prev].load() == SIZE_MAX) && scheduled[prev].load() == s) {
        /* prepare to return */
        /* note that scheduled[k] is *not* bumped here */
        locking::unlock(m + prev);
        /* we emulate the equivalent of ::complete(), and terminate */
        locking::lock(m + k);
        status.catchup(completed[k], s, k);
        active[k]--;
        locking::signal_broadcast(bored + k);
        locking::unlock(m + k);
        return nullptr;
    }
    // ASSERT(scheduled[k] < a + completed[prev]);
    scheduled[k].increment();
    const size_t slot = s & (SIZE_BUF_REL - 1);
    earlyparsed_relation_ptr rel = rels[slot];
    status.update_shouldbealreadyok(s, k-1);
    locking::unlock(m + prev);
    return rel;
}
/*}}}*/

/*{{{ ::complete() (generic)*/
template<typename locking, int n>
void
inflight_rels_buffer<locking, n>::complete(int k,
        earlyparsed_relation_srcptr rel)
{
    // coverity[result_independent_of_operands]
    ASSERT(active[k] <= locking::max_supported_concurrent);
    const int slot = rel - (earlyparsed_relation_srcptr) rels.get();

    locking::lock(m + k);

    size_t my_absolute_index;
    if (locking::max_supported_concurrent == 1) {       /* static check */
        my_absolute_index = completed[k].load();
    } else {
        /* recover the integer relation number being currently processed from
         * the one modulo SIZE_BUF_REL.
         *
         * We have (using ck = completed[k]):
         *          ck <= zs < ck + N
         *          ck <= s+xN < ck + N <= s+(x+1)N
         *          xN < ck-s + N <= (x+1) N
         *
         */
        const size_t c = completed[k].load();
        my_absolute_index = slot;
        my_absolute_index += ((c - slot + SIZE_BUF_REL - 1) & -SIZE_BUF_REL);
    }

    /* morally, this is completed[k]++ */
    status.catchup_until_mine_completed(completed[k], my_absolute_index, k);
    locking::signal_broadcast(bored + k);
    locking::unlock(m + k);
}
/*}}}*/

/*}}}*/

/* }}} */

/************************************************************************/

/* {{{ early parsing routines, for filling the earlyparsed_relation
 * structures. */

/* malloc()'s are avoided as long as there are less than NB_PRIMES_OPT in
 * the relation
 */
void realloc_buffer_primes(earlyparsed_relation_ptr buf)
{
    if (buf->nb_alloc == NB_PRIMES_OPT) {
	buf->nb_alloc += buf->nb_alloc >> 1;
	prime_t *p = buf->primes;
	buf->primes = (prime_t*) malloc(buf->nb_alloc * sizeof(prime_t));
	if (!buf->primes) {
            fprintf(stderr, "malloc failure: %s\n", __func__);
            abort();
        }
	memcpy(buf->primes, p, NB_PRIMES_OPT * sizeof(prime_t));
    } else {
	buf->nb_alloc += buf->nb_alloc >> 1;
	checked_realloc(buf->primes, buf->nb_alloc);
	if (!buf->primes) {
            fprintf(stderr, "malloc failure: %s\n", __func__);
            abort();
        }
    }
#if DEBUG >= 2
    fprintf(stderr, "realloc_buffer_primes: num=%" PRIu64 " nb_alloc=%u\n",
	    buf->num, buf->nb_alloc);
#endif
}

#define PARSER_ASSERT_ALWAYS(got, expect, sline, ptr) do {		\
    if (UNLIKELY((got)!=(expect))) {					\
        fprintf(stderr, "Parse error in %s at %s:%d\n"			\
                "Expected character '%c', got '%c'"			\
                " after reading %zd bytes from:\n"			\
                "%s\n",							\
                __func__,__FILE__,__LINE__,				\
                expect, got, (ptr) - (sline), (sline));			\
        abort();							\
    }									\
} while (0)

/* decimal a,b is now an exception, and only for the first step of dup2
 * (when primes are not yet renumbered)
 *
 * "normal" processing of a,b, i.e. everywhere except for the very first
 * stage of dup2, involves reading a,b in hex.
 */
static inline int earlyparser_inner_read_ab_withbase(ringbuf_ptr r, const char ** pp, earlyparsed_relation_ptr rel, const uint64_t base)
{
    const char * p = *pp;
    int c;
    uint64_t v,w;
    RINGBUF_GET_ONE_BYTE(c, r, p);
    int negative = 0;
    if (c == '-') {
        negative = 1;
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    for (w = 0; (v = ugly[c]) < base;) {
        w = w * base + v;
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    PARSER_ASSERT_ALWAYS(c, ',', *pp, p);
    if (negative)
        rel->a = -runtime_numeric_cast<int64_t>(w);
    else
        rel->a = runtime_numeric_cast<int64_t>(w);
    RINGBUF_GET_ONE_BYTE(c, r, p);
    for (w = 0; (v = ugly[c]) < base;) {
        w = w * base + v;
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    *pp = p;
    rel->b = w;
    return c;
}
static int earlyparser_inner_read_ab_decimal(ringbuf_ptr r, const char ** pp, earlyparsed_relation_ptr rel)
{
    return earlyparser_inner_read_ab_withbase(r, pp, rel, 10);
}

static int earlyparser_inner_read_ab_hexa(ringbuf_ptr r, const char ** pp, earlyparsed_relation_ptr rel)
{
    return earlyparser_inner_read_ab_withbase(r, pp, rel, 16);
}

static int earlyparser_inner_read_prime(ringbuf_ptr r, const char ** pp, uint64_t * pr)
{
    uint64_t v,w;
    int c;
    const char * p = *pp;
#define BASE 16         /* primes are always written in hexa */
    /* copy-paste code blob above */
    RINGBUF_GET_ONE_BYTE(c, r, p);
    for (w = 0; (v = ugly[c]) < BASE;) {
        w = w * BASE + v;       /* *16 ought to be optimized */
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
#undef BASE
    *pr = w;
    *pp = p;
    return c;
}

static int earlyparser_inner_read_mpz(ringbuf_ptr r, const char ** pp, mpz_t z)
{
    uint64_t v;
    int c;
    const char * p = *pp;
#define BASE 10
    /* copy-paste code blob above */
    RINGBUF_GET_ONE_BYTE(c, r, p);
    for (mpz_set_ui(z, 0); (v = ugly[c]) < BASE;) {
        mpz_mul_ui(z, z, 10);
        mpz_add_ui(z, z, v);
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
#undef BASE
    *pp = p;
    return c;
}

static int earlyparser_inner_skip_ab(ringbuf_ptr r, const char ** pp)
{
    const char * p = *pp;
    int c = 0;
    for( ; (c != ':') ; ) {
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    *pp = p;
    return c;
}

static int
earlyparser_inner_read_active_sides(ringbuf_ptr r, const char ** pp, earlyparsed_relation_ptr rel)
{
    const char * p = *pp;
    if (*p == '@') {
        p++;
        uint64_t v,w;
        int c;
#define BASE 10
        /* copy-paste code blob above */
        RINGBUF_GET_ONE_BYTE(c, r, p);
        for (w = 0; (v = ugly[c]) < BASE;) {
            w = w * BASE + v;       /* *16 ought to be optimized */
            RINGBUF_GET_ONE_BYTE(c, r, p);
        }
        rel->active_sides[0] = runtime_numeric_cast<int>(w);
        PARSER_ASSERT_ALWAYS(c, ',', r->rhead, p);
        p++;
        RINGBUF_GET_ONE_BYTE(c, r, p);
        for (w = 0; (v = ugly[c]) < BASE;) {
            w = w * BASE + v;       /* *16 ought to be optimized */
            RINGBUF_GET_ONE_BYTE(c, r, p);
        }
        rel->active_sides[1] = runtime_numeric_cast<int>(w);
#undef BASE
    } else {
        rel->active_sides[0] = 0;
        rel->active_sides[1] = 1;
    }
    *pp = p;
    return 1;
}


struct prime_t_cmp{
    static int cmp(prime_t const * a, prime_t const * b)
    {
        int r = (a->side > b->side) - (b->side > a->side);
        if (r) return r;
        r = (a->p > b->p) - (b->p > a->p);
        return r;
    }
    bool operator()(prime_t const & a, prime_t const & b) {
        return cmp(&a, &b) < 0;
    }
};

struct prime_t_cmp_indices{
    static int cmp(prime_t const * a, prime_t const * b)
    {
        int const r = (a->h > b->h) - (b->h > a->h);
        return r;
    }
    bool operator()(prime_t const & a, prime_t const & b) {
        return cmp(&a, &b) < 0;
    }
};

/* This earlyparser is the pass which is used to perform the renumbering
 * of the relation files. This has some implications.
 *  - a,b are in decimal (for factoring -- for FFS it's hex anyway)
 *  - the primes need not be sorted (XXX update: now they are, so we may
 *    save time here).
 *  - we have 'nb_polys' input sides, but a flat, unique side on output.
 *  - at some point in the process we must compute r=a/b mod p, and this
 *    is going to end up in the renumber table. We use the .side field in the
 *    prime_t structure to pass information to the routine which does this.
 */
static inline int earlyparser_abp_withbase(earlyparsed_relation_ptr rel, ringbuf_ptr r, uint64_t base);
static int earlyparser_abp_decimal(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_abp_withbase(rel, r, 10);
}
static int earlyparser_abp_hexa(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_abp_withbase(rel, r, 16);
}

static
unsigned int sort_and_compress_rel_primes(prime_t * primes, unsigned int n)
{
    prime_t_cmp P;
    /* sort ; note that we're sorting correctly w.r.t the side as
     * well. We could of course exploit the fact that the sides
     * themselves are always in order, but the benefit is likely to
     * be small. Anyway this branch is not critical, as we prefer
     * to have las create sorted files */
    std::sort(primes, primes + n, P);
    /* compress. idiomatic albeit subtle loop. */
    unsigned int i,j;
    prime_t * qq = primes;
    for(i = j = 0 ; i < n ; j++) {
        qq[j] = qq[i];
        for(i++ ; i < n && !P(qq[j], qq[i]) ; i++) {
            qq[j].e += qq[i].e;
        }
    }
    return j;
}

static
unsigned int sort_and_compress_rel_indices(prime_t * primes, unsigned int n)
{
    prime_t_cmp_indices P;
    /* sort ; note that we're sorting correctly w.r.t the side as
     * well. We could of course exploit the fact that the sides
     * themselves are always in order, but the benefit is likely to
     * be small. Anyway this branch is not critical, as we prefer
     * to have las create sorted files */
    std::sort(primes, primes + n, P);
    /* compress. idiomatic albeit subtle loop. */
    unsigned int i,j;
    prime_t * qq = primes;
    for(i = j = 0 ; i < n ; j++) {
        qq[j] = qq[i];
        for(i++ ; i < n && !P(qq[j], qq[i]) ; i++) {
            qq[j].e += qq[i].e;
        }
    }
    return j;
}

static inline int earlyparser_abp_withbase(earlyparsed_relation_ptr rel, ringbuf_ptr r, uint64_t base)
{
    const char * p = r->rhead;

    int c = earlyparser_inner_read_ab_withbase(r, &p, rel, base);

    earlyparser_inner_read_active_sides(r, &p, rel);

    unsigned int n = 0;

    uint64_t last_prime = 0;
    int sorted = 1;
    int side = -1;

    for( ; ; ) {
        uint64_t pr;
        if (c == '\n') break;
        else if (c == ':') { last_prime = 0; side++; }
        else PARSER_ASSERT_ALWAYS(c, ',', r->rhead, p);
        c = earlyparser_inner_read_prime(r, &p, &pr);
	// it can be that c = ':', during the descent or more often in MNFS
	// and this implies pr = 0 which is used to detect this case
	// FIXME: do we need a specific ASSERT in the normal case?
	if(pr == 0)
	    continue;
        // not enforcing at the moment.
        // ASSERT_ALWAYS(pr >= last_prime);        /* relations must be sorted */
        sorted = sorted && pr >= last_prime;
        if (n && pr == rel->primes[n-1].p) {
            rel->primes[n-1].e++;
        } else {
            if (rel->nb_alloc == n) realloc_buffer_primes(rel);
            // rel->primes[n++] = (prime_t) { .h = (index_t) side,.p = (p_r_values_t) pr,.e = 1};
            rel->primes[n].side = side;
            rel->primes[n].p = (p_r_values_t) pr;
            rel->primes[n].e = 1;
            n++;
        }
        last_prime = pr;
    }
    if (!sorted)
        n = sort_and_compress_rel_primes(rel->primes, n);
    rel->nb = n;

    return 1;
}



static int
earlyparser_ab(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    earlyparser_inner_read_ab_hexa(r, &p, rel);
    earlyparser_inner_read_active_sides(r, &p, rel);

    return 1;
}

static int
earlyparser_line(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    const char *p = r->rhead;

    int i = 0;
    int n = 0;
    /* the total number of commas an colons in a relation with k textual
     * sides is
     * 1 comma + k * (1 colon + max(nprimes-1,0) commas)
     *
     * If we assume that no side has zero prime, the following invariant holds:
     *      ncolon+ncommas = nprimes+1, whence nprimes = ncolons+ncommas-1.
     * With MNFS, we cannot assume anymore that no side has zero prime. But we
     * can still have this invariant if we consider any number of ... colons as
     * only one colon.
     */
    int prev_is_colon = 0;
    for(int c = 0 ; ; )
    {
      if (c == '\n')
        break;
      else if (c == ',')
        n++;
      else if (c == ':')
      {
        if (prev_is_colon == 0)
          n++;
        prev_is_colon = 1;
      }
      else
        prev_is_colon = 0;
      RINGBUF_GET_ONE_BYTE(c, r, p);
      rel->line[i++] = c;
      ASSERT_ALWAYS(i < RELATION_MAX_BYTES);
    }
    rel->line[i++] = '\0';
    ASSERT_ALWAYS(i < RELATION_MAX_BYTES);
    rel->nb = n - 1;    /* see explanation above */

    return 1;
}

/* e.g. for dup1 for FFS with -abhexa command-line flag */
static int
earlyparser_abline_hexa(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    earlyparser_inner_read_ab_hexa(r, &p, rel);
    return earlyparser_line(rel, r);
}

/* e.g. for dup1 */
static int
earlyparser_abline_decimal(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    earlyparser_inner_read_ab_decimal(r, &p, rel);
    return earlyparser_line(rel, r);
}


/* Note: for these routine, the sorting of the primes is not considered, and at
 * least for the 1st pass of purge, so far we've been using this code on
 * unsorted relations.
 */
static int
earlyparser_index_maybeabhexa(earlyparsed_relation_ptr rel, ringbuf_ptr r,
        int parseab, int parsesm, int sort)
{
    const char *p = r->rhead;

    if (parseab) {
        earlyparser_inner_read_ab_hexa(r, &p, rel);
    } else {
        earlyparser_inner_skip_ab(r, &p);
    }
    
    unsigned int n = 0;
    int is_sorted = 1;

    char const next_delim = parsesm ? ':' : '\n';
    int c = '\0';
    for( ; ; ) {
        uint64_t pr;
        int sgn = 1;
        if (c == next_delim) break;
        /* In all cases, we want to stop processing at newlines, in case
         * we expect to find SM info but in fact it's missing */
        if (c == '\n') break;
        if (p[0] == '-') {
        //if (c == '-') {
            sgn = -1;
            RINGBUF_GET_ONE_BYTE(c, r, p);
        }
        c = earlyparser_inner_read_prime(r, &p, &pr);
        if (sort && is_sorted && n && pr < rel->primes[n-1].h)
            is_sorted = 0;
        if (n && pr == rel->primes[n-1].h) {
            rel->primes[n-1].e += sgn;
        } else {
            if (rel->nb_alloc == n) realloc_buffer_primes(rel);
            rel->primes[n].h = (index_t) pr;
            rel->primes[n].p = 0;
            rel->primes[n].e = sgn;
            n++;
        }
    }
    if (sort && !is_sorted)
        n = sort_and_compress_rel_indices(rel->primes, n);
    rel->nb = n;
    if (parsesm) {
        /* We need some more stuff. */
        for(rel->sm_size = 0 ; c != '\n' ; rel->sm_size++) {
            if (rel->sm_size >= rel->sm_alloc) {
                checked_realloc(rel->sm, rel->sm_alloc + 1);
                mpz_init(rel->sm[rel->sm_size]);
                rel->sm_alloc++;
            }
            c = earlyparser_inner_read_mpz(r, &p, rel->sm[rel->sm_size]);
        }
    }

    return 1;
}

static int
earlyparser_index(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa(rel, r, 0, 0, 0);
}

/* merge wants sorted indices */
static int
earlyparser_index_sorted(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa(rel, r, 0, 0, 1);
}

static int
earlyparser_indexline(earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    earlyparser_index_maybeabhexa(rel, r, 0, 0, 0);
    r->rhead = p; // rewind ringbuf
    return earlyparser_line(rel, r);
}

static int
earlyparser_abindex_hexa (earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa(rel, r, 1, 0, 0);
}

static int
earlyparser_abindex_hexa_sm (earlyparsed_relation_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa(rel, r, 1, 1, 0);
}


/*}}}*/

/************************************************************************/

/*{{{ filter_rels producer thread */

static void * filter_rels_producer_thread(
    ringbuf_ptr r,
    std::vector<std::string> const & input_files,
    timingstats_dict_ptr stats)
{
    for(auto const & filename : input_files) {
        int status;
        /* We expect all the "filenames" to have been returned by
         * prepare_grouped_command_lines, thus in fact be commands to be
         * passed through popen()
         */
        FILE * f = cado_popen(filename.c_str(), "r");
        ssize_t const rc = ringbuf_feed_stream(r, f);
        int const saved_errno = errno;
#ifdef  HAVE_GETRUSAGE
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
        struct rusage rus;
        status = cado_pclose2(f, &rus);
        if (stats) timingstats_dict_add(stats, "feed-in", &rus);
#else
        status = cado_pclose(f);
#endif
        if (rc < 0) errno = saved_errno;

        if (rc < 0 || status == -1
#if defined(WIFEXITED) && defined(WEXITSTATUS)
            || !WIFEXITED(status) || WEXITSTATUS(status) != 0
#endif
           ) {
            fprintf(stderr,
                    "%s: load error (%s) while %s\n%s\n",
                    __func__,
                    strerror(errno),
                    rc < 0 ? "reading from" : "closing",
                    filename.c_str());
            abort();
        }
    }
    ringbuf_mark_done(r);
    if (stats) timingstats_dict_add_mythread(stats, "producer");
    /*
    double thread_times[2];
    thread_seconds_user_sys(thread_times);
    fprintf(stderr, "Producer thread ends after having spent %.2fs+%.2fs on cpu\n",
            thread_times[0],
            thread_times[1]);
            */
    return nullptr;
}
/*}}}*/

/*{{{ filter_rels consumer thread */
template<typename inflight_t>
static void filter_rels_consumer_thread(
    void *(*callback_fct) (void *, earlyparsed_relation_ptr),
    void * callback_arg,
    int k,
    inflight_t * inflight,
    timingstats_dict_ptr stats)
{
    inflight->enter(k);
    earlyparsed_relation_ptr slot;
    for( ; (slot = inflight->schedule(k)) != nullptr ; ) {
        (*callback_fct)(callback_arg, slot);
        inflight->complete(k, slot);
    }

    inflight->leave(k);
    /*
       double thread_times[2];
       thread_seconds_user_sys(thread_times);
       fprintf(stderr, "Consumer thread (level %d) ends after having spent %.2fs+%.2fs on cpu\n", k, thread_times[0], thread_times[1]);
       */
    if (stats) timingstats_dict_add_mythread(stats, "consumer");
}

/*}}}*/


/*
 * the earlyparse_needed_data is a bitwise OR of the EARLYPARSE_NEED_*
 * constants defined in filter_io.h ; this bitmask defines which function
 * is used here in the consumer thread to achieve the parsing of the ring
 * buffer data, prior to shipping the data to the callback thread. These
 * fields end up in the earlyparse_relation structure (which, by
 * definition, is always somewhat incomplete). These structures form the
 * meat of the inflight buffer, and their allocation is made according to
 * what earlyparse_needed_data asks for.
 *
 * if the bit_vector_srcptr argument [active] is non-NULL, then only
 * relations marked with 1 in this bit vector are processed, while the
 * others are skipped.
 *
 * returns the number of relations filtered
 */

/* see non-templated filter_rels2 below to see how this template is
 * instantiated */
template<typename inflight_t>
static uint64_t filter_rels2_inner(std::vector<std::string> const & input_files,
        filter_rels_description * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    stats_data_t infostats;  /* for displaying progress */
    uint64_t nrels = 0, nactive = 0;
    size_t nB = 0;

    /* {{{ setup and start the producer thread (for the first pipe) */
    auto commands = prepare_grouped_command_lines(input_files);

    ringbuf rb;
    ringbuf_init(rb, 0);

    std::thread producer(filter_rels_producer_thread, rb, commands, stats);

    /* }}} */

    /* {{{ configure the (limited) parsing we will do on the relations
     * read from the files. */
    int (*earlyparser)(earlyparsed_relation_ptr rel, ringbuf_ptr r);
#define _(X) EARLYPARSE_NEED_ ## X      /* convenience */
    switch(earlyparse_needed_data) {
        case _(AB_DECIMAL)|_(LINE):
            /* e.g. for dup1 */
            earlyparser = earlyparser_abline_decimal;
            break;
        case _(AB_HEXA)|_(LINE):
            /* e.g. for dup1 -- ffs reaches here via a command-line flag
             * -abhexa. */
            earlyparser = earlyparser_abline_hexa;
            break;
        case _(AB_HEXA):
            /* e.g. for dup2 (for renumbered files) */
            earlyparser = earlyparser_ab;       /* in hex ! */
            break;
        
        /* dup2/pass2 decides between the two settings here by
         * differenciation of the binaries (dup2-ffs versus dup2)
         */
        case _(AB_DECIMAL)|_(PRIMES):
            /* dup2/pass2 */
            earlyparser = earlyparser_abp_decimal;
            break;
        case _(AB_HEXA)|_(PRIMES):
            /* dup2/pass2 */
            earlyparser = earlyparser_abp_hexa;
            break;

        case _(INDEX)|_(SORTED):
            earlyparser = earlyparser_index_sorted;
            break;
        case _(INDEX):
            /* all binaries after dup2 which do not need a,b*/
            earlyparser = earlyparser_index;
            break;
        case _(INDEX) | _(AB_HEXA):
            /* e.g. reconstructlog */
            earlyparser = earlyparser_abindex_hexa;
            break;
        case _(INDEX) | _(AB_HEXA) | _(SM):
            /* e.g. reconstructlog */
            earlyparser = earlyparser_abindex_hexa_sm;
            break;
        case _(LINE):
            /* e.g. for purge/2 */
            earlyparser = earlyparser_line;
            break;
        case _(LINE) | _(INDEX):
            earlyparser = earlyparser_indexline;
            break;
        default:
            fprintf(stderr, "Unexpected bitmask in %s, please fix\n", __func__);
            /* Fixing is not hard. Just needs writing another parsing
             * function above */
            abort();
    }
#undef _
    /* }}} */

    int n;      /* number of levels of the pipe */
    int ncons = 0;      /* total number of consumers (levels >=1) */
    for(n = 0 ; desc[n].f ; n++) ncons += desc[n].n;
    n++;        /* match with the "number of levels" counter
                   in the inflight_rels_buffer structure definition. */
    ASSERT_ALWAYS(n >= 2);

    /* now prepare the inflight buffer, and the appropriate threads */
    inflight_t inflight_obj(ncons + 1);
    inflight_t * inflight = &inflight_obj;

    /* {{{ complement the inflight rels allocation depending on our
     * parsing needs */
#define _(X) EARLYPARSE_NEED_ ## X      /* convenience */
    if (earlyparse_needed_data & (_(PRIMES) | _(INDEX))) {
        for (unsigned int i = 0 ; i < SIZE_BUF_REL; i++) {
            inflight->rels[i]->primes = inflight->rels[i]->primes_data;
            inflight->rels[i]->nb_alloc = NB_PRIMES_OPT;
        }
    }
    if (earlyparse_needed_data & _(LINE)) {
        for (unsigned int i = 0 ; i < SIZE_BUF_REL; i++) {
            inflight->rels[i]->line = (char*) malloc(RELATION_MAX_BYTES);
	    if (!inflight->rels[i]->line)
	      {
		fprintf (stderr, "Cannot allocate memory\n");
		abort ();
	      }
        }
    }
#undef _
    /* }}} */

    /* {{{ setup and start all the consumer threads (at all levels) */
    /* these first few linse are also found in the non-templated function
     * which instantiates and calls us.
     */

    std::vector<std::thread> consumers;
    consumers.reserve(ncons);
    for(int j = 0, k = 1 ; j < ncons ; j+=desc[k-1].n, k++) {
        ASSERT_ALWAYS(k < n);
        for(int i = 0 ; i < desc[k-1].n ; i++) {
            consumers.emplace_back(filter_rels_consumer_thread<inflight_t>, 
                    desc[k-1].f,
                    desc[k-1].arg,
                    k,
                    inflight,
                    stats);
        }
    }

    /* }}} */

    /* {{{ main loop */
    /* will print report at 2^10, 2^11, ... 2^23 read rels and every 2^23 rels
     * after that */
    if (!active)
      stats_init (infostats, stdout, &nrels, 23, "Read", "relations", "",
                  "rels");
    else
      stats_init (infostats, stdout, &nrels, 23, "Read", "active rels",
                  "read rels", "rels");

    inflight->enter(0);
    for(size_t avail_seen = 0 ; ; ) {
        pthread_mutex_lock(rb->mx);
        while(rb->avail_to_read == avail_seen && !rb->done) {
            pthread_cond_wait(rb->bored, rb->mx);
        }
        avail_seen = rb->avail_to_read; /* must be before mutex unlock ! */
        int const done = rb->done;            /* must be before mutex unlock ! */
        pthread_mutex_unlock(rb->mx);
        if (avail_seen == 0 && done) {
            /* end of producer1 is with rb->done = 1 -- which is
             * compatible with bytes still being in the pipe ! */
            break;
        }
        /* We may have one or several lines which have just been
         * produced. As long as we succeed reading complete lines, we
         * consume them, and feed the second pipe.
         */
        int nl;
        for(size_t avail_offset = 0; avail_offset < avail_seen && (nl = ringbuf_strchr(rb, '\n', 0)) > 0 ; ) {
            if (*rb->rhead != '#') {
                uint64_t const relnum = nrels++;
                if (!active || bit_vector_getbit(active, relnum)) {
                    earlyparsed_relation_ptr slot = inflight->schedule(0);
                    slot->num = relnum;
                    (*earlyparser)(slot, rb);
                    inflight->complete(0, slot);
                    nactive++;
                }
            }
            /* skip the newline byte as well */
            nl++;
            nB += nl;
            ringbuf_skip_get(rb, nl);
            avail_seen -= nl;
            avail_offset += nl;
            if (stats_test_progress(infostats))
            {
              if (!active)
                stats_print_progress (infostats, nrels, 0, nB, 0);
              else
                stats_print_progress (infostats, nactive, nrels, nB, 0);
            }
        }
    }
    inflight->drain();
    inflight->leave(0);
    /*}}}*/

    /* {{{ join all threads */
    for(auto & t : consumers)
        t.join();
    producer.join();
    /*}}}*/

    /* NOTE: the inflight dtor is called automatically */

    if (!active)
      stats_print_progress (infostats, nrels, 0, nB, 1);
    else
      stats_print_progress (infostats, nactive, nrels, nB, 1);

    /* clean producer stuff */
    ringbuf_clear(rb);

    return nactive;
}

uint64_t filter_rels2(char const ** input_files,
        struct filter_rels_description * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    std::vector<std::string> stl_input_files;
    for(char const ** x = input_files ; *x ; x++) {
        stl_input_files.emplace_back(*x);
    }
    return filter_rels2(stl_input_files, desc, earlyparse_needed_data, active, stats);
}



uint64_t filter_rels2(std::vector<std::string> const & input_files,
        filter_rels_description * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    int multi = 0;
    int n;      /* number of levels of the pipe */
    // int ncons = 0;      /* total number of consumers (levels >=1) */

    for(auto const & f : input_files) {
      if (f == "-") {
        fprintf (stderr, "Error, using - to read from standard input does "
                         "not work.\nPlease use named pipes instead.\n");
        abort ();
      }
    }

    for(n = 0 ; desc[n].f ; n++) {
        // ncons += desc[n].n;
        multi += desc[n].n > 1;
    }
    n++;        /* match with the "number of levels" counter
                   in the inflight_rels_buffer structure definition. */
    ASSERT_ALWAYS(n >= 2);
    /* Currently we only have had use for n==2 or n==3 */

    if (n == 2 && !multi && !filter_rels_force_posix_threads) {
        using inflight_t = inflight_rels_buffer<ifb_locking_lightweight, 2>;
        return filter_rels2_inner<inflight_t>(input_files, desc, earlyparse_needed_data, active, stats);
    } else if (n == 2) {
        using inflight_t = inflight_rels_buffer<ifb_locking_posix, 2>;
        return filter_rels2_inner<inflight_t>(input_files, desc, earlyparse_needed_data, active, stats);
    } else if (n == 3) {
        using inflight_t = inflight_rels_buffer<ifb_locking_posix, 3>;
        return filter_rels2_inner<inflight_t>(input_files, desc, earlyparse_needed_data, active, stats);
    } else {
        fprintf(stderr, "filter_rels2 is not explicitly configured (yet) to support deeper pipes\n");
        abort();
    }
}

uint64_t filter_rels(char const ** input_files,
        filter_rels_callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    struct filter_rels_description desc[2] = {
        { f, arg, 1, }, { nullptr, nullptr, 0, },
    };
    std::vector<std::string> stl_input_files;
    for(char const ** x = input_files ; *x ; x++) {
        stl_input_files.emplace_back(*x);
    }
    return filter_rels2(stl_input_files, desc, earlyparse_needed_data, active, stats);
}
