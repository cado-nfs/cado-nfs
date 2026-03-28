#ifndef CADO_FILTER_IO_HPP
#define CADO_FILTER_IO_HPP

#include <cstdint>
#include <ctime>

#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <barrier>

#include "cado_type_traits.hpp"

#include <gmp.h>        // for mpz_t
                        //
#include "bit_vector.h" // for bit_vector_srcptr
#include "timing.h"     // for timingstats_dict_ptr
#include "typedefs.h"   // for prime_t, weight_t
#include "macros.h"

#define MAX_FILES 1000000

#define RELATION_MAX_BYTES 4096

/* Size of relations buffer between parsing & processing.
 * CAREFUL! SIZE_BUF_REL must be greater (at least double) than (1<<(NNFR+1)).
 * Stores the sentences precomputed but not inserted.
 * About 64K sentences for the optimal.
 */
// #define SIZE_BUF_REL MAX((1<<(NFR+NNFR+1+(NFR==0))),(1<<6))
#define SIZE_BUF_REL (1U << 15U)

/* we want the programs to specify in a completely explicit way whether
 * they want stuff in base 10 or 16 */
#define EARLYPARSE_NEED_AB_DECIMAL 1U
#define EARLYPARSE_NEED_AB_HEXA 2U
#define EARLYPARSE_NEED_LINE 4U
/* for reading ideals (e.g. output of las) */
#define EARLYPARSE_NEED_PRIMES 8U
/* for reading index (i.e. renumber ideal) */
#define EARLYPARSE_NEED_INDEX 16U
#define EARLYPARSE_NEED_SM 32U
#define EARLYPARSE_NEED_SORTED 64U
#define EARLYPARSE_NEED_INDEX_SORTED                                           \
    (EARLYPARSE_NEED_INDEX | EARLYPARSE_NEED_SORTED)

/* Initial size of primes_data array in earlyparsed_relation_s,
   If more than NB_PRIMES_OPT is needed (should be rare), *primes is
   allocated
*/
#define NB_PRIMES_OPT 31

/* This field does not contain a proper relation, but only something
 * which has undergone quite limited parsing within a thread whose job is
 * to read data fast, and not to bother about the fine points of the
 * relation semantics. Because of this, there is no unique definition of
 * which are the valid fields below. This depends on how this field is
 * meant to be used downstream. Depending on the earlyparse_needed_data
 * bitmask argument fed to filter_relations, we may fill only some fields.
 * Which fields are filled is controlled by which of the
 * EARLYPARSE_NEED_* appear in the earlyparse_needed_data bitmask. The
 * callback thread function given to process_rels is then called for each
 * such "relation"
 */
struct earlyparsed_relation_s {
    int64_t a;
    uint64_t b;
    int active_sides[2];
    prime_t * primes; /*if nb_alloc <= NB_PRIME_OPT, primes == primes_data*/
    prime_t primes_data[NB_PRIMES_OPT];
    weight_t nb;       /* number of primes */
    weight_t nb_alloc; /* allocated space for primes
                        * (if > NB_PRIMES_OPT: indirect addressing,
                        * otherwise primes == primes_data) */
    /* nb_above_min_index is counted only when ->primes is needed anyway,
     * so we defer it to the callback function instead.
     */
    // weight_t nb_above_min_index; /* nb of primes above min_index, must be
    // <=nb */
    uint64_t num; /* (absolute) relation number */
    char * line; /* If not NULL, contains the relation with a '\n' at the end */
    mpz_t * sm;
    int sm_size;
    int sm_alloc;
};
typedef struct earlyparsed_relation_s earlyparsed_relation[1];
typedef struct earlyparsed_relation_s * earlyparsed_relation_ptr;
typedef const struct earlyparsed_relation_s * earlyparsed_relation_srcptr;

struct earlyparsed_relation_mpz_s {
    mpz_t a;
    mpz_t b;
    int active_sides[2];
    prime_t * primes; /*if nb_alloc <= NB_PRIME_OPT, primes == primes_data*/
    prime_t primes_data[NB_PRIMES_OPT];
    weight_t nb;       /* number of primes */
    weight_t nb_alloc; /* allocated space for primes
                        * (if > NB_PRIMES_OPT: indirect addressing,
                        * otherwise primes == primes_data) */
    /* nb_above_min_index is counted only when ->primes is needed anyway,
     * so we defer it to the callback function instead.
     */
    // weight_t nb_above_min_index; /* nb of primes above min_index, must be
    // <=nb */
    uint64_t num; /* (absolute) relation number */
    char * line; /* If not NULL, contains the relation with a '\n' at the end */
    mpz_t * sm;
    int sm_size;
    int sm_alloc;
};
typedef struct earlyparsed_relation_mpz_s earlyparsed_relation_mpz[1];
typedef struct earlyparsed_relation_mpz_s * earlyparsed_relation_mpz_ptr;
typedef const struct earlyparsed_relation_mpz_s * earlyparsed_relation_mpz_srcptr;

static unsigned char const ugly[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   255, 255,
    255, 255, 255, 255, 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 10,  11,  12,  13,  14,  15,  255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255};

void realloc_buffer_primes_c(earlyparsed_relation_ptr buf);

extern int filter_rels_force_posix_threads;
/*
 * A pointer to such a structure must be provided to filter_rels, and
 * termination is indicated by f==NULL. The function specified in the
 * k-th member (starting from 1) in this array describes the operation
 * performed at level k in the process, level 0 being the implicit
 * production of relation from the input files.  Alongside with the
 * relation on which the thread is allowed to work, all threads at level
 * k receive the void* argument specified in the array member. n
 * specifies the number of worker threads to be used for each step.
 */

struct filter_rels_description {
    void * (*f)(void *, earlyparsed_relation_ptr);
    void * arg;
    int n;
};

struct filter_rels_mpz_description {
    void * (*f)(void *, earlyparsed_relation_mpz_ptr);
    void * arg;
    int n;
};

extern uint64_t filter_rels2(char const ** input_files,
                             struct filter_rels_description * desc,
                             int earlyparse_needed_data,
                             bit_vector_srcptr active, timingstats_dict_ptr);

typedef void * (*filter_rels_callback_t)(void *, earlyparsed_relation_ptr);

extern uint64_t filter_rels(char const ** input_files, filter_rels_callback_t f,
                            void * arg, int earlyparse_needed_data,
                            bit_vector_srcptr active,
                            timingstats_dict_ptr stats);

typedef void * (*filter_rels_mpz_callback_t)(void *, earlyparsed_relation_mpz_ptr);
extern uint64_t filter_rels_mpz(
        char const ** input_files,
        filter_rels_mpz_callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats);

namespace cado::filter_io_details {

template<typename T>
concept filter_io_config = std::is_empty_v<T>
    && requires { typename T::rel_t; }
    && requires { typename T::rel_ptr; }
    && requires { typename T::rel_srcptr; }
    && requires { typename T::description_t; }
    && requires { typename T::callback_t; }
    && requires { T::ab_init; }
    && requires { T::ab_clear; };

    /* {{{ concept is_locking_layer_v */
    template<typename T>
        concept is_locking_layer_v = std::is_empty_v<T>
        && requires { typename T::template critical_datatype<int>; }
    && requires { typename T::lock_t; }
    && requires { typename T::cond_t; }
    && requires { T::lock; }
    && requires { T::unlock; }
    && requires { T::wait; }
    && requires { T::signal; }
    && requires { T::signal_broadcast; }
    ; /* }}} */

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
    static_assert(cado::filter_io_details::is_locking_layer_v<ifb_locking_posix>);
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
    static void wait(cond_t *, lock_t *) {
        /* {{{ define NANOSLEEP */
        /* The realistic minimal non-CPU waiting with nanosleep is about
         * 10 to 40 microseconds (1<<13 for nanosleep).  But all the I/O
         * between the threads have been buffered, and a thread does a
         * nanosleep only if its buffer is empty.  So I use here ~2ms
         * (1<<21) to optimize CPU scheduler.  Max pause is about 4 to
         * 8ms (1<<22, 1<<23); above that, the program is slowed down.
         */
#ifndef HAVE_NANOSLEEP
#ifdef HAVE_USLEEP
        usleep((unsigned long) (1<<21 / 1000UL));
#else
        sleep(0);
#endif
#else
        struct timespec wait_classical = { .tv_sec=0, .tv_nsec=1<<21 }; // …
        nanosleep(&wait_classical, nullptr); // W: no header providing "nan…
#endif
        /*}}}*/
    }
    static void signal(cond_t *) {}
    static void signal_broadcast(cond_t *) {}
    static int isposix() { return 0; }
};
static_assert(cado::filter_io_details::is_locking_layer_v<ifb_locking_lightweight>);
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

/*{{{ inflight buffer. See filter_io.tex for documentation. */

/* {{{ inflight_rels_buffer: n-level buffer, with underyling locking
 * mechanism specified by the template class.  */
template<typename locking, int n, cado::filter_io_details::filter_io_config config>
struct inflight_rels_buffer {
    using locking_layer = locking;
    std::barrier<cado::nop_function> sync_point;
    using cfg = config;
    using relation_type = typename cfg::rel_t;
    using csize_t = locking_layer::template critical_datatype<size_t>::t;
    using lock_t = locking_layer::lock_t;
    using cond_t = locking_layer::cond_t;

    std::unique_ptr<typename cfg::rel_t[]> rels;        /* always malloc()-ed to SIZE_BUF_REL,
                                           which is a power of two */
    /* invariant:
     * scheduled_0 >= ... >= completed_{n-1} >= scheduled_0 - SIZE_BUF_REL
     */
    csize_t completed[n];
    csize_t scheduled[n];
    status_table<locking_layer> status;
    lock_t m[n];
    cond_t bored[n];
    int active[n] = { 0, } ;     /* number of active threads */

    explicit inflight_rels_buffer(int nthreads_total)
        : sync_point(nthreads_total)
        , rels(std::make_unique<relation_type[]>(SIZE_BUF_REL))
    {
        for(int i = 0 ; i < n ; i++) {
            completed[i] = 0;
            scheduled[i] = 0;
        }
    }

    ~inflight_rels_buffer();

    inflight_rels_buffer(inflight_rels_buffer const&) = delete;
    inflight_rels_buffer(inflight_rels_buffer &&) = delete;
    inflight_rels_buffer& operator=(inflight_rels_buffer const&) = delete;
    inflight_rels_buffer& operator=(inflight_rels_buffer &&) = delete;

    void drain()
        /* This belongs to the buffer closing process.  The out condition of this
         * call is that all X(k) for k>0 terminate.  This call (as well as
         * init/clear) must be called on the producer side (step 0) (in a
         * multi-producer context, only one thread is entitled to call this) */
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

    typename cfg::rel_ptr schedule(int k) /* {{{ */
        /* Schedule a new relation slot for processing at level k.
         *
         * This call may block until a relation is processed by level k-1
         * (or, if k==0, until a slot is made available in the relation
         * buffer).
         *
         * The relation is free for use by the current (consumer) thread
         * until it calls inflight_rels_buffer_completed. When the
         * producing stream ends, this function returns NULL. */
    {
        size_t const prev = k ? (k-1) : (n-1);
        // coverity[result_independent_of_operands]
        ASSERT(active[k] <= locking::max_supported_concurrent);
        size_t s;
        size_t const a = k ? 0 : SIZE_BUF_REL;
        /* in 1-thread scenario, scheduled[k] == completed[k] */
        locking::lock(m + prev);
        if constexpr (locking::max_supported_concurrent == 1) {
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
        typename cfg::rel_ptr rel = rels[slot];
        status.update_shouldbealreadyok(s, k-1);
        locking::unlock(m + prev);
        return rel;
    }
    /*}}}*/
    void complete(int k, typename cfg::rel_srcptr rel) /* {{{ */
    {
        // coverity[result_independent_of_operands]
        ASSERT(active[k] <= locking::max_supported_concurrent);
        const int slot = rel - (typename cfg::rel_srcptr) rels.get();

        locking::lock(m + k);

        size_t my_absolute_index;
        if constexpr (locking::max_supported_concurrent == 1) {
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

/*{{{ ::~inflight_rels_buffer() */
/* must be called on producer side */
template<typename locking, int n, cado::filter_io_details::filter_io_config config>
inflight_rels_buffer<locking, n, config>::~inflight_rels_buffer()
{
    for(int i = 0 ; i < n ; i++) {
        ASSERT_ALWAYS_NOTHROW(active[i] == 0);
    }
    for(size_t i = 0 ; i < SIZE_BUF_REL ; i++) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Waddress"
        if constexpr (cfg::ab_clear != nullptr) {
            cfg::ab_clear(rels[i]->a);
            cfg::ab_clear(rels[i]->b);
        }
#pragma GCC diagnostic pop
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

/*}}}*/

/* }}} */

} /* namespace cado::filter_io_details */

struct filter_io_default_cfg
{
    using rel_t = earlyparsed_relation;
    using rel_ptr = earlyparsed_relation_ptr;
    using rel_srcptr = earlyparsed_relation_srcptr;
    using description_t = filter_rels_description;
    using callback_t = filter_rels_callback_t;
    static constexpr void (*ab_init)(uint64_t) = nullptr;
    static constexpr void (*ab_clear)(uint64_t) = nullptr;
};
static_assert(cado::filter_io_details::filter_io_config<filter_io_default_cfg>);

struct filter_io_large_ab_cfg
{
    using rel_t = earlyparsed_relation_mpz;
    using rel_ptr = earlyparsed_relation_mpz_ptr;
    using rel_srcptr = earlyparsed_relation_mpz_srcptr;
    using description_t = filter_rels_mpz_description;
    using callback_t = filter_rels_mpz_callback_t;
    static constexpr void (*ab_init)(mpz_ptr) = &mpz_init;
    static constexpr void (*ab_clear)(mpz_ptr) = &mpz_clear;
};
static_assert(cado::filter_io_details::filter_io_config<filter_io_large_ab_cfg>);

template<cado::filter_io_details::filter_io_config cfg>
void realloc_buffer_primes(typename cfg::rel_ptr buf);

template<cado::filter_io_details::filter_io_config cfg>
extern uint64_t filter_rels2(
        std::vector<std::string> const & input_files,
        typename cfg::description_t * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats);

template<cado::filter_io_details::filter_io_config cfg>
uint64_t
filter_rels(
        std::vector<std::string> const & input_files,
        typename cfg::callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    typename cfg::description_t desc[2] = {
        { f, arg, 1, }, { nullptr, nullptr, 0, },
    };
    return filter_rels2<cfg>(input_files, desc, earlyparse_needed_data, active, stats);
}

#endif /* CADO_FILTER_IO_HPP */
