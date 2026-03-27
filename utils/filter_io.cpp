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

#include "arithxx/u64arith.h"          // for constepxr u64arith_clz
#include "barrier.h"                   // for barrier_destroy, barrier_init
#include "bit_vector.h"
#include "cado_popen.h"                // for cado_pclose2, cado_popen
#include "cxx_mpz.hpp"
#include "filter_io.hpp"
#include "gzip.h"                      // prepare_grouped_command_lines
#include "macros.h"                    // for ASSERT_ALWAYS, ASSERT, UNLIKELY
#include "ringbuf.hpp"                   // for ringbuf_s, ringbuf_ptr, RINGBU...
#include "stats.h"                     // stats_data_t
#include "portability.h" // sleep // IWYU pragma: keep
#include "runtime_numeric_cast.hpp"
#include "timing.h"
#include "typedefs.h"
#include "utils_cxx.hpp"
#include "cado_compile_time_hacks.hpp"

/* This is a configuration variable which may be set by the caller (it's
 * possible to bind it to a command-line argument)
 */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
int filter_rels_force_posix_threads = 0;

/************************************************************************/

/* {{{ early parsing routines, for filling the earlyparsed_relation
 * structures. */

/* malloc()'s are avoided as long as there are less than NB_PRIMES_OPT in
 * the relation
 */
template<cado::filter_io_details::filter_io_config cfg>
void realloc_buffer_primes(typename cfg::rel_ptr buf)
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

void realloc_buffer_primes_c(earlyparsed_relation_ptr buf)
{
    realloc_buffer_primes<filter_io_default_cfg>(buf);
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
template<uint64_t base>
static inline int
earlyparser_inner_read_ab(
        ringbuf_ptr r,
        const char ** pp,
        earlyparsed_relation_ptr rel)
{
    static_assert(1u <= base && base <= 16u, "base should be in [1, 16]");
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

template<uint64_t base>
static inline int
earlyparser_inner_read_ab(
        ringbuf_ptr r,
        const char ** pp,
        earlyparsed_relation_mpz_ptr rel)
{
    using cado_math_aux::log2_ct;
    static_assert(1u <= base && base <= 16u, "base should be in [1, 16]");
    const char * p = *pp;
    int c;
    unsigned long v;
    RINGBUF_GET_ONE_BYTE(c, r, p);
    int negative = 0;
    if (c == '-') {
        negative = 1;
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    mpz_set_ui(rel->a, 0);
    for (; (v = ugly[c]) < base;) {
        if constexpr ((base & (base - 1)) == 0) {
            mpz_mul_2exp(rel->a, rel->a, log2_ct(base));
        } else {
            mpz_mul_ui(rel->a, rel->a, base);
        }
        mpz_add_ui(rel->a, rel->a, v);
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    PARSER_ASSERT_ALWAYS(c, ',', *pp, p);
    if (negative)
        mpz_neg(rel->a, rel->a);
    RINGBUF_GET_ONE_BYTE(c, r, p);
    mpz_set_ui(rel->b, 0);
    for (; (v = ugly[c]) < base;) {
        if constexpr ((base & (base - 1)) == 0) {
            mpz_mul_2exp(rel->b, rel->b, log2_ct(base));
        } else {
            mpz_mul_ui(rel->b, rel->b, base);
        }
        mpz_add_ui(rel->b, rel->b, v);
        RINGBUF_GET_ONE_BYTE(c, r, p);
    }
    *pp = p;
    return c;
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

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_inner_read_active_sides(int c, ringbuf_ptr r, const char ** pp, typename cfg::rel_ptr rel)
{
    const char * p = *pp;
    if (c == '@') {
        uint64_t v,w;
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
    return c;
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

template<cado::filter_io_details::filter_io_config cfg, uint64_t base>
static inline int earlyparser_abp(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;

    int c = earlyparser_inner_read_ab<base>(r, &p, rel);
    c = earlyparser_inner_read_active_sides<cfg>(c, r, &p, rel);

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
            if (rel->nb_alloc == n) realloc_buffer_primes<cfg>(rel);
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


template<cado::filter_io_details::filter_io_config cfg, uint64_t base>
static int
earlyparser_ab(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    int c = earlyparser_inner_read_ab<base>(r, &p, rel);
    earlyparser_inner_read_active_sides<cfg>(c, r, &p, rel);

    return 1;
}

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_line(typename cfg::rel_ptr rel, ringbuf_ptr r)
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

/* e.g. for dup1 */
template<cado::filter_io_details::filter_io_config cfg, uint64_t base>
static int
earlyparser_abline(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    int c = earlyparser_inner_read_ab<base>(r, &p, rel);
    earlyparser_inner_read_active_sides<cfg>(c, r, &p, rel);
    return earlyparser_line<cfg>(rel, r);
}

/* Note: for these routine, the sorting of the primes is not considered, and at
 * least for the 1st pass of purge, so far we've been using this code on
 * unsorted relations.
 */
template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_index_maybeabhexa(
        typename cfg::rel_ptr rel,
        ringbuf_ptr r,
        int parseab,
        int parsesm,
        int sort)
{
    const char *p = r->rhead;

    if (parseab) {
        earlyparser_inner_read_ab<16u>(r, &p, rel);
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
            if (rel->nb_alloc == n) realloc_buffer_primes<cfg>(rel);
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

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_index(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa<cfg>(rel, r, 0, 0, 0);
}

/* merge wants sorted indices */
template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_index_sorted(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa<cfg>(rel, r, 0, 0, 1);
}

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_indexline(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    const char * p = r->rhead;
    earlyparser_index_maybeabhexa<cfg>(rel, r, 0, 0, 0);
    r->rhead = p; // rewind ringbuf
    return earlyparser_line<cfg>(rel, r);
}

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_abindex_hexa(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa<cfg>(rel, r, 1, 0, 0);
}

template<cado::filter_io_details::filter_io_config cfg>
static int
earlyparser_abindex_hexa_sm(typename cfg::rel_ptr rel, ringbuf_ptr r)
{
    return earlyparser_index_maybeabhexa<cfg>(rel, r, 1, 1, 0);
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
    void *(*callback_fct) (void *, typename inflight_t::cfg::rel_ptr),
    void * callback_arg,
    int k,
    inflight_t * inflight,
    timingstats_dict_ptr stats)
{
    inflight->enter(k);
    typename inflight_t::cfg::rel_ptr slot;
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
 * constants defined in filter_io.hpp ; this bitmask defines which function
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
static uint64_t filter_rels2_inner(
        std::vector<std::string> const & input_files,
        typename inflight_t::cfg::description_t * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    using cfg = inflight_t::cfg;
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
    int (*earlyparser)(typename cfg::rel_ptr rel, ringbuf_ptr r);
#define _(X) EARLYPARSE_NEED_ ## X      /* convenience */
    switch(earlyparse_needed_data) {
        case _(AB_DECIMAL)|_(LINE):
            /* e.g. for dup1 */
            earlyparser = earlyparser_abline<cfg, 10u>;
            break;
        case _(AB_HEXA)|_(LINE):
            /* e.g. for dup1 -- ffs reaches here via a command-line flag
             * -abhexa. */
            earlyparser = earlyparser_abline<cfg, 16u>;
            break;
        case _(AB_HEXA):
            /* e.g. for dup2 (for renumbered files) */
            earlyparser = earlyparser_ab<cfg, 16u>;
            break;
        case _(AB_DECIMAL)|_(PRIMES):
            /* dup2/pass2 */
            earlyparser = earlyparser_abp<cfg, 10u>;
            break;
        case _(AB_HEXA)|_(PRIMES):
            /* dup2/pass2 */
            earlyparser = earlyparser_abp<cfg, 16u>;
            break;

        case _(INDEX)|_(SORTED):
            earlyparser = earlyparser_index_sorted<cfg>;
            break;
        case _(INDEX):
            /* all binaries after dup2 which do not need a,b*/
            earlyparser = earlyparser_index<cfg>;
            break;
        case _(INDEX) | _(AB_HEXA):
            /* e.g. reconstructlog */
            earlyparser = earlyparser_abindex_hexa<cfg>;
            break;
        case _(INDEX) | _(AB_HEXA) | _(SM):
            /* e.g. reconstructlog */
            earlyparser = earlyparser_abindex_hexa_sm<cfg>;
            break;
        case _(LINE):
            /* e.g. for purge/2 */
            earlyparser = earlyparser_line<cfg>;
            break;
        case _(LINE) | _(INDEX):
            earlyparser = earlyparser_indexline<cfg>;
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
                    typename inflight_t::cfg::rel_ptr slot = inflight->schedule(0);
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
    return filter_rels2<filter_io_default_cfg>(stl_input_files, desc, earlyparse_needed_data, active, stats);
}


template<cado::filter_io_details::filter_io_config cfg>
uint64_t filter_rels2(std::vector<std::string> const & input_files,
        typename cfg::description_t * desc,
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
        using inflight_t = cado::filter_io_details::inflight_rels_buffer<cado::filter_io_details::ifb_locking_lightweight, 2, cfg>;
        return filter_rels2_inner<inflight_t>(input_files, desc, earlyparse_needed_data, active, stats);
    } else if (n == 2) {
        using inflight_t = cado::filter_io_details::inflight_rels_buffer<cado::filter_io_details::ifb_locking_posix, 2, cfg>;
        return filter_rels2_inner<inflight_t>(input_files, desc, earlyparse_needed_data, active, stats);
    } else if (n == 3) {
        using inflight_t = cado::filter_io_details::inflight_rels_buffer<cado::filter_io_details::ifb_locking_posix, 3, cfg>;
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
    return filter_rels2<filter_io_default_cfg>(stl_input_files, desc, earlyparse_needed_data, active, stats);
}

uint64_t filter_rels_mpz(
        char const ** input_files,
        filter_rels_mpz_callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    struct filter_rels_mpz_description desc[2] = {
        { f, arg, 1, }, { nullptr, nullptr, 0, },
    };
    std::vector<std::string> stl_input_files;
    for(char const ** x = input_files ; *x ; x++) {
        stl_input_files.emplace_back(*x);
    }
    return filter_rels2<filter_io_large_ab_cfg>(stl_input_files, desc, earlyparse_needed_data, active, stats);
}
