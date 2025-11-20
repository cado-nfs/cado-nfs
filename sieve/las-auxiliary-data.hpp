#ifndef CADO_LAS_AUXILIARY_DATA_HPP
#define CADO_LAS_AUXILIARY_DATA_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <utility>
#include <vector>

#include "las-info.hpp"
#include "las-report-stats.hpp"
#include "las-where-am-i-proxy.hpp"
#include "lock_guarded_container.hpp"
#include "threadpool.hpp"
#include "tdict.hpp"
#include "timing.h"

struct las_output; // IWYU pragma: keep
struct special_q; // IWYU pragma: keep


/* Compute a checksum over the bucket region.
 *
 * We import the bucket region into an mpz_t and take it modulo
 * checksum_prime. The checksums for different bucket regions are added up,
 * modulo checksum_prime. This makes the combined checksum independent of
 * the order in which buckets are processed, but it is dependent on size of
 * the bucket region. Note that the selection of the sieve region, i.e., of J
 * depends somewhat on the number of threads, as we want an equal number of
 * bucket regions per thread. Thus the checksums are not necessarily
 * clonable between runs with different numbers of threads.
 */

struct report_and_timer {
    las_report rep;
    timetree_t timer;
    std::mutex mm;
};

class sieve_checksum {
  static const unsigned int checksum_prime = 4294967291U; /* < 2^32 */
  unsigned int checksum = 0;
  void update(unsigned int);

  public:
  unsigned int get_checksum() const { return checksum; }

  /* Combine two checksums */ 
  void update(sieve_checksum const & other) {
    /* Simply (checksum+checksum2) % checksum_prime, but using
       ularith_addmod_ul_ul() to handle sums >= 2^32 correctly. */
    this->update(other.checksum);
  }
  /* Update checksum with the pointed-to data */
  void update(const unsigned char *, size_t);
};

/* This structure is here to gather mostly timing and bookkeeping
 * information for all threads that work on a given special-q, and more
 * precisely one one "attempt", meaning one that does not change the
 * bkmult. If bkmult changes, we get a new structure, because we want to
 * collect timings differently.
 * There's one important exception: the already_printed_for_q data member
 * is a reference to an object that is persistent across all attempts for
 * a given q.
 *
 * This structure lives through a shared_ptr, and so does the
 * already_printed_for_q member we're pointing to. This makes it possible
 * for objects to persist in memory while another special-q is being
 * processed.
 */
class nfs_aux {/*{{{*/
    /* okay, it's hidden. We *only* use it in the dtor, where we decide
     * whether we print some stuff or not. But beyond that, really, it
     * seems wrong to keep a tie to las_info here (well, actually
     * anywhere, to be honest -- but for the core algorithmic data,
     * that's slightly more entangled).
     */
    las_info const & las;
    public:

    /* This lives somewhere in special_q_task_collection. We don't have
     * ownership, as the special_q_task_collection is persistent anyway.
     * Note that &doing might be dynamic_cast-able to
     * special_q_task_simple or special_q_task_tree, and we do make use
     * of this.
     */
    special_q_task & doing;
    
    /* we rarely have ownership, if ever, of course. In the typical case,
     * there just one output file and that's it.
     * However in client-server file we may want several output files. In
     * that case, the caller that reads the todo list will create a
     * dedicated output file, and its dtor will do the final work,
     * whenever it gets called.
     */
    std::shared_ptr<las_output> output_p;

    using abpair_t = std::pair< int64_t, uint64_t>;

    struct abpair_hash_t {
        unsigned long operator()(abpair_t const& o) const {
            return 314159265358979323UL * o.first + 271828182845904523UL + o.second;
        }
    };

    using rel_hash_t = lock_guarded_container<std::unordered_set<abpair_t, abpair_hash_t>>;

    std::shared_ptr<rel_hash_t> rel_hash_p;
    rel_hash_t & get_rel_hash() { return * rel_hash_p ; }

    where_am_I w;

    report_and_timer rt;

    /* The creator scope must fill these fields by hand. The final report
     * will be collated into these two. A reasonable way to go is to
     * first set a default destination that is essentially a trash can,
     * and modify it to something that is a better-defined destination
     * once we're sure we won't have exceptions.
     *
     * The boolean value "complete" is typically used to distinguish
     * between the case where this destination is a trash can versus when
     * we successfully run all the sieving without any need for
     *
     * See las_subjob()
     */
    report_and_timer * dest_rt = nullptr;
    bool complete = false;

    /* This gets completed somewhat late */
    std::vector<sieve_checksum> checksum_post_sieve;

    struct thread_data {/*{{{*/
        //nfs_aux & common;
        /* each thread has its own, and we'll summarize at the end */
        las_report rep;
        timetree_t timer;
        std::vector<sieve_checksum> checksum_post_sieve;
        where_am_I w;
        //thread_data(nfs_aux & t) : common(t) {}
        void update_checksums(int side, const unsigned char *data, const size_t len) {
            if (!data)
                return;
            /* It's simpler to auto-vivify */
            for( ; checksum_post_sieve.size() <= (size_t) side ; )
                checksum_post_sieve.emplace_back();
            checksum_post_sieve[side].update(data, len);
        }
    };/*}}}*/

    std::vector<thread_data> th;

    timetree_t & get_timer(worker_thread * worker) {
        return worker->is_synchronous() ? rt.timer : th[worker->rank()].timer;
    }

    double qt0;
    double wct_qt0;

    nfs_aux(las_info const & las,
            special_q_task & doing,
            std::shared_ptr<rel_hash_t> & rel_hash_p,
            int nthreads)
        : las(las)   /* shame... */
        , doing(doing)
        , rel_hash_p(rel_hash_p)
        , checksum_post_sieve(las.cpoly->nb_polys)
        , th(nthreads)
          //, thread_data(*this))
        , qt0(seconds())
        , wct_qt0(wct_seconds())
    {
    }


    /* This dtor does stuff! It collates all the auxiliary data for the
     * different threads, and does some final printing.
     */
    ~nfs_aux();
};/*}}}*/

#ifndef DISABLE_TIMINGS

extern tdict::slot_parametric tdict_slot_for_side;
extern tdict::slot tdict_slot_for_alloc_buckets;
extern tdict::slot tdict_slot_for_threads;
extern tdict::slot tdict_slot_for_fibt;

#define ENTER_THREAD_TIMER(timer)       \
    ACTIVATE_TIMER_IF_NOT_RUNNING(timer);                           \
    const typename std::remove_reference<decltype(timer)>::type::accounting_child UNIQUE_ID(dummy)(timer, tdict_slot_for_threads)

#define ENTER_THREAD_FUZZY_TIMER(T, U)       \
    ACTIVATE_TIMER_IF_NOT_RUNNING(T);                           \
    const tdict::tie_timer<typename TIMER_TYPE_(T)::timer_type, fast_timetree_t::timer_type> U(T, tdict_slot_for_threads)

#define MARK_TIMER_FOR_SIDE(timer, side)       \
    const typename std::remove_reference<decltype(timer)>::type::accounting_child UNIQUE_ID(dummy)(timer, tdict_slot_for_side(side))

#else /* DISABLE_TIMINGS */

#define ENTER_THREAD_TIMER(timer) timer.nop() /**/
#define MARK_TIMER_FOR_SIDE(timer, side) timer.nop() /**/
#define ENTER_THREAD_FUZZY_TIMER(T, U) T.nop(); tdict::tie_timer U

#endif  /* DISABLE_TIMINGS */

#endif	/* CADO_LAS_AUXILIARY_DATA_HPP */
