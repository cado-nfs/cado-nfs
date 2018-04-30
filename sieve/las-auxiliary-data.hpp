#ifndef LAS_AUXILIARY_DATA_HPP_
#define LAS_AUXILIARY_DATA_HPP_

#include <unordered_set>
#include "las-threads.hpp"
#include "las-todo-entry.hpp"
#include "las-report-stats.hpp"
#include "utils/tdict.hpp"
#include "utils/timing.h"

typedef std::pair< int64_t, uint64_t> abpair_t;

struct abpair_hash_t {
    inline unsigned long operator()(abpair_t const& o) const {
        return 314159265358979323UL * o.first + 271828182845904523UL + o.second;
    }
};

/* Compute a checksum over the bucket region.

   We import the bucket region into an mpz_t and take it modulo
   checksum_prime. The checksums for different bucket regions are added up,
   modulo checksum_prime. This makes the combined checksum independent of
   the order in which buckets are processed, but it is dependent on size of
   the bucket region. Note that the selection of the sieve region, i.e., of J
   depends somewhat on the number of threads, as we want an equal number of
   bucket regions per thread. Thus the checksums are not necessarily
   clonable between runs with different numbers of threads. */

class sieve_checksum {
  static const unsigned int checksum_prime = 4294967291u; /* < 2^32 */
  unsigned int checksum;
  void update(const unsigned int);

  public:
  sieve_checksum() : checksum(0) {}
  unsigned int get_checksum() {return checksum;}

  /* Combine two checksums */ 
  void update(const sieve_checksum &other) {
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
 * Someday, most of this will be done with shared_ptr's. Both the
 * nfs_aux structure, and the already_printed_for_q member.
 */
class nfs_aux {/*{{{*/
    /* okay, it's hidden. We *only* use it in the dtor, where we decide
     * whether we print some stuff or not. But beyond that, really, it
     * seems wrong to keep a tie to las_info here (well, actually
     * anywhere, to be honest -- but for the core algorithmic data, that'
     * slightly more entangled).
     */
    las_info const & las;
    public:
    las_todo_entry doing;

    /* Will be converted to a shared_ptr. */
    std::unordered_set<abpair_t, abpair_hash_t>& already_printed_for_q;

    /* These two are initialized by the caller, and the caller itself
     * will collate them with the global counters.
     *
     * Unfortunately, this is incompatible with asynchronous processing.
     * We would then need to arrange them in lists (and no other
     * containers), and have the exception handlers splice them to a
     * "botched" list, before final synchronization.
     */
    las_report & rep;
    timetree_t & timer_special_q;
    where_am_I w;

    /* This boolean is set to true when we successfully run all the
     * sieving without any need for reallocation */
    bool complete = false;

#if 0
    /* It _seems_ that it's ok to have sieve_info as a reference, even
     * though asynchronous special-q's might lead us to continue
     * computing with a special-q beyond a point where the sieve_info
     * structure, that is fetched from las.sievers, might have been taken
     * over by the next special q.
     * This is because we're only doing ECM, past this point.
     *
     * This being said, it's probably cleaner to not have this in this
     * structure.
     */
    sieve_info &si;
#endif

    std::array<sieve_checksum,2> checksum_post_sieve;

    /* This is meant to replace at least part
     * nfs_work::thread_data. Maybe not everything.
     */

    struct thread_data {/*{{{*/
        nfs_aux & common;
        /* each thread has its own, and we'll summarize at the end */
        las_report rep;
        timetree_t timer;
        std::array<sieve_checksum,2> checksum_post_sieve;
        where_am_I w;
        thread_data(nfs_aux & t) : common(t) {}
        void update_checksums(nfs_work::thread_data & tws);

    };/*}}}*/

    std::vector<thread_data> th;

    double qt0;

    nfs_aux(las_info const & las, las_todo_entry const & doing, std::unordered_set<abpair_t, abpair_hash_t>& a, las_report & rep, timetree_t & t, int nthreads)
        :
            las(las),   /* shame... */
            doing(doing),
            already_printed_for_q(a),
            rep(rep),
            timer_special_q(t),
            th(nthreads, thread_data(*this))
    {
        qt0 = seconds();
        ASSERT_ALWAYS(!timer_special_q.running());
        timer_special_q.start();
    }


    /* This dtor does stuff! It collates all the auxiliary data for the
     * different threads, and does some final printing.
     */
    ~nfs_aux();
};/*}}}*/


#endif	/* LAS_AUXILIARY_DATA_HPP_ */
