#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cinttypes>                    // for PRIi64
#include <cstdarg>                      // for va_arg, va_list // IWYU pragma: keep
#include <cstdint>                      // for uint64_t, uint8_t, int64_t
#include <cstdio>                       // for fprintf, size_t, FILE, asprintf
#include <cstdlib>                      // for free
#include <algorithm>                    // for is_sorted, sort
#include <array>                        // for array
#include <initializer_list>             // for initializer_list
#include <memory>                       // for allocator_traits<>::value_type
#include <vector>                       // for vector<>::iterator, vector, swap

#include "las-smallsieve.hpp"           // for resieve_small_bucket_region

#include "macros.h"                     // for ASSERT, ASSERT_ALWAYS, MAYBE_...

#include "bucket-push-update.hpp"       // for bucket_single::push_update
#include "bucket.hpp"                   // for bucket_update_t, bucket_primes_t
#include "fb-types.h"                   // for fbprime_t, FBROOT_FORMAT, FBP...
#include "fb.hpp"                       // for fb_entry_general, fb_factorba...
#include "gcd.h"       // for gcd_ul
#include "las-arith.hpp"                // for invmod_32
#include "las-config.h"                 // for LOG_BUCKET_REGION
#include "las-where-am-i.hpp"           // for where_am_I, WHERE_AM_I_UPDATE
#include "las-forwardtypes.hpp"         // for spos_t, long_spos_t
#include "las-fbroot-qlattice.hpp"      // for fb_root_in_qlattice
#include "las-qlattice.hpp"             // for qlattice...
#include "las-sieve2357.hpp"            // for sieve2357base::prime_t, sieve...
#include "las-smallsieve-glue.hpp"      // for small_sieve, small_sieve::super
#include "las-smallsieve-lowlevel.hpp"  // for SMALLSIEVE_COMMON_DEFS
#include "las-smallsieve-types.hpp"     // for ssp_t, small_sieve_data_t
#include "las-todo-entry.hpp"           // for las_todo_entry
#include "las-where-am-i-proxy.hpp"          // for where_am_I
#include "portability.h"  // asprintf // IWYU pragma: keep
#include "verbose.h"


/* small sieve and resieving */

/* {{{ Some documentation first.
 *
 * Small primes or powers of small primes p^k with projective root.
 * These hit at 
 *   i*v == j*u (mod p^k) 
 * for some u,v in Z, but gcd(v, p^k) > 1.
 * We may assume gcd(u,p)==1, or we divide the entire equation by p.
 * XXX [ET]: we should also assume that v is a prime power, and that u
 * XXX [ET]: is within [0..p^k/v-1[ ; 
 * We store g = gcd(v, p^k), q = p^k / g, and U = u * (v/g)^(-1) (mod q).
 * XXX [ET]: which would then imply g==v, q=p^k/v, and U=u

 * Then we have
 *   i*v == j*u (mod p^k)  <==>  i == (j/g)*U (mod q)
 * with g|j.
 * 
 * In other words, we can sieve this projective prime (power) much like a 
 * normal prime (power) q with root U, except that after sieving a line 
 * we don't advance by one line, but by g lines.
 * The case where g = q^k and thus q = 1 can be sieved more efficiently,
 * of course, since every entry in each g-th line will be hit, so that
 * the sieving should use long word transfers.

 * Just like for normal primes, the ssdpos value points at the first
 * position to sieve relative to the start of the current sieve region.

 * Within a line that starts at index line_start, for array element of 
 * index x, we have x - line_start = i+I/2. 
 * We skip j=0, as it contains only the single possible relation 
 * (i,j) = (1,0). 
 * For j=1*g, we want i=U (mod q), so x - line_start == I/2+U (mod q),
 * so we initialise 
 *   ssdpos = I*g + (I/2 + U) % q
 * to get the first array index in line j=g, 
 * then within a line sieve ssdpos + t*q < I, t in N,
 * and update 
 *   ssdpos = (ssdpos - line_start + U) % q + line_start + g*I 
 * to get the first position to sieve in the next suitable line.
 * }}} */


/* {{{ Some code for information purposes only */

void ssp_simple_t::print(FILE *f) const
{
    fprintf(f, "# p = %" FBPRIME_FORMAT
            ", r = %" FBROOT_FORMAT
            ", logp = %hhu",
            p,
            r,
            logp);
}

void ssp_t::print(FILE *f) const
{
    if (!is_proj()) {
        fprintf(f, "# p = %" FBPRIME_FORMAT
                ", r = %" FBROOT_FORMAT
                ", logp = %hhu",
                get_p(),
                get_r(),
                logp);
        if (is_pow2()) fprintf(f, " (power of 2)");
        if (is_pattern_sieved()) fprintf(f, " (pattern-sieved)");
        if (is_discarded_proj()) fprintf(f, "(discarded because of projective root)");
        if (is_discarded_sublat()) fprintf(f, "(discarded because not compatible with sub lattices)");
    } else {
        fprintf(f, "# q = %" FBPRIME_FORMAT ", g = %" FBROOT_FORMAT ", U = %" FBPRIME_FORMAT ", logp = %hhu",
            get_q(), get_g(), get_U(), logp);
        if (is_pow2()) fprintf(f, " (power of 2)");
        if (is_pattern_sieved()) fprintf(f, " (pattern-sieved)");
        fprintf(f, " (projective root)");
        if (is_discarded_proj()) fprintf(f, "(discarded because of projective root)");
        if (is_discarded_sublat()) fprintf(f, "(discarded because not compatible with sub lattices)");
    }
}

/* Uses va_list argument for use with verbose_output_vfprint() */
int
small_sieve_dump(FILE *f, const char *header, va_list va)
{
    const small_sieve_data_t * p_ssd = va_arg(va, const small_sieve_data_t *);

    fprintf(f, "%s", header);
    for(auto const & x : p_ssd->ssps) {
        x.print(f);
        fprintf(f, "\n");
    }
    for (auto const & x : p_ssd->ssp) {
        x.print(f);
        fprintf(f, "\n");
    }
    return 1;
}

static void small_sieve_print_contents(const char * prefix, small_sieve_data_t const & ssd)
{
    int nnice=ssd.ssps.size();
    int nproj=0;
    int npow2=0;
    int npattern=0;
    int ndiscard=0;
    for(auto const & ssp : ssd.ssp) {
        nproj += ssp.is_proj();
        npow2 += ssp.is_pow2();
        npattern += ssp.is_pattern_sieved();
        ndiscard += ssp.is_discarded();
        nnice += ssp.is_nice();
    }

    verbose_output_start_batch();
    verbose_output_print(0, 3, "# %s: %d nice primes", prefix, nnice);
    /* Primes may be both even and projective... */
    if (npow2) verbose_output_print(0, 3, ", %d powers of 2", npow2);
    if (npattern) verbose_output_print(0, 3, ", %d pattern-sieved", npattern);
    if (nproj) verbose_output_print(0, 3, ", and %d projective primes", nproj);
    verbose_output_print(0, 3, ".");
    if (ndiscard) verbose_output_print(0, 3, " %d discarded.", ndiscard);
    verbose_output_print(0, 3, "\n");
    /* With -v -v -v, dump all the small sieve data */
    verbose_output_vfprint (0, 4, small_sieve_dump, "# Dump of small sieve data:\n", &ssd);
    verbose_output_end_batch();
}


void small_sieve_info(const char * what, int side, small_sieve_data_t const & r)
{
    char * tmp;
    int rc = asprintf(&tmp, "%s(side %d)", what, side);
    ASSERT_ALWAYS(rc >= 0);
    small_sieve_print_contents(tmp, r);
    free(tmp);
}

/* }}} */

/* {{{ Sieve initialization / clearing : first the easy ones */
void small_sieve_clear(small_sieve_data_t & ssd)
{
    ssd.ssps.clear();
    ssd.ssp.clear();
    ssd.ssdpos_many.clear();
    ssd.ssdpos_many_next.clear();
    ssd.offsets.clear();
}

/* }}} */

/* {{{ Sieve initialization: now the real stuff */

typedef sieve2357base::preferred_simd_type preferred_simd_type;
#if GNUC_VERSION_ATLEAST(6,1,0)
/* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69884 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
typedef sieve2357<preferred_simd_type, uint8_t> preferred_sieve2357;
#if GNUC_VERSION_ATLEAST(6,1,0)
#pragma GCC diagnostic pop
#endif

ssp_t::ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, bool proj) /*{{{*/
: ssp_t(_p, _r, _logp) /* First, initialize everything as if proj=false */
{
    if (_p % 2 == 0) {
        set_pow2();
    }
    if (proj) {
        set_proj();
        unsigned int v = r; /* have consistent notations */
        unsigned int g = gcd_ul(p, v);
        fbprime_t q = p / g;
        set_q(q);
        set_g(g);
        if (q == 1) {
            ASSERT(v == 0); /* We know p|v, make sure v=0 */
            set_U(0);
        } else {
            int rc;
            uint32_t U = v / g; /* coprime to q */
            ASSERT(gcd_ul(U, p) == 1);
            rc = invmod_32(&U, q);
            ASSERT_ALWAYS(rc != 0);
            set_U(U);
        }
        if (preferred_sieve2357::can_sieve(this->get_q())) {
            set_pattern_sieved();
        }
    } else {
        if (preferred_sieve2357::can_sieve(this->get_p())) {
            set_pattern_sieved();
        }
    }
}/*}}}*/

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// ssdpos points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION
// It could actually be larger than 32 bits when I > 16.

/* Functor for sorting ssp_t in the order in which sieve2357 expects them */
struct order_ssp_t {
    bool operator()(const ssp_t &ssp1, const ssp_t &ssp2) const {
        fbprime_t q1 = ssp1.is_proj() ? ssp1.get_q() : ssp1.get_p();
        fbprime_t q2 = ssp2.is_proj() ? ssp2.get_q() : ssp2.get_p();
        return sieve2357base::order_lt(q1, q2);
    }
};
static order_ssp_t order_ssp;

void small_sieve_init(small_sieve_data_t & ssd,
                      std::vector<fb_entry_general> const & resieved,
                      std::vector<fb_entry_general> const & rest,
                      int logI,
                      int side,
                      fb_factorbase::key_type const & fbK,
                      qlattice_basis const & Q,
                      double scale)
{
    const unsigned int thresh = fbK.thresholds[0];
    const int verbose = 0;
    where_am_I w MAYBE_UNUSED;

    ssd.fbK = fbK;

    // This zeroes out all vectors, but keeps storage around nevertheless
    small_sieve_clear(ssd);

    ssd.ssps.reserve(resieved.size() + rest.size());
    if (logI < LOG_BUCKET_REGION)
        ssd.offsets.reserve(resieved.size() + rest.size());
    if (logI > LOG_BUCKET_REGION)
        ssd.offsets.reserve(resieved.size() + rest.size());

    // Do a pass on fb and projective primes, to fill in the data
    // while we have any regular primes or projective primes < thresh left
    
    // The processing of bucket region by nb_threads is interleaved.
    // It means that the positions for the small sieve must jump
    // over the (nb_threads - 1) regions after each region.
    // For typical primes, this jump is easily precomputed and goes into
    // the ssp struct.
    
    // If we are doing sublattices modulo m, then we jump virtually m
    // times faster.
    unsigned int sublatm = Q.sublat.m;

    for (auto const & c : { &resieved, &rest }) {
        if (c == &rest)
            ssd.resieve_end_offset = ssd.ssps.size();
        for (auto const & e : *c) {
            /* p=pp^k, the prime or prime power in this entry, and pp is prime */
            const fbprime_t p = e.q, pp = e.p;
            WHERE_AM_I_UPDATE(w, p, p);

            ASSERT_ALWAYS(p <= thresh);
            if (p > thresh) {
                continue;
            }

            for (int nr = 0; nr < e.nr_roots; nr++) {
                const fb_general_root &root = e.roots[nr];
                fb_root_p1 Rab { root.r, root.proj };

                if (!Q.doing.is_coprime_to(pp)) {
                    continue;
                }

                if (sublatm) {
                    // In sublat mode, disable pattern sieving and primes
                    // dividing m. (pp is the prime, here)
                    //
                    // FIXME. ok, they're certainly not "nice", but we should
                    // sieve them nonetheless.
                    if (pp == 3 || (sublatm % pp) == 0) {
                        continue;
                    }
                }

                const unsigned char logp = fb_log_delta (pp, root.exp, root.oldexp, scale);

                WHERE_AM_I_UPDATE(w, r, Rab.to_old_format(p));
                auto Rq = fb_root_in_qlattice(p, Rab, e.invq, Q);
                fbroot_t r_q = Rq.r;
                bool is_proj_in_ij = Rq.is_projective();
                /* If this root is somehow interesting (projective in (a,b) or
                   in (i,j) plane), print a message */
                if (verbose && (Rab.is_projective() || is_proj_in_ij))
                    verbose_output_print(0, 1, "# small_sieve_init: side %d, prime %"
                            FBPRIME_FORMAT " root %s%" FBROOT_FORMAT " (logp %hhu) "
                            " -> %s%" FBROOT_FORMAT "\n", side, p,
                            Rab.is_projective() ? "1/" : "", Rab.r % p, logp,
                            is_proj_in_ij ? "1/" : "", r_q);

                ssp_t new_ssp(p, r_q, logp, is_proj_in_ij);

                if (p != pp)
                    new_ssp.set_pow(pp);

                /* pattern-sieved primes go to ssp */
                if (new_ssp.is_proj()) {
#if 0
                    if (new_ssp.get_g() >= J) {
                        /* some projective primes never hit (number of lines
                         * to skip is >= J). We're tempted to remove them,
                         * but:
                         *  - we lose hits to (+-1,0) this way (the two
                         *    locations are equal up to sign, but we should
                         *    sieve one of them!) -- see bug 21505.
                         *  - the cost of having them in the list is
                         *    ridiculously small anyway.
                         *
                         * this being said, we should probably deal with
                         * projective primes a bit differently, because the
                         * arithmetic constraints we're imposing on the pos
                         * fields are different from what we have in the
                         * normal case (see the computation of gI), and it
                         * isn't neat.
                         *
                         * (on top of all that, notice that J is no
                         * longer something we have access to from this
                         * point of the code)
                         */
                        if (verbose) {
                            verbose_output_print(0, 1,
                                    "# small_sieve_init: not adding projective prime"
                                    " (1:%" FBROOT_FORMAT ") mod %" FBPRIME_FORMAT ")"
                                    " to small sieve  because g=%d >= J = %d\n",
                                    r_q-p, p, new_ssp.get_g(), J);
                        }
                        continue;
                    }
#endif
                    /* projective primes of all sorts go to ssp anyway */
                    ssd.ssp.push_back(new_ssp);
                } else if (new_ssp.is_pow2() || new_ssp.is_pattern_sieved()) {
                    ssd.ssp.push_back(new_ssp);
                } else if (new_ssp.is_nice()) {
                    ssd.ssps.push_back(new_ssp);
                    /* We're only pushing the offsets for the ssps anyway.  */
                    if (logI < LOG_BUCKET_REGION) {
                        /* The "up" offset is equal to r whenever logI >=
                         * LOG_BUCKET_REGION, and then we'll use r to compute
                         * the start positions.
                         */
                        int v = LOG_BUCKET_REGION - logI;
                        /* faster than a modular reduction as long as v<=10,
                         * which will be our use case anyway.
                         */
                        fbprime_t offset = r_q;
                        for(int k = 0 ; k < v ; k++) {
                            offset = offset + offset;
                            if (offset > p) offset -= p;
                        }
                        ASSERT(offset == (r_q << v) % p);
                        ssd.offsets.push_back(offset);
                    }
                    if (logI > LOG_BUCKET_REGION) {
                        /* The "right" offset is needed only in this case.
                         *
                         * It is such that (B+c, 0) is in L_p.
                         *
                         * Since LOG_BUCKET_REGION is a priori something as
                         * large as 16 bits, and p might also be of roughly
                         * the same size, maybe even more (up to logI), we
                         * must pay attention to overflows.
                         */
                        unsigned long offset = ((int64_t) (p - 1) << LOG_BUCKET_REGION) % p;
                        ssd.offsets.push_back(offset);
                    }
                } else {
                    ASSERT_ALWAYS(0);
                }
            }
        }
    }

    /* arrange so that the small_sieve() ctor is happy */
    /* I _think_ that normally, if the new code does its job correctly,
     * then this should be already sorted */
    ASSERT(std::is_sorted(ssd.ssps.begin(), ssd.ssps.begin() + ssd.resieve_end_offset));
    ASSERT(std::is_sorted(ssd.ssps.begin() + ssd.resieve_end_offset, ssd.ssps.end()));
    // std::sort(ssd.ssps.begin(), ssd.ssps.begin() + ssd.resieve_end_offset);
    // std::sort(ssd.ssps.begin() + ssd.resieve_end_offset, ssd.ssps.end());

    /* Sort general ssp vector in the order in which sieve2357::sieve expects
       them. small_sieve::do_pattern_sieve may drop some of these entries but
       preserves the ordering. */
    std::sort(ssd.ssp.begin(), ssd.ssp.end(), order_ssp);
}
/* }}} */

/* {{{ Creation of the ssdpos tables */

/*{{{ doc */
/* The places to be sieved are governed by the shape of the underlying
 * lattice of points.
 *
 * In the generic case, the lattice has determinant p (which may be a
 * prime power), and its basis in the (i,j) plane may be written as: (p
 * 0) (r 1) so that the lattice points are exactly the ones which satisfy
 * the equations i-r*j=0 mod p.  This means, in particular, that there
 * are hits on every line.
 *
 * With the extra assumption that p<=I and p<=BUCKET_REGION, this implies
 * that we'll have a starting point on the very first line of each bucket
 * region. However, versatility of the code may lead us to be slightly
 * tolerant about this assumption, and we'd like to see whether it is
 * possible to deal with larger p.
 *
 * The non-generic case is "projective". The lattice may be written as (q
 * 0) (U g) with q*g=p (which may be a prime power pp^k, and pp|g).
 *
 * What follows can be seen as a generalization of the generic case,
 * since it adapts by taking q=p and g=1.
 *
 * Note that when g>1, we no longer have a hit on every line. So the
 * first step is to try and see on which line we'll have a hit. Then the
 * (i,j) points will have to satisfy the equation i-(j/g)*U=0 mod q.
 *
 *
 * First task: the position of the first hit
 * -----------------------------------------
 *
 * We have a starting (current) bucket region with begin coordinates
 * (i0,j0). We want to find the first (i,j), >= (i0,j0) in processing
 * order, where p hits. Note that when this means j>j0, we may have i<i0.
 *
 * Specific for affine primes: in fact, because of the way we do the
 * sieve, we're rather going to compute the first i in the line j=j0,
 * relative to the current starting position i0, so that we have a hit.
 * The adjustments we make in the course of the computation will allow us
 * to catch up with the alignment constraints.
 *
 * Second task: actual sieve
 * -------------------------
 *
 * We have the first hit in ssdpos. We increment position by q (p in the
 * generic case), and update each point. The only difficult questions
 * are: - how do we deal with q even, given that we want to make sure we
 * don't sieve the (even,even) positions ?  - how do we deal, within the
 * sieve, with p reaching end of line ?  - when the position goes beyond
 * BUCKET_REGION, we need to move on to the next prime. We'll store this
 * "first off" position for further dealings with p. Do we need an
 * adjustment.  - the timing for end-of-bucket and end-of-line
 * adjustments is subtle and requires some care.
 *
 * Third task: reposition the sieve to jump over bucket regions
 * ------------------------------------------------------------
 *
 * In multithreaded mode, bucket regions are processed in an interleaved
 * fashion (it's not the only possible way to proceed). This means that
 * with N threads working, the position values per prime must be adjusted
 * after processing so that the N-1 bucket regions which have been sieved
 * by other threads are skipped.
 *
 * This adjustment may be made in a different manner depending on whether
 * we consider affine or projective primes.
 *
 *
 * end-of-bucket and end-of-line adjustments
 * -----------------------------------------
 *
 * This is a tricky part. When several bucket regions are in a line, once
 * we've marked hits in a bucket region and the sieve position has grown
 * past the end-of-bucket, the only needed adjustment is to subtract the
 * bucket region size from the position value: this will properly
 * represent the position, relative to the next bucket (if buckets are
 * processed in order). We do the same when there are no hits, so the
 * case "first hit is in the next bucket" does not need to be treated in
 * any specific way.
 *
 * On the other hand, at the end of a line, it's a bit different.
 * (i1-1)+j*I and i0+(j+1)*I are represented by consecutive bytes in
 * memory, so the coordinates (i1,j) and (i0,j+1) are equivalent
 * memory-wise. They're not equivalent with respect to the sieve
 * constraint that (X(i,j)=i-jr = 0 mod p). So for the "next position",
 * we have X(i,j+1)=-r if we don't pay attention. We fix that by adding r
 * to the sieved locations for line j+1, incorporating a reduction mod p
 * of the first sieved location.  This must be done for every line.
 *
 *
 */
/*}}}*/

/* Only compute the initial ssdpos fields. */
void small_sieve_start(std::vector<spos_t> & ssdpos,
        small_sieve_data_t & ssd,
        unsigned int first_region_index,
        int logI, sublat_t const & sl)
{
    /* We want to compute the index of the "next" hit, counted from the
     * starting offset of the "current" bucket region at (i0,j0). The
     * next hit means that it has to be the first in the current bucket
     * region, or in the first bucket region that comes after it in
     * processing order.
     *
     * In processing order, the memory location for cell (i,j), relative
     * to (i0,j0), is (i-i0) + (j-j0)*I -- note that (i-i0) may be
     * negative if j>j0.
     *
     */
    /* We store start positions for the simple case only.
     * For all other cases (which are rare enough
     * -- typically at most 20, counting powers of two and such), we
     *  compute the starting point from within sieve_small_bucket_region
     *  for each bucket region.
     */
    small_sieve_base C(logI, first_region_index, sl);
    ssdpos.clear();
    ssdpos.reserve(ssd.ssps.size());
    for (ssp_simple_t const & ssp : ssd.ssps) {
        ssdpos.push_back(C.first_position_ordinary_prime(ssp));
    }
}

void small_sieve_activate_many_start_positions(small_sieve_data_t & ssd)
{
    std::swap(ssd.ssdpos_many, ssd.ssdpos_many_next);
}

void small_sieve_prepare_many_start_positions(
        small_sieve_data_t & ssd,
        unsigned int first_region_index,
        int nregions,
        int logI, sublat_t const & sl)
{
    /* We're going to stage the next batch of init values in
     * ssdpos_many_next, while the init values in ssdpos_many are
     * potentially still in use, *except for the few last ones*, because by
     * design we always compute the "next" batch of init values as well.
     * Here, we mean that the init positions for a full row of bucket
     * regions above the bucket regions that are currently considered are
     * always computed. This means 1<<(max(0,logI-logB)) bucket regions.
     *
     * Therefore it is legitimate to steal these last elements, they were
     * constructed exactly for us !
     */

    auto & res(ssd.ssdpos_many_next);
    res.clear();

    int logB = LOG_BUCKET_REGION;
    int v = logI - logB;
    int w = (v > 0) ? (1 << v) : 1;

    res.assign(nregions + w, std::vector<spos_t>(ssd.ssps.size(), 0));

    int k;

    if (ssd.ssdpos_many.empty()) {
        small_sieve_start(res.front(), ssd, first_region_index, logI, sl);
        k = 0;
        /* must complete the first row */
        for(int i = 1 ; i < w ; i++) {
            for(size_t s = 0 ; s < ssd.ssps.size(); ++s) {
                fbprime_t x = res[k+i-1][s];
                x = x + ssd.offsets[s];
                if (x >= ssd.ssps[s].get_p()) x -= ssd.ssps[s].get_p();
                res[k+i][s] = x;
            }
        }
        k = w;
    } else {
        ASSERT_ALWAYS(ssd.ssdpos_many.size() >= (size_t) w);
        for(int i = 0 ; i < w ; ++i)
            std::swap(ssd.ssdpos_many[ssd.ssdpos_many.size()-w+i], res[i]);
    }
    ASSERT(res.front().size() == ssd.ssps.size());

    if (logI >= logB) {
        ASSERT(nregions % w == 0);
        for(int k = w; k < nregions + w ; k += w) {
            ASSERT(res[k].size() == ssd.ssps.size());
            /* Would it be possible to do all this with SIMD
             * instructions? It seems fairly likely, in fact.
             */

            /* infer from previous row of bucket regions.  */
            for(size_t s = 0 ; s < ssd.ssps.size(); ++s) {
                fbprime_t x = res[k-w][s];
                x = x + ssd.ssps[s].get_r();
                if (x >= ssd.ssps[s].get_p()) x -= ssd.ssps[s].get_p();
                res[k][s] = x;
            }

            /* Note that when logI == logB, we have w==1, so that ssd.offsets
             * is not needed at all.
             */
            for(int i = 1 ; i < w ; i++) {
                /* complete this row */
                for(size_t s = 0 ; s < ssd.ssps.size(); ++s) {
                    fbprime_t x = res[k+i-1][s];
                    x = x + ssd.offsets[s];
                    if (x >= ssd.ssps[s].get_p()) x -= ssd.ssps[s].get_p();
                    res[k+i][s] = x;
                }
            }
        }
    } else {
        for(int k = w; k < nregions + w ; k+= w) {
            /* infer from previous bucket region  */
            ASSERT(res[k].size() == ssd.ssps.size());
            for(size_t s = 0 ; s < ssd.ssps.size(); ++s) {
                fbprime_t x = res[k-w][s];
                /* ssd.offsets, not just ssd.ssps[s].get_r() */
                x = x + ssd.offsets[s];
                if (x >= ssd.ssps[s].get_p()) x -= ssd.ssps[s].get_p();
                res[k][s] = x;
            }
        }
    }
#ifndef NDEBUG
    for(int k = 0; k < nregions + w ; k++) {
        ASSERT(res[k].size() == ssd.ssps.size());
        int N = first_region_index + k;
        small_sieve_base C(logI, N, sl);
        for(size_t s = 0 ; s < ssd.ssps.size(); ++s) {
            ASSERT(res[k][s] == C.first_position_ordinary_prime(ssd.ssps[s]));
        }
    }
#endif
}

/* }}} */

void small_sieve::handle_projective_prime(ssp_t const & ssp, where_am_I & w MAYBE_UNUSED) {/*{{{*/
    /* This code also covers projective powers of 2 */
    const fbprime_t q = ssp.get_q();
    const fbprime_t g = ssp.get_g();
    const size_t gI = (size_t)ssp.get_g() << logI;
    const fbprime_t U = ssp.get_U();
    const fbprime_t p MAYBE_UNUSED = g * q;
    WHERE_AM_I_UPDATE(w, p, p);
    const unsigned char logp = ssp.logp;
    /* Sieve the projective primes. We have
     *         p^index | fij(i,j)
     * for i,j such that
     *         i * g == j * U (mod p^index)
     * where g = p^l and gcd(U, p) = 1.
     * This hits only for g|j, then j = j' * g, and
     *         i == j' * U (mod p^(index-l)).
     * In every g-th line, we sieve the entries with
     *         i == (j/g)*U (mod q).
     * In ssd we have stored g, q = p^(index-l), U, and ssdpos so
     * that S + ssdpos is the next sieve entry that needs to be
     * sieved.  So if S + ssdpos is in the current bucket region,
     * we update all  S + ssdpos + n*q  where ssdpos + n*q < I,
     * then set ssdpos = ((ssdpos % I) + U) % q) + I * g.  */
    if (!test_divisibility && ssp.get_q() == 1)
    {
        /* q = 1, therefore U = 0, and we sieve all entries in lines
           with g|j, beginning with the line starting at S[ssdpos] */
        unsigned long logps;
        long_spos_t pos = super::first_position_projective_prime(ssp);

        // The following is for the case where p divides the norm
        // at the position (i,j) = (1,0).
        if (UNLIKELY(super::has_origin && pos == (long_spos_t) gI)) {
#ifdef TRACE_K
            if (trace_on_spot_Nx(w->N, (1-i0))) {
                WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                unsigned int v = logp;
                sieve_increase_logging(S + w->x, v, w);
            }
#endif
            S[1 - i0] += logp;
        }
        // The event SSP_DISCARD might have occurred due to
        // the first row to sieve being larger than J. The row
        // number 0 must still be sieved in that case, but once
        // it's done, we can indeed skip the next part of
        // sieving.
        if (ssp.is_discarded_proj())
            return;
        ASSERT (ssp.get_U() == 0);
        ASSERT (pos % F() == 0);
        ASSERT (F() % (4 * sizeof (unsigned long)) == 0);
        for (size_t x = 0; x < sizeof (unsigned long); x++)
            ((unsigned char *)&logps)[x] = logp;
        unsigned int j = j0 + (pos >> logI);
        for( ; j < j1 ; j += g) {
            /* our loop is over line fragments that have a hit,
             * and by the condition q=1 above we'll sieve them
             * completely */
            unsigned long *S_ptr = (unsigned long *) (S + pos);
            unsigned long *S_end = (unsigned long *) (S + pos + F());
            unsigned long logps2 = logps;
            if (!(j&1)) {
                /* j is even. We update only odd i-coordinates */
                /* Use array indexing to avoid endianness issues. */
                for (size_t x = 0; x < sizeof (unsigned long); x += 2)
                    ((unsigned char *)&logps2)[x] = 0;
            }
#ifdef TRACE_K
            if (trace_on_range_Nx(w->N, pos, pos + F())) {
                WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                unsigned int x = trace_Nx.x;
                unsigned int index = x % I();
                unsigned int v = (((unsigned char *)(&logps2))[index%sizeof(unsigned long)]);
                sieve_increase_logging(S + w->x, v, w);
            }
#endif
            for( ; S_ptr < S_end ; S_ptr += 4) {
                S_ptr[0] += logps2;
                S_ptr[1] += logps2;
                S_ptr[2] += logps2;
                S_ptr[3] += logps2;
            }
            pos += gI;
        }
#if 0
        ssdpos[index] = pos - (1U << LOG_BUCKET_REGION);
#endif
    } else {
        /* q > 1, more general sieving code. */
        spos_t pos = super::first_position_projective_prime(ssp);
        const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
        unsigned char * S_ptr = S;
        S_ptr += pos - (pos & (I() - 1U));
        unsigned int j = j0 + (pos >> logI);
        pos = pos & (I() - 1U);
        ASSERT (U < q);
        for( ; j < j1 ; S_ptr += gI, j+=g) {
            WHERE_AM_I_UPDATE(w, j, j - j0);
            unsigned int q_or_2q = q;
            int i = pos;
            if (!(j&1)) {
                /* j even: sieve only odd i values */
                if (i % 2 == 0) /* Make i odd */
                    i += q;
                q_or_2q = evenq;
            }
            if ((i|j) & 1) { /* odd j, or not i and q both even */
                for ( ; i < F(); i += q_or_2q) {
                    WHERE_AM_I_UPDATE(w, x, (S_ptr - S) + i);
                    sieve_increase (S_ptr + i, logp, w);
                }
            }

            pos += U;
            if (pos >= (spos_t) q) pos -= q;
        }
#if 0
        ssdpos[index] = linestart + lineoffset - (1U << LOG_BUCKET_REGION);
#endif
    }
}/*}}}*/

/* {{{ Pattern-sieve primes with the is_pattern_sieved flag */
void small_sieve::do_pattern_sieve(where_am_I & w MAYBE_UNUSED)
{
    sieve2357base::prime_t psp[ not_nice_primes.size() + 1];

    unsigned int j = j0;

/* This can be enabled or disabled. If enabled, only location i = 1 will be
   updated in line j = 0. If disabled, the more general pattern-sieving
   code below will be used for line 0, too, but that will hit all locations
   with odd i. */
    if (skip_line_jj0 && j == 0 && super::sublatj0 == 0) {
        const int verbose = 0;
        WHERE_AM_I_UPDATE(w, j, 0);
        /* If sublattice, does this sublattice contain ii = 1 ?
           If we sieve fragments of a line, does this fragment contain
           i = 1? We assume that a fragment contains i = 1 iff it contains
           the origin */
        if ((super::sublatm == 1 || super::sublati0 == 1) &&
            super::has_origin) {
            for (auto const & ssp : not_nice_primes) {
                /* Primes that are not pattern-sieved are handled elsewhere */
                if (!ssp.is_pattern_sieved())
                    continue;
                /* Does this prime hit location ii = 1 ? In line jj = 0, all
                   prime (powers) q hit at ii = 0 (where norm = 0) and also
                   at ii = +-q, +-2q, +-3q, ... Thus the only case that hits
                   ii = 1 is a projective root with q = 1. */
                if (ssp.is_proj() && ssp.get_q() == 1) {
                    if (verbose) {
                        ssp.print(stdout);
                        printf(" hits at ii=1, jj=0\n");
                    }
                    WHERE_AM_I_UPDATE(w, p, ssp.get_g());
                    WHERE_AM_I_UPDATE(w, r, 0);
                    const unsigned int x = -super::i0 + 1;
                    ASSERT (x < (1U << LOG_BUCKET_REGION));
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_increase(S + x, ssp.logp, w);
                }
            }
        }
        j++;
    }

    for ( ; j < j1; j++) {
        size_t i = 0;   // Info: this is not an abscissa, here, but a plain counter
        const unsigned int dj = j - j0;
        const unsigned int jj = j * super::sublatm + super::sublatj0;
        const size_t x0 = (size_t) dj << logI;
        int skip_mod_2 = 0;
        WHERE_AM_I_UPDATE(w, j, dj);
#ifdef TRACE_K
        unsigned char orig_Sx = 0;
        if (trace_on_range_Nx(super::N, x0, x0 + super::F())) {
            orig_Sx = S[trace_Nx.x];
        }
#endif
        if (jj % 2 == 0) {
            if (super::sublatm == 3 && super::sublati0 == 1) {
                /* Skip odd indices which correspond to even ii */
                skip_mod_2 = 2;
            } else {
                /* Skip even indices which correspond to even ii */
                skip_mod_2 = 1;
            }
        }
        for (auto const & ssp : not_nice_primes) {
            ASSERT_ALWAYS(i < not_nice_primes.size());
            if (!ssp.is_pattern_sieved()) {
                /* Nothing to do here */
                continue;
            } else if (ssp.is_proj()) {
                if (jj % ssp.get_g() != 0) {
                    /* This projective root does not hit in this line */
                    continue;
                }
                WHERE_AM_I_UPDATE(w, p, ssp.get_q() * ssp.get_g());
                WHERE_AM_I_UPDATE(w, r, ssp.get_U());
                const fbprime_t pos = super::first_position_in_line(ssp, dj);
                psp[i++] = {ssp.get_q(), pos, ssp.logp};
            } else {
                WHERE_AM_I_UPDATE(w, p, ssp.get_p());
# if 0 /* what do about this one? */
                if (ssp.is_pow2() && sieve_only_odd)
                    continue;
#endif
                const fbprime_t pos = super::first_position_in_line(ssp, dj);
                psp[i++] = {ssp.get_p(), pos, ssp.logp};
            }
#ifdef TRACE_K
            if (trace_on_range_Nx(super::N, x0, x0 + super::F())) {
                /* We are in the correct line (fragment). */
                const fbprime_t q = psp[i - 1].q, pos = psp[i - 1].idx;
                const unsigned char logp = psp[i - 1].logp;
                if (0) {
                    printf("# Pattern sieve side %i, line %u (N=%d, x0=%zu, trace_Nx.x=%u): Adding psp[%zu] = {%" FBPRIME_FORMAT", %" FBPRIME_FORMAT", %hhu}, from  ",
                        w->side, jj, super::N, x0, trace_Nx.x, i - 1, q, pos, logp);
                    ssp.print(stdout);
                    printf("\n");
                }
                if ((trace_Nx.x - x0) % q == pos &&
                    (skip_mod_2 == 0 || trace_Nx.x % 2 == (unsigned) skip_mod_2 % 2)) {
                    WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                    sieve_increase(S + trace_Nx.x, logp, w);
                }
            }
#endif
        }
#ifdef TRACE_K
        unsigned char new_Sx = 0;
        if (trace_on_range_Nx(super::N, x0, x0 + super::F())) {
            new_Sx = S[trace_Nx.x];
            S[trace_Nx.x] = orig_Sx;
        }
#endif
        psp[i++] = {0, 0, 0};
        preferred_sieve2357::sieve(
            (preferred_simd_type *) (S + x0), F(), psp,
            skip_mod_2, sieve2357base::update_add, w);
#ifdef TRACE_K
        if (trace_on_range_Nx(super::N, x0, x0 + super::F())) {
            ASSERT_ALWAYS(new_Sx == S[trace_Nx.x]);
        }
#endif
    }
}
/* }}} */


/* {{{ Normal small sieve */
// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, unsigned int N,
                               small_sieve_data_t const & ssd,
                               std::vector<spos_t> const & ssdpos,
                               int logI, sublat_t const & sl,
                               where_am_I & w)
{
    small_sieve SS(ssdpos, ssd.ssps, ssd.ssp, S, logI, N, sl);
    SS.do_pattern_sieve(w);
    SS.exceptional_sieve(w);
    SS.normal_sieve(w);
}
/*}}}*/

/* {{{ resieving. Different interface, since it plays with buckets as well.
 */


/* SMALLSIEVE_COMMON_DEFS() is defined in
 * sieve/las-smallsieve-lowlevel.hpp because the testing framework uses
 * it.
 */
/* Sieve small primes (p < I) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd. */
void
resieve_small_bucket_region (bucket_primes_t *BP,
        unsigned char *S, unsigned int N,
        small_sieve_data_t & ssd,
        std::vector<spos_t> const & ssdpos,
        int logI, sublat_t const & sl,
        where_am_I & w MAYBE_UNUSED)
{
    SMALLSIEVE_COMMON_DEFS();
    small_sieve_base C(logI, N, sl);

    unsigned char *S_ptr;
    const int resieve_very_verbose = 0;

    unsigned int i_compens_sublat = sublati0 & 1;

    // Odd/even property of j is the same as for j+2, even with
    // sublat, unless sublat.m is even, which is not handled right
    // now. Same for i.
    ASSERT_ALWAYS(!sublatm || ((sublatm & 1) == 1));

    for(size_t index = 0 ; index < ssd.resieve_end_offset ; index++) {
        auto const & ssps(ssd.ssps[index]);

        const fbprime_t p = ssps.get_p();
        fbprime_t r = ssps.get_r();
        WHERE_AM_I_UPDATE(w, p, p);
        int pos = ssdpos[index];
        S_ptr = S;
        ASSERT(pos < (spos_t) p);
        /* for j even, we sieve only odd index. This translates into loops
         * which look as follows:
         *
         * j even: (sieve only odd index)
         *   for(index = pos + (p & -!(pos&1)) ; index < I ; index += p+p)
         *   (where (p & -!(pos&1)) is 0 if pos is odd, and p otherwise)
         * j odd: (sieve all values of index)
         *   for(index = pos                  ; index < I ; index += p)
         *
         * we may merge the two by setting q=p&-!((j&1)^row0_is_oddj)
         *
         * which, when (j+row0_is_oddj) is even, is p, and is 0
         * otherwise.
         *
         * In turn, since q changes for each j, 1 xor within the loop
         * is enough to make it alternate between 0 and p, once the
         * starting value is correct.
         *
         * TODO: ok, this is nice and good, but:
         *
         *  - I haven't seen this win for simple small sieve, so I
         *    doubt it's a good idea. Relying on the branch predictor
         *    or the compiler does not seem to be so stupid after all.
         *  - as present, the behaviour is obviously buggy for
         *    sublatm even.
         *  - we really want to have the same structure both for
         *    small sieve and resieving.
         */
        bool row0_even = (((j0&sublatm)+sublatj0) & 1) == 0;
        unsigned int q = row0_even ? p : 0;
        int i = 0; /* placate gcc */
        for (unsigned int j = j0; j < j1; j ++) {
            WHERE_AM_I_UPDATE(w, j, j);
            for (i = pos + (q& -!((pos+i_compens_sublat)&1)) ; i < C.F() ; i += p+q) {
                if (LIKELY(S_ptr[i] == 255)) continue;
                bucket_update_t<1, primehint_t> prime;
                unsigned int x = ((size_t) (j-j0) << logI) + i;
                if (resieve_very_verbose) {
                    verbose_output_print(0, 1, "resieve_small_bucket_region: root %"
                            FBROOT_FORMAT ",%" FBPRIME_FORMAT" divides at x = "
                            "%d = %u * %u + %d\n",
                            p, r, x, j, I, i);
                }
                prime.p = p;
                prime.x = x;
                ASSERT(prime.p >= ssd.fbK.td_thresh);
                BP->push_update(prime);
            }
            pos += r;
            if (pos >= (spos_t) p) pos -= (spos_t) p;

            S_ptr += I;
            q ^= p;
        }
    }

    for(auto const & ssp : ssd.ssp) {
        /* "not nice" cases are either projective or power of two. we
         * obviously won't resieve powers of two, so we're bound to deal
         * with only projective primes here.
         */
        ASSERT(ssp.is_pow2() || ssp.is_proj() || ssp.is_pattern_sieved());

        /* FIXME: I should not have to do this test */
        if (ssp.is_pow() || ssp.is_pow2())
            continue;

        /* TODO: it doesn't seem very smart to resieve projective primes
         */
        if (ssp.is_proj()) {
            const fbprime_t g = ssp.get_g();
            const fbprime_t q = ssp.get_q();
            const fbprime_t p = g * q;

            /* the code below definitely does not deal with projective
             * primes that do not sieve full lines.  */
            if (q > 1) continue;
            if (p < ssd.fbK.td_thresh) continue;

            WHERE_AM_I_UPDATE(w, p, p);

            const uint64_t gI = (uint64_t)g << logI;

            /* Test every p-th line, starting at S[ssdpos] */
            long_spos_t pos = C.first_position_projective_prime(ssp);
            // This block is for the case where p divides at (1,0).
            if (UNLIKELY(has_origin && pos == (long_spos_t) gI)) {
                bucket_update_t<1, primehint_t> prime;
                prime.p = p;
                prime.x = 1 - i0;
                ASSERT(prime.p >= ssd.fbK.td_thresh);
                BP->push_update(prime);
            }
            // Same as in sieving: we discard after checking for row 0.
            if (ssp.is_discarded_proj())
                continue;

            /* make sure ssdpos points at start of line or region when
             * we're sieving whole lines. */
            ASSERT(q > 1 || !(pos % (1 << MIN(logI, LOG_BUCKET_REGION))));

            if (resieve_very_verbose) {
                verbose_output_print(0, 1, "# resieving projective prime %" FBPRIME_FORMAT
                        ", i0 = %" PRIi64 "\n", q, i0 + pos);
            }
            if (pos >> LOG_BUCKET_REGION)
                continue;
            unsigned int j = j0 + (pos >> logI);
            for (; j < j1; j += g) {
                unsigned char *S_ptr = S + pos;
                /* FIXME: sublat */
                if (!(j & 1)) {
                    /* Even j: test only odd ii-coordinates */
                    for (int ii = 1; ii < C.F(); ii += 2) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region even j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        BP->push_update(prime);
                    }
                } else {
                    /* Odd j: test all ii-coordinates */
                    for (int ii = 0; ii < C.F(); ii++) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region odd j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        BP->push_update(prime);
                    }
                }
                pos += gI;
            }
#if 0
            ssdpos[index] = pos - bucket_region;
            if (resieve_very_verbose) {
                verbose_output_print(0, 1, "# resieving: new pos = %" PRIu64
                        ", bucket_region = %d, "
                        "new ssdpos = %" PRIu64 "\n",
                        pos, bucket_region, ssdpos[index]);
            }
#endif
        }
    }
}
/* }}} */
