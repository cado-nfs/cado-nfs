#ifndef CADO_LAS_SMALLSIEVE_GLUE_HPP
#define CADO_LAS_SMALLSIEVE_GLUE_HPP

#include <cstdio>
#include <climits>
#include <cstdint>

#include <vector>
#include <list>

#ifdef HAVE_SSE41
#include <x86intrin.h>
#endif

#include "las-forwardtypes.hpp"         // spos_t
#include "fb-types.hpp"
#include "gcd.h"
#include "las-smallsieve-lowlevel.hpp"
#include "las-smallsieve.hpp"

#include "macros.h"

/* generated code for small sieve critical routine. */

struct list_nil {};

template<typename T, typename U, int b, typename F> struct choice_list_car {};

struct small_sieve_best_code_choices {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
// && !defined(TRACK_CODE_PATH)

/* The selection below is really scrapping everything we've done to
 * design a code path by bits and pieces, and strives to use almost a
 * one-size-fits-all approach (which admittedly also has its advantages).
 *
 * currently it does seem to be the case that this one-size-fits-all
 * choice wins over the tailored-fit choice above (at least on my
 * laptop), reaching a 12% win as per the report of test-smallsieve -A 31
 * -I 16 -B 16 --only-complete-functions (commit faa0396)
 */

    typedef choice_list_car<
                assembly_generic_oldloop,
                assembly_generic_oldloop,
                5, /* use the two functions above for primes that have a
                      max number of hits on odd lines that is at least 32 */
            choice_list_car<
                assembly3, /* 4 to 8 hits on even lines */
                assembly4, /* 16 to 32 hits on odd lines */
                4, /* use the two functions above for primes that have a
                      max number of hits on odd lines that is at least 16 */
            choice_list_car<
                assembly2, /* 2 to 4 hits on even lines */
                assembly3, /* 4 to 8 hits on odd lines */
                3, /* use the two functions above for primes that have a
                      max number of hits on odd lines that is at least 8 */
            choice_list_car<
                assembly1x, /* 0 to 2 hits on even lines */
                assembly2x, /* 0 to 4 hits on odd lines */
                0, /* use the two functions above for the remaining primes */
            list_nil
            >>>> type;
#else
    typedef choice_list_car<
                manual_oldloop,
                manual_oldloop,
                0,
            list_nil
            > type;
#endif
};

/*{{{ small sieve classes */
/* So many things are used in common for many small sieve routines that
 * it makes sense to gather them in a common object */
struct small_sieve_base {/*{{{*/
    int min_logI_logB;
    int logI;
    int N;
    unsigned int log_regions_per_line;
    unsigned int region_rank_in_line;
    bool last_region_in_line;
    unsigned int j0;
    unsigned int j1;
    int i0;
    // int i1;
    /* those are (1,0,0) in the standard case */
    sublat_runtime_t sublat;
    bool has_origin;
    inline int F() const { return 1 << min_logI_logB; }
    inline int I() const { return 1 << logI; }
    static const bool skip_line_jj0 = false;
    /* returned by first_position_power_of_two() when the prime power has
     * no hit at all in the residue class being sieved, which can happen
     * once the sublattice modulus is even */
    static constexpr spos_t NO_POSITION = INT32_MIN;
    small_sieve_base(int logI, int N, sublat_runtime_t const & sublat)/*{{{*/
        : logI(logI)
        , N(N)
        , sublat(sublat)
    {
        unsigned int log_lines_per_region;
        //min_logI_logB         = is_fragment::test(LOG_BUCKET_REGION, logI, MIN(LOG_BUCKET_REGION, logI));
        min_logI_logB         = std::min(LOG_BUCKET_REGION, logI);
        log_lines_per_region  = (LOG_BUCKET_REGION-min_logI_logB);
        log_regions_per_line  = (logI-min_logI_logB);
        unsigned int regions_per_line      = (1<<log_regions_per_line);
        region_rank_in_line   = (N&(regions_per_line-1));
        // last_region_in_line   = (region_rank_in_line==(regions_per_line-1));
        j0    = ((N>>log_regions_per_line)<<log_lines_per_region);
        j1    = (j0+(1<<log_lines_per_region));
        // I     = (1<<logI);
        i0    = ((region_rank_in_line<<LOG_BUCKET_REGION)-(1 << (logI-1)));
        // i1    = (i0+(1<<min_logI_logB));

        bool has_haxis = !j0;
        bool has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);
        has_origin = has_haxis && has_vaxis;
        ASSERT(sublat.m == 1 || sublat.m == 2 || sublat.m == 3 || sublat.m == 6);
    }/*}}}*/

    /* {{{ parity of the real (i,j) coordinates
     *
     * The sieve position x on row dj of this bucket region has real
     * coordinates
     *
     *      ii = sublat.m * (i0 + x)  + sublat.i0
     *      jj = sublat.m * (j0 + dj) + sublat.j0
     *
     * and the small sieve never has to touch a position where ii and jj
     * are both even, since gcd(ii,jj) >= 2 there. Note that i0 is always
     * even, so it never contributes to the parity of ii.
     *
     * When sublat.m is odd -- which includes the ordinary case
     * sublat.m == 1 -- ii and jj have the parities of x + sublat.i0 and
     * dj + sublat.j0. So the parity of a row flips as dj grows, and on an
     * even row every other x must be skipped.
     *
     * When sublat.m is even, both parities are *constant* over the whole
     * sublattice: ii has the parity of sublat.i0, jj that of sublat.j0.
     * The caller never asks for the class where both are even (las.cpp
     * skips i_cong == j_cong == 0), so at least one of the two is odd
     * everywhere in the class, and there is nothing to skip on any row.
     *
     * That is why an even modulus -- in particular m == 2, the classical
     * "oddness type" decomposition used by the Franke-Kleinjung sievers
     * -- is the easy case rather than a harder one: the alternating
     * logic simply does not apply.
     */
    bool has_even_sublatm() const { return (sublat.m & 1) == 0; }

    /* Does row dj == 0 of this region need the skip-every-other-x
     * treatment, and does that property flip from one row to the next? */
    bool row0_needs_parity_skip() const {
        if (has_even_sublatm()) {
            ASSERT(((sublat.i0 | sublat.j0) & 1) != 0);
            return false;
        }
        return ((j0 + sublat.j0) & 1) == 0;
    }
    bool parity_skip_alternates() const { return !has_even_sublatm(); }

    /* On a row whose jj is even, which x have ii even and must therefore
     * be skipped: 0 for none at all, 1 for the even x, 2 for the odd x. */
    int parity_skip_class() const {
        if (has_even_sublatm()) return 0;
        return (sublat.i0 & 1) ? 2 : 1;
    }
    /* }}} */

    /* {{{ prime powers that are not coprime to the sublattice modulus
     *
     * Write g = gcd(q, sublat.m). The congruence
     *
     *      ii == r * jj  (mod q),   ii = m*x + sublat.i0, jj = m*y + sublat.j0
     *
     * becomes, in sublattice coordinates,
     *
     *      m*x == m*r*y + c  (mod q),      c = r*sublat.j0 - sublat.i0
     *
     * which has a solution only when g divides c, and which then reads
     *
     *      x == r*y + (c/g) * (m/g)^-1   (mod q/g).
     *
     * So a prime power q behaves, in sublattice coordinates, exactly like
     * the prime power q/g -- same root, shifted offset -- and either it is
     * absent from this residue class altogether, or it hits with stride
     * q/g rather than q. Verified exhaustively for m in {2,3,6} by
     * bench/cado/verify-sublat-powers.py in the siever-archaeology repo.
     *
     * For g = 1, which is every case except a prime dividing the modulus,
     * all of this collapses to the identity and nothing changes.
     *
     * Concretely for m = 2 and q = 2^k: the class (1,0) -- i odd, j even --
     * never has any power of two dividing, while the classes (0,1) and
     * (1,1) have stride 2^(k-1) for the roots of the right parity.
     */
    fbprime_t sublat_reduced_modulus(fbprime_t q) const {
        return q / (fbprime_t) gcd_ul(q, (unsigned long) sublat.m);
    }

    /* Like fix_sublat_i(), but works when q and sublat.m share a factor,
     * and says so when the congruence has no solution in this class
     * (in which case *res is untouched). */
    bool fix_sublat_i_checked(int64_t ii, fbprime_t q, int & res) const {
        unsigned long const g = gcd_ul(q, (unsigned long) sublat.m);
        if (g == 1) {
            res = fix_sublat_i(ii, q);
            return true;
        }
        /* Adding q to ii cannot change ii mod g, so either it is already
         * right or there is nothing to be done. */
        int64_t d = ii - (int64_t) sublat.i0;
        if (d % (int64_t) g)
            return false;
        res = fix_sublat_i(ii, q);
        return true;
    }
    /* }}} */

    /* {{{ entries that are sieved row by row under sublattices
     *
     * This covers the entries that small_sieve_init() marks with
     * set_sublat_rowwise(): prime powers whose prime divides the
     * sublattice modulus (and which are not flat over the class), and
     * the projective entries that are not pattern-sieved.
     *
     * The condition for a hit is
     *
     *      ii == T  (mod q)
     *
     * where T = r*jj and q = p for an affine entry, and T = (jj/g)*U for
     * a projective one, which furthermore needs g | jj. With
     * ii = m*x + sublat.i0, this becomes
     *
     *      m*x == T - sublat.i0  (mod q)
     *
     * which, with d = gcd(m, q), is solvable iff d divides the right
     * hand side, and then reads
     *
     *      x == ((T - sublat.i0)/d) * (m/d)^-1  (mod q/d).
     *
     * (m/d is coprime to q/d because m is squarefree.) This is the
     * statement of sublat_reduced_modulus() above, with the row
     * condition of projective entries on top.
     *
     * There are few such entries, so we simply do this row by row. We
     * do not bother with skipping the positions where ii and jj are both
     * even (the survivor search skips them anyway).
     *
     * f is called with the row index relative to j0, the first position
     * in that row, and the stride.
     */
    template<typename F>
    void foreach_sublat_rowwise_row(ssp_t const & ssp, F && f) const
    {
        bool const proj = ssp.is_proj();
        uint64_t const q = proj ? ssp.get_q() : ssp.get_p();
        uint64_t const g = proj ? ssp.get_g() : 1;
        uint64_t const U = proj ? ssp.get_U() : ssp.get_r();
        uint64_t const m = sublat.m;
        uint64_t const d = gcd_ul(q, m);
        uint64_t const s = q / d;
        uint64_t const md = m / d;
        /* (m/d)^-1 mod s. m/d is 1, 2 or 3, so one of s+1 and 2*s+1 is
         * a multiple of it. */
        uint64_t minv = 0;
        for (uint64_t t = 0; t < md; t++) {
            if ((t * s + 1) % md == 0) {
                minv = ((t * s + 1) / md) % s;
                break;
            }
        }
        uint64_t const si0 = sublat.i0 % q;
        for (unsigned int j = j0; j < j1; j++) {
            uint64_t const jj = uint64_t(j) * m + sublat.j0;
            if (jj % g)
                continue;
            uint64_t const c = ((jj / g) % q * U + q - si0) % q;
            if (c % d)
                continue;
            int64_t const x = int64_t((c / d) * minv % s);
            int64_t pos = (x - i0) % int64_t(s);
            if (pos < 0)
                pos += s;
            f(j - j0, pos, (int64_t) s);
        }
    }
    /* }}} */

    /* Returns (ii - sublat.i0 + k*q) / sublat.m with k >= 0 minimal so that
       result is an integer */
    int fix_sublat_i(int64_t ii, const fbprime_t q) const {
        if (sublat.m == 1) {
            ASSERT(ii <= INT_MAX);
            return (int) ii;
        } else if (sublat.m == 2) {
            ASSERT(q % 2 != 0 || ii % 2 == sublat.i0);
            for( ; ii % 2 != sublat.i0 ; ii += q);
            ii = (ii - sublat.i0) / 2;
            ASSERT(ii <= INT_MAX);
            return ii;
        } else if (sublat.m == 3) {
            ASSERT(q % 3 != 0 || ii % 3 == sublat.i0);
            for( ; ii % 3 != sublat.i0 ; ii += q);
            ii = (ii - sublat.i0) / 3;
            ASSERT(ii <= INT_MAX);
            return ii;
        } else if (sublat.m == 6) {
            ASSERT(q % 2 != 0 || ii % 2 == sublat.i0);
            ASSERT(q % 3 != 0 || ii % 3 == sublat.i0);
            for( ; ii % 6 != sublat.i0 ; ii += q);
            ii = (ii - sublat.i0) / 6;
            ASSERT(ii <= INT_MAX);
            return ii;
        } else {
            ASSERT_ALWAYS(false);
        };
    }

    spos_t first_position_ordinary_prime(ssp_simple_t const & ssp) const/*{{{*/
    {
        /* equation here: i-r*j = 0 mod p */

        /* Expanded wrt first row j0 and sublats: smallest pos0 s.t.:
         *  sublat.m*(i0+pos0)+sublat.i0 - r*(j0*sublat.m+sublat.j0) = -k*p
         *  (p is coprime to sublat.m).
         *
         * In non-sublat mode, this means:  (i0+pos0)-r*j0 = 0 mod p
         * so that pos0 = (r*j0 - i0) mod p where i0 is traditionally
         * negative (that helps for sign stuff, of course).
         *
         * So S = sublat.m*pos0 must be such that:
         *
         * S = r*(j0*sublat.m+sublat.j0) - (sublat.m*i0+sublat.i0) + k*p
         *
         * pos0 being an integer, of course. So there's only one
         * congruence class for k mod sublat.m such that the above
         * happens. We start by computing it.
         *
         *  k:=((sublat.i0 - r * sublat.j0) * invp_mod_m) mod sublat.m;
         *
         * then we deduce pos0 as:
         *
         *  pos0:=j0 * r - i0 + (((r*sublat.j0-sublat.i0)+k*p) div
         *  sublat.m)
         *
         * and finally we take pos0 mod p.
         */
        /* Ouch. That can't work if we have largeish small-sieved primes!
         */
        // spos_t x = (spos_t)(j0) * (spos_t)ssp.get_r() - i0;
        int64_t x = (int64_t)(j0) * (int64_t)ssp.get_r() - i0;
        if (sublat.m > 1) {
            ASSERT(ssp.get_p() % sublat.m);
            /* alternative code. not clear it's better.
               ASSERT(sublat.m == 2 || sublat.m == 3 || sublat.m == 6);
               uint64_t invp_mod_m = p;
               int k = ((sublat.i0 - r * sublat.j0) * invp_mod_m) % sublat.m;
               x += ((r * sublat.j0 - sublat.i0) + k * p) / sublat.m;
               */
            spos_t y = ssp.get_r() * sublat.j0;
            x += fix_sublat_i(y, ssp.get_p());
        }
        x = x % (int64_t)ssp.get_p();
        /* As long as i0 <= 0, which holds in the normal case where
         * logI <= LOG_BUCKET_REGION, we can be sure that x >= 0, so
         * we have no issue with the sign.  However, when i0 is a
         * positive number, some extra care is needed.
         */
        if (x < 0) x += ssp.get_p();
        return (spos_t) x;
    }/*}}}*/

    /* This return value is typically logI bits larger than for ordinary
     * primes, so it makes sense to return it as a 64-bit integer.
     *
     * (more precisely, if p is, say, a 25-bit prime, and logI is 20,
     * then we return 45 bits. note however that projective primes are
     * not typically that big -- but undoubtedly it may happen).
     *
     * Quite probably we would have an interest in doing those projective
     * primes another way, and not deal with artificially shifted values.
     */
    int64_t first_position_projective_prime(ssp_t const & ssp) const/*{{{*/
    {
        /* equation here: i == (j/g)*U (mod q) */
        /* we're super-non-critical, here.
         *
         * Note that the code below could actually cover the generic
         * case as well, if we adapt by letting g=1, q=p, U=r.
         */

        /*
         * Because of the sublat feature, (i0,j0) actually correspond
         * to "expanded" indices computed as:
         *      ii = i0 * sublat.m + sublat.i0
         *      jj = j0 * sublat.m + sublat.j0
         */

        unsigned int jj = j0*sublat.m + sublat.j0;

        /* First question is to find the next multiple of g above jj.
         * This is done as follows */

        /* next multiple of g above j0 */
        unsigned int jjmod = jj % ssp.get_g();
        if (jjmod) jj += ssp.get_g();
        jj -= jjmod;

        /* All of the fragment above trivializes to a no-op in the
         * non-projective case (which has g=1)
         */

        /* Now we'd like to avoid row number 0 (so jj == 0). */
        /* The sieving code may do some special stuff to fill the
         * position (1,0). At least at some point it did */
        if (skip_line_jj0 && jj == 0) jj += ssp.get_g();

        // In sublat mode, we also need jj congruent to sublat.j0 mod m.
        // XXX A very nasty situation: when ssp.g and sublat.m are
        // not coprime, we may very well have an infinite loop here
        // (say we want to sieve only even lines and sublat.j0=1,
        // sublat.m=1. For the moment, we do SSP_DISCARD_SUBLAT
        // whenever p and sublat.m have a common divisor. If we want
        // to be finer grain, we need to selectively discard some
        // small sieved primes depending on the sublattice we're
        // considering (or but ssdpos to ULONG_MAX ?).
        if (sublat.m > 1) {
            for( ; jj % sublat.m != sublat.j0 ; jj += ssp.get_g());
        }
        // Find the corresponding i
        int ii = int(jj/ssp.get_g())*int(ssp.get_U());

        ii = fix_sublat_i(ii, ssp.get_q());
        // In the sublat mode, switch back to reduced convention
        // (exact divisions)
        if (sublat.m > 1) {
            jj = (jj-sublat.j0) / sublat.m;
        }

        /* At this point we know that the point (ii,jj) is one of the
         * points to be sieved, and jj-j0 perhaps is within reach for
         * this bucket region (not always). How do we adjust with
         * respect to this region's i0 ?
         */

        /* Given the value of ii, we need to know the next multiple
         * of q which is relative to the current i0. So something
         * congruent to ii mod q, but little above i0.
         */
        /* when skipping to the next line, we want to adjust with
         * respect to the very first possible i0, so virtually with
         * i0=-I/2.
         */
        int i0ref = (j0 == jj) ? i0 : (-I()/2);
        int64_t x = (ii-i0ref) % (uint64_t)ssp.get_q();
        if (x < 0) x += ssp.get_q();
        if (jj > j0) {
            x -= region_rank_in_line << LOG_BUCKET_REGION;
            x += (((uint64_t)(jj - j0))<<logI);
        }
        return x;
    }/*}}}*/

    spos_t first_position_power_of_two(ssp_t const & ssp) const/*{{{*/
    {
        /* equation here: i-r*j = 0 mod p, p a power of 2. */
        /* only difference with ordinary case is that we want to
         * sieve only odd lines.
         */
        /* For powers of 2, we sieve only odd lines (*) and 
         * ssdpos needs to point at line j=1. We assume
         * that in this case (I/2) % p == 0
         * (*) for lines with j even, we have a root mod the prime
         * power for i-j*r multiple of our power of 2, which means
         * i even too. Thus a useless report.
         */
        /* With an even sublattice modulus the parity of the real j is
         * fixed over the whole class, so either every row is odd and
         * there is no "next odd line" to look for, or no row is and this
         * power of two has no hit here at all. See the discussion of
         * sublat_reduced_modulus() above. */
        unsigned int j = j0;
        uint64_t jj = j*sublat.m + sublat.j0;
        // uint64_t ii = i0*sublat.m + sublat.i0;
        /* next odd line */
        // This was: jj |= 1;
        int i0ref = i0;
        if ((jj & 1) == 0) {
            if (has_even_sublatm())
                return NO_POSITION;
            jj += sublat.m;
            i0ref = (-I()/2);
        }
        fbprime_t const q = sublat_reduced_modulus(ssp.get_p());
        int xi;
        if (!fix_sublat_i_checked((spos_t)jj * (spos_t)ssp.get_r(),
                                  ssp.get_p(), xi))
            return NO_POSITION;
        spos_t x = xi;
        x = (x - i0ref) & (q - 1);
        if (x < 0) x += q;
        /* our target is position x in the bucket region which starts
         * at coordinates (i0ref, jj). How far is that from us ?
         */
        // The condition is: if ((jj / sublat.m) > j)
        // but we separate the case sublat.m=1 to avoid the div in the
        // general case.
        if (((sublat.m == 1) && (jj > j)) ||
                ((jj / sublat.m) > j)) {
            x -= region_rank_in_line << LOG_BUCKET_REGION;
            x += (((uint64_t)((jj / sublat.m) - j0))<<logI);
            /* For the case of several bucket regions per line, it's
             * clear that this position will be outside the current
             * bucket region. Note that some special care is needed
             * when we want to incrementally update the position, but we
             * rarely (if ever) want to do that anyway.
             */
        }
        return x;
    }/*}}}*/

    /* This function returns the first location hit, relative to the start
       of the current *line fragment*. dj is the line number being sieved
       relative to the bucket region start. If a bucket region contains only
       one line, or only a fragment of a line, then dj must be 0. If ssp is
       a projective root, then the line index jj being sieved must be
       divisible by ssp.g. Line jj=0 is not treated any different than non-
       zero lines. For powers of 2, even lines are not skipped by this
       method - the caller needs to do it, if desired. */
    fbroot_t first_position_in_line(ssp_t const & ssp, unsigned int dj = 0) const
    {
        const unsigned int j = j0 + dj;
        const unsigned int jj = j * sublat.m + sublat.j0;

        if (ssp.is_proj()) {
            ASSERT (jj % ssp.get_g() == 0);

            const int64_t ii = int64_t(jj / ssp.get_g()) * int64_t(ssp.get_U());
            const int i = fix_sublat_i(ii, ssp.get_q());
            int64_t x = (i - i0) % (uint64_t)ssp.get_q();
            if (x < 0) x += ssp.get_q();
            ASSERT_ALWAYS(0 <= x && x < int64_t(ssp.get_q()));
            return (fbroot_t) x;
        } else {
            /* Affine roots hit in every line, but not necessarily in every line fragment */
            const int64_t ii = int64_t(jj) * int64_t(ssp.get_r());
            const int i = fix_sublat_i(ii, ssp.get_p());
            int64_t x = (i - i0) % (uint64_t)ssp.get_p();
            if (x < 0) x += ssp.get_p();
            ASSERT_ALWAYS(0 <= x && x < int64_t(ssp.get_p()));
            return (fbroot_t) x;
        }
    }

};/*}}}*/

struct small_sieve : public small_sieve_base {/*{{{*/
    typedef small_sieve_base super;
    std::vector<int> const & positions;
    std::vector<ssp_simple_t> const& primes;
    size_t sorted_limit;
    std::list<size_t> sorted_subranges;
    std::vector<ssp_t> const& not_nice_primes;
    unsigned char * S;

    static const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    static const int test_divisibility = 0; /* very slow, but nice for debugging */

    /* This "index" field is actually the loop counter that is shared by
     * the various instantiations that are triggered by
     * handle_nice_primes_meta_loop<>.
     */
    size_t index = 0;

    small_sieve(std::vector<int> const& positions,
            std::vector<ssp_simple_t> const & primes,
            std::vector<ssp_t> const & not_nice_primes,
            unsigned char*S, int logI,
            unsigned int N,
            sublat_runtime_t const & sublat
            )
        : super(logI, N, sublat),
            positions(positions),
            primes(primes), 
            not_nice_primes(not_nice_primes),
            S(S) {
                auto s = begin(primes);
                auto s0 = s;
                auto e = end(primes);
                for( ; e > s ; ) {
                    if ((e-s) < 32) {
                        /* we don't want to bother adding an extra
                         * control loop for a small bunch of primes 
                         */
                        break;
                    }
                    int c;
                    for(c = 1 ; c < (e-s) && !(s[c] < s[c-1]) ; c++);
                    if (c <= 16) {
                        if (s != s0) {
                            /* s != s0 : cf mail ET->{PG,PZ} 201906071825
                             */
                            fprintf(stderr, "warning, the prime list looks really ugly\n");
                        }
                    }
                    /*
                    fprintf(stderr, "ssp entries [%zd..%zd[ (out of %zu) are sorted\n",
                            s-s0, s+c-s0, primes.size());
                            */
                    sorted_subranges.push_back((s+c) - s0);
                    s += c;
                }
                /*
                fprintf(stderr, "ssp entries [%zd..%zd[ (tail) do not need to be sorted\n",
                        s-s0, primes.size());
                        */
            }

    bool finished() const { return index == primes.size(); }
    bool finished_sorted_prefix() const { return index == sorted_limit; }

    /* most fields from the parent class must be redeclared here if we
     * want to use them transparently -- it's only a shorthand
     * convenience. We don't re-expose everything though, since some are
     * only seldom used.
     */
    using super::logI;
    using super::i0;
    using super::j0;
    using super::j1;
    using super::sublat;
    using super::I;
    using super::F;
    using super::skip_line_jj0;

    void handle_projective_prime(ssp_t const & ssp, where_am_I & w);
    void handle_sublat_rowwise(ssp_t const & ssp, where_am_I & w MAYBE_UNUSED) {/*{{{*/
        WHERE_AM_I_UPDATE(w, p, ssp.is_proj() ? ssp.get_g() * ssp.get_q() : ssp.get_p());
        const unsigned char logp = ssp.logp;
        super::foreach_sublat_rowwise_row(ssp,
                [&](unsigned int dj, int64_t pos, int64_t s) {
                    unsigned char * S_ptr = S + ((size_t) dj << logI);
                    for ( ; pos < F() ; pos += s) {
                        WHERE_AM_I_UPDATE(w, x, ((size_t) dj << logI) + pos);
                        sieve_increase (S_ptr + pos, logp, w);
                    }
                });
    }/*}}}*/
    void handle_power_of_2(ssp_t const & ssp, where_am_I & w MAYBE_UNUSED) {/*{{{*/
        /* Powers of 2 are treated separately */
        /* Don't sieve powers of 2 again that were pattern-sieved */
        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        WHERE_AM_I_UPDATE(w, p, p);

        if (ssp.is_pattern_sieved())
            return;

        const unsigned char logp = ssp.logp;
        unsigned char *S_ptr = S;

        int pos = super::first_position_power_of_two(ssp);
        if (pos == super::NO_POSITION)
            return;

        /* In sublattice coordinates this power of two hits with stride
         * q = p/gcd(p,sublat.m), not p. With an even modulus the parity of
         * the real j is constant, so every sublattice row is sieved and
         * consecutive rows are one apart rather than two; the position
         * then advances by r per row instead of 2r. See
         * sublat_reduced_modulus(). */
        const fbprime_t q = super::sublat_reduced_modulus(p);
        const unsigned int dj = super::has_even_sublatm() ? 1 : 2;

        unsigned int j = j0;
        /* Our encoding is that when the first row is even, the
         * pos we keep in the ssdpos array deliberately includes an
         * additional quantity I to mean ``next row''. Compensate
         * that, move on to an odd row, and proceed.
         *
         * With several bucket regions in a line, we adjust by only
         * (i1-i0).  Note that several bucket regions in a line also
         * means only one one line per bucket... So here we'll end up
         * subtracting I in several independent chunks, while I know
         * in advance that I could have done that subtraction all in
         * one go. But really, who cares, at the end of the day.
         */

        if ((j*super::sublat.m + super::sublat.j0) % 2 == 0) {
            /* cannot happen for an even modulus: first_position_power_of_two
             * returned NO_POSITION in that case */
            ASSERT(!super::has_even_sublatm());
            ASSERT(pos >= F());
            pos -= F(); S_ptr += F();
            j++;
        }
        if (j < j1) pos &= (q-1);
        for( ; j < j1 ; j += dj) {
            for (int i = pos; i < F(); i += q) {
                WHERE_AM_I_UPDATE(w, x, ((size_t) (j-j0) << logI) + i);
                sieve_increase (S_ptr + i, logp, w);
            }
            pos = (pos + r * dj) & (q - 1);
            for (unsigned int k = 0; k < dj; k++)
                S_ptr += I();
        }
#if 0
        /* see above. Because we do j+=2, we have either j==j1 or
         * j==j1+1 */
        if (j > j1) pos += I;
        ssdpos[index] = pos;
#endif
    }/*}}}*/
    template<typename even_code, typename odd_code, int bits_off>
    bool handle_nice_prime(ssp_simple_t const & ssp, spos_t pos, where_am_I & w) {/*{{{*/
        const fbprime_t p = ssp.get_p();
        WHERE_AM_I_UPDATE(w, p, p);
        if (bits_off != 0 && (p >> (super::min_logI_logB + 1 - bits_off))) {
            /* time to move on to the next bit size; */
            return false;
        }

        const fbprime_t r = ssp.get_r();
        WHERE_AM_I_UPDATE(w, r, r);
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        unsigned char * S1 = S + F();

        /* we sieve over the area [S0..S0+(i1-i0)] (F() is (i1-i0)),
         * which may actually be just a fragment of a line. After
         * that, if (i1-i0) is different from I, we'll break anyway.
         * So whether we add I or (i1-i0) to S0 does not matter much.
         */

        bool even = super::row0_needs_parity_skip();
        const bool alternates = super::parity_skip_alternates();

        for(unsigned int j = j0; j < j1; j++) {
            WHERE_AM_I_UPDATE(w, j, j - j0);
            if (skip_line_jj0 && sublat.j0 == 0 && j0 == 0 && j == 0) {
                /* A nice prime p hits in line j=0 only in locations
                   where p|i, so no need to sieve those. */
            } else if (even) {
                /* for j even, we sieve only odd pi, so step = 2p. */
                {
                    spos_t xpos = ((super::sublat.i0 + pos) & 1) ? pos : (pos+p);
                    even_code()(S0, S1, S0 - S, xpos, p+p, logp, w);
                }
            } else {
                odd_code()(S0, S1, S0 - S, pos, p, logp, w);
            }
            S0 += I();
            S1 += I();
            pos += r; if (pos >= (spos_t) p) pos -= p;
            even ^= alternates;
        }
        return true;
    }/*}}}*/

    /* This function is responsible for small-sieving primes of a
     * specific bit size */
    template<typename even_code, typename odd_code, int bits_off>
        void handle_nice_primes(where_am_I & w MAYBE_UNUSED) /* {{{ */
        {
            /* here, we can sieve for primes p < 2 * F() / 2^bits_off,
             * (where F() is i1-i0 = 2^min(logI, logB)).
             *
             * meaning that the number of hits in a line is at least
             * floor(F() / p) = 2^(bits_off-1)
             *
             * Furthermore, if p >= 2 * F() / 2^(bits_off+1), we can also
             * say that the number of hits is at most 2^bits_off
             */

#ifdef HAVE_SSE41
            const bool row0_even = super::row0_needs_parity_skip();
            const bool alternates = super::parity_skip_alternates();

            for( ; index + 3 < sorted_limit ; index+=4) {
                /* find 4 index values with no special prime */
                ssp_simple_t const & ssp0(primes[index]);
                spos_t pos0 = (positions[index]);
                ssp_simple_t const & ssp1(primes[index+1]);
                spos_t pos1 = (positions[index+1]);
                ssp_simple_t const & ssp2(primes[index+2]);
                spos_t pos2 = (positions[index+2]);
                ssp_simple_t const & ssp3(primes[index+3]);
                spos_t pos3 = (positions[index+3]);


                const fbprime_t p0 = ssp0.get_p();
                const fbprime_t p1 = ssp1.get_p();
                const fbprime_t p2 = ssp2.get_p();
                const fbprime_t p3 = ssp3.get_p();

                if (bits_off != 0 && (p0 >> (super::min_logI_logB + 1 - bits_off))) {
                    /* time to move on to the next bit size; */
                    return;
                }
                if (bits_off != 0 && (p3 >> (super::min_logI_logB + 1 - bits_off))) {
                    break;
                }

                static_assert(std::is_same<spos_t, int32_t>::value, "spos_t must be an int32_t for the sse2 code to be valid");

                __m128i p = _mm_setr_epi32(p0, p1, p2, p3);
                __m128i r = _mm_setr_epi32(ssp0.get_r(), ssp1.get_r(), ssp2.get_r(), ssp3.get_r());
                __m128i pos = _mm_setr_epi32(pos0, pos1, pos2, pos3);
                const unsigned char logp0 = ssp0.logp;
                const unsigned char logp1 = ssp1.logp;
                const unsigned char logp2 = ssp2.logp;
                const unsigned char logp3 = ssp3.logp;

                unsigned char * S0 = S;
                unsigned char * S1 = S + F();

                /* we sieve over the area [S0..S0+(i1-i0)] (F() is (i1-i0)),
                 * which may actually be just a fragment of a line. After
                 * that, if (i1-i0) is different from I, we'll break anyway.
                 * So whether we add I or (i1-i0) to S0 does not matter much.
                 */
                __m128i ones = _mm_set1_epi32(1);
                __m128i sublati0 = _mm_set1_epi32(super::sublat.i0);
                bool even = row0_even;
                for(unsigned int j = j0 ; j < j1; j++) {
                    WHERE_AM_I_UPDATE(w, j, j - j0);
                    /* the if() branch here that the compiler and/or the
                     * branch predictor are smart enough to make its code
                     * reduce to almost zero */
                    if (skip_line_jj0 && sublat.j0 == 0 && j0 == 0 && j == 0) {
                        /* A nice prime p hits in line j=0 only in locations
                           where p|i, so no need to sieve those. */
                    } else if (even) {
                        /* for j even, we sieve only odd pi, so step = 2p. */
                        __m128i xpos =
                            _mm_add_epi32(pos,
                                    _mm_and_si128(p,
                                        _mm_cmplt_epi32(
                                            _mm_and_si128(
                                                _mm_add_epi32(sublati0, pos),
                                                ones),
                                            ones)
                                        )
                                    )
                            ;
                        WHERE_AM_I_UPDATE(w, p, p0);
                        even_code()(S0, S1, S0 - S, _mm_extract_epi32(xpos, 0), p0+p0, logp0, w);
                        WHERE_AM_I_UPDATE(w, p, p1);
                        even_code()(S0, S1, S0 - S, _mm_extract_epi32(xpos, 1), p1+p1, logp1, w);
                        WHERE_AM_I_UPDATE(w, p, p2);
                        even_code()(S0, S1, S0 - S, _mm_extract_epi32(xpos, 2), p2+p2, logp2, w);
                        WHERE_AM_I_UPDATE(w, p, p3);
                        even_code()(S0, S1, S0 - S, _mm_extract_epi32(xpos, 3), p3+p3, logp3, w);
                    } else {
                        WHERE_AM_I_UPDATE(w, p, p0);
                        odd_code()(S0, S1, S0 - S, _mm_extract_epi32(pos, 0), p0, logp0, w);
                        WHERE_AM_I_UPDATE(w, p, p1);
                        odd_code()(S0, S1, S0 - S, _mm_extract_epi32(pos, 1), p1, logp1, w);
                        WHERE_AM_I_UPDATE(w, p, p2);
                        odd_code()(S0, S1, S0 - S, _mm_extract_epi32(pos, 2), p2, logp2, w);
                        WHERE_AM_I_UPDATE(w, p, p3);
                        odd_code()(S0, S1, S0 - S, _mm_extract_epi32(pos, 3), p3, logp3, w);
                    }
                    S0 += I();
                    S1 += I();
                    pos = _mm_add_epi32(pos, r);
                    pos = _mm_sub_epi32(pos, _mm_andnot_si128(_mm_cmplt_epi32(pos, p), p));
                    even ^= alternates;
                }
            }
#endif

            for( ; index < sorted_limit ; index++) {
                ssp_simple_t const & ssp(primes[index]);
                spos_t pos = positions[index];
                if (!handle_nice_prime<even_code, odd_code, bits_off>(ssp, pos, w))
                    /* time to move on to the next bit size; */
                    return;
            }

        }/*}}}*/
    private:
    /* {{{ template machinery to make one single function out of several */
    /* we'll now craft all the specific handle_nice_primes_meta_loop
     * functions into one big function. Because we're playing tricks with
     * types and lists of types and such, we need to work with partial
     * specializations at the class level, which is admittedly messy. */
    template<typename T, int max_bits_off = INT_MAX>
        struct handle_nice_primes_meta_loop
        {
            void operator()(small_sieve & SS, where_am_I &) {
                /* default, should be at end of list. We require that we
                 * are done processing, at this point. */
                ASSERT_ALWAYS(SS.finished_sorted_prefix());
            }
        };

    /* optimization: do not split into pieces when we have several times
     * the same code anyway. */
    template<typename E0, typename O0, int b0, int b1, typename T, int bn>
        struct handle_nice_primes_meta_loop<choice_list_car<E0,O0,b0,
                     choice_list_car<E0,O0,b1,
                    T>>, bn>
        {
            static_assert(b0 > b1, "choice list is in wrong order");
            void operator()(small_sieve & SS, where_am_I & w) {
                /*
                SS.handle_nice_primes<E0, O0, b1>(w);
                handle_nice_primes_meta_loop<T>()(SS, w);
                */
                handle_nice_primes_meta_loop<choice_list_car<E0,O0,b1, T>, bn>()(SS, w);
            }
        };
    template<typename E0, typename O0, int b0, int bn>
        struct is_compatible_for_range {
            static_assert(bn > b0, "choice list is in wrong order");
            static const int value =
                E0::template is_compatible<bn-1>::value && 
                O0::template is_compatible<bn>::value && 
                is_compatible_for_range<E0, O0, b0, bn-1>::value;
        };
    template<typename E0, typename O0, int b0>
        struct is_compatible_for_range<E0, O0, b0, b0> {
            static const int value =
                E0::template is_compatible<b0-1>::value && 
                O0::template is_compatible<b0>::value;
        };
    template<typename E0, typename O0, int b0>
        struct is_compatible_for_range<E0, O0, b0, INT_MAX> {
            static const int value =
                E0::template is_compatible<INT_MAX>::value && 
                O0::template is_compatible<INT_MAX>::value &&
                is_compatible_for_range<E0, O0, b0, b0 + 5>::value;
        };

    template<typename E0, typename O0, int b0, typename T, int bn>
        struct handle_nice_primes_meta_loop<choice_list_car<E0,O0,b0,T>, bn>
        {
            static_assert(is_compatible_for_range<E0, O0, b0, bn>::value, "Cannot use these two code fragments for primes of the current size");
            void operator()(small_sieve & SS, where_am_I & w) {
                SS.handle_nice_primes<E0, O0, b0>(w);
                handle_nice_primes_meta_loop<T, b0-1>()(SS, w);
            }
        };
    /* }}} */

    public:
    void do_pattern_sieve(where_am_I &);

    void normal_sieve(where_am_I & w) {
        for(size_t s : sorted_subranges) {
            sorted_limit = s;
            /* This function will eventually call handle_nice_primes on
             * sub-ranges of the set of small primes. */
            typedef small_sieve_best_code_choices::type choices;
            small_sieve::handle_nice_primes_meta_loop<choices>()(*this, w);
        }

        /* This is for the tail of the list. We typically have powers,
         * here. These are ordinary, nice, simple prime powers, but the
         * only catch is that these don't get resieved (because the prime
         * itself was already divided out, either via trial division or
         * earlier resieving. By handling them here, we benefit from the
         * ssdpos table. */
        for( ; index < primes.size() ; index++) {
            ssp_simple_t const & ssp(primes[index]);
            spos_t pos = positions[index];
            WHERE_AM_I_UPDATE(w, p, ssp.get_p());
            typedef default_smallsieve_inner_loop even_code;
            typedef default_smallsieve_inner_loop odd_code;
            handle_nice_prime<even_code, odd_code, 0>(ssp, pos, w);
        }
    }

    void exceptional_sieve(where_am_I & w) {
        /* a priori we'll never have "nice" primes here, but we're not
         * forced to rule it out completely, given that we have the code
         * available at hand. The only glitch is that we're storing the
         * start positions *only* for the primes in the ssps array. */

        for(auto const & ssp : not_nice_primes) {
            if (ssp.is_sublat_rowwise()) {
                handle_sublat_rowwise(ssp, w);
            } else if (ssp.is_pattern_sieved()) {
                /* This ssp is pattern-sieved, nothing to do here */
            } else if (ssp.is_proj()) {
                handle_projective_prime(ssp, w);
            } else if (ssp.is_pow2()) {
                handle_power_of_2(ssp, w);
            } else {
                /* I don't think we can end up here.  */
                ASSERT_ALWAYS(0);
            }
        }
    }

};/*}}}*/

/*}}}*/


#endif	/* CADO_LAS_SMALLSIEVE_GLUE_HPP */
