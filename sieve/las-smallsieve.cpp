#include "cado.h"
#include <pthread.h>
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#include "las-config.h"
#include "las-smallsieve.hpp"
#include "las-debug.hpp"
#include "las-qlattice.hpp"
#include "misc.h"
#include "portability.h"
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
    fprintf(f, "# p = %" FBPRIME_FORMAT ", r = %" FBROOT_FORMAT ", offset = %" FBPRIME_FORMAT ", logp = %hhu",
        this->p, this->r, this->offset, this->logp);
}

void ssp_t::print(FILE *f) const
{
    fprintf(f, "# p = %" FBPRIME_FORMAT ", r = %" FBROOT_FORMAT ", offset = %" FBPRIME_FORMAT ", logp = %hhu",
        this->p, this->r, this->offset, this->logp);
    if (this->is_pow2()) fprintf(f, " (power of 2)");
    if (this->is_proj()) fprintf(f, " (projective root)");
    if (this->is_discarded_proj()) fprintf(f, "(discarded) because of projective root");
    if (this->is_discarded_sublat()) fprintf(f, "(discarded because not compatible with sub lattices)");
}

int
small_sieve_dump(FILE *f, const char *header, va_list va)
{
    const small_sieve_data_t * ssd = va_arg(va, const small_sieve_data_t *);

    fprintf(f, "%s", header);
    for (int i = 0; i < ssd->nb_ssps; i++) {
        ssd->ssps[i].print(f);
        fprintf(f, "\n");
    }
    for (int i = 0; i < ssd->nb_ssp; i++) {
        ssd->ssp[i].print(f);
        fprintf(f, "\n");
    }
    return 1;
}

static void small_sieve_print_contents(const char * prefix, small_sieve_data_t * ssd)
{
    int nice=ssd->nb_ssps;
    int nproj=0;
    int npow2=0;
    int ndiscard=0;
    for(int i = 0; i < ssd->nb_ssp; i++) {
        const ssp_t &ssp = ssd->ssp[i];
        nproj += ssp.is_proj();
        npow2 += ssp.is_pow2();
        ndiscard += ssp.is_discarded();
        nice += ssp.is_nice();
    }

    verbose_output_start_batch();
    verbose_output_print(0, 2, "# %s: %d nice primes", prefix, nice);
    /* Primes may be both even and projective... */
    if (npow2) verbose_output_print(0, 2, ", %d powers of 2", npow2);
    if (nproj) verbose_output_print(0, 2, ", and %d projective primes", nproj);
    verbose_output_print(0, 2, ".");
    if (ndiscard) verbose_output_print(0, 2, " %d discarded.", ndiscard);
    verbose_output_print(0, 2, "\n");
    /* With -v -v -v, dump all the small sieve data */
    verbose_output_vfprint (0, 4, small_sieve_dump, "# Dump of small sieve data:\n", ssd);
    verbose_output_end_batch();
}


void small_sieve_info(const char * what, int side, small_sieve_data_t * r)
{
    char * tmp;
    int rc = asprintf(&tmp, "%s(side %d)", what, side);
    ASSERT_ALWAYS(rc >= 0);
    small_sieve_print_contents(tmp, r);
    free(tmp);
}

/* }}} */

/* {{{ Sieve initialization / clearing : first the easy ones */
void small_sieve_clear(small_sieve_data_t * ssd)
{
    free(ssd->ssps); ssd->ssps = NULL;
    free(ssd->ssp); ssd->ssp = NULL;
}

void small_sieve_extract_interval(small_sieve_data_t * r, small_sieve_data_t * s, int bounds[2])
{
    memset(r, 0, sizeof(small_sieve_data_t));
    r->nb_ssp = bounds[1] - bounds[0];
    r->ssp = (ssp_t *) malloc (r->nb_ssp * sizeof (ssp_t));
    FATAL_ERROR_CHECK(r->nb_ssp > 0 && r->ssp == NULL, "malloc failed");

    memcpy(r->ssp, s->ssp + bounds[0], r->nb_ssp * sizeof (ssp_t));
}

/* }}} */

/* {{{ Sieve initialization: now the real stuff */

void ssp_t::init_proj(fbprime_t p, fbprime_t r, unsigned char _logp, unsigned int skip MAYBE_UNUSED)/*{{{*/
{
    set_proj();
    logp = _logp;
    unsigned int v = r; /* have consistent notations */
    unsigned int g = gcd_ul(p, v);
    fbprime_t q = p / g;
    set_q(q);
    set_g(g);
    if (q == 1) {
        ASSERT(r == 0);
        set_U(0);
    } else {
        int rc;
        uint32_t U = v / g; /* coprime to q */
        rc = invmod_32(&U, q);
        ASSERT_ALWAYS(rc != 0);
        set_U(U);
    }
}/*}}}*/

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// ssdpos points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION
// It could actually be larger than 32 bits when I > 16.

void small_sieve_init(small_sieve_data_t *ssd, unsigned int interleaving,
                      const std::vector<fb_general_entry> *fb,
                      sieve_info const & si, const int side)
{
    const unsigned int thresh = si.conf.bucket_thresh;
    const int verbose = 0;
    where_am_I w MAYBE_UNUSED;

    size_t size = si.sides[side].fb_parts_x->rest[1];

    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->ssp = (ssp_t *) malloc(size * sizeof(ssp_t));
    FATAL_ERROR_CHECK(ssd->ssp == NULL, "malloc failed");
    memset(ssd->ssp, 0, size * sizeof(ssp_t));

    ssd->ssps = (ssp_simple_t *) malloc(size * sizeof(ssp_simple_t));
    FATAL_ERROR_CHECK(ssd->ssps == NULL, "malloc failed");
    memset(ssd->ssps, 0, size * sizeof(ssp_simple_t));

    // Do another pass on fb and projective primes, to fill in the data
    // while we have any regular primes or projective primes < thresh left
    ssp_t * tail = ssd->ssp;
    ssp_simple_t * tails MAYBE_UNUSED = ssd->ssps;

    // The processing of bucket region by nb_threads is interleaved.
    // It means that the positions for the small sieve must jump
    // over the (nb_threads - 1) regions after each region.
    // For typical primes, this jump is easily precomputed and goes into
    // the ssp struct.
    
    // If we are doing sublattices modulo m, then we jump virutally m
    // times faster.
    unsigned int sublatm = si.conf.sublat.m;
    const unsigned int skiprows = ((interleaving-1) << LOG_BUCKET_REGION) >> si.conf.logI_adjusted;

    for (std::vector<fb_general_entry>::const_iterator iter = fb->begin() ; iter != fb->end() ; iter++) {
        /* p=pp^k, the prime or prime power in this entry, and pp is prime */
        const fbprime_t p = iter->q, pp = iter->p;
        WHERE_AM_I_UPDATE(w, p, p);

        if (p > thresh) {
            continue;
        }

        for (int nr = 0; nr < iter->nr_roots; nr++) {
            const fb_general_root &root = iter->roots[nr];
            /* Convert into old format for projective roots by adding p if projective.
               FIXME, this sucks. */
            const fbroot_t r = root.r + (root.proj ? p : 0);
            /* skip q: mark it in place of p to be recognized quickly.
               TODO: maybe use the SSP_DISCARD mechanism ?
               */
            if (!si.doing.is_coprime_to(pp)) {
                tail->set_p(mpz_get_ui(si.doing.p));
                tail++;
                continue;
            }

            if (sublatm) {
                // In sublat mode, disable pattern sieving and primes
                // dividing m. (pp is the prime, here)
                if (pp == 2 || pp == 3 || (sublatm % pp) == 0) {
                    tail->set_discarded_sublat();
                }
            }

            if ((p & 1)==0) tail->set_pow2();

            const unsigned char logp = fb_log_delta (pp, root.exp, root.oldexp,
                si.sides[side].lognorms->scale);
            
            WHERE_AM_I_UPDATE(w, r, r);
            const fbroot_t r_q = fb_root_in_qlattice(p, r, iter->invq, si.qbasis);
            /* If this root is somehow interesting (projective in (a,b) or
               in (i,j) plane), print a message */
            if (verbose && (r > p || r_q >= p))
                verbose_output_print(0, 1, "# small_sieve_init: side %d, prime %"
                        FBPRIME_FORMAT " root %s%" FBROOT_FORMAT " -> %s%" 
                        FBROOT_FORMAT "\n", side, p, 
                        r >= p ? "1/" : "", r % p,
                        r_q >= p ? "1/" : "", r_q % p);

            /* Handle projective roots */
            if (r_q >= p && !tail->is_discarded_sublat()) {
                /* Compute the init data in any case, since the gcd
                 * dominates (and anyway we won't be doing this very
                 * often). */
                tail->init_proj(p, r_q - p, logp, skiprows);
                /* If g exceeds J, then the only reached locations in the
                 * sieving area will be on line (j=0), thus (1,0) only since
                 * the other are equivalent.
                 */
                if (tail->get_g() >= si.J) {
                    if (verbose) {
                        verbose_output_print(0, 1,
                                "# small_sieve_init: not adding projective prime"
                                " (1:%" FBROOT_FORMAT ") mod %" FBPRIME_FORMAT ")"
                                " to small sieve  because g=%d >= si.J = %d\n",
                                r_q-p, p, tail->get_g(), si.J);
                    }
                    tail->set_discarded();
                }
            } else if (!tail->is_discarded()) {
                tail->init(p, r_q, logp, skiprows);
            }

            tail++;
        }
    }
    ssd->nb_ssp = size;
}
/* }}} */

/* {{{ Creation of the ssdpos tables */
int64_t * small_sieve_copy_start(int64_t * base, int bounds[2])
{
    int64_t *res = (int64_t *)malloc((bounds[1] - bounds[0]) * sizeof(int64_t));
    memcpy(res, base + bounds[0], (bounds[1] - bounds[0]) * sizeof(int64_t));
    return res;
}

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


/* So many things are used in common for many small sieve routines that
 * it makes sense to gather them in a common object */
struct small_sieve_context {
    int logI;
    int N;
    unsigned int log_lines_per_region;
    unsigned int log_regions_per_line;
    unsigned int regions_per_line;
    unsigned int region_rank_in_line;
    bool last_region_in_line;
    unsigned int j0;
    unsigned int j1;
    int I;
    int i0;
    int i1;
    /* those are (1,0,0) in the standard case */
    int sublatm;
    unsigned int sublati0;
    unsigned int sublatj0;
    int row0_is_oddj;
    bool has_haxis;
    bool has_vaxis;
    bool has_origin;
    small_sieve_context(int logI, int N, sublat_t const & sublat)
        : logI(logI), N(N)
    {
        log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);
        log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);
        regions_per_line = 1 << log_regions_per_line;           
        region_rank_in_line = N & (regions_per_line - 1);       
        last_region_in_line = region_rank_in_line == (regions_per_line - 1); 
        j0 = (N >> log_regions_per_line) << log_lines_per_region;    
        j1 = j0 + (1 << log_lines_per_region);    
        I = 1 << logI;
        i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          
        i1 = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     
        /* those are (1,0,0) in the standard case */
        sublatm = sublat.m ? sublat.m : 1; 
        sublati0 = sublat.i0;       
        sublatj0 = sublat.j0;       
        row0_is_oddj = (j0*sublatm + sublatj0) & 1;
        has_haxis = !j0;                                               
        has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);   
        has_origin = has_haxis && has_vaxis;              
    }

    int64_t first_position_ordinary_prime(const ssp_t * ssp, unsigned int dj = 0)
    {
        /* equation here: i-r*j = 0 mod p */

        /* Expanded wrt first row j0 and sublats: smallest pos0 s.t.:
         *  sublatm*(i0+pos0)+sublati0 - r*(j0*sublatm+sublatj0) = -k*p
         *  (p is coprime to sublatm).
         *
         * In non-sublat mode, this means:  (i0+pos0)-r*j0 = 0 mod p
         * so that pos0 = (r*j0 - i0) mod p where i0 is traditionally
         * negative (that helps for sign stuff, of course).
         *
         * So S = sublatm*pos0 must be such that:
         *
         * S = r*(j0*sublatm+sublatj0) - (sublatm*i0+sublati0) + k*p
         *
         * pos0 being an integer, of course. So there's only one
         * congruence class for k mod sublatm such that the above
         * happens. We start by computing it.
         *
         *  k:=((sublati0 - r * sublatj0) * invp_mod_m) mod sublatm;
         *
         * then we deduce pos0 as:
         *
         *  pos0:=j0 * r - i0 + (((r*sublatj0-sublati0)+k*p) div sublatm)
         *
         * and finally we take pos0 mod p.
         */
        int64_t x = (int64_t)(j0 + dj) * (int64_t)ssp->get_r() - i0;
        if (sublatm > 1) {
            ASSERT(ssp->get_p() % sublatm);
            /* alternative code. not clear it's better.
               ASSERT(sublatm == 2 || sublatm == 3 || sublatm == 6);
               uint64_t invp_mod_m = p;
               int k = ((sublati0 - r * sublatj0) * invp_mod_m) % sublatm;
               x += ((r * sublatj0 - sublati0) + k * p) / sublatm;
               */
            int64_t y = ssp->get_r() * sublatj0;
            for( ; y % sublatm != sublati0 ; y += ssp->get_p());
            x += (y - sublati0) / sublatm;
        }
        x= x % (int64_t)ssp->get_p();
        /* As long as i0 <= 0, which holds in the normal case where
         * logI <= LOG_BUCKET_REGION, we can be sure that x >= 0, so
         * we have no issue with the sign.  However, when i0 is a
         * positive number, some extra care is needed.
         */
        if (x < 0) x += ssp->get_p();
        return x;
    }

    /* This return value is typically logI bits larger than for ordinary
     * primes, so it makes sense to return it as a 64-bit integer.
     */
    int64_t first_position_projective_prime(const ssp_t * ssp)
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
         *      ii = i0 * sublatm + sublati0
         *      jj = j0 * sublatm + sublatj0
         */

        uint64_t jj = j0*sublatm + sublatj0;

        /* First question is to find the next multiple of g above jj.
         * This is done as follows */

        /* next multiple of g above j0 */
        uint64_t jjmod = jj % ssp->get_g();
        if (jjmod) jj += ssp->get_g();
        jj -= jjmod;

        /* All of the fragment above trivializes to a no-op in the
         * non-projective case (which has g=1)
         */

        /* Now we'd like to avoid row number 0 (so jj == 0). */
        /* The sieving code may do some special stuff to fill the
         * position (1,0). At least at some point it did */
        if (jj == 0) jj += ssp->get_g();

        // In sublat mode, we also need jj congruent to sublatj0 mod m.
        // XXX A very nasty situation: when ssp->g and sublatm are
        // not coprime, we may very well have an infinite loop here
        // (say we want to sieve only even lines and sublatj0=1,
        // sublatm=1. For the moment, we do SSP_DISCARD_SUBLAT
        // whenever p and sublatm have a common divisor. If we want
        // to be finer grain, we need to selectively discard some
        // small sieved primes depending on the sublattice we're
        // considering (or but ssdpos to ULONG_MAX ?).
        if (sublatm > 1) {
            for( ; jj % sublatm != sublatj0 ; jj += ssp->get_g());
        }
        // Find the corresponding i
        uint64_t ii = uint64_t(jj/ssp->get_g())*uint64_t(ssp->get_U());

        if (sublatm > 1) {
            for( ; ii % sublatm != sublati0 ; ii += ssp->get_q());
        }
        // In the sublat mode, switch back to reduced convention
        // (exact divisions)
        if (sublatm > 1) {
            jj = (jj-sublatj0) / sublatm;
            ii = (ii-sublati0) / sublatm;
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
        int i0ref = (j0 == jj) ? i0 : (-I/2);
        int64_t x = (ii-i0ref) % (uint64_t)ssp->get_q();
        if (x < 0) x += ssp->get_q();
        if (jj > j0) {
            x -= region_rank_in_line << LOG_BUCKET_REGION;
            x += (((uint64_t)(jj - j0))<<logI);
        }
        return x;
    }

    int64_t first_position_power_of_two(const ssp_t * ssp, unsigned int dj = 0)
    {
        /* equation here: i-r*j = 0 mod p, p a power of 2. */
        /* only difference with ordinary case is that we want to
         * sieve only odd lines.
         */
        // FIXME: this is probably broken with sublattices (well,
        // let's give it a try).
        /* For powers of 2, we sieve only odd lines (*) and 
         * ssdpos needs to point at line j=1. We assume
         * that in this case (si.I/2) % p == 0
         * (*) for lines with j even, we have a root mod the prime
         * power for i-j*r multiple of our power of 2, which means
         * i even too. Thus a useless report.
         */
        unsigned int j = j0 + dj;
        uint64_t jj = j*sublatm + sublatj0;
        // uint64_t ii = i0*sublatm + sublati0;
        /* next odd line */
        jj |= 1;
        int i0ref = (j == jj) ? i0 : (-I/2);
        int64_t x = (int64_t)jj * (int64_t)ssp->get_r();
        if (sublatm > 1) {
            for( ; x % sublatm != sublati0 ; x += ssp->get_p());
            x = (x - sublati0) / sublatm;
        }
        x = (x - i0ref) & (ssp->get_p() - 1);
        if (x < 0) x += ssp->get_p();
        /* our target is position x in the bucket region which starts
         * at coordinates (i0ref, jj). How far is that from us ?
         */
        if (jj > j) {
            x -= region_rank_in_line << LOG_BUCKET_REGION;
            x += (((uint64_t)(jj - j0))<<logI);
            /* For the case of several bucket regions per line, it's
             * clear that this position will be outside the current
             * bucket region. Note that some special care is needed
             * when we want to incrementally update the position, but we
             * rarely (if ever) want to do that anyway.
             */
        }
        return x;
    }
};

/* This is copied from LOGNORM_FILL_COMMON_DEFS in las-norms.cpp ; from
 * logI, N, and LOG_BUCKET_REGION, define the integers i0, i1, j0, j1,
 * and I.
 */
#define SMALLSIEVE_COMMON_DEFS()                                         \
    const unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);\
    const unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);\
    const unsigned int regions_per_line = 1 << log_regions_per_line;           \
    const unsigned int region_rank_in_line = N & (regions_per_line - 1);       \
    const bool last_region_in_line MAYBE_UNUSED = region_rank_in_line == (regions_per_line - 1); \
    const unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    \
    const unsigned int j1 MAYBE_UNUSED = j0 + (1 << log_lines_per_region);    \
    const int I = 1 << logI;                                            \
    const int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          \
    const int i1 MAYBE_UNUSED = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     \
    /* those are (1,0,0) in the standard case */                        \
    const int sublatm MAYBE_UNUSED = si.conf.sublat.m ? si.conf.sublat.m : 1; \
    const unsigned int sublati0 MAYBE_UNUSED = si.conf.sublat.i0;       \
    const unsigned int sublatj0 MAYBE_UNUSED = si.conf.sublat.j0;       \
    const int row0_is_oddj MAYBE_UNUSED = (j0*sublatm + sublatj0) & 1;  \
    bool has_haxis = !j0;                                               \
    bool has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);   \
    bool has_origin MAYBE_UNUSED = has_haxis && has_vaxis;              \
    do {} while (0)

/* Only compute the initial ssdpos fields. */
int64_t *small_sieve_start(small_sieve_data_t *ssd,
        unsigned int first_region_index,
        sieve_info const & si)
{
    small_sieve_context C(si.conf.logI_adjusted, first_region_index, si.conf.sublat);
    int64_t * ssdpos = (int64_t *) malloc(ssd->nb_ssp * sizeof(int64_t));

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
    for(int index = 0 ; index < ssd->nb_ssp ; index++) {
        /* generic case only. For all other cases (which are rare enough
         * -- typically at most 20, counting powers of two and such), we
         *  compute the starting point from within
         *  sieve_small_bucket_region for each bucket region.
         */
        const ssp_t & ssp = ssd->ssp[index];
        if (ssp.is_discarded()) {
            /* Do nothing for discarded entries */
        } else if (ssp.is_nice()) {
            ssdpos[index] = C.first_position_ordinary_prime(&ssp);
        } else if (ssp.is_proj()) {
            /* This also handles powers of 2 with projective root */
            ssdpos[index] = C.first_position_projective_prime(&ssp);
        } else if (ssp.is_pow2()) {
            /* Powers of 2 with affine root */
            ssdpos[index] = C.first_position_power_of_two(&ssp);
        } else {
            abort(); /* How did we get here? */
        }
    }
    return ssdpos;
}
/* }}} */

/* This adds extra logging for pattern sieving. Very slow.
 */
#define xxxUGLY_DEBUGGING

/* #define exactly one of these */
// #define xxxSMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
// #define xxxSMALLSIEVE_CRITICAL_MANUAL_UNROLL
// #define SMALLSIEVE_CRITICAL_PLAIN

/* on an rsa155 test with I=14, the UGLY_ASSEMBLY version gains
 * practically nothing over the plain C version. The MANUAL_UNROLL
 * version is slower.
 *
 * On an rsa220 test with I=16, the UGLY_ASSEMBLY version gains about 3%
 * on the small sieve part strictly speaking, and then less than 1%
 * overall.
 */
#if 1
#if defined( HAVE_SSE2 ) && !defined( TRACK_CODE_PATH ) /* x86 optimized code */
#define SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
#else
#define SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#endif
#else
#warning "using plain C code for critical part of small sieve"
#define SMALLSIEVE_CRITICAL_PLAIN
#endif


/* {{{ critical part of the small sieve */
#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
      /* Some defines for the critical part of small sieve.
         0. The C code and the asm X86 code have the same algorithm.
         Read first the C code to understand easily the asm code.
         1. If there are less than 12 "T" in the line, the goal is to do
         only one jump (pipe-line breaks) and of course the instructions
         minimal number.
         2. 12 "T" seems the best in the critical loop. Before, gcc tries
         to optimize in a bad way the loop. For gcc generated code, the
         best is here a systematic code (12*2 instructions), like the C code.
         3. The asm X86 optimization uses addb logp,(pi+p_or_2p*[0,1,2])
         & three_p_or_2p = 3 * p_or_2p; the lea (add simulation) & real
         add interlace seems a bit interesting.
         So, the loop is smaller & faster (19 instructions versus 27 for
         gcc best X86 asm).
         4. Of course, the gain between the 2 versions is light, because
         the main problem is the access time of the L0 cache: read + write
         with sieve_increase(pi,logp,w), or *pi += logp in fact.
      */
#define U1                                                              \
        "addb %4,(%1)\n"                /* sieve_increase(pi,logp,w) */ \
        "addb %4,(%1,%3,1)\n"   /* sieve_increase(p_or_2p+pi,logp,w) */ \
        "addb %4,(%1,%3,2)\n" /* sieve_increase(p_or_2p*2+pi,logp,w) */ \
        "lea (%1,%2,1),%1\n"                    /* pi += p_or_2p * 3 */
#define U2                                                              \
        "cmp %5, %1\n"                    /* if (pi >= S1) break; */    \
        "jae 2f\n"                                                      \
        "addb %4,(%1)\n"                /* sieve_increase(pi,logp,w) */ \
        "lea (%1,%3,1),%1\n"                        /* pi += p_or_2p */
#define U do {                                                          \
        unsigned char *pi_end;                                          \
        size_t three_p_or_2p;                                           \
        __asm__ __volatile__ (                                          \
        "lea (%3,%3,2), %2\n"         /* three_p_or_2p = p_or_2p * 3 */ \
        "lea (%1,%2,4), %0\n"            /* pi_end = pi + p_or_2p*12 */ \
        "cmp %5, %0\n"                /* if (pi_end > S1) no loop */    \
        "jbe 0f\n"                                                      \
        "1:\n"                                                          \
        U2 U2 U2 U2 U2 U2 U2 U2 U2 U2 U2                                \
        "cmp %5, %1\n"                                                  \
        "jae 2f\n"                                                      \
        "addb %4,(%1)\n"                                                \
        "jmp 2f\n"                                                      \
        ".balign 8\n 0:\n"                              /* Main loop */ \
        U1 U1           /* sieve_increase(p_or_2p*[0..11]+pi,logp,w) */ \
        U1 U1                                    /* pi += p_or_2p*12 */ \
        "lea (%1,%2,4), %0\n"    /* if (pi+p_or_2p*12 > S1) break */    \
        "cmp %5, %0\n"                                                  \
        "jbe 0b\n"                                                      \
        "jmp 1b\n"                                                      \
        ".balign 8\n 2:\n"                                              \
        : "=&r"(pi_end), "+&r"(pi), "=&r"(three_p_or_2p)                \
        : "r"(p_or_2p), "q"(logp), "r"(S1) : "cc");                     \
        pi = pi_end;                                                    \
      } while (0)
#endif

#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#define T do {                                                          \
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);                              \
    sieve_increase (pi, logp, w); pi += p_or_2p;                        \
} while(0)
#endif

inline size_t sieve_full_line(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
    U;
#endif
#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
    while (UNLIKELY(pi + p_or_2p * 12 <= S1))
    { T; T; T; T; T; T; T; T; T; T; T; T; }
    do {
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
    } while (0);
#endif
#ifdef SMALLSIEVE_CRITICAL_PLAIN
    for ( ; pi <= S1 ; pi += p_or_2p) {
        WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);
        sieve_increase (pi, logp, w);
    }
#endif
    return pi - S1;
}

#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
#undef U1
#undef U2
#undef U
#undef T
#endif
#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#undef T
#endif

/* }}} */

/* {{{ Normal small sieve */

// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, int N,
                               small_sieve_data_t * ssd, int64_t * ssdpos,
                               sieve_info & si, int side,
                               int interleaving MAYBE_UNUSED,
                               where_am_I & w MAYBE_UNUSED)
{
    int logI = si.conf.logI_adjusted;
    SMALLSIEVE_COMMON_DEFS();
    small_sieve_context C(logI, N, si.conf.sublat);

    const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    const int test_divisibility = 0; /* very slow, but nice for debugging */

    /* Pattern-sieve powers of 2 up to 2 * sizeof(long). {{{
     * TODO: use SSE2 (well, cost does not seem to be significant anyway).
     */
    WHERE_AM_I_UPDATE(w, p, 2);

    /* First collect updates for powers of two in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (unsigned int j = j0; j < j1; j++) {
        WHERE_AM_I_UPDATE(w, j, j - j0);
        unsigned long pattern[2] = {0,0};

        /* Prepare the pattern */

        /* Process only powers of 2. For now we do a full pass over the ssp_t
           entries. TODO: split off powers of 2 into their own array. */
        const bool verbose_pattern_2 = false;
        for(int index = 0 ; index < ssd->nb_ssp ; index++) {
            const ssp_t &ssp = ssd->ssp[index];
            if (!ssp.is_pow2())
                continue;
            if (ssp.is_proj()) {
                /* Do nothing. It's a projective power of 2, * but for the
                   moment these are not pattern-sieved. TODO: pattern sieve
                   them in those lines where they hit. */
                continue;
            }
            /* Powers of 2 greater than the pattern size get sieved below */
            if (ssp.get_p() > 2*sizeof(unsigned long))
                continue;
            /* *Affine* powers of two are relevant only for odd lines anyway:
             * indeed, if i-j*r=0 mod 2^k and j even, then i even too, so
             * useless.
             */
            if (j&1) {
                ssp_t * ssp = &(ssd->ssp[index]);
                const fbprime_t p = ssp->get_p();
                if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
                    continue;
                }

                if (verbose_pattern_2) {
                    FILE *out;
                    verbose_output_start_batch();
                    for (size_t out_i=0; (out=verbose_output_get(TRACE_CHANNEL, 0, out_i)) != NULL; out_i++) {
                        fprintf(out, "# Adding ");
                        ssp->print(out);
                        fprintf(out, " to sieving pattern for powers of 2\n");
                    }
                    verbose_output_end_batch();
                }

                int pos = C.first_position_power_of_two(ssp, j - j0);

                if (pos < (i1 - i0)) {
                    // At the start of the bucket region, we are not sure
                    // that pos is reduced mod p when we enter here.
                    // (see comment at the end of small_sieve_skip_stride()).
                    // Hence, we have to do this reduction here:
                    pos &= (p-1);
                    ASSERT (pos < (int) p);
                    ASSERT (j % 2);
                    for (int x = pos; x < (int) pattern2_size; x += p) {
                        if (verbose_pattern_2)
                            verbose_output_print(TRACE_CHANNEL, 0,
                                "# Hits at pattern[%d]\n", x);
                        ((unsigned char *)pattern)[x] += ssp->logp;
                    }
#ifdef UGLY_DEBUGGING
                    for (int x = pos; x < (i1 - i0) ; x+= p) {
                        WHERE_AM_I_UPDATE(w, x, (w.j << logI) + x);
                        sieve_increase(S + x, ssp->logp, w);
                        /* cancel the above action */
                        S[x] += ssp->logp;
                    }
#endif
#if 0
                    /* Skip two lines above, since we sieve only odd lines.
                     * Even lines would correspond to useless reports.
                     */
                    pos = ((pos + 2 * ssp->get_r()) & (p - 1)) + (2 << logI);
#endif
                }
#if 0
                /* In this loop, ssdpos gets updated to the first 
                   index to sieve relative to the start of the next line, 
                   but after all lines of this bucket region are processed, 
                   it will point to the first position to sieve relative
                   to the start of the next bucket region, as required */
                ssdpos[index] = pos - (i1 - i0);
#endif
            }
        }

        /* Apply the pattern */
        if (verbose_pattern_2)
            verbose_output_print(TRACE_CHANNEL, 0,
                "# Sieving pattern for powers of 2 is: %lx, %lx\n",
            pattern[0], pattern[1]);
        if (pattern[0] || pattern[1]) {
            unsigned long *S_ptr = (unsigned long *) (S + ((size_t) (j-j0) << logI));
            const unsigned long *S_end = (unsigned long *)(S + ((size_t) (j-j0) << logI) + (i1-i0));

#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w.N, w.j*I, w.j*I+I)) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I;
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))&1)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_increase_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            ASSERT((S_end - S_ptr) % 4 == 0);

            for ( ; S_ptr < S_end ; S_ptr += 4) {
                S_ptr[0] += pattern[0];
                S_ptr[1] += pattern[1];
                S_ptr[2] += pattern[0];
                S_ptr[3] += pattern[1];
            }
        }
    }
    /* }}} */

    /* {{{ Pattern-sieve 3 */
    /* For the time being, it's really 3. But the only thing we care about is
     * the hit rate to be 1 every 3rd, on lines for which we hit. So the code
     * below could be improved to handle more cases, although presently it
     * does not. */
    WHERE_AM_I_UPDATE(w, p, 3);
    /* First collect updates for powers of three in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (unsigned int j = j0; j < j1; j++)
    {
        WHERE_AM_I_UPDATE(w, j, j - j0);
        unsigned long pattern[3];

        pattern[0] = pattern[1] = pattern[2] = 0UL;

        int * interval = si.sides[side].fb_parts_x->pow3;
        for(int index = interval[0] ; index < interval[1] ; index++) {
            ssp_t * ssp = &(ssd->ssp[index]);
            if (ssp->is_nice()) {
                ASSERT_ALWAYS(ssp->get_p() == 3);
                const fbprime_t p = 3;
                WHERE_AM_I_UPDATE(w, p, p);

                int pos = C.first_position_ordinary_prime(ssp, j - j0);

                ASSERT (pos < (int) p);
                for (unsigned int x = pos; x < 3 * sizeof(unsigned long); x += p)
                    ((unsigned char *)pattern)[x] += ssp->logp;
            } else {
                /* index points to a power of 3, and we have an exceptional
                 * event. Sure it can neither be SSP_END nor SSP_POW2.
                 * It's thus almost surely SSP_PROJ, although we could
                 * conceivably have SSP_DISCARD as well
                 */
                /* We should / could do something, anyway. Given that at
                 * this point, we have only 3 ulongs for the pattern,
                 * we're certain that a projective prime is trivial*/
                /* TODO */
            }
        }

        if (!(j&1)) {
            /* We have an even j. There, we must not sieve even i's either ! */
            for (unsigned int x = 0; x < 3 * sizeof(unsigned long); x += 2)
                ((unsigned char *)pattern)[x] = 0;
        }

        /* We want to test if there is a non-zero entry in the pattern
         * within the first 6 entries (we sieve mod 3 and mod 2).
         * For compatibility with 32-bit machines, we test the first two
         * pattern unsigned longs */
        if (pattern[0] || pattern[1]) {
            unsigned long *S_ptr = (unsigned long *) (S + ((size_t) (j-j0) << logI));
            const unsigned long *S_end = (unsigned long *)(S + ((size_t) (j-j0) << logI) + (i1-i0)) - 2;

#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w.N, w.j*I, w.j*I+I)) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I;
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))%3)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_increase_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            for( ; S_ptr < S_end ; S_ptr += 3) {
                S_ptr[0] += pattern[0];
                S_ptr[1] += pattern[1];
                S_ptr[2] += pattern[2];
            }

            S_end += 2;
            if (S_ptr < S_end)
                *(S_ptr++) += pattern[0];
            if (S_ptr < S_end)
                *(S_ptr) += pattern[1];
        }
    }
    /* }}} */

    // sieve with everyone, since pattern-sieving may miss some of the
    // small primes.

    /* use size_t for k, i and twop to speed up in x86 with LEA asm instructions */
    for(size_t index = 0 ; index < (size_t) ssd->nb_ssp ; index++) {
        ssp_t * ssp = &(ssd->ssp[index]);
        if(ssp->is_nice()) {
            const fbprime_t p = ssp->get_p();
            const fbprime_t r = ssp->get_r();

            if (mpz_cmp_ui(si.qbasis.q, p) == 0) continue;

            /* Don't sieve 3 again as it was pattern-sieved -- unless
             * it's projective, see TODO above. */
            if (p == 3) continue;

            WHERE_AM_I_UPDATE(w, p, p);
            const unsigned char logp = ssp->logp;
            unsigned char *S0 = S;
            int pos = ssdpos[index];
#ifdef TRACE_K
            /* we're keeping track of it with ssdpos, but just in case,
             * make sure we get it right ! */
            ASSERT_ALWAYS(pos == C.first_position_ordinary_prime (ssp));
#endif
            ASSERT(pos < (int) p);
            unsigned int i_compens_sublat = si.conf.sublat.i0 & 1;
            unsigned int j = j0;

            /* we sieve over the area [S0..S0+(i1-i0)], which may
             * actually be just a fragment of a line. After that, if
             * (i1-i0) is different from I, we'll break anyway. So
             * whether we add I or (i1-i0) to S0 does not matter much.
             */
            size_t overrun = 0; /* tame gcc */
            if (j & 1) {
                WHERE_AM_I_UPDATE(w, j, j - j0);
                overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                        pos, p, logp, w);
                S0 += I;
                j++;
                pos += r; if (pos >= (int) p) pos -= p;
            }
            for( ; j < j1 ; ) {
                /* for j even, we sieve only odd pi, so step = 2p. */
                int xpos = ((i_compens_sublat + pos) & 1) ? pos : (pos+p);
                overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                        xpos, p+p, logp, w);
                S0 += I;
                pos += r; if (pos >= (int) p) pos -= p; 
                if (++j >= j1) break;

                /* now j odd again */
                WHERE_AM_I_UPDATE(w, j, j - j0);
                overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                        pos, p, logp, w);
                S0 += I;
                pos += r; if (pos >= (int) p) pos -= p;
                ++j;
            }

            if (logI > LOG_BUCKET_REGION) {
                /* quick notes for incremental adjustment in case I>B (B =
                 * LOG_BUCKET_REGION).
                 *
                 * Let q = 2^(I-B).
                 * Let N = a*q+b, and N'=N+interleaving=a'*q+b' ; N' is the
                 * next bucket region we'll handle.
                 *
                 * Let interleaving = u*q+v
                 *
                 * The row increase is dj = (N' div q) - (N div q) = a'-a
                 * The fragment increase is di = (N' mod q) - (N mod q) = b'-b
                 *
                 * Of course we have -q < b'-b < q
                 *

                 * dj can be written as (N'-N-(b'-b)) div q, which is an
                 * exact division. we rewrite that as:
                 *
                 * dj = u + (v - (b'-b)) div q
                 *
                 * where the division is again exact. Now (v-(b'-b))
                 * satisfies:
                 * -q < v-(b'-b) < 2*q-1
                 *
                 * so that the quotient may only be 0 or 1.
                 *
                 * It is 1 if and only if v >= q + b'-b, which sounds like a
                 * reasonable thing to check.
                 */
                int N1 = N + interleaving;
                int Q = logI - LOG_BUCKET_REGION;
                int dj = (N1>>Q) - j0;
                int di = (N1&((1<<Q)-1)) - (N&((1<<Q)-1));
                /* Note that B_mod p is not reduced. It may be <0, and
                 * may also be >= p if we sieved with 2p because of even j
                 */
                int B_mod_p = overrun - ssdpos[index];
                /* We may avoid some of the cost for the modular
                 * reduction, here:
                 *
                 * dj is either always the same thing, or that same thing
                 * + 1.
                 *
                 * di is within a small interval (albeit a centered one).
                 * 
                 * So it seems feasible to get by with a fixed number of
                 * conditional subtractions.
                 */
                pos = (ssdpos[index] + B_mod_p * di + dj * r) % p;
                if (pos < 0) pos += p;
#ifndef NDEBUG
                small_sieve_context C1(logI, N1, si.conf.sublat);
                ASSERT(pos == C1.first_position_ordinary_prime(ssp));
#endif
                ssdpos[index] = pos;
            } else {
                /* skip stride */
                pos += ssp->get_offset();
                if (pos >= (int) p) pos -= p;
                ssdpos[index] = pos;
            }
        } else {

            // if (ssp->is_discarded_proj()) continue;
            if (ssp->is_proj()) {
                /* This code also covers projective powers of 2 */
                ssp_t * ssp = &(ssd->ssp[index]);
                const fbprime_t q = ssp->get_q();
                const fbprime_t g = ssp->get_g();
                const size_t gI = (size_t)ssp->get_g() << logI;
                const fbprime_t U = ssp->get_U();
                const fbprime_t p MAYBE_UNUSED = g * q;
                WHERE_AM_I_UPDATE(w, p, p);
                const unsigned char logp = ssp->logp;
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
                if (!test_divisibility && ssp->get_q() == 1)
                {
                    /* q = 1, therefore U = 0, and we sieve all entries in lines
                       with g|j, beginning with the line starting at S[ssdpos] */
                    unsigned long logps;
                    int64_t pos = C.first_position_projective_prime(ssp);

                    // The following is for the case where p divides the norm
                    // at the position (i,j) = (1,0).
                    if (UNLIKELY(has_origin && pos == (int64_t) gI)) {
#ifdef TRACE_K
                        if (trace_on_spot_Nx(w.N, (1-i0))) {
                            WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                            unsigned int v = logp;
                            sieve_increase_logging(S + w.x, v, w);
                        }
#endif
                        S[1 - i0] += logp;
                    }
                    // The event SSP_DISCARD might have occurred due to
                    // the first row to sieve being larger than J. The row
                    // number 0 must still be sieved in that case, but once
                    // it's done, we can indeed skip the next part of
                    // sieving.
                    if (ssp->is_discarded_proj())
                        continue;
                    ASSERT (ssp->get_U() == 0);
                    ASSERT (pos % (i1 - i0) == 0);
                    ASSERT ((i1 - i0) % (4 * sizeof (unsigned long)) == 0);
                    for (size_t x = 0; x < sizeof (unsigned long); x++)
                        ((unsigned char *)&logps)[x] = logp;
                    unsigned int j = j0 + (pos >> logI);
                    for( ; j < j1 ; ) {
                        /* our loop is over line fragments that have a hit,
                         * and by the condition q=1 above we'll sieve them
                         * completely */
                        unsigned long *S_ptr = (unsigned long *) (S + pos);
                        unsigned long *S_end = (unsigned long *) (S + pos + i1 - i0);
                        unsigned long logps2 = logps;
                        if (!(j&1)) {
                            /* j is even. We update only odd i-coordinates */
                            /* Use array indexing to avoid endianness issues. */
                            for (size_t x = 0; x < sizeof (unsigned long); x += 2)
                                ((unsigned char *)&logps2)[x] = 0;
                        }
#ifdef TRACE_K
                        if (trace_on_range_Nx(w.N, pos, pos + (i1 - i0))) {
                            WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                            unsigned int x = trace_Nx.x;
                            unsigned int index = x % I;
                            unsigned int v = (((unsigned char *)(&logps2))[index%sizeof(unsigned long)]);
                            sieve_increase_logging(S + w.x, v, w);
                        }
#endif
                        for( ; S_ptr < S_end ; S_ptr += 4) {
                            S_ptr[0] += logps2;
                            S_ptr[1] += logps2;
                            S_ptr[2] += logps2;
                            S_ptr[3] += logps2;
                        }
                        pos += gI;
                        j += g;
                    }
#if 0
                    ssdpos[index] = pos - (1U << LOG_BUCKET_REGION);
#endif
                } else {
                    /* q > 1, more general sieving code. */
                    int64_t pos = C.first_position_projective_prime(ssp);
                    const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
                    unsigned char * S_ptr = S;
                    S_ptr += pos - (pos & (I - 1U));
                    unsigned int j = j0 + (pos >> logI);
                    pos = pos & (I - 1U);
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
                            for ( ; i < (i1 - i0); i += q_or_2q) {
                                WHERE_AM_I_UPDATE(w, x, (S_ptr - S) + i);
                                sieve_increase (S_ptr + i, logp, w);
                            }
                        }

                        pos += U;
                        if (pos >= (int) q) pos -= q;
                    }
#if 0
                    ssdpos[index] = linestart + lineoffset - (1U << LOG_BUCKET_REGION);
#endif
                }
            } else if (ssp->is_pow2()) {
                /* Powers of 2 are treated separately */
                /* Don't sieve powers of 2 again that were pattern-sieved */
                ssp_t * ssp = &(ssd->ssp[index]);
                const fbprime_t p = ssp->get_p();
                const fbprime_t r = ssp->get_r();
                WHERE_AM_I_UPDATE(w, p, p);

                if (p <= pattern2_size)
                    continue;

                const unsigned char logp = ssp->logp;
                unsigned char *S_ptr = S;

                int pos = C.first_position_power_of_two(ssp);

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

                if (j % 2 == 0) {
                    ASSERT(pos >= (i1 - i0));
                    pos -= (i1 - i0); S_ptr += (i1 - i0);
                    j++;
                }
                if (j < j1) pos &= (p-1);
                for( ; j < j1 ; j+= 2) {
                    for (int i = pos; i < (i1 - i0); i += p) {
                        WHERE_AM_I_UPDATE(w, x, ((size_t) (j-j0) << logI) + i);
                        sieve_increase (S_ptr + i, logp, w);
                    }
                    // odd lines only.
                    pos = (pos + (r << 1)) & (p - 1);
                    S_ptr += I;
                    S_ptr += I;
                }
#if 0
                /* see above. Because we do j+=2, we have either j==j1 or
                 * j==j1+1 */
                if (j > j1) pos += I;
                ssdpos[index] = pos;
#endif
            }
        }
    }
}
/* }}} */

/* {{{ resieving. Different interface, since it plays with buckets as well.
 */

/* Sieve small primes (p < I) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd. */
void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t *ssd, int64_t * ssdpos,
        sieve_info const & si, int interleaving MAYBE_UNUSED, where_am_I & w MAYBE_UNUSED)
{
    int logI = si.conf.logI_adjusted;
    SMALLSIEVE_COMMON_DEFS();
    small_sieve_context C(logI, N, si.conf.sublat);

    unsigned char *S_ptr;
    const int resieve_very_verbose = 0, resieve_very_verbose_bad = 0;

    unsigned int i_compens_sublat = sublati0 & 1;

    if (sublatm>1) {
        // Odd/even property of j is the same as for j+2, even with
        // sublat, unless sublat.m is even, which is not handled right
        // now. Same for i.
        ASSERT_ALWAYS((sublatm & 1) == 1);
    }

    for(int index = 0 ; index < ssd->nb_ssp ; index++) {
        ssp_t * ssp = &(ssd->ssp[index]);
        if (ssp->is_pow2())
            continue;
        if (ssp->is_nice()) {
            if (ssp->is_discarded())
                continue;
            const fbprime_t p = ssp->get_p();
            if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
                continue;
            }
            fbprime_t r = ssp->get_r();
            WHERE_AM_I_UPDATE(w, p, p);
            int pos = ssdpos[index];
            S_ptr = S;
            ASSERT(pos < (int) p);
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
             */
            unsigned int q = p&-!row0_is_oddj;
            for (unsigned int j = j0; j < j1; j ++) {
                WHERE_AM_I_UPDATE(w, j, j);
                for (int i = pos + (q& -!((pos+i_compens_sublat)&1)) ; i < (i1-i0); i += p+q) {
                    if (LIKELY(S_ptr[i] == 255)) continue;
                    bucket_update_t<1, primehint_t> prime;
                    unsigned int x = ((size_t) (j-j0) << logI) + i;
                    if (resieve_very_verbose) {
                        verbose_output_print(0, 1, "resieve_small_bucket_region: root %"
                                FBROOT_FORMAT ",%d divides at x = "
                                "%d = %u * %u + %d\n",
                                p, r, x, j, I, i);
                    }
                    prime.p = p;
                    prime.x = x;
                    ASSERT(prime.p >= si.conf.td_thresh);
                    BP->push_update(prime);
                }
                pos += r;
                if (pos >= (int) p) pos -= (int) p;

                S_ptr += I;
                q ^= p;
            }
            pos += ssp->get_offset();
            if (pos >= (int) p) pos -= p;

            ssdpos[index] = pos;
        } else if (ssp->is_proj()) {
            const fbprime_t g = ssp->get_g();
            const uint64_t gI = (uint64_t)g << logI;

            WHERE_AM_I_UPDATE(w, p, g * ssp->get_q());

            /* Test every p-th line, starting at S[ssdpos] */
            int64_t pos = C.first_position_projective_prime(ssp);
            // This block is for the case where p divides at (1,0).
            if (UNLIKELY(has_origin && pos == (int64_t) gI)) {
                bucket_update_t<1, primehint_t> prime;
                prime.p = g;
                prime.x = 1 - i0;
                ASSERT(prime.p >= si.conf.td_thresh);
                BP->push_update(prime);
            }
            // Same as in sieving: we discard after checking for row 0.
            if (ssp->is_discarded_proj())
                continue;

            /* make sure ssdpos points at start of line or region */
            ASSERT(!(pos % (1 << MIN(logI, LOG_BUCKET_REGION))));

            if (resieve_very_verbose_bad) {
                verbose_output_print(0, 1, "# resieving projective prime %" FBPRIME_FORMAT
                        ", i0 = %" PRIu64 "\n", g, i0 + pos);
            }
            if (pos >> LOG_BUCKET_REGION)
                continue;
            unsigned int j = j0 + (pos >> logI);
            for (; j < j1; j += g) {
                unsigned char *S_ptr = S + pos;
                if (!(j & 1)) {
                    /* Even j: test only odd ii-coordinates */
                    for (int ii = 1; ii < (i1 - i0); ii += 2) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose_bad) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region even j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        ASSERT(prime.p >= si.conf.td_thresh);
                        BP->push_update(prime);
                    }
                } else {
                    /* Odd j: test all ii-coordinates */
                    for (int ii = 0; ii < (i1 - i0); ii++) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose_bad) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region odd j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        ASSERT(prime.p >= si.conf.td_thresh);
                        BP->push_update(prime);
                    }
                }
                pos += gI;
            }
#if 0
            ssdpos[index] = pos - bucket_region;
            if (resieve_very_verbose_bad) {
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
