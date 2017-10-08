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


static const int bucket_region = 1 << LOG_BUCKET_REGION;

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
 * In other words, we can sieve this bad prime (power) much like a 
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


int
small_sieve_dump(FILE *f, const char *header, va_list va)
{
    const small_sieve_data_t * ssd = va_arg(va, const small_sieve_data_t *);
    const ssp_marker_t * next_marker = ssd->markers;
    const ssp_t *ssp = ssd->ssp;
    const unsigned char *logp = ssd->logp;

    fprintf(f, "%s", header);
    for (int i = 0; !(next_marker->index == i && next_marker->event == SSP_END); i++) {
        fprintf(f, "# p = %" FBPRIME_FORMAT ", r = %" FBROOT_FORMAT ", offset = %" FBPRIME_FORMAT ", logp = %hhu",
            ssp[i].p, ssp[i].r, ssp[i].offset, logp[i]);
        for ( ; next_marker->index == i; next_marker++) {
            switch (next_marker->event) {
                case SSP_POW2 : fprintf(f, " (power of 2)"); break;
                case SSP_PROJ : fprintf(f, " (projective root)"); break;
                case SSP_DISCARD : fprintf(f, "(discarded)"); break;
            }
        }
        fprintf(f, "\n");
    }
    return 1;
}

static void small_sieve_print_contents(const char * prefix, small_sieve_data_t * ssd)
{
    ssp_marker_t * next_marker = ssd->markers;
    int nice=ssd->nb_ssp;
    int nproj=0;
    int npow2=0;
    int ndiscard=0;
    for( ; next_marker->event != SSP_END ; next_marker++) {
        unsigned int event = next_marker->event;
        nproj += ((event & SSP_PROJ) != 0);
        npow2 += ((event & SSP_POW2) != 0);
        ndiscard += ((event & SSP_DISCARD) != 0);
        nice -= (nproj || npow2 || ndiscard) != 0;
    }
    ASSERT_ALWAYS(next_marker->index == ssd->nb_ssp);

    verbose_output_print(0, 2, "# %s: %d nice primes", prefix, nice);
    /* Primes may be both even and projective... */
    if (npow2) verbose_output_print(0, 2, ", %d powers of 2", npow2);
    if (nproj) verbose_output_print(0, 2, ", and %d projective primes", nproj);
    verbose_output_print(0, 2, ".");
    if (ndiscard) verbose_output_print(0, 2, " %d discarded.", ndiscard);
    verbose_output_print(0, 2, "\n");
    /* With -v -v -v, dump all the small sieve data */
    verbose_output_vfprint (0, 4, small_sieve_dump, "# Dump of small sieve data:\n", ssd);
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

/* This macro is used for sieve initialization */
void push_ssp_marker(small_sieve_data_t *ssd, int &nmarkers, const int __index, const int __event)
{
    ssd->markers = (ssp_marker_t *) realloc(ssd->markers, (nmarkers + 1) * sizeof(ssp_marker_t));
    ASSERT_ALWAYS(ssd->markers != NULL);
    ssd->markers[nmarkers].index = __index;
    ssd->markers[nmarkers].event = __event;
    nmarkers++;
}

/* {{{ Sieve initialization / clearing : first the easy ones */
void small_sieve_clear(small_sieve_data_t * ssd)
{
    free(ssd->ssp); ssd->ssp = NULL;
    free(ssd->logp); ssd->logp = NULL;
    free(ssd->markers); ssd->markers = NULL;
}

void small_sieve_extract_interval(small_sieve_data_t * r, small_sieve_data_t * s, int bounds[2])
{
    memset(r, 0, sizeof(small_sieve_data_t));
    r->nb_ssp = bounds[1] - bounds[0];
    r->ssp = (ssp_t *) malloc (r->nb_ssp * sizeof (ssp_t));
    FATAL_ERROR_CHECK(r->nb_ssp > 0 && r->ssp == NULL, "malloc failed");
    r->logp = (unsigned char *) malloc (r->nb_ssp);
    FATAL_ERROR_CHECK(r->nb_ssp > 0 && r->logp == NULL, "malloc failed");
    r->markers = NULL;
    int r_nmarkers = 0;

    memcpy(r->ssp, s->ssp + bounds[0], r->nb_ssp * sizeof (ssp_t));
    memcpy(r->logp, s->logp + bounds[0], r->nb_ssp);

    ssp_marker_t * next_marker = s->markers;
    for( ; next_marker->index < bounds[0] ; next_marker++);
    for( ; next_marker->index < bounds[1] ; next_marker++) {
        int fence = next_marker->index;
        unsigned int event = next_marker->event;
        push_ssp_marker(r, r_nmarkers, fence - bounds[0], event);
    }
    push_ssp_marker(r, r_nmarkers, bounds[1] - bounds[0], SSP_END);
}

/* }}} */

/* {{{ Sieve initialization: now the real stuff */

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// ssdpos points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION
// It could actually be larger than 32 bits when I > 16.

/* Initialization procedures for the ssp data */

static inline void ssp_init_oa(ssp_t * tail, fbprime_t p, fbprime_t r, unsigned int skip, where_am_I & w MAYBE_UNUSED)/*{{{*/
{
    tail->p = p;
    tail->r = r;
    tail->offset = (r * skip) % p;
}/*}}}*/

static inline void ssp_init_op(ssp_bad_t * tail, fbprime_t p, fbprime_t r, unsigned int skip MAYBE_UNUSED, where_am_I & w MAYBE_UNUSED)/*{{{*/
{
    unsigned int v = r; /* have consistent notations */
    unsigned int g = gcd_ul(p, v);
    fbprime_t q = p / g;
    tail->g = g;
    tail->q = q;
    if (q == 1) {
        ASSERT(r == 0);
        tail->U = 0;
    } else {
        int rc;
        uint32_t U = v / g; /* coprime to q */
        rc = invmod_32(&U, q);
        ASSERT_ALWAYS(rc != 0);
        tail->U = U;
    }
}/*}}}*/

void small_sieve_init(small_sieve_data_t *ssd, las_info & las,
                      const std::vector<fb_general_entry> *fb,
                      sieve_info const & si, const int side)
{
    const unsigned int thresh = si.conf.bucket_thresh;
    const int verbose = 0;
    where_am_I w;

    size_t size = si.sides[side].fb_parts_x->rest[1];

    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->ssp = (ssp_t *) malloc(size * sizeof(ssp_t));
    FATAL_ERROR_CHECK(ssd->ssp == NULL, "malloc failed");
    ssd->markers = NULL;
    int nmarkers = 0;
    ssd->logp = (unsigned char *) malloc(size);
    FATAL_ERROR_CHECK(ssd->logp == NULL, "malloc failed");
    // Do another pass on fb and badprimes, to fill in the data
    // while we have any regular primes or bad primes < thresh left
    ssp_t * tail = ssd->ssp;

    unsigned int index = 0;

    // The processing of bucket region by nb_threads is interleaved.
    // It means that the positions for the small sieve must jump
    // over the (nb_threads - 1) regions after each region.
    // For typical primes, this jump is easily precomputed and goes into
    // the ssp struct.
    
    // If we are doing sublattices modulo m, then we jump virutally m
    // times faster.
    unsigned int sublatm = si.conf.sublat.m;
    const unsigned int skiprows = (bucket_region >> si.conf.logI_adjusted)*(las.nb_threads-1);

    for (std::vector<fb_general_entry>::const_iterator iter = fb->begin() ; iter != fb->end() && index < size ; iter++) {
        /* p=pp^k, the prime or prime power in this entry, and pp is prime */
        const fbprime_t p = iter->q, pp = iter->p;
        WHERE_AM_I_UPDATE(w, p, p);

        if (p > thresh) {
            continue;
        }

        const double log_scale = si.sides[side].lognorms->scale;

        for (int nr = 0; nr < iter->nr_roots; nr++, index++) {
            const fb_general_root *root = &(iter->roots[nr]);
            /* Convert into old format for projective roots by adding p if projective */
            const fbroot_t r = root->r + (root->proj ? p : 0);
            /* skip q: mark it in place of p to be recognized quickly.
               TODO: maybe use the SSP_DISCARD mechanism ?
               */
            if (!si.doing.is_coprime_to(pp)) {
                tail->p = mpz_get_ui(si.doing.p);
                tail++;
                continue;
            }
            unsigned int event = 0;

            if (sublatm) {
                // In sublat mode, disable pattern sieving and primes
                // dividing m.
                if (pp == 2 || pp == 3 || (sublatm % pp) == 0) {
                    event |= SSP_DISCARD;
                }
            }

            if ((p & 1)==0) event |= SSP_POW2;

            /* p may already have occurred before, and was taken into account
               to the power 'oldexp', with a log contribution of 'old_log'. We
               now want to take into account the extra contribution of going
               from p^oldexp to p^exp. */
            const double old_log = (root->oldexp == 0) ? 0. :
                fb_log (fb_pow (pp, root->oldexp), log_scale, 0.);
            ssd->logp[index] = fb_log (fb_pow (pp, root->exp), log_scale, - old_log);
            
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
            if (r_q >= p && !(event & SSP_DISCARD)) {
                /* Compute the init data in any case, since the gcd
                 * dominates (and anyway we won't be doing this very
                 * often). */
                event |= SSP_PROJ;
                ssp_bad_t * ssp = (ssp_bad_t *) tail;
                ssp_init_op(ssp, p, r_q - p, skiprows, w);
                /* If g exceeds J, then the only reached locations in the
                 * sieving area will be on line (j=0), thus (1,0) only since
                 * the other are equivalent.
                 */
                if (ssp->g >= si.J) {
                    if (verbose) {
                        verbose_output_print(0, 1,
                                "# small_sieve_init: not adding bad prime"
                                " (1:%" FBROOT_FORMAT ") mod %" FBPRIME_FORMAT ")"
                                " to small sieve  because g=%d >= si.J = %d\n",
                                r_q-p, p, ssp->g, si.J);
                    }
                    event |= SSP_DISCARD;
                }
            } else if (!(event & SSP_DISCARD)) {
                ssp_init_oa(tail, p, r_q, skiprows, w);
            }

            tail++;
            if (event)
                push_ssp_marker(ssd, nmarkers, index, event);
        }
    }
    push_ssp_marker(ssd, nmarkers, index, SSP_END);
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

/* Only compute the initial ssdpos fields. */
int64_t *small_sieve_start(small_sieve_data_t *ssd, unsigned int j0,
        sieve_info const & si)
{
    ssp_marker_t * next_marker = ssd->markers;
    int64_t * ssdpos = (int64_t *) malloc(ssd->nb_ssp * sizeof(int64_t));
    const uint64_t sublatm = si.conf.sublat.m;
    const uint64_t sublati0 = si.conf.sublat.i0;
    const uint64_t sublatj0 = si.conf.sublat.j0;

    for(int i = 0 ; i < ssd->nb_ssp ; i++) {
        int fence;
        unsigned int event;
        event = next_marker->event;
        fence = next_marker->index;
        next_marker++;
        for( ; i < fence ; i++) {
            ssp_t * ssp = &(ssd->ssp[i]);
            if (!sublatm) {
                uint64_t compensate = si.I / 2;
                compensate += (uint64_t)j0 * (uint64_t)ssp->r;
                ssdpos[i] = compensate % (uint64_t)ssp->p;
            } else {
                uint64_t jj = j0*sublatm + sublatj0;
                uint64_t ii = (jj*uint64_t(ssp->r)) % uint64_t(ssp->p);
                while ((ii % sublatm) != sublati0) {
                    ii += ssp->p;
                }
                ii = (ii-sublati0) / sublatm; // exact division
                uint64_t compensate = (si.I / 2) + ii;
                ssdpos[i] = compensate % (uint64_t)ssp->p;
            }
        }
        if (event & SSP_DISCARD) continue;
        if (event & SSP_END) break;
        // Remark: even in the case of discard, we want to precompute the
        // data (we are in the projective case), for the (rare) case
        // where we have to update the (i,j)=(1,0) position.
        if (event & SSP_PROJ) {
            uint64_t jj0;
            if (!sublatm) {
                jj0 = j0;
            } else {
                jj0 = j0*sublatm + sublatj0;
            }
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[i]);
            /* Compute the next multiple of g above j0 */
            uint64_t j1 = jj0 - (jj0 % ssp->g);
            // FIXME: is it (j1<jj0) the right condition ?
            if (j1<jj0) { /* most often j1 is < jj0 -- in this case,
                         the j1 we're looking for needs +g */
                j1 += ssp->g;
            }
            ASSERT(j1 >= jj0);
            ASSERT(j1 % ssp->g == 0);
            /* Now we'd like to avoid row number 0 (so j1 == 0).  */
            /* The position (1,0) is updated with a special code in the */
            /* sieving code below */
            if (j1 == 0) {
                j1 += ssp->g;
            }
            // In sublat mode, we also need j1 congruent to sublatj0 mod m.
            if (sublatm) {
                while ((j1 % sublatm) != sublatj0) {
                    j1 += ssp->g;
                }
            }
            // Find the corresponding i
            uint64_t ii = uint64_t(j1/ssp->g)*uint64_t(ssp->U);
            if (sublatm) {
                while ((ii % sublatm) != sublati0) {
                    ii += ssp->q;
                }
            }
            // In the sublat mode, switch back to reduced convention
            // (exact divisions)
            if (sublatm) {
                j1 = (j1-sublatj0) / sublatm;
                ii = (ii-sublati0) / sublatm;
            }
            uint64_t compensate = si.I / 2;
            compensate += ii;
            ssdpos[i] = (((uint64_t)(j1 - j0))<<si.conf.logI_adjusted)
                + compensate % (uint64_t)ssp->q;
        } else if (event & SSP_POW2) {
            // FIXME: this is probably broken with sublattices.
            /* For powers of 2, we sieve only odd lines (*) and 
             * ssdpos needs to point at line j=1. We assume
             * that in this case (si.I/2) % p == 0
             * (*) for lines with j even, we have a root mod the prime
             * power for i-j*r multiple of our power of 2, which means
             * i even too. Thus a useless report.
             */
            ssp_t * ssp = &(ssd->ssp[i]);
            /* Note that j0 may perfectly be odd, in the case I==16 ! */
            unsigned int j1 = j0 | 1;
            unsigned int compensate = si.I / 2;
            compensate += j1 * ssp->r;
            ssdpos[i] = (j1-j0) * si.I + compensate % ssp->p;
        }
    }
    return ssdpos;
}
/* }}} */



/* {{{ Skip stride -- used in multithreaded context only */
void small_sieve_skip_stride(small_sieve_data_t *ssd, int64_t * ssdpos,
        unsigned int skip, sieve_info const & si)
{
    if (skip == 0) return;

    ssp_marker_t * next_marker = ssd->markers;
    unsigned int skipI = skip << si.conf.logI_adjusted;

    for(int i = 0 ; i < ssd->nb_ssp ; i++) {
        int fence;
        unsigned int event;
        event = next_marker->event;
        fence = next_marker->index;
        next_marker++;
        for( ; i < fence ; i++) {
            ssp_t * ssp = &(ssd->ssp[i]);
            ssdpos[i] += ssp->offset;
            if (ssdpos[i] >= (int) ssp->p)
                ssdpos[i] -= ssp->p;
        }
        if (event & SSP_DISCARD) continue;
        if (event & SSP_END) break;
        if (event & SSP_PROJ) {
            // FIXME: this is probably broken with sublattices.
            /* Don't bother. Pay attention to the fact that we have offsets to
             * the (current) bucket base. */
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[i]);
            uint64_t x = ssdpos[i];
            const unsigned int I = 1U << si.conf.logI_adjusted;
            unsigned int imask = I-1;
            unsigned int j = x >> si.conf.logI_adjusted;
            if (j >= skip) {
                /* The ``ssdpos'' is still ahead of us, so there's
                 * no adjustment to make */
              x -= skipI;
            } else {
                /* We've hit something in this bucket, but the
                 * ssdpos field lands in the blank space between
                 * this bucket and the next one to be handled. So we must
                 * advance: add g to j enough times so that j>=skip.
                 * Which means j+g*ceil((skip-j)/g)
                 */
                uint64_t i = x & imask;
                uint64_t jI = x - i;
                uint64_t nskip = iceildiv(skip-j, ssp->g);
                jI = jI
                    + (((uint64_t)(nskip * ssp->g - skip)) << si.conf.logI_adjusted);
                i = (i + nskip * (uint64_t)ssp->U) % ssp->q;
                x = jI + i;
            }
            ssdpos[i] = x;
        } else if (event & SSP_POW2) {
            // FIXME: this is probably broken with sublattices.
            ssp_t * ssp = &(ssd->ssp[i]);
            ssdpos[i] += ssp->offset;
            // If there is only one line per bucket region, we should
            // skip the even lines. So we add I to ensure that we start
            // at the next odd line.
            const unsigned long nj = bucket_region >> si.conf.logI_adjusted;
            if (nj == 1)
                ssdpos[i] += bucket_region;
            /* Pay attention to the fact that at the moment, ssdpos
             * may still point to the _second line_ in the area. So we
             * must not cancel the high bits.
             */
            // ssdpos[i] &= ssp->p - 1;

            // However, when p == I, we might get something >= 2*I, which
            // is just an artifact.
            // See bug #18814
            const unsigned int I = 1U << si.conf.logI_adjusted;
            if ((ssp->p == I) && ((unsigned int)ssdpos[i] >= 2*I)) {
                ssdpos[i] -= I;
            }
        }
    }
}
/* }}} */

/* This adds extra logging for pattern sieving. Very slow.
 */
#define xxxUGLY_DEBUGGING

/* {{{ Normal small sieve */

// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, int N,
			       small_sieve_data_t * ssd, int64_t * ssdpos,
                               sieve_info & si, int side,
			       where_am_I & w MAYBE_UNUSED)
{
    const uint32_t I = si.I;
    const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    unsigned long j;
    const int test_divisibility = 0; /* very slow, but nice for debugging */
    const unsigned long nj = bucket_region >> si.conf.logI_adjusted; /* Nr. of lines 
                                                           per bucket region */
    /* In order to check whether a j coordinate is even, we need to take
     * into account the bucket number, especially in case buckets are as
     * large as the sieve region. The row number corresponding to a given
     * i0 is i0/I, but we also need to add bucket_nr*bucket_size/I to
     * this, which is what this flag is for.
     * Sublat must also be taken into account.
     */
    int row0_is_oddj;
    if (si.conf.sublat.m == 0) {
        row0_is_oddj = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted)) & 1;
    } else {
        int row0 = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted));
        row0_is_oddj = (row0 * si.conf.sublat.m + si.conf.sublat.j0) & 1;
        // Odd/even property of j is the same as for j+2, even with
        // sublat, unless sublat.m is even, which is not handled right
        // now. Same for i.
        ASSERT_ALWAYS((si.conf.sublat.m & 1) == 1);
    }

    /* Handle powers of 2 up to 2 * sizeof(long) separately. 
     * TODO: use SSE2 */
    WHERE_AM_I_UPDATE(w, p, 2);
    /* First collect updates for powers of two in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (j = 0; j < nj; j++)
    {
        WHERE_AM_I_UPDATE(w, j, j);
        unsigned long pattern[2] = {0,0};

        /* Prepare the pattern */

        ssp_marker_t * next_marker = ssd->markers;
        int fence = -1;
        unsigned int event = 0;
        int * interval = si.sides[side].fb_parts_x->pow2;
        for(int n = interval[0] ; n < interval[1] ; n++) {
            for( ; fence < n || event == SSP_POW2 ; next_marker++) {
                event = next_marker->event;
                fence = next_marker->index;
            }
            if (n < fence) {
                const fbprime_t p = ssd->ssp[n].p;
                if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
                    continue;
                }
                uint64_t i0 = ssdpos[n];
                if (i0 < I) {
                    // At the start of the bucket region, we are not sure
                    // that i0 is reduced mod p when we enter here.
                    // (see comment at the end of small_sieve_skip_stride()).
                    // Hence, we have to do this reduction here:
                    i0 &= (p-1);
                    ASSERT (i0 < p);
                    ASSERT ((nj * N + j) % 2 == 1);
                    for (unsigned int i = i0; i < pattern2_size; i += p)
                        ((unsigned char *)pattern)[i] += ssd->logp[n];
#ifdef UGLY_DEBUGGING
                    for (unsigned int j = i0; j < I ; j+= p) {
                        WHERE_AM_I_UPDATE(w, x, (w.j << si.conf.logI_adjusted) + j);
                        sieve_increase(S + j, ssd->logp[n], w);
                        /* cancel the above action */
                        S[j] += ssd->logp[n];
                    }
#endif
                    /* Skip two lines above, since we sieve only odd lines.
                     * Even lines would correspond to useless reports.
                     */
                    i0 = ((i0 + 2 * ssd->ssp[n].r) & (p - 1))
                      + (2 << si.conf.logI_adjusted);
                }
                /* In this loop, ssdpos gets updated to the first 
                   index to sieve relative to the start of the next line, 
                   but after all lines of this bucket region are processed, 
                   it will point to the first position to sieve relative
                   to the start of the next bucket region, as required */
                ssdpos[n] = i0 - I;
            } else {
                /* nothing. It's a (presumably) projective power of 2,
                 * but for the moment these are not pattern-sieved. */
            }
        }

        /* Apply the pattern */
        if (pattern[0] || pattern[1]) {
          unsigned long *S_ptr = (unsigned long *) (S + (j << si.conf.logI_adjusted));
          const unsigned long *end = (unsigned long *)(S + (j << si.conf.logI_adjusted) + I);

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

            while (S_ptr < end)
            {
                *(S_ptr) += pattern[0];
                *(S_ptr + 1) += pattern[1];
                *(S_ptr + 2) += pattern[0];
                *(S_ptr + 3) += pattern[1];
                S_ptr += 4;
            }
        }
    }


    /* Handle 3 */
    /* For the time being, it's really 3. But the only thing we care about is
     * the hit rate to be 1 every 3rd, on lines for which we hit. So the code
     * below could be improved to handle more cases, although presently it
     * does not. */
    WHERE_AM_I_UPDATE(w, p, 3);
    /* First collect updates for powers of three in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (j = 0; j < nj; j++)
    {
        WHERE_AM_I_UPDATE(w, j, j);
        unsigned long pattern[3];

        pattern[0] = pattern[1] = pattern[2] = 0UL;

        ssp_marker_t * next_marker = ssd->markers;
        int fence = -1;
        // unsigned int event = 0;
        int * interval = si.sides[side].fb_parts_x->pow3;
        for(int n = interval[0] ; n < interval[1] ; n++) {
            for( ; fence < n ; next_marker++) {
                // event = next_marker->event;
                fence = next_marker->index;
            }
            if (n < fence) { /* a nice prime */
                ASSERT_ALWAYS(ssd->ssp[n].p == 3);
                const fbprime_t p = 3;
                WHERE_AM_I_UPDATE(w, p, p);
                unsigned int i0 = ssdpos[n];
                unsigned int i;
                ASSERT (i0 < p);
                for (i = i0; i < 3 * sizeof(unsigned long); i += p)
                    ((unsigned char *)pattern)[i] += ssd->logp[n];
                i0 += ssd->ssp[n].r;
                if (i0 >= p)
                    i0 -= p;
                ssdpos[n] = i0;
            } else {
                /* n points to a power of 3, and we have an exceptional
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
        
        if (!((j ^ row0_is_oddj)&1)) {
            /* We have an even j. There, we must not sieve even i's either !
             */
            for (unsigned int i = 0; i < 3 * sizeof(unsigned long); i += 2)
                ((unsigned char *)pattern)[i] = 0;
        }

        /* We want to test if there is a non-zero entry in the pattern
         * within the first 6 entries (we sieve mod 3 and mod 2).
         * For compatibility with 32-bit machines, we test the first two
         * pattern unsigned longs */
        if (pattern[0] || pattern[1]) {
          unsigned long *S_ptr = (unsigned long *) (S + (j << si.conf.logI_adjusted));
          const unsigned long *end = (unsigned long *)(S + (j << si.conf.logI_adjusted) + I) - 2;
            
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

            while (S_ptr < end)
            {
                *(S_ptr) += pattern[0];
                *(S_ptr + 1) += pattern[1];
                *(S_ptr + 2) += pattern[2];
                S_ptr += 3;
            }

            end += 2;
            if (S_ptr < end)
                *(S_ptr++) += pattern[0];
            if (S_ptr < end)
                *(S_ptr) += pattern[1];
        }
    }
    
    ssp_marker_t * next_marker = ssd->markers;
    
    // sieve with everyone, since pattern-sieving may miss some of the
    // small primes.
    
    /* use size_t for k, i and twop to speed up in x86 with LEA asm instructions */
    for(size_t k = 0 ; k < (size_t) ssd->nb_ssp ; k++) {
      int fence = next_marker->index;
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
#if defined( HAVE_SSE2 ) && !defined( TRACK_CODE_PATH ) /* x86 optimized code */
#define T
#define U1								\
        "addb %4,(%1)\n"                /* sieve_increase(pi,logp,w) */	\
	"addb %4,(%1,%3,1)\n"   /* sieve_increase(p_or_2p+pi,logp,w) */ \
	"addb %4,(%1,%3,2)\n" /* sieve_increase(p_or_2p*2+pi,logp,w) */ \
	"lea (%1,%2,1),%1\n"                    /* pi += p_or_2p * 3 */
#define U2								\
        "cmp %5, %1\n"		          /* if (pi >= S_ptr) break; */	\
	"jae 2f\n"							\
        "addb %4,(%1)\n"		/* sieve_increase(pi,logp,w) */	\
        "lea (%1,%3,1),%1\n"	                    /* pi += p_or_2p */		       
#define U do {								\
	unsigned char *pi_end;						\
	size_t three_p_or_2p;						\
	__asm__ __volatile__ (						\
        "lea (%3,%3,2), %2\n"         /* three_p_or_2p = p_or_2p * 3 */ \
	"lea (%1,%2,4), %0\n"            /* pi_end = pi + p_or_2p*12 */ \
	"cmp %5, %0\n"	              /* if (pi_end > S_ptr) no loop */	\
	"jbe 0f\n"							\
	"1:\n"								\
        U2 U2 U2 U2 U2 U2 U2 U2 U2 U2 U2				\
	"cmp %5, %1\n"							\
	"jae 2f\n"							\
        "addb %4,(%1)\n"						\
	"jmp 2f\n"							\
	".balign 8\n 0:\n"                              /* Main loop */	\
	U1 U1	        /* sieve_increase(p_or_2p*[0..11]+pi,logp,w) */ \
        U1 U1		                         /* pi += p_or_2p*12 */	\
        "lea (%1,%2,4), %0\n"	 /* if (pi+p_or_2p*12 > S_ptr) break */	\
	"cmp %5, %0\n"							\
	"jbe 0b\n"							\
	"jmp 1b\n"							\
	".balign 8\n 2:\n"						\
	: "=&r"(pi_end), "+&r"(pi), "=&r"(three_p_or_2p)		\
	: "r"(p_or_2p), "q"(logp), "r"(S_ptr) : "cc");			\
      } while (0)
#else
#define T do {							\
	WHERE_AM_I_UPDATE(w, x, j * I + pi + I - S_ptr);	\
	sieve_increase (pi, logp, w); pi += p_or_2p;		\
      } while(0)
#define U1
#define U2
#define U do {								\
	while (UNLIKELY(pi + p_or_2p * 12 <= S_ptr))			\
	  { T; T; T; T; T; T; T; T; T; T; T; T; }			\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
	if (pi >= S_ptr) break; T; if (pi >= S_ptr) break; T;		\
      } while (0)
#endif
      for( ; k < (size_t) fence ; k++) {
	const fbprime_t p = ssd->ssp[k].p;
        if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
            continue;
        }
	/* Don't sieve 3 again as it was pattern-sieved -- unless
	 * it's projective, but in this branch we have no projective
	 * primes. */
	if (p == 3) continue;
	const fbprime_t r =  ssd->ssp[k].r;
	WHERE_AM_I_UPDATE(w, p, p);
	const unsigned char logp = ssd->logp[k];
	unsigned char *S_ptr = S;
	size_t p_or_2p = p;
	unsigned int i0 = ssdpos[k]; ASSERT(i0 < p);
        unsigned int i_compens_sublat = 0;
        if (si.conf.sublat.m != 0) {
            i_compens_sublat = si.conf.sublat.i0 & 1;
        }
	j = 0;
	if (row0_is_oddj) goto j_odd;
	do {
	  unsigned char *pi;
	  WHERE_AM_I_UPDATE(w, j, j);
	  pi = S_ptr + i0;
	  S_ptr += I;
	  /* for j even, we sieve only odd pi, so step = 2p. */
	  p_or_2p += p_or_2p;
          if (!((i_compens_sublat + i0) & 1)) {
              pi += p;
          }
	  U;
	  i0 += r; if (i0 >= p) i0 -= p; 
	  /* Next line */
	  if (++j >= nj) break;
	  p_or_2p >>= 1;
	j_odd:
	  WHERE_AM_I_UPDATE(w, j, j);
	  pi = S_ptr + i0;
	  S_ptr += I;
	  U;
	  i0 += r; if (i0 >= p) i0 -= p;
	} while (++j < nj);
	ssdpos[k] = i0;
      }
#undef T
#undef U1
#undef U2
#undef U
        unsigned int event = (next_marker++)->event;
        if (event & SSP_DISCARD) {
            continue;
        }
        if (event & SSP_END) {
            ASSERT_ALWAYS(fence == ssd->nb_ssp);
            break;
        }
        if (event & SSP_PROJ) {
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[k]);
            const fbprime_t g = ssp->g;
            const uint64_t gI = (uint64_t)ssp->g << si.conf.logI_adjusted;
            const fbprime_t q = ssp->q;
            const fbprime_t U = ssp->U;
            const fbprime_t p MAYBE_UNUSED = g * q;
            WHERE_AM_I_UPDATE(w, p, p);
            const unsigned char logp = ssd->logp[k];
            /* Sieve the bad primes. We have p^k | fij(i,j) for i,j such
             * that i * g == j * U (mod p^k) where g = p^l and gcd(U, p)
             * = 1.  This hits only for g|j, then j = j' * g, and i == j'
             * * U (mod p^(k-l)).  In every g-th line, we sieve the
             * entries with i == (j/g)*U (mod q).  In ssd we have stored
             * g, q = p^(k-l), U, and ssdpos so that S +
             * ssdpos is the next sieve entry that needs to be
             * sieved.  So if S + ssdpos is in the current bucket
             * region, we update all  S + ssdpos + n*q  where
             * ssdpos + n*q < I, then set ssdpos =
             * ((ssdpos % I) + U) % q) + I * g.  */
            if (!test_divisibility && ssp->q == 1)
            {
                /* q = 1, therefore U = 0, and we sieve all entries in lines
                   with g|j, beginning with the line starting at S[ssdpos] */
                unsigned long logps;
                uint64_t i0 = ssdpos[k];
                // The following is for the case where p divides the norm
                // at the position (i,j) = (1,0).
                if (UNLIKELY(N == 0 && i0 == gI)) {
#ifdef TRACE_K
                    if (trace_on_spot_Nx(w.N, 1+(I>>1))) {
                        WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                        unsigned int v = logp;
                        sieve_increase_logging(S + w.x, v, w);
                    }
#endif
                    S[1+(I>>1)] += logp;
                }
                // The event SSP_DISCARD might have occurred due to
                // the first row to sieve being larger than J. The row
                // number 0 must still be sieved in that case, but once
                // it's done, we can indeed skip the next part of
                // sieving.
                if (event & SSP_DISCARD)
                    continue;
                ASSERT (ssp->U == 0);
                ASSERT (i0 % I == 0);
                ASSERT (I % (4 * sizeof (unsigned long)) == 0);
                for (j = 0; j < sizeof (unsigned long); j++)
                    ((unsigned char *)&logps)[j] = logp;
                while (i0 < (unsigned int) bucket_region)
                {
                    unsigned long *S_ptr = (unsigned long *) (S + i0);
                    unsigned long *end = S_ptr + I / sizeof (unsigned long);
                    unsigned long logps2 = logps;
                    /* Check whether the j coordinate is even. */
                    if (((i0 & I) == 0) ^ row0_is_oddj) {
                        /* Yes, j is even. We update only odd i-coordinates */
                        /* Use array indexing to avoid endianness issues. */
                        for (j = 0; j < sizeof (unsigned long); j += 2)
                            ((unsigned char *)&logps2)[j] = 0;
                    }
#ifdef TRACE_K
                    if (trace_on_range_Nx(w.N, i0, i0 + I)) {
                        WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                        unsigned int x = trace_Nx.x;
                        unsigned int k = x % I;
                        unsigned int v = (((unsigned char *)(&logps2))[k%sizeof(unsigned long)]);
                        sieve_increase_logging(S + w.x, v, w);
                    }
#endif
                    while (S_ptr < end)
                    {
                        *(S_ptr) += logps2;
                        *(S_ptr + 1) += logps2;
                        *(S_ptr + 2) += logps2;
                        *(S_ptr + 3) += logps2;
                        S_ptr += 4;
                    }
                    i0 += gI;
                }
                ssdpos[k] = i0 - (1U << LOG_BUCKET_REGION);
            } else {
                /* q > 1, more general sieving code. */
                const uint64_t i0 = ssdpos[k];
                const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
                size_t lineoffset = i0 & (I - 1U),
                       linestart = i0 - lineoffset;
                ASSERT (U < q);
                while (linestart < (1U << LOG_BUCKET_REGION))
                {
                    WHERE_AM_I_UPDATE(w, j, linestart / I);
                    size_t i = lineoffset;
                    if (((linestart & I) == 0) ^ row0_is_oddj) /* Is j even? */
                    {
                        /* Yes, sieve only odd i values */
                        if (i % 2 == 0) /* Make i odd */
                            i += q;
                        if (i % 2 == 1) /* If not both i,q are even */
                            for ( ; i < I; i += evenq)
                            {
                                WHERE_AM_I_UPDATE(w, x, linestart + i);
                                sieve_increase (S + linestart + i, logp, w);
                            }
                    }
                    else
                    {
                        for ( ; i < I; i += q)
                        {
                            WHERE_AM_I_UPDATE(w, x, linestart + i);
                            sieve_increase (S + linestart + i, logp, w);
                        }
                    }

                    linestart += gI;
                    lineoffset += U;
                    if (lineoffset >= q)
                        lineoffset -= q;
                }
                ssdpos[k] = linestart + lineoffset - 
                    (1U << LOG_BUCKET_REGION);
            }
        } else if (event & SSP_POW2) {
            /* Powers of 2 are treated separately */
            /* Don't sieve powers of 2 again that were pattern-sieved */
            ssp_t * ssp = &(ssd->ssp[k]);
            const fbprime_t p = ssp->p;
            const fbprime_t r = ssp->r;
            WHERE_AM_I_UPDATE(w, p, p);

            if (p <= pattern2_size)
                continue;

            const unsigned char logp = ssd->logp[k];
            unsigned char *S_ptr = S;
            unsigned int linestart = 0;

            unsigned int i0 = ssdpos[k];

            j=0;
            /* Our encoding is that when the first row is even, the
             * i0 we keep in the ssdpos array deliberately includes an
             * additional quantity I to mean ``next row''. Compensate
             * that, move on to an odd row, and proceed. */
            if ((nj * N + j) % 2 == 0) {
                ASSERT(i0 >= I);
                i0 -= I; linestart += I; S_ptr += I;
                j++;
            }
            if (j < nj) i0 &= (p-1);
            for( ; j < nj ; j+= 2) {
                for (unsigned int i = i0; i < I; i += p) {
                    WHERE_AM_I_UPDATE(w, x, (j << si.conf.logI_adjusted) + i);
                    sieve_increase (S_ptr + i, logp, w);
                }
                // odd lines only.
                i0 = (i0 + (r << 1)) & (p - 1);
                linestart += I; S_ptr += I;
                linestart += I; S_ptr += I;
            }
            if (j > nj) i0 += I;
            ssdpos[k] = i0;
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
        sieve_info const & si, where_am_I & w MAYBE_UNUSED)
{
    const uint32_t I = si.I;
    unsigned char *S_ptr;
    unsigned long j, nj;
    const int resieve_very_verbose = 0, resieve_very_verbose_bad = 0;
    /* See comment above about the variable of the same name */
    int row0_is_oddj;
    unsigned int i_compens_sublat = 0;
    if (si.conf.sublat.m == 0) {
        row0_is_oddj = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted)) & 1;
    } else {
        int row0 = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted));
        row0_is_oddj = (row0 * si.conf.sublat.m + si.conf.sublat.j0) & 1;
        i_compens_sublat = si.conf.sublat.i0 & 1;
        // Odd/even property of j is the same as for j+2, even with
        // sublat, unless sublat.m is even, which is not handled right
        // now. Same for i.
        ASSERT_ALWAYS((si.conf.sublat.m & 1) == 1);
    }

    nj = (bucket_region >> si.conf.logI_adjusted);

    ssp_marker_t * next_marker = ssd->markers;

    for(int i = 0 ; i < ssd->nb_ssp ; i++) {
        int fence;
        unsigned int event;
        event = next_marker->event;
        fence = next_marker->index;
        next_marker++;
        for( ; i < fence ; i++) {
            ssp_t * ssp = &(ssd->ssp[i]);
            const fbprime_t p = ssp->p;
            if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
                continue;
            }
            fbprime_t r = ssp->r;
            WHERE_AM_I_UPDATE(w, p, p);
            unsigned int i0 = ssdpos[i];
            S_ptr = S;
            ASSERT(i0 < p);
            /* for j even, we sieve only odd i. This translates into loops
             * which look as follows:
             *
             * j even: (sieve only odd i)
             *   for(i = i0 + (p & -!(i0&1)) ; i < I ; i += p+p)
             *   (where (p & -!(i0&1)) is 0 if i0 is odd, and p otherwise)
             * j odd: (sieve all values of i)
             *   for(i = i0                  ; i < I ; i += p)
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
            for (j = 0; j < nj; j ++) {
                WHERE_AM_I_UPDATE(w, j, j);
                for (unsigned int i = i0 + (q& -!((i0+i_compens_sublat)&1)) ; i < I; i += p+q) {
                    if (LIKELY(S_ptr[i] == 255)) continue;
                    bucket_update_t<1, primehint_t> prime;
                    unsigned int x = (j << (si.conf.logI_adjusted)) + i;
                    if (resieve_very_verbose) {
                        verbose_output_print(0, 1, "resieve_small_bucket_region: root %"
                                FBROOT_FORMAT ",%d divides at x = "
                                "%d = %lu * %u + %d\n",
                                p, r, x, j, 1 << si.conf.logI_adjusted, i);
                    }
                    prime.p = p;
                    prime.x = x;
                    ASSERT(prime.p >= si.conf.td_thresh);
                    BP->push_update(prime);
                }
                i0 += r;
                if (i0 >= p)
                    i0 -= p;
                S_ptr += I;
                q ^= p;
            }
            ssdpos[i] = i0;
        }
        if (event == SSP_END) {
            break;
        }
        if (event & SSP_PROJ) {
            ssp_bad_t * ssp = (ssp_bad_t * ) &(ssd->ssp[i]);
            const fbprime_t g = ssp->g;
            const uint64_t gI = (uint64_t)g << si.conf.logI_adjusted;

            WHERE_AM_I_UPDATE(w, p, g * ssp->q);

            /* Test every p-th line, starting at S[ssdpos] */
            uint64_t i0 = ssdpos[i];
            // This block is for the case where p divides at (1,0).
            if (UNLIKELY(N == 0 && i0 == gI)) {
                bucket_update_t<1, primehint_t> prime;
                prime.p = g;
                prime.x = 1+(I>>1);
                ASSERT(prime.p >= si.conf.td_thresh);
                BP->push_update(prime);
            }
            // Same as in sieving: we discard after checking for row 0.
            if (event & SSP_DISCARD)
                continue;

            unsigned int ii;
            ASSERT (i0 % I == 0); /* make sure ssdpos points at start
                                     of line */
            if (resieve_very_verbose_bad) {
                verbose_output_print(0, 1, "# resieving bad prime %" FBPRIME_FORMAT
                        ", i0 = %" PRIu64 "\n", g, i0);
            }
            while (i0 < (unsigned int) bucket_region) {
                unsigned char *S_ptr = S + i0;
                if (((i0 & I) == 0) ^ row0_is_oddj) { /* Even j coordinate? */
                    /* Yes, test only odd ii-coordinates */
                    for (ii = 1; ii < I; ii += 2) {
                        if (S_ptr[ii] != 255) {
                            bucket_update_t<1, primehint_t> prime;
                            const unsigned int x = i0 + ii;
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
                    }
                } else {
                    /* No, test all ii-coordinates */
                    for (ii = 0; ii < I; ii++) {
                        if (S_ptr[ii] != 255) {
                            bucket_update_t<1, primehint_t> prime;
                            const unsigned int x = i0 + ii;
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
                }
                i0 += gI;
            }
            ssdpos[i] = i0 - bucket_region;
            if (resieve_very_verbose_bad) {
                verbose_output_print(0, 1, "# resieving: new i0 = %" PRIu64
                        ", bucket_region = %d, "
                        "new ssdpos = %" PRIu64 "\n",
                        i0, bucket_region, ssdpos[i]);
            }
        }
    }
}
/* }}} */
