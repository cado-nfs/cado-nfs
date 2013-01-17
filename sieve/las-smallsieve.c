#include "cado.h"
#include <pthread.h>
#include "las-config.h"
#include "las-smallsieve.h"
#include "las-types.h"
#include "las-debug.h"
#include "las-qlattice.h"
#include "basicnt.h"
#include "misc.h"

/* It's defined as a global variable in las.c */
extern pthread_mutex_t io_mutex;

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
 * In other words, we can sieve this bad prime (powers) much like a 
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

static void small_sieve_print_contents(FILE * f, const char * prefix, small_sieve_data_t * ssd)
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

    fprintf (f, "# %s: %d nice primes", prefix, nice);
    /* Primes may be both even and projective... */
    if (npow2) fprintf (f, ", %d powers of 2", npow2);
    if (nproj) fprintf (f, ", and %d projective primes", nproj);
    fprintf (f, ".");
    if (ndiscard) fprintf (f, " %d discarded.", ndiscard);
    fprintf (f, "\n");
}


void small_sieve_info(sieve_info_srcptr si, const char * what, int side, small_sieve_data_t * r)
{
    if (!si->verbose) return;
    char * tmp;
    int rc = asprintf(&tmp, "%s(%s side)", what, sidenames[side]);
    ASSERT_ALWAYS(rc >= 0);
    small_sieve_print_contents(si->output, tmp, r);
    free(tmp);
}

/* }}} */

/* This macro is used for sieve initialization */

#define PUSH_SSP_MARKER(ssd, nmarkers, __index, __event) do {		\
    ssd->markers = (ssp_marker_t *) realloc(ssd->markers,               \
            (nmarkers + 1) * sizeof(ssp_marker_t));                     \
    ssd->markers[nmarkers].index = __index;				\
    ssd->markers[nmarkers].event = __event;				\
    nmarkers++;								\
} while (0)

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
    FATAL_ERROR_CHECK(r->ssp == NULL, "malloc failed");
    r->logp = malloc (r->nb_ssp);
    FATAL_ERROR_CHECK(r->logp == NULL, "malloc failed");
    r->markers = NULL;
    int r_nmarkers = 0;

    memcpy(r->ssp, s->ssp + bounds[0], r->nb_ssp * sizeof (ssp_t));
    memcpy(r->logp, s->logp + bounds[0], r->nb_ssp);

    ssp_marker_t * next_marker = s->markers;
    for( ; next_marker->index < bounds[0] ; next_marker++);
    for( ; next_marker->index < bounds[1] ; next_marker++) {
        int fence = next_marker->index;
        unsigned int event = next_marker->event;
        PUSH_SSP_MARKER(r, r_nmarkers, fence - bounds[0], event);
    }
    PUSH_SSP_MARKER(r, r_nmarkers, bounds[1] - bounds[0], SSP_END);
}

/* }}} */

/* {{{ Sieve initialization: now the real stuff */

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// ssdpos points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION

/* Initialization procedures for the ssp data */

static inline void ssp_init_oa(ssp_t * tail, fbprime_t p, fbprime_t r, unsigned int skip, where_am_I_ptr w MAYBE_UNUSED)/*{{{*/
{
    tail->p = p;
    tail->r = r;
    tail->offset = (r * skip) % p;
}/*}}}*/

static inline void ssp_init_op(ssp_bad_t * tail, fbprime_t p, fbprime_t r, unsigned int skip MAYBE_UNUSED, where_am_I_ptr w MAYBE_UNUSED)/*{{{*/
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
        uint64_t U = v / g; /* coprime to q */
        rc = invmod(&U, q);
        ASSERT_ALWAYS(rc != 0);
        tail->U = U;
    }
}/*}}}*/

void small_sieve_init(small_sieve_data_t *ssd, const factorbase_degn_t *fb,
                      sieve_info_srcptr si, int side)
{
    const factorbase_degn_t *fb_sav = fb;
    int size = 0;
    const unsigned int thresh = si->bucket_thresh;
    const int verbose = 0;
    const int do_bad_primes = 1;
    where_am_I w;

    // Count prime ideals of factor base primes p < thresh
    while (fb->p != FB_END && fb->p < thresh) {
        size += fb->nr_roots;
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fb = fb_sav;

    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->ssp = (ssp_t *) malloc(size * sizeof(ssp_t));
    FATAL_ERROR_CHECK(ssd->ssp == NULL, "malloc failed");
    ssd->markers = NULL;
    int nmarkers = 0;
    ssd->logp = (unsigned char *) malloc(size);
    // Do another pass on fb and badprimes, to fill in the data
    // while we have any regular primes or bad primes < thresh left
    ssp_t * tail = ssd->ssp;

    int index;

    // The processing of bucket region by nb_threads is interleaved.
    // It means that the positions for the small sieve must jump
    // over the (nb_threads - 1) regions after each region.
    // For typical primes, this jump is easily precomputed and goes into
    // the ssp struct.
    
    unsigned int skiprows = (bucket_region >> si->logI)*(si->nb_threads-1);
    for (index = 0 ; fb->p != FB_END && fb->p < thresh ; fb = fb_next(fb)) {
        const fbprime_t p = fb->p;
        WHERE_AM_I_UPDATE(w, p, p);

        int nr;
        fbprime_t r;

        for (nr = 0; nr < fb->nr_roots; nr++, index++) {
            unsigned int event = 0;
            if ((fb->p&1)==0) event |= SSP_POW2;
            ssd->logp[index] = fb->plog;
            WHERE_AM_I_UPDATE(w, r, fb->roots[nr]);
            r = fb_root_in_qlattice_generic (p, fb->roots[nr], fb->invp, si);
            /* If this root is somehow interesting (projective in (a,b) or
               in (i,j) plane), print a message */
            if (verbose && (fb->roots[nr] >= p || r >= p))
                fprintf (si->output, "# small_sieve_init: %s side, prime " 
                        FBPRIME_FORMAT " root " FBPRIME_FORMAT " -> " 
                        FBPRIME_FORMAT "\n", sidenames[side], p, fb->roots[nr], r);

            /* Handle projective roots */
            if (r >= p) {
                /* Compute the init data in any case, since the gcd
                 * dominates (and anyway we won't be doing this very
                 * often). */
                event |= SSP_PROJ;
                ssp_bad_t * ssp = (ssp_bad_t *) tail;
                ssp_init_op(ssp, p, r - p, skiprows, w);
                /* If g exceeds J, then the only reached locations in the
                 * sieving area will be on line (j=0), thus (1,0) only since
                 * the other are equivalent.
                 */
                if (!do_bad_primes) {
                    if (verbose) {
                        fprintf (si->output,
                                "# small_sieve_init: not adding bad prime"
                                " (1:"FBPRIME_FORMAT") mod "FBPRIME_FORMAT")"
                                " to small sieve because do_bad_primes = 0\n",
                                r-p, p);
                    }
                    event |= SSP_DISCARD;
                } else if (ssp->g >= si->J) {
                    if (verbose) {
                        fprintf (si->output,
                                "# small_sieve_init: not adding bad prime"
                                " (1:"FBPRIME_FORMAT") mod "FBPRIME_FORMAT")"
                                " to small sieve  because g=%d >= si->J = %d\n",
                                r-p, p, ssp->g, si->J);
                    }
                    event |= SSP_DISCARD;
                }
            } else {
                ssp_init_oa(tail, p, r, skiprows, w);
            }
            tail++;
            if (event)
                PUSH_SSP_MARKER(ssd, nmarkers, index, event);
        }
    }
    PUSH_SSP_MARKER(ssd, nmarkers, index, SSP_END);
    ssd->nb_ssp = size;
}
/* }}} */

#undef PUSH_SSP_MARKER

/* {{{ Creation of the ssdpos tables */
int * small_sieve_copy_start(int * base, int bounds[2])
{
    int * res = malloc((bounds[1] - bounds[0]) * sizeof(int));
    memcpy(res, base + bounds[0], (bounds[1] - bounds[0]) * sizeof(int));
    return res;
}
/* Only compute the initial ssdpos fields. */
int * small_sieve_start(small_sieve_data_t *ssd, unsigned int j0, sieve_info_srcptr si)
{
    ssp_marker_t * next_marker = ssd->markers;
    int * ssdpos = (int *) malloc(ssd->nb_ssp * sizeof(int));

    for(int i = 0 ; i < ssd->nb_ssp ; i++) {
        int fence;
        unsigned int event;
        event = next_marker->event;
        fence = next_marker->index;
        next_marker++;
        for( ; i < fence ; i++) {
            ssp_t * ssp = &(ssd->ssp[i]);
            unsigned int compensate = si->I / 2;
            compensate += j0 * ssp->r;
            ssdpos[i] = compensate %ssp->p;
        }
        if (event & SSP_DISCARD) continue;
        if (event & SSP_END) break;
        if (event & SSP_PROJ) {
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[i]);
            /* Compute the next multiple of g above j0 */
            unsigned int j1 = j0 - (j0 % ssp->g);
            unsigned int compensate = si->I / 2;
            if (j0) { /* most often j1 is < j0 -- in this case,
                         the j1 we're looking for needs +g */
                j1 += ssp->g;
            }
            ASSERT(j1 >= j0);
            ASSERT(j1 % ssp->g == 0);
            /* Now we'd like to avoid row number 0 (so j1 == 0).  */
            /* Note that we avoid it entirely -- we could fathom sieving
             * (1,0), but it's probably not really worth it */
            if (j1 == 0) {
                j1 += ssp->g;
            }
            compensate += j1 * ssp->U;
            ssdpos[i] = (j1-j0) * si->I + compensate % ssp->q;
        } else if (event & SSP_POW2) {
            /* For powers of 2, we sieve only odd lines (*) and 
             * ssdpos needs to point at line j=1. We assume
             * that in this case (si->I/2) % p == 0
             * (*) for lines with j even, we have a root mod the prime
             * power for i-j*r multiple of our power of 2, which means
             * i even too. Thus a useless report.
             */
            ssp_t * ssp = &(ssd->ssp[i]);
            /* Note that j0 may perfectly be odd, in the case I==16 ! */
            unsigned int j1 = j0 | 1;
            unsigned int compensate = si->I / 2;
            compensate += j1 * ssp->r;
            ssdpos[i] = (j1-j0) * si->I + compensate % ssp->p;
        }
    }
    return ssdpos;
}
/* }}} */



/* {{{ Skip stride -- used in multithreaded context only */
void small_sieve_skip_stride(small_sieve_data_t *ssd, int * ssdpos, unsigned int skip, sieve_info_srcptr si)
{
    if (skip == 0) return;

    ssp_marker_t * next_marker = ssd->markers;

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
            /* Don't bother. Pay attention to the fact that we have offsets to
             * the (current) bucket base. */
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[i]);
            unsigned int x = ssdpos[i];
            const unsigned int I = 1U << si->logI;
            unsigned int imask = I-1;
            unsigned int j = x >> si->logI;
            if (j >= skip) {
                /* The ``ssdpos'' is still ahead of us, so there's
                 * no adjustment to make */
                x -= skip*I;
            } else {
                /* We've hit something in this bucket, but the
                 * ssdpos field lands in the blank space between
                 * this bucket and the next one to be handled. So we must
                 * advance: add g to j enough times so that j>=skip.
                 * Which means j+g*ceil((skip-j)/g)
                 */
                unsigned int i = x & imask;
                unsigned int jI = x - i;
                unsigned int nskip = iceildiv(skip-j, ssp->g);
                jI = jI + (nskip * ssp->g - skip) * I;
                i = (i + nskip * ssp->U) % ssp->q;
                x = jI + i;
            }
            ssdpos[i] = x;
        } else if (event & SSP_POW2) {
            ssp_t * ssp = &(ssd->ssp[i]);
            ssdpos[i] += ssp->offset;
            /* Pay attention to the fact that the moment, ssdpos
             * may still point to the _second line_ in the area. So we
             * must not cancel the high bits
             */
            // ssdpos[i] &= ssp->p - 1;
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
			       small_sieve_data_t * ssd, int * ssdpos, sieve_info_ptr si,
                               int side,
			       where_am_I_ptr w MAYBE_UNUSED)
{
    const uint32_t I = si->I;
    const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    unsigned long j;
    const int test_divisibility = 0; /* very slow, but nice for debugging */
    const unsigned long nj = bucket_region >> si->logI; /* Nr. of lines 
                                                           per bucket region */
    /* In order to check whether a j coordinate is even, we need to take
     * into account the bucket number, especially in case buckets are as
     * large as the sieve region. The row number corresponding to a given
     * i0 is i0/I, but we also need to add bucket_nr*bucket_size/I to
     * this, which is what this flag is for.
     */
    const int row0_is_oddj = (N << (LOG_BUCKET_REGION - si->logI)) & 1;


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
        int * interval = si->sides[side]->fb_parts_x->pow2;
        for(int n = interval[0] ; n < interval[1] ; n++) {
            for( ; fence < n || event == SSP_POW2 ; next_marker++) {
                event = next_marker->event;
                fence = next_marker->index;
            }
            if (n < fence) {
                const fbprime_t p = ssd->ssp[n].p;
                unsigned int i0 = ssdpos[n];
                if (i0 < I) {
                    ASSERT (i0 < p);
                    ASSERT ((nj * N + j) % 2 == 1);
                    for (unsigned int i = i0; i < pattern2_size; i += p)
                        ((unsigned char *)pattern)[i] += ssd->logp[n];
#ifdef UGLY_DEBUGGING
                    for (unsigned int j = i0; j < I ; j+= p) {
                        WHERE_AM_I_UPDATE(w, x, (w->j << si->logI) + j);
                        sieve_decrease(S + j, ssd->logp[n], w);
                        /* cancel the above action */
                        S[j] += ssd->logp[n];
                    }
#endif
                    /* Skip two lines above, since we sieve only odd lines.
                     * Even lines would correspond to useless reports.
                     */
                    i0 = ((i0 + 2 * ssd->ssp[n].r) & (p - 1)) + 2 * I;
                }
                /* In this loop, ssdpos gets updated to the first 
                   index to sieve relative to the start of the next line, 
                   but after all lines of this bucket region are processed, 
                   it will point the the first position to sieve relative  
                   to the start of the next bucket region, as required */
                ssdpos[n] = i0 - I;
            } else {
                /* nothing. It's a (presumably) projective power of 2,
                 * but for the moment these are not pattern-sieved. */
            }
        }

        /* Apply the pattern */
        if (pattern[0] || pattern[1]) {
            unsigned long *S_ptr = (unsigned long *) (S + j * I);
            const unsigned long *end = (unsigned long *)(S + j * I + I);

#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w->N, w->j*I, w->j*I+I)) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I;
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))&1)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_decrease_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            while (S_ptr < end)
            {
                *(S_ptr) -= pattern[0];
                *(S_ptr + 1) -= pattern[1];
                *(S_ptr + 2) -= pattern[0];
                *(S_ptr + 3) -= pattern[1];
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
        int * interval = si->sides[side]->fb_parts_x->pow3;
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

        if (pattern[0]) {
            unsigned long *S_ptr = (unsigned long *) (S + j * I);
            const unsigned long *end = (unsigned long *)(S + j * I + I) - 2;
            
#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w->N, w->j*I, w->j*I+I)) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I;
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))%3)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_decrease_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            while (S_ptr < end)
            {
                *(S_ptr) -= pattern[0];
                *(S_ptr + 1) -= pattern[1];
                *(S_ptr + 2) -= pattern[2];
                S_ptr += 3;
            }

            end += 2;
            if (S_ptr < end)
                *(S_ptr++) -= pattern[0];
            if (S_ptr < end)
                *(S_ptr) -= pattern[1];
        }
    }

    ssp_marker_t * next_marker = ssd->markers;

    // sieve with everyone, since pattern-sieving may miss some of the
    // small primes.

    for(int i = 0 ; i < ssd->nb_ssp ; i++) {
        int fence;
        unsigned int event;
        event = next_marker->event;
        fence = next_marker->index;
        next_marker++;
        for( ; i < fence ; i++) {
            ssp_t * ssp = &(ssd->ssp[i]);
            const fbprime_t p = ssp->p;
            const fbprime_t r = ssp->r;
            WHERE_AM_I_UPDATE(w, p, p);
            const unsigned char logp = ssd->logp[i];
            unsigned char *S_ptr = S;
            fbprime_t twop;
            unsigned int linestart = 0;
            /* Always S_ptr = S + linestart. S_ptr is used for the actual array
               updates, linestart keeps track of position relative to start of
               bucket region and is used only for computing i,j-coordinates
               in overflow and divisibility checking, and relation tracing. */

            unsigned int i0 = ssdpos[i];

            /* Don't sieve 3 again as it was pattern-sieved -- unless
             * it's projective, but in this branch we have no projective
             * primes. */
            if (p == 3)
                continue;

            ASSERT(i0 < p);
            for (j = 0; j < nj; j++) {
                WHERE_AM_I_UPDATE(w, j, j);
                twop = p;
                unsigned int i = i0;
                if ((((nj & N) ^ j) & 1) == 0) { /* (nj*N+j)%2 */
                    /* for j even, we sieve only odd i. */
                    twop += p;
                    i += (i0 & 1) ? 0 : p;
                }
                for ( ; i < I; i += twop) {
                    WHERE_AM_I_UPDATE(w, x, j * I + i);
                    sieve_decrease (S_ptr + i, logp, w);
                }
                i0 += r;
                if (i0 >= p)
                    i0 -= p;
                S_ptr += I;
                linestart += I;
            }
            ssdpos[i] = i0;
        }
        if (event == SSP_END) {
            ASSERT_ALWAYS(fence == ssd->nb_ssp);
            break;
        }
        if (event & SSP_DISCARD)
            continue;
        if (event & SSP_PROJ) {
            ssp_bad_t * ssp = (ssp_bad_t *) &(ssd->ssp[i]);
            const fbprime_t g = ssp->g;
            const fbprime_t q = ssp->q;
            const fbprime_t U = ssp->U;
            const fbprime_t p MAYBE_UNUSED = g * q;
            WHERE_AM_I_UPDATE(w, p, p);
            const unsigned char logp = ssd->logp[i];
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
                unsigned int i0 = ssdpos[i];
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
                    if (trace_on_range_Nx(w->N, i0, i0 + I)) {
                        WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                        unsigned int x = trace_Nx.x;
                        unsigned int k = x % I;
                        unsigned int v = (((unsigned char *)(&logps2))[k%sizeof(unsigned long)]);
                        sieve_decrease_logging(S + w->x, v, w);
                    }
#endif
                    while (S_ptr < end)
                    {
                        *(S_ptr) -= logps2;
                        *(S_ptr + 1) -= logps2;
                        *(S_ptr + 2) -= logps2;
                        *(S_ptr + 3) -= logps2;
                        S_ptr += 4;
                    }
                    i0 += ssp->g * I;
                }
                ssdpos[i] = i0 - (1U << LOG_BUCKET_REGION);
            } else {
                /* q > 1, more general sieving code. */
                const unsigned int i0 = ssdpos[i];
                const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
                unsigned int lineoffset = i0 & (I - 1U),
                             linestart = i0 - lineoffset;
                ASSERT (U < q);
                while (linestart < (1U << LOG_BUCKET_REGION))
                {
                    WHERE_AM_I_UPDATE(w, j, linestart / I);
                    unsigned int i = lineoffset;
                    if (((linestart & I) == 0) ^ row0_is_oddj) /* Is j even? */
                    {
                        /* Yes, sieve only odd i values */
                        if (i % 2 == 0) /* Make i odd */
                            i += q;
                        if (i % 2 == 1) /* If not both i,q are even */
                            for ( ; i < I; i += evenq)
                            {
                                WHERE_AM_I_UPDATE(w, x, linestart + i);
                                sieve_decrease (S + linestart + i, logp, w);
                            }
                    }
                    else
                    {
                        for ( ; i < I; i += q)
                        {
                            WHERE_AM_I_UPDATE(w, x, linestart + i);
                            sieve_decrease (S + linestart + i, logp, w);
                        }
                    }

                    linestart += g * I;
                    lineoffset += U;
                    if (lineoffset >= q)
                        lineoffset -= q;
                }
                ssdpos[i] = linestart + lineoffset - 
                    (1U << LOG_BUCKET_REGION);
            }
        } else if (event & SSP_POW2) {
            /* Powers of 2 are treated separately */
            /* Don't sieve powers of 2 again that were pattern-sieved */
            ssp_t * ssp = &(ssd->ssp[i]);
            const fbprime_t p = ssp->p;
            const fbprime_t r = ssp->r;
            WHERE_AM_I_UPDATE(w, p, p);

            if (p <= pattern2_size)
                continue;

            const unsigned char logp = ssd->logp[i];
            unsigned char *S_ptr = S;
            unsigned int linestart = 0;

            unsigned int i0 = ssdpos[i];
            for (j = 0; j < nj; j++) {
                WHERE_AM_I_UPDATE(w, j, j);
                if (i0 < I) {
                    ASSERT(i0 < p);
                    ASSERT ((nj * N + j) % 2 == 1);
                    for (unsigned int i = i0; i < I; i += p) {
                        WHERE_AM_I_UPDATE(w, x, j * I + i);
                        sieve_decrease (S_ptr + i, logp, w);
                    }
                    // odd lines only.
                    i0 = ((i0 + 2 * r) & (p - 1)) + 2 * I;
                }
                i0 -= I;
                linestart += I;
                S_ptr += I;
            }
            ssdpos[i] = i0;
        }
    }
}
/* }}} */

/* {{{ resieving. Different interface, since it plays with buckets as well.
 */

/* Sieve small primes (p < I, p not in trialdiv_primes list) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd.
   Primes in trialdiv_primes must be in increasing order. */
void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t *ssd, int * ssdpos,
        sieve_info_srcptr si, where_am_I_ptr w MAYBE_UNUSED)
{
    const uint32_t I = si->I;
    unsigned char *S_ptr;
    unsigned long j, nj;
    const int resieve_very_verbose = 0, resieve_very_verbose_bad = 0;
    /* See comment above about the variable of the same name */
    const int row0_is_oddj = (N << (LOG_BUCKET_REGION - si->logI)) & 1;

    nj = (bucket_region >> si->logI);

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
                for (unsigned int i = i0 + (q& -!(i0&1)) ; i < I; i += p+q) {
                    if (S_ptr[i] == 255) continue;
                    bucket_prime_t prime;
                    unsigned int x = (j << (si->logI)) + i;
                    if (resieve_very_verbose) {
                        pthread_mutex_lock(&io_mutex);
                        fprintf (stderr, "resieve_small_bucket_region: root "
                                FBPRIME_FORMAT ",%d divides at x = "
                                "%d = %lu * %u + %d\n",
                                p, r, x, j, 1 << si->logI, i);
                        pthread_mutex_unlock(&io_mutex);
                    }
                    prime.p = p;
                    prime.x = x;
                    ASSERT(prime.p >= si->td_thresh);
                    push_bucket_prime (BP, prime);
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
        if (event == SSP_DISCARD) {
            continue;
        }
        if (event == SSP_PROJ) {
            ssp_bad_t * ssp = (ssp_bad_t * ) &(ssd->ssp[i]);
            const fbprime_t g = ssp->g;

            WHERE_AM_I_UPDATE(w, p, ssp->g * ssp->q);

            /* Test every p-th line, starting at S[ssdpos] */
            unsigned int i0 = ssdpos[i];
            unsigned int ii;
            ASSERT (i0 % I == 0); /* make sure ssdpos points at start
                                     of line */
            if (resieve_very_verbose_bad) {
                pthread_mutex_lock(&io_mutex);
                fprintf (stderr, "# resieving bad prime " FBPRIME_FORMAT
                        ", i0 = %u\n", g, i0);
                pthread_mutex_unlock(&io_mutex);
            }
            while (i0 < (unsigned int) bucket_region) {
                unsigned char *S_ptr = S + i0;
                if ((i0 >> si->logI) % 2 == 0) { /* Even j coordinate? */
                    /* Yes, test only odd ii-coordinates */
                    for (ii = 1; ii < I; ii += 2) {
                        if (S_ptr[ii] != 255) {
                            bucket_prime_t prime;
                            const unsigned int x = i0 + ii;
                            if (resieve_very_verbose_bad) {
                                pthread_mutex_lock(&io_mutex);
                                fprintf (stderr, "resieve_small_bucket_region even j: root "
                                        FBPRIME_FORMAT ",inf divides at x = %u\n",
                                        g, x);
                                pthread_mutex_unlock(&io_mutex);
                            }
                            prime.p = g;
                            prime.x = x;
                            ASSERT(prime.p >= si->td_thresh);
                            push_bucket_prime (BP, prime);
                        }
                    }
                } else {
                    /* No, test all ii-coordinates */
                    for (ii = 0; ii < I; ii++) {
                        if (S_ptr[ii] != 255) {
                            bucket_prime_t prime;
                            const unsigned int x = i0 + ii;
                            if (resieve_very_verbose_bad) {
                                pthread_mutex_lock(&io_mutex);
                                fprintf (stderr, "resieve_small_bucket_region odd j: root "
                                        FBPRIME_FORMAT ",inf divides at x = %u\n",
                                        g, x);
                                pthread_mutex_unlock(&io_mutex);
                            }
                            prime.p = g;
                            prime.x = x;
                            ASSERT(prime.p >= si->td_thresh);
                            push_bucket_prime (BP, prime);
                        }
                    }
                }
                i0 += g * I;
            }
            ssdpos[i] = i0 - bucket_region;
            if (resieve_very_verbose_bad) {
                pthread_mutex_lock(&io_mutex);
                fprintf (stderr, "# resieving: new i0 = %u, bucket_region = %d, "
                        "new ssdpos = %d\n",
                        i0, bucket_region, ssdpos[i]);
                pthread_mutex_unlock(&io_mutex);
            }
        }
    }
}
/* }}} */
