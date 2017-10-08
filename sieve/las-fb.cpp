#include "cado.h"
#include <stdarg.h>
#include <gmp.h>
#include "las-types.hpp"
#include "fb.hpp"

/*  Factor base handling */
/* {{{ sieve_info::{init,clear,share}_factor_bases */
void sieve_info::init_factor_bases(param_list_ptr pl)
{
    fb_factorbase *fb[2];

    for(int side = 0 ; side < 2 ; side++) {
        fbprime_t bk_thresh = conf.bucket_thresh;
        fbprime_t bk_thresh1 = conf.bucket_thresh1;
        const fbprime_t fbb = conf.sides[side].lim;
        if (!fbb) {
            fb[side] = NULL;
            continue;
        }
        const fbprime_t powlim = conf.sides[side].powlim;
        if (bk_thresh > fbb) {
            bk_thresh = fbb;
        }
        if (bk_thresh1 == 0 || bk_thresh1 > fbb) {
            bk_thresh1 = fbb;
        }
        const fbprime_t thresholds[4] = {bk_thresh, bk_thresh1, fbb, fbb};
        const bool only_general[4]={true, false, false, false};
        sides[side].fb = std::make_shared<fb_factorbase>(thresholds, powlim, only_general);
        fb[side] = sides[side].fb.get();
    }

    const char * fbcfilename = param_list_lookup_string(pl, "fbc");

    if (fbcfilename != NULL) {
      /* Try to read the factor base cache file. If that fails, because
         the file does not exist or is not compatible with our parameters,
         it will be written after we generate the factor bases. */
      verbose_output_print(0, 1, "# Mapping memory image of factor base from file %s\n",
                           fbcfilename);
      if (fb_mmap_fbc(fb, fbcfilename)) {
        verbose_output_print(0, 1, "# Finished mapping memory image of factor base\n");
        return;
      } else {
        verbose_output_print(0, 1, "# Could not map memory image of factor base\n");
      }
    }

    for(int side = 0 ; side < 2 ; side++) {
        cxx_mpz_poly pol;

        if (!conf.sides[side].lim) continue;

        mpz_poly_set(pol, cpoly->pols[side]);

        if (pol->deg > 1) {
            double tfb = seconds ();
            double tfb_wct = wct_seconds ();
            char fbparamname[4];
            std::string polystring = pol.print_poly("x");
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# Reading side-%d factor base from disk"
                    " for polynomial f%d(x) = %s\n",
                    side, side, polystring.c_str());
            snprintf(fbparamname, sizeof(fbparamname), "fb%d", side);
            const char * fbfilename = param_list_lookup_string(pl, fbparamname);
            if (!fbfilename) {
                fprintf(stderr, "Error: factor base file for side %d is not given\n", side);
                exit(EXIT_FAILURE);
            }
            verbose_output_print(0, 1, "# Reading side-%d factor base from %s\n", side, fbfilename);
            if (!fb[side]->read(fbfilename))
                exit(EXIT_FAILURE);
            tfb = seconds () - tfb;
            tfb_wct = wct_seconds () - tfb_wct;
            verbose_output_print(0, 1,
                    "# Reading side-%d factor base of %zuMb took %1.1fs (%1.1fs real)\n",
                    side, fb[side]->size() >> 20, tfb, tfb_wct);
        } else {
            double tfb = seconds ();
            double tfb_wct = wct_seconds ();
            /* note: we parse again the -t option here -- it gets parsed
             * in the las_info ctor too */
            int nb_threads = 1;		/* default value */
            param_list_parse_int(pl, "t", &nb_threads);
            fb[side]->make_linear_threadpool (pol->coeff, nb_threads);
            tfb = seconds () - tfb;
            tfb_wct = wct_seconds() - tfb_wct;
            verbose_output_print(0, 1,
                    "# Creating side-%d rational factor base of %zuMb took %1.1fs (%1.1fs real)\n",
                    side, fb[side]->size() >> 20, tfb, tfb_wct);
        }
    }

    if (fbcfilename != NULL) {
        verbose_output_print(0, 1, "# Writing memory image of factor base to file %s\n", fbcfilename);
        fb_dump_fbc(fb, fbcfilename);
        verbose_output_print(0, 1, "# Finished writing memory image of factor base\n");
    }

    /* Note that max_bucket_fill_ratio and friends are set from within
     * print_fb_statistics, which is a bit ugly.
     */
}

void sieve_info::share_factor_bases(sieve_info& other)
{
    for(int side = 0 ; side < 2 ; side++) {
        ASSERT_ALWAYS(conf.sides[side].lim == other.conf.sides[side].lim);
        ASSERT_ALWAYS(conf.sides[side].powlim == other.conf.sides[side].powlim);
        sieve_info::side_info & sis(sides[side]);
        sieve_info::side_info & sis0(other.sides[side]);
        sis.fb = sis0.fb;
        for (int i = 0; i < FB_MAX_PARTS; i++)
            sis.max_bucket_fill_ratio[i] = sis0.max_bucket_fill_ratio[i];
    }
}

/* }}} */

/*  reordering of the small factor base
 *
 * We split the small factor base in several non-overlapping, contiguous
 * zones:
 *
 *      - powers of 2 (up until the pattern sieve limit)
 *      - powers of 3 (up until the pattern sieve limit)
 *      - trialdiv primes (not powers)
 *      - resieved primes
 *      (- powers of trialdiv primes)
 *      - rest.
 *
 * Problem: bad primes may in fact be pattern sieved, and we might want
 * to pattern-sieve more than just the ``do it always'' cases where p is
 * below the pattern sieve limit.
 *
 * The answer to this is that such primes are expected to be very very
 * rare, so we don't really bother. If we were to do something, we could
 * imagine setting up a schedule list for projective primes -- e.g. a
 * priority queue. But it feels way overkill.
 *
 * Note that the pre-treatment (splitting the factor base in chunks) can
 * be done once and for all.
 */

static size_t count_roots(const std::vector<fb_general_entry> &v)
{
    size_t count = 0;
    for(std::vector<fb_general_entry>::const_iterator it = v.begin();
        it != v.end(); it++) {
        count += it->nr_roots;
    }
    return count;
}

/* {{{ sieve_info::{init,clear}_fb_smallsieved */
void sieve_info::init_fb_smallsieved(int side)
{
    sieve_info & si(*this);

    if (!si.sides[side].fb) return;

    /* We go through all the primes in FB part 0 and sort them into one of
       5 vectors, and some we discard */

    /* alloc the 6 vectors */
    enum {POW2, POW3, TD, RS, REST, SKIPPED};
    std::vector<fb_general_entry> *pieces = new std::vector<fb_general_entry>[6];

    fb_part *small_part = si.sides[side].fb->get_part(0);
    ASSERT_ALWAYS(small_part->is_only_general());
    const fb_vector<fb_general_entry> *small_entries = small_part->get_general_vector();

    fbprime_t plim = si.conf.bucket_thresh;
    fbprime_t costlim = si.conf.td_thresh;

    const size_t pattern2_size = sizeof(unsigned long) * 2;
    for (fb_vector<fb_general_entry>::const_iterator it = small_entries->begin();
         it != small_entries->end();
         it++)
    {
        /* The extra conditions on powers of 2 and 3 are related to how
         * pattern-sieving is done.
         */
        if (it->q < si.conf.skipped) {
            pieces[SKIPPED].push_back(*it);
        } else if (it->p == 2 && it->q <= pattern2_size) {
            pieces[POW2].push_back(*it);
        } else if (it->q == 3) { /* Currently only q=3 is pattern sieved */
            pieces[POW3].push_back(*it);
        } else if (it->k != 1) {
            /* prime powers always go into "rest" */
            pieces[REST].push_back(*it);
        } else if (it->q <= plim && it->q <= costlim * it->nr_roots) {
            pieces[TD].push_back(*it);
        } else if (it->q <= plim) {
            pieces[RS].push_back(*it);
        } else {
            abort();
        }
    }
    /* Concatenate the 6 vectors into one, and store the beginning and ending
       index of each part in fb_parts_x */
    std::vector<fb_general_entry> *s = new std::vector<fb_general_entry>;
    /* FIXME: hack to be able to access the struct fb_parts_x entries via
       an index */
    typedef int interval_t[2];
    interval_t *parts_as_array = &si.sides[side].fb_parts_x->pow2;
    for (int i = 0; i < 6; i++) {
        parts_as_array[i][0] = count_roots(*s);
        std::sort(pieces[i].begin(), pieces[i].end());
        s->insert(s->end(), pieces[i].begin(), pieces[i].end());
        parts_as_array[i][1] = count_roots(*s);
    }
    delete[] pieces;
    si.sides[side].fb_smallsieved = std::shared_ptr<std::vector<fb_general_entry> >(s);
}

/* }}} */

/* {{{ Print some statistics about the factor bases
 * This also fills the field si.sides[*]->max_bucket_fill_ratio, which
 * is used to verify that per-thread allocation for buckets is
 * sufficient.
 * */
void sieve_info::print_fb_statistics(int side)
{
    side_info & s(sides[side]);
    if (!s.fb) return;
    for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++)
    {
        size_t nr_primes, nr_roots;
        double weight;
        s.fb->get_part(i_part)->count_entries(&nr_primes, &nr_roots, &weight);
        if (nr_primes != 0 || weight != 0.) {
            verbose_output_print(0, 1, "# Number of primes in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_primes);
            verbose_output_print(0, 1, "# Number of prime ideals in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_roots);
            verbose_output_print(0, 1, "# Weight of primes in side-%d factor base part %d = %0.5g\n",
                    side, i_part, weight);
        }
        s.max_bucket_fill_ratio[i_part] = weight * 1.07;
    }
}
/*}}}*/



