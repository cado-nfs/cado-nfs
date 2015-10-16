#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h> /* for PRIx64 macro and strtoumax */
#include <cmath>   // for ceiling, floor in cfrac
#include <ctype.h>
#include <float.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <stdarg.h> /* Required so that GMP defines gmp_vfprintf() */
#include <algorithm>
#include <vector>
#include "threadpool.h"
#include "fb.h"
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "ecm/facul.h"
#include "bucket.h"
#include "trialdiv.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "las-debug.h"
#include "las-duplicate.h"
#include "las-report-stats.h"
#include "las-norms.h"
#include "las-unsieve.h"
#include "las-arith.h"
#include "las-plattice.h"
#include "las-qlattice.h"
#include "las-smallsieve.h"
#include "las-descent-trees.h"
#include "las-cofactor.h"
#include "las-fill-in-buckets.h"
#include "las-threads.h"
#include "las-todo.h"
#include "memusage.h"
#ifdef  DLP_DESCENT
#include "las-dlog-base.h"
#endif

#ifdef HAVE_SSE41
/* #define SSE_SURVIVOR_SEARCH 1 */
#include <smmintrin.h>
#endif

// #define HILIGHT_START   "\e[01;31m"
// #define HILIGHT_END   "\e[00;30m"

#define HILIGHT_START   ""
#define HILIGHT_END   ""

int recursive_descent = 0;
int prepend_relation_time = 0;
int exit_after_rel_found = 0;

double general_grace_time_ratio = DESCENT_DEFAULT_GRACE_TIME_RATIO;


double tt_qstart;

/* {{{ for cofactorization statistics */

int cof_stats = 0; /* write stats file: necessary to generate our strategy.*/
FILE *cof_stats_file;
uint32_t **cof_call; /* cof_call[r][a] is the number of calls of the
                        cofactorization routine with a cofactor of r bits on
                        the rational side, and a bits on the algebraic side */
uint32_t **cof_succ; /* cof_succ[r][a] is the corresponding number of
                        successes, i.e., of call that lead to a relation */
/* }}} */

/*****************************/

/* siever_config stuff */

// FIXME: MNFS
void siever_config_display(siever_config_srcptr sc)/*{{{*/
{
    verbose_output_print(0, 1, "# Sieving parameters for q~2^%d on the %s side\n",
            sc->bitsize, sidenames[sc->side]);
    /* Strive to keep these output lines untouched */
    verbose_output_print(0, 1,
	    "# Sieving parameters: lim0=%lu lim1=%lu lpb0=%d lpb1=%d\n",
	    sc->sides[0]->lim, sc->sides[1]->lim,
	    sc->sides[0]->lpb, sc->sides[1]->lpb);
    verbose_output_print(0, 1,
	    "#                     mfb0=%d mfb1=%d\n",
	    sc->sides[0]->mfb, sc->sides[1]->mfb);
    if (sc->sides[0]->lambda != 0 || sc->sides[1]->lambda != 0) {
        verbose_output_print(0, 1,
                "#                     lambda0=%1.1f lambda1=%1.1f\n",
            sc->sides[0]->lambda, sc->sides[1]->lambda);
    }
}/*}}}*/

int siever_config_cmp(siever_config_srcptr a, siever_config_srcptr b)/*{{{*/
{
    return memcmp(a, b, sizeof(siever_config));
}/*}}}*/

/* siever_config_match only tells whether this config entry is one we
 * want to use for a given problem size/side */
int siever_config_match(siever_config_srcptr a, siever_config_srcptr b)/*{{{*/
{
    return a->side == b->side && a->bitsize == b->bitsize;
}/*}}}*/



/* sieve_info stuff */
static void sieve_info_init_trialdiv(sieve_info_ptr si, int side)/*{{{*/
{
    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    /* XXX This function consider the full contents of the list
     * si->sides[side]->fb to be trialdiv primes, therefore that
     * sieve_info_split_bucket_fb_for_threads must have already been
     * called, so that the only primes which are still in s->fb are the
     * trial-divided ones */
    sieve_side_info_ptr s = si->sides[side];
    unsigned long pmax = MIN((unsigned long) si->conf->bucket_thresh,
                             trialdiv_get_max_p());
    std::vector<unsigned long> *trialdiv_primes = new std::vector<unsigned long>;
    s->fb->extract_bycost(*trialdiv_primes, pmax, si->conf->td_thresh);
    std::sort(trialdiv_primes->begin(), trialdiv_primes->end());
    size_t n = trialdiv_primes->size();
    int skip2 = n > 0 && (*trialdiv_primes)[0] == 2;
    s->trialdiv_data = trialdiv_init (&trialdiv_primes->front() + skip2, n - skip2);
    delete trialdiv_primes;
}/*}}}*/

static void sieve_info_clear_trialdiv(sieve_info_ptr si, int side)/*{{{*/
{
        trialdiv_clear (si->sides[side]->trialdiv_data);
}/*}}}*/

/* {{{ Factor base handling */
/* {{{ Initialize the factor bases */
void sieve_info_init_factor_bases(las_info_ptr las, sieve_info_ptr si, param_list pl)
{
    double tfb;
    const char * fbcfilename = param_list_lookup_string(pl, "fbc");

    for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
        mpz_poly_ptr pol = las->cpoly->pols[side];
        sieve_side_info_ptr sis = si->sides[side];

        const fbprime_t bk_thresh = si->conf->bucket_thresh;
        fbprime_t bk_thresh1 = si->conf->bucket_thresh1;
        const fbprime_t fbb = si->conf->sides[side]->lim;
        const fbprime_t powlim = si->conf->sides[side]->powlim;
        if (bk_thresh > fbb) {
            fprintf(stderr, "Error: lim is too small compared to bk_thresh\n");
            ASSERT_ALWAYS(0);
        }
        if (bk_thresh1 == 0 || bk_thresh1 > fbb) {
            bk_thresh1 = fbb;
        }
        const fbprime_t thresholds[4] = {bk_thresh, bk_thresh1, fbb, fbb};
        const bool only_general[4]={true, false, false, false};
        sis->fb = new fb_factorbase(thresholds, powlim, only_general);

        if (fbcfilename != NULL) {
            /* Try to read the factor base cache file. If that fails, because
               the file does not exist or is not compatible with our parameters,
               it will be written after we generate the factor bases. */
            verbose_output_print(0, 1, "# Mapping memory image of factor base from file %s\n",
                   fbcfilename);
            if (sis->fb->mmap_fbc(fbcfilename)) {
                verbose_output_print(0, 1, "# Finished mapping memory image of factor base\n");
                continue;
            } else {
                verbose_output_print(0, 1, "# Could not map memory image of factor base\n");
            }
        }

        if (pol->deg > 1) {
            tfb = seconds ();
            char fbparamname[4];
            snprintf(fbparamname, sizeof(fbparamname), "fb%d", side);
            const char * fbfilename = param_list_lookup_string(pl, fbparamname);
            if (!fbfilename) {
                fprintf(stderr, "Error: factor base file for algebraic side %d is not given\n", side);
                exit(EXIT_FAILURE);
            }
            verbose_output_print(0, 1, "# Reading %s factor base from %s\n", sidenames[side], fbfilename);
            tfb = seconds () - tfb;
            if (!sis->fb->read(fbfilename))
                exit(EXIT_FAILURE);
            verbose_output_print(0, 1,
                    "# Reading %s factor base of %zuMb took %1.1fs\n",
                    sidenames[side],
                    sis->fb->size() >> 20, tfb);
        } else {
            tfb = seconds ();
            double tfb_wct = wct_seconds ();
            sis->fb->make_linear_threadpool ((const mpz_t *) pol->coeff,
                    las->nb_threads);
            tfb = seconds () - tfb;
            tfb_wct = wct_seconds() - tfb_wct;
            verbose_output_print(0, 1,
                    "# Creating rational factor base of %zuMb took %1.1fs (%1.1fs real)\n",
                    sis->fb->size() >> 20, tfb, tfb_wct);
        }

        if (fbcfilename != NULL) {
            verbose_output_print(0, 1, "# Writing memory image of factor base to file %s\n", fbcfilename);
            sis->fb->dump_fbc(fbcfilename);
            verbose_output_print(0, 1, "# Finished writing memory image of factor base\n");
        }
    }
}
/*}}}*/

/* {{{ reordering of the small factor base
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

static size_t
count_roots(const std::vector<fb_general_entry> &v)
{
    size_t count = 0;
    for(std::vector<fb_general_entry>::const_iterator it = v.begin();
        it != v.end(); it++) {
        count += it->nr_roots;
    }
    return count;
}

void reorder_fb(sieve_info_ptr si, int side)
{
    /* We go through all the primes in FB part 0 and sort them into one of
       5 vectors, and some we discard */

    /* alloc the 5 vectors */
    enum {POW2, POW3, TD, RS, REST};
    std::vector<fb_general_entry> *pieces = new std::vector<fb_general_entry>[5];

    fb_part *small_part = si->sides[side]->fb->get_part(0);
    ASSERT_ALWAYS(small_part->is_only_general());
    const std::vector<fb_general_entry> *small_entries = small_part->get_general_vector();

    fbprime_t plim = si->conf->bucket_thresh;
    fbprime_t costlim = si->conf->td_thresh;

    const size_t pattern2_size = sizeof(unsigned long) * 2;
    for(std::vector<fb_general_entry>::const_iterator it = small_entries->begin();
        it != small_entries->end();
        it++)
    {
        /* The extra conditions on powers of 2 and 3 are related to how
         * pattern-sieving is done.
         */
        if (it->p == 2 && it->q <= pattern2_size) {
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
    /* Concatenate the 5 vectors into one, and store the beginning and ending
       index of each part in fb_parts_x */
    std::vector<fb_general_entry> *s = new std::vector<fb_general_entry>;
    /* FIXME: hack to be able to access the struct fb_parts_x entries via
       an index */
    typedef int interval_t[2];
    interval_t *parts_as_array = &si->sides[side]->fb_parts_x->pow2;
    for (int i = 0; i < 5; i++) {
        parts_as_array[i][0] = count_roots(*s);
        std::sort(pieces[i].begin(), pieces[i].end());
        s->insert(s->end(), pieces[i].begin(), pieces[i].end());
        parts_as_array[i][1] = count_roots(*s);
    }
    delete[] pieces;
    si->sides[side]->fb_smallsieved = s;
}



/* {{{ Print some statistics about the factor bases
 * This also fills the field si->sides[*]->max_bucket_fill_ratio, which
 * is used to verify that per-thread allocation for buckets is
 * sufficient.
 * */
void sieve_info_print_fb_statistics(las_info_ptr las MAYBE_UNUSED, sieve_info_ptr si, int side)
{
    sieve_side_info_ptr s = si->sides[side];

    for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++)
    {
        size_t nr_primes, nr_roots;
        double weight;
        s->fb->get_part(i_part)->count_entries(&nr_primes, &nr_roots, &weight);
        if (nr_primes != 0 || weight != 0.) {
            verbose_output_print(0, 1, "# Number of primes in %s factor base part %d = %zu\n",
                    sidenames[side], i_part, nr_primes);
            verbose_output_print(0, 1, "# Number of prime ideals in %s factor base part %d = %zu\n",
                    sidenames[side], i_part, nr_roots);
            verbose_output_print(0, 1, "# Weight of primes in %s factor base part %d = %0.5g\n",
                    sidenames[side], i_part, weight);
        }
        s->max_bucket_fill_ratio[i_part] = weight * 1.07;
    }
}
/*}}}*/

/*}}}*/
/* }}} */

/* {{{ sieve_info_init_... */
static void
sieve_info_init_from_siever_config(las_info_ptr las, sieve_info_ptr si, siever_config_srcptr sc, param_list pl)
{
    memset(si, 0, sizeof(sieve_info));
    /* This is a kludge, really. We don't get the chance to call the ctor
     * and dtor of sieve_info properly here. This ought to be fixed...
     * Someday */
    mpz_init(si->qbasis.q);

    si->cpoly = las->cpoly;
    memcpy(si->conf, sc, sizeof(siever_config));
    ASSERT_ALWAYS(sc->logI > 0);
    si->I = 1 << sc->logI;
    si->J = 1 << (sc->logI - 1);

    /* Allocate memory for transformed polynomials */
    sieve_info_init_norm_data(si);

    /* This function in itself is too expensive if called often.
     * In the descent, where we initialize a sieve_info for each pair
     * (bitsize_of_q, side), we would call it dozens of time.
     * For the moment, we will assume that all along the hint file, the
     * values of I, lim0, lim1 are constant, so that the factor bases
     * are exactly the same and can be shared.
     */
    if (las->sievers == si) {
        verbose_output_print(0, 1, "# bucket_region = %" PRIu64 "\n",
                BUCKET_REGION);
        sieve_info_init_factor_bases(las, si, pl);
        for (int side = 0; side < si->cpoly->nb_polys; side++) {
            sieve_info_print_fb_statistics(las, si, side);
        }
    } else {
        // We are in descent mode, it seems, so let's not duplicate the
        // factor base data.
        // A few sanity checks, first.
        ASSERT_ALWAYS(las->sievers->conf->logI == si->conf->logI);
        ASSERT_ALWAYS(las->sievers->conf->bucket_thresh == si->conf->bucket_thresh);
        ASSERT_ALWAYS(las->sievers->conf->bucket_thresh1 == si->conf->bucket_thresh1);
        // Then, copy relevant data from the first sieve_info
        verbose_output_print(0, 1, "# Do not regenerate factor base data: copy it from first siever\n");
        for (int side = 0; side < si->cpoly->nb_polys; side++) {
            ASSERT_ALWAYS(las->sievers->conf->sides[side]->lim == si->conf->sides[side]->lim);
            ASSERT_ALWAYS(las->sievers->conf->sides[side]->powlim == si->conf->sides[side]->powlim);
            sieve_side_info_ptr sis = si->sides[side];
            sieve_side_info_ptr sis0 = las->sievers->sides[side];
            sis->fb = sis0->fb;
            /* important ! Otherwise we'll have 0, and badness will
             * occur. */
            for (int i = 0; i < FB_MAX_PARTS; i++)
                sis->max_bucket_fill_ratio[i] = sis0->max_bucket_fill_ratio[i];
        }
    }

    // Now that fb have been initialized, we can set the toplevel.
    si->toplevel = MAX(si->sides[0]->fb->get_toplevel(),
            si->sides[1]->fb->get_toplevel());

    /* Initialize the number of buckets */

    /* If LOG_BUCKET_REGION == sc->logI, then one bucket (whose size is the
     * L1 cache size) is actually one line. This changes some assumptions
     * in sieve_small_bucket_region and resieve_small_bucket_region, where
     * we want to differentiate on the parity on j.
     */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= sc->logI);

    /* set the maximal value of the number of buckets (might be less
       for a given special-q if J is smaller) */
    uint32_t XX[FB_MAX_PARTS] = { 0, NB_BUCKETS_2, NB_BUCKETS_3, 0};
    uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
    XX[si->toplevel] = 1 +
      ((((uint64_t)si->J) << si->conf->logI) - UINT64_C(1)) / BRS[si->toplevel];
    for (int i = si->toplevel+1; i < FB_MAX_PARTS; ++i)
        XX[i] = 0;
    // For small Jmax, the number of buckets at toplevel-1 could also
    // be less than the maximum allowed.
    if (si->toplevel > 1 && XX[si->toplevel] == 1) {
        XX[si->toplevel-1] = 1 +
            ((((uint64_t)si->J) << si->conf->logI) - UINT64_C(1))
            / BRS[si->toplevel-1];
        ASSERT_ALWAYS(XX[si->toplevel-1] != 1);
    }
    for (int i = 0; i < FB_MAX_PARTS; ++i) {
        si->nb_buckets_max[i] = XX[i];
        si->nb_buckets[i] = XX[i];
    }

    si->j_div = init_j_div(si->J);
    si->us = init_unsieve_data(si->I);
    si->doing = NULL;


    /* TODO: We may also build a strategy book, given that several
     * strategies will be similar. Presently we spend some time creating
     * each of them for the descent case.
     */
    const char *cofactfilename = param_list_lookup_string (pl, "file-cofact");
    /*
     * We create our strategy book from the file given by
     * 'file-cofact'. Otherwise, we use a default strategy given by
     * the function facul_make_default_strategy ().
     */

    FILE* file = NULL;
    if (cofactfilename != NULL) /* a file was given */
      file = fopen (cofactfilename, "r");
    double time_strat = seconds();
    si->strategies = facul_make_strategies (sc->sides[0]->lim,
					    sc->sides[0]->lpb,
					    sc->sides[0]->mfb,
					    sc->sides[1]->lim,
					    sc->sides[1]->lpb,
					    sc->sides[1]->mfb,
					    sc->sides[0]->ncurves,
					    sc->sides[1]->ncurves,
					    file, 0);
    verbose_output_print(0, 1, "# Building/reading strategies took %1.1fs\n",
            seconds() - time_strat);

    if (si->strategies == NULL)
      {
	fprintf (stderr, "impossible to read %s\n",
		 cofactfilename);
	abort ();
      }
    if (file != NULL)
      fclose (file);

    for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
	/* init_norms (si, side); */ /* only depends on scale, logmax, lognorm_table */
	sieve_info_init_trialdiv(si, side); /* Init refactoring stuff */
	mpz_init (si->BB[side]);
	mpz_init (si->BBB[side]);
	mpz_init (si->BBBB[side]);
	unsigned long lim = si->conf->sides[side]->lim;
	mpz_ui_pow_ui (si->BB[side], lim, 2);
	mpz_mul_ui (si->BBB[side], si->BB[side], lim);
	mpz_mul_ui (si->BBBB[side], si->BBB[side], lim);

        // This variable is obsolete (but las-duplicates.cpp still reads
        // it) FIXME
        si->sides[side]->strategy = NULL;

        reorder_fb(si, side);
        verbose_output_print(0, 2, "# small %s factor base", sidenames[side]);
        size_t nr_roots;
        si->sides[side]->fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
    }
}
/* }}} */

/* This is an ownership transfer from the current head of the todo list
 * into si->doing.  The field todo->next is not accessed, since it does
 * not make sense for si->doing, which is only a single item. The head of
 * the todo list is pruned */
void sieve_info_pick_todo_item(sieve_info_ptr si, las_todo_stack * todo)
{
    delete si->doing;
    si->doing = todo->top();
    todo->pop();
    ASSERT_ALWAYS(mpz_poly_is_root(si->cpoly->pols[si->doing->side],
                  si->doing->r, si->doing->p));
    /* sanity check */
    if (!mpz_probab_prime_p(si->doing->p, 1)) {
        verbose_output_vfprint(1, 0, gmp_vfprintf, "Error, %Zd is not prime\n",
                               si->doing->p);
        exit(1);
    }
    ASSERT_ALWAYS(si->conf->side == si->doing->side);
}

static void sieve_info_update (sieve_info_ptr si, int nb_threads,
    const size_t nr_workspaces)/*{{{*/
{
  /* essentially update the fij polynomials and J value */
  sieve_info_update_norm_data(si, nb_threads);

  /* update number of buckets at toplevel */
  uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
  si->nb_buckets[si->toplevel] = 1 +
      ((((uint64_t)si->J) << si->conf->logI) - UINT64_C(1)) / 
      BRS[si->toplevel];
  // maybe there is only 1 bucket at toplevel and less than 256 at
  // toplevel-1, due to a tiny J.
  if (si->toplevel > 1 && si->nb_buckets[si->toplevel] == 1) {
        si->nb_buckets[si->toplevel-1] = 1 +
            ((((uint64_t)si->J) << si->conf->logI) - UINT64_C(1)) /
            BRS[si->toplevel-1]; 
        // we forbid skipping two levels.
        ASSERT_ALWAYS(si->nb_buckets[si->toplevel-1] != 1);
  }

  /* Update the slices of the factor base according to new log base */
  for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
      /* The safety factor controls by how much a single slice should fill a
         bucket array at most, i.e., with .5, a single slice should never fill
         a bucket array more than half-way. */
      const double safety_factor = .5;
      sieve_side_info_ptr sis = si->sides[side];
      double max_weight[FB_MAX_PARTS];
      for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++) {
        max_weight[i_part] = sis->max_bucket_fill_ratio[i_part] / nr_workspaces
            * safety_factor;
      }
      sis->fb->make_slices(sis->scale * LOG_SCALE, max_weight);
  }

}/*}}}*/

static void sieve_info_clear (las_info_ptr las, sieve_info_ptr si)/*{{{*/
{
    clear_unsieve_data(si->us);
    si->us = NULL;
    clear_j_div(si->j_div);
    si->j_div = NULL;

    for(int s = 0 ; s < si->cpoly->nb_polys ; s++) {
        if (si->sides[s]->strategy != NULL)
	    facul_clear_strategy (si->sides[s]->strategy);
	si->sides[s]->strategy = NULL;
	sieve_info_clear_trialdiv(si, s);
        delete si->sides[s]->fb_smallsieved;
        si->sides[s]->fb_smallsieved = NULL;
        if (si == las->sievers) {
            delete si->sides[s]->fb;
        }

        mpz_clear (si->BB[s]);
        mpz_clear (si->BBB[s]);
        mpz_clear (si->BBBB[s]);
    }
    delete si->doing;
    facul_clear_strategies (si->strategies);
    si->strategies = NULL;
    sieve_info_clear_norm_data(si);
    /* This is a kludge, really. We don't get the chance to call the ctor
     * and dtor of sieve_info properly here. This ought to be fixed...
     * Someday */
    mpz_clear(si->qbasis.q);
}/*}}}*/

/* las_info stuff */

#ifdef  DLP_DESCENT
static void las_info_init_hint_table(las_info_ptr las, param_list pl)/*{{{*/
{
    const char * filename = param_list_lookup_string(pl, "descent-hint-table");
    if (filename == NULL) return;
    char line[1024];
    unsigned int hint_alloc = 0;
    unsigned int hint_size = 0;
    FILE * f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        /* There's no point in proceeding, since it would really change
         * the behaviour of the program to do so */
        exit(1);
    }
    las->max_hint_bitsize[0] = 0;
    las->max_hint_bitsize[1] = 0;
    for(;;) {
        char * x = fgets(line, sizeof(line), f);
        double t;
        unsigned long z;
        /* Tolerate comments and blank lines */
        if (x == NULL) break;
        for( ; *x && isspace(*x) ; x++) ;
        if (*x == '#') continue;
        if (!*x) continue;

        /* We have a new entry to parse */
        if (hint_size >= hint_alloc) {
            hint_alloc = 2 * hint_alloc + 8;
            las->hint_table = (descent_hint *) realloc(las->hint_table, hint_alloc * sizeof(descent_hint));
        }

        descent_hint_ptr h = las->hint_table[hint_size++];
        siever_config_ptr sc = h->conf;

        z = strtoul(x, &x, 10);
        ASSERT_ALWAYS(z > 0);
        sc->bitsize = z;
        switch(*x++) {
            case '@' :
                       sc->side = strtoul(x, &x, 0);
                       ASSERT_ALWAYS(sc->side < 2);
                       break;
            default:
                       fprintf(stderr, "%s: parse error at %s\n", filename, line);
                       exit(1);
        }
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h->expected_time = t;
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h->expected_success = t;
        for( ; *x && isspace(*x) ; x++) ;
        for( ; *x && !isdigit(*x) ; x++) ;
        z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
        sc->logI = z;
        
        for(int s = 0 ; s < 2 ; s++) {
            for( ; *x && isspace(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc->sides[s]->lim = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc->sides[s]->lpb = z;
            /* recognize this as a double. If it's < 10, we'll consider
             * this means lambda */
            {
                for( ; *x && !isdigit(*x) ; x++) ;
                double t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                if (t < 10) {
                    sc->sides[s]->lambda = t;
                    sc->sides[s]->mfb = t * sc->sides[s]->lpb;
                    /* Then no "lambda" is allowed */
                    continue;
                } else {
                    sc->sides[s]->mfb = t;
                }
            }
            if (*x == ',') {
                for( ; *x && !isdigit(*x) ; x++) ;
                t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                sc->sides[s]->lambda = t;
            } else {
                /* this means "automatic" */
                sc->sides[s]->lambda = 0;
            }
        }
        for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));

        if (sc->bitsize > las->max_hint_bitsize[sc->side])
            las->max_hint_bitsize[sc->side] = sc->bitsize;
        // Copy default value for non-given parameters
        // sc->skewness = las->default_config->skewness;
        sc->bucket_thresh = las->config_base->bucket_thresh;
        sc->bucket_thresh1 = las->config_base->bucket_thresh1;
        sc->td_thresh = las->config_base->td_thresh;
        sc->unsieve_thresh = las->config_base->unsieve_thresh;
        for(int side = 0 ; side < las->cpoly->nb_polys ; side++) {
            sc->sides[side]->powlim = las->config_base->sides[side]->powlim;
            sc->sides[side]->ncurves = las->config_base->sides[side]->ncurves;
        }
    }
    if (las->hint_table == NULL) {
        fprintf(stderr, "%s: no data ??\n", filename);
        exit(1);
    }
    if (hint_size >= hint_alloc) {
        hint_alloc = 2 * hint_alloc + 8;
        las->hint_table = (descent_hint *) realloc(las->hint_table, hint_alloc * sizeof(descent_hint));
    }
    las->hint_table[hint_size++]->conf->bitsize = 0;
    las->hint_table = (descent_hint *) realloc(las->hint_table, hint_size * sizeof(descent_hint));

    /* Allocate the quick lookup tables */

    for(int s = 0 ; s < 2 ; s++) {
        unsigned int n = las->max_hint_bitsize[s] + 1;
        las->hint_lookups[s] = (int *)malloc(n * sizeof(int));
        for(unsigned int i = 0 ; i < n ; i++) {
            las->hint_lookups[s][i] = -1;
        }
    }
    for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
        siever_config_ptr sc = h->conf;
        int d = las->hint_lookups[sc->side][sc->bitsize];
        if (d >= 0) {
            fprintf(stderr, "Error: two hints found for %d%c\n",
                    sc->bitsize, sidenames[sc->side][0]);
            exit(1);
        }
        las->hint_lookups[sc->side][sc->bitsize] = h - las->hint_table[0];

        /* We must make sure that the default factor base bounds are
         * larger than what we have for all the hint cases, so that it
         * remains reasonable to base our work on the larger factor base
         * (thus doing incomplete sieving).
         */
        /* But now, the .poly file does not contain lim data anymore.
         * This is up to the user to do this check, anyway, because
         * it is not sure that makefb was run with the lim given in the
         * poly file.
         */
        /*
        for(int s = 0 ; s < 2 ; s++) {
            ASSERT_ALWAYS(sc->sides[s]->lim <= las->cpoly->pols[s]->lim);
        }
        */
    }

    fclose(f);
}/*}}}*/
#endif /* DLP_DESCENT */

static void las_info_clear_hint_table(las_info_ptr las)/*{{{*/
{
    if (las->hint_table == NULL) return;
    for(int s = 0 ; s < 2 ; s++) {
        free(las->hint_lookups[s]);
        las->hint_lookups[s] = NULL;
        las->max_hint_bitsize[s] = 0;
    }
    free(las->hint_table);
    las->hint_table = NULL;
}/*}}}*/

static void las_info_init(las_info_ptr las, param_list pl)/*{{{*/
{
    memset(las, 0, sizeof(las_info));
    /* {{{ Parse and install general operational flags */
    las->outputname = param_list_lookup_string(pl, "out");
    /* Init output file */
    las->output = stdout;
    if (las->outputname) {
	if (!(las->output = fopen_maybe_compressed(las->outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", las->outputname);
	    exit(EXIT_FAILURE);
	}
    }
    setvbuf(las->output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */

    las->galois = param_list_lookup_string(pl, "galois");

    las->verbose = param_list_parse_switch(pl, "-v");

    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, las->output, las->verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, las->output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }

#ifdef TRACE_K
    const char *trace_file_name = param_list_lookup_string(pl, "traceout");
    FILE *trace_file = stderr;
    if (trace_file_name != NULL) {
        trace_file = fopen(trace_file_name, "w");
        DIE_ERRNO_DIAG(trace_file == NULL, "fopen", trace_file_name);
    }
    verbose_output_add(TRACE_CHANNEL, trace_file, 1);
#endif

    verbose_interpret_parameters(pl);
    param_list_print_command_line(las->output, pl);
    las_display_config_flags();

    las->suppress_duplicates = param_list_parse_switch(pl, "-dup");
    las->nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "t", &las->nb_threads);
    if (las->nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    
    /* }}} */
    /* {{{ Parse polynomial */
    cado_poly_init(las->cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }

    if (!cado_poly_read(las->cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

#ifdef BATCH
    if ((tmp = param_list_lookup_string (pl, "batch0")) == NULL)
      {
        fprintf (stderr, "Error: -batch0 is missing\n");
        param_list_print_usage (pl, NULL, stderr);
        exit (EXIT_FAILURE);
      }
#endif

    // sc->skewness = las->cpoly->skew;
    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    param_list_parse_double(pl, "skew", &(las->cpoly->skew));

    if (las->cpoly->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	cado_poly_clear(las->cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    /* }}} */

    /* {{{ Parse default siever config (fill all possible fields) */
    {
        siever_config_ptr sc = las->config_base;
        memset(sc, 0, sizeof(siever_config));


        /* The default config is not necessarily a complete bit of
         * information.
         *
         * The field with the bitsize of q is here filled with q0 as found in
         * the command line. Note though that this is only one choice among
         * several possible (q1, or just no default).
         * For the descent, we do not intend to have a default config, thus
         * specifying q0 makes no sense. Likewise, for file-based todo lists,
         * we have no default either, and no siever is configured to provide
         * this ``default'' behaviour (yet, the default bits here are used to
         * pre-fill the config data later on).
         */
        mpz_t q0;
        mpz_init(q0);
        if (param_list_parse_mpz(pl, "q0", q0)) {
            sc->bitsize = mpz_sizeinbase(q0, 2);
        }
        mpz_clear(q0);
        sc->side = 1; // Legacy.
        if (!param_list_parse_int(pl, "sqside", &sc->side)) {
            verbose_output_print(0, 1, "# Warning: sqside not given, "
                    "assuming side 1 for backward compatibility.\n");
        }
        param_list_parse_double(pl, "lambda0", &(sc->sides[0]->lambda));
        param_list_parse_double(pl, "lambda1", &(sc->sides[1]->lambda));
        param_list_parse_double(pl, "lambda2", &(sc->sides[2]->lambda));
        int complete = 1;
        complete &= param_list_parse_int  (pl, "I",    &(sc->logI));
        complete &= param_list_parse_ulong(pl, "lim0", &(sc->sides[0]->lim));
        complete &= param_list_parse_ulong(pl, "lim1", &(sc->sides[1]->lim));
#ifdef DLP_MNFSL
        complete &= param_list_parse_ulong(pl, "lim2", &(sc->sides[2]->lim));
#endif
        if (!complete) {
            /* ok. Now in fact, for the moment we really need these to be
             * specified, because the call to "new fb_interface" of
             * course depends on the factor base limits. For the very
             * reason that presently, we want these to be common across
             * several siever config values in the hint table, we cannot
             * leave default_config half-baked.
             *
             * since bucket_thresh depends on I, we need I too.
             */
            fprintf(stderr, "Error: as long as per-qrange factor bases are not fully supported, we need to know at least the I and lim[01] fields\n");
            exit(EXIT_FAILURE);
        }

        complete &= param_list_parse_int(pl, "lpb0",  &(sc->sides[0]->lpb));
        complete &= param_list_parse_int(pl, "mfb0",  &(sc->sides[0]->mfb));
        complete &= param_list_parse_int(pl, "lpb1",  &(sc->sides[1]->lpb));
        complete &= param_list_parse_int(pl, "mfb1",  &(sc->sides[1]->mfb));
#ifdef DLP_MNFSL
        complete &= param_list_parse_int(pl, "lpb2",  &(sc->sides[2]->lpb));
        complete &= param_list_parse_int(pl, "mfb2",  &(sc->sides[2]->mfb));
#endif
        if (!complete) {
            verbose_output_print(0, 1, "# default siever configuration is incomplete ; required parameters are I, lim[01], lpb[01], mfb[01]\n");

        }


        /* Parse optional siever configuration parameters */
        sc->td_thresh = 1024;	/* default value */
        param_list_parse_uint(pl, "tdthresh", &(sc->td_thresh));

        sc->unsieve_thresh = 100;
        if (param_list_parse_uint(pl, "unsievethresh", &(sc->unsieve_thresh))) {
            verbose_output_print(0, 1, "# Un-sieving primes > %u\n",
                    sc->unsieve_thresh);
        }

        sc->bucket_thresh = 1 << sc->logI;	/* default value */
        sc->bucket_thresh1 = 0;	/* default value */
        /* overrides default only if parameter is given */
        param_list_parse_ulong(pl, "bkthresh", &(sc->bucket_thresh));
        param_list_parse_ulong(pl, "bkthresh1", &(sc->bucket_thresh1));

        const char *powlim_params[NB_POLYS_MAX] = {"powlim0", "powlim1", "powlim2"};
        for (int side = 0; side < las->cpoly->nb_polys; side++) {
            if (!param_list_parse_ulong(pl, powlim_params[side],
                        &sc->sides[side]->powlim)) {
                sc->sides[side]->powlim = sc->bucket_thresh - 1;
                verbose_output_print(0, 1,
                        "# Using default value of %lu for -%s\n",
                        sc->sides[side]->powlim, powlim_params[side]);
            }
        }

        const char *ncurves_params[NB_POLYS_MAX] = {"ncurves0", "ncurves1", "ncurves2"};
        for (int side = 0; side < las->cpoly->nb_polys; side++)
            if (!param_list_parse_int(pl, ncurves_params[side],
                        &sc->sides[side]->ncurves))
                sc->sides[side]->ncurves = -1;

        if (complete)
            las->default_config = sc;
    }
    /* }}} */

    if (!las->default_config &&
            !param_list_lookup_string(pl, "descent-hint-table")) {
        fprintf(stderr,
                "Error: no default config set, and no hint table either\n");
        exit(EXIT_FAILURE);
    }

    //{{{ parse option stats_cofact
    const char *statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (statsfilename != NULL) { /* a file was given */
        siever_config_srcptr sc = las->default_config;
        if (sc == NULL) {
            fprintf(stderr, "Error: option stats-cofact works only "
                    "with a default config\n");
            exit(EXIT_FAILURE);
        } else if (param_list_lookup_string(pl, "descent-hint-table")) {
            verbose_output_print(0, 1, "# Warning: option stats-cofact "
                    "only applies to the default siever config\n");
        }

        cof_stats_file = fopen (statsfilename, "w");
        if (cof_stats_file == NULL)
        {
            fprintf (stderr, "Error, cannot create file %s\n",
                    statsfilename);
            exit (EXIT_FAILURE);
        }
        cof_stats = 1;
        //allocate cof_call and cof_succ
        int mfb0 = sc->sides[0]->mfb;
        int mfb1 = sc->sides[1]->mfb;
        cof_call = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        cof_succ = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        for (int i = 0; i <= mfb0; i++) {
            cof_call[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            cof_succ[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            for (int j = 0; j <= mfb1; j++)
                cof_call[i][j] = cof_succ[i][j] = 0;
        }
    }
    //}}}
    
    /* {{{ Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    las->todo = new las_todo_stack;
    mpz_init(las->todo_q0);
    mpz_init(las->todo_q1);
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        las->todo_list_fd = fopen(filename, "r");
        if (las->todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            cado_poly_clear(las->cpoly);
            param_list_clear(pl);
            exit(EXIT_FAILURE);
        }
    }
    /* }}} */
#ifdef  DLP_DESCENT
    las_info_init_hint_table(las, pl);
#endif
    /* Allocate room for only one sieve_info */
    las->sievers = (sieve_info_ptr) malloc(sizeof(sieve_info));
    memset(las->sievers, 0, sizeof(sieve_info));

    las->nq_max = UINT_MAX;

    gmp_randinit_default(las->rstate);

    if (param_list_parse_uint(pl, "random-sample", &las->nq_max)) {
        las->random_sampling = 1;
        unsigned long seed = 0;
        if (param_list_parse_ulong(pl, "seed", &seed)) {
            gmp_randseed_ui(las->rstate, seed);
        }
    } else {
        if (param_list_parse_uint(pl, "nq", &las->nq_max)) {
            if (param_list_lookup_string(pl, "q1") || param_list_lookup_string(pl, "rho")) {
                fprintf(stderr, "Error: argument -nq is incompatible with -q1 or -rho\n");
                exit(EXIT_FAILURE);
            }
        }
    }

#ifdef  DLP_DESCENT
    las->dlog_base = new las_dlog_base();
    las->dlog_base->lookup_parameters(pl);
    las->dlog_base->read();
#endif

#ifdef BATCH
    cofac_list_init (las->L);
#endif

}/*}}}*/

void las_todo_pop(las_todo_stack * stack)/*{{{*/
{
    delete stack->top();
    stack->pop();
}

int las_todo_pop_closing_brace(las_todo_stack * stack)
{
    if (stack->top()->side >= 0)
        return 0;
    las_todo_pop(stack);
    return 1;
}
/*}}}*/

void las_info_clear(las_info_ptr las)/*{{{*/
{
    gmp_randclear(las->rstate);

    las_info_clear_hint_table(las);

    for(sieve_info_ptr si = las->sievers ; si->conf->bitsize ; si++) {
        sieve_info_clear(las, si);
    }
    free(las->sievers);
    if (las->outputname) {
        fclose_maybe_compressed(las->output, las->outputname);
    }
    verbose_output_clear();
    mpz_clear(las->todo_q0);
    mpz_clear(las->todo_q1);
    if (las->todo_list_fd)
        fclose(las->todo_list_fd);
    cado_poly_clear(las->cpoly);
    // The stack might not be empty, with the option -exit-early
    while (!las->todo->empty()) {
        las_todo_pop(las->todo);
    }
    delete las->todo;
#ifdef DLP_DESCENT
    delete las->dlog_base;
#endif

#ifdef BATCH
    cofac_list_clear (las->L);
#endif
}/*}}}*/

/* Look for an existing sieve_info in las->sievers with configuration matching
   that in sc; if none exists, create one. */
sieve_info_ptr get_sieve_info_from_config(las_info_ptr las, siever_config_srcptr sc, param_list pl)/*{{{*/
{
    int n = 0;
    sieve_info_ptr si;
    for(si = las->sievers ; si->conf->bitsize ; si++, n++) {
        if (siever_config_cmp(si->conf, sc) == 0)
            break;
    }
    if (si->conf->bitsize)
        return si;
    /* We've hit the end marker. Need to add a new config. */
    las->sievers = (sieve_info_ptr) realloc(las->sievers, (n+2) * sizeof(sieve_info));
    si = las->sievers + n;
    verbose_output_print(0, 1, "# Creating new sieve configuration for q~2^%d on the %s side\n",
            sc->bitsize, sidenames[sc->side]);
    sieve_info_init_from_siever_config(las, si, sc, pl);
    memset(si + 1, 0, sizeof(sieve_info));
    siever_config_display(sc);
    return si;
}/*}}}*/

void las_todo_push_withdepth(las_todo_stack * stack, mpz_srcptr p, mpz_srcptr r, int side, int depth, int iteration = 0)/*{{{*/
{
    stack->push(new las_todo_entry(p, r, side, depth, iteration));
}
void las_todo_push(las_todo_stack * stack, mpz_srcptr p, mpz_srcptr r, int side)
{
    las_todo_push_withdepth(stack, p, r, side, 0);
}
void las_todo_push_closing_brace(las_todo_stack * stack, int depth)
{
    stack->push(new las_todo_entry(-1, depth));
}

/*}}}*/


/* Put in r the smallest legitimate special-q value that it at least
   s + diff (note that if s+diff is already legitimate, then r = s+diff
   will result. */
static void
next_legitimate_specialq(mpz_t r, const mpz_t s, const unsigned long diff)
{
    mpz_add_ui(r, s, diff);
    /* At some point in the future, we might want to allow prime-power or 
       composite special-q here. */
    /* mpz_nextprime() returns a prime *greater than* its input argument,
       which we don't always want, so we subtract 1 first. */
    mpz_sub_ui(r, r, 1);
    mpz_nextprime(r, r);
}

static void
parse_command_line_q0_q1(las_todo_stack *stack, mpz_ptr q0, mpz_ptr q1, param_list pl, const int qside)
{
    ASSERT_ALWAYS(param_list_parse_mpz(pl, "q0", q0));
    if (param_list_parse_mpz(pl, "q1", q1)) {
        next_legitimate_specialq(q0, q0, 0);
        return;
    }

    /* We don't have -q1. If we have -rho, we sieve only <q0, rho>. If we
       don't have -rho, we sieve only q0, but all roots of it. If -q0 does
       not give a legitimate special-q value, advance to the next legitimate
       one and print a warning. */
    mpz_t t;
    mpz_init_set(t, q0);
    next_legitimate_specialq(q0, q0, 0);
    /*
    if (mpz_cmp(t, q0) != 0)
        verbose_output_vfprint(1, 0, gmp_vfprintf, "Warning: fixing q=%Zd to next prime q=%Zd\n", t, q0);
        */

    mpz_set(q1, q0);
    if (param_list_parse_mpz(pl, "rho", t)) {
        las_todo_push(stack, q0, t, qside);
        /* Set empty interval [q0 + 1, q0] as special-q interval */
        mpz_add_ui (q0, q0, 1);
    } else {
        /* Special-q are chosen from [q, q]. Nothing more to do here. */
    }
    mpz_clear(t);
}

/* Galois automorphisms
   autom2.1: 1/x
   autom2.2: -x
   autom3.1: 1-1/x
   autom3.2: -1-1/x
   autom4.1: -(x+1)/(x-1)
   autom6.1: -(2x+1)/(x-1)
*/
static int
skip_galois_roots(const int orig_nroots, const mpz_t q, mpz_t *roots,
		  const char *galois_autom)
{
    int nroots = orig_nroots;
    if(nroots == 0)
	return 0;
    int ord = 0;
    int A, B, C, D; // x -> (A*x+B)/(C*x+D)
    if(strcmp(galois_autom, "autom2.1") == 0 
       || strcmp(galois_autom, "1/y") == 0){
	ord = 2; A = 0; B = 1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom2.2") == 0
	    || strcmp(galois_autom, "_y") == 0){
	ord = 2; A = -1; B = 0; C = 0; D = 1;
    }
    else if(strcmp(galois_autom, "autom3.1") == 0){
	// 1-1/x = (x-1)/x
	ord = 3; A = 1; B = -1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom3.2") == 0){
	// -1-1/x = (-x-1)/x
	ord = 3; A = -1; B = -1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom4.1") == 0){
	// -(x+1)/(x-1)
	ord = 4; A= -1; B = -1; C = 1; D = -1;
    }
    else if(strcmp(galois_autom, "autom6.1") == 0){
	// -(2*x+1)/(x-1)
	ord = 6; A = -2; B = -1; C = 1; D = -1;
    }
    else{
	fprintf(stderr, "Unknown automorphism: %s\n", galois_autom);
	ASSERT_ALWAYS(0);
    }
    if (nroots % ord) {
        fprintf(stderr, "Number of roots modulo q is not divisible by %d. Don't know how to interpret -galois.\n", ord);
        ASSERT_ALWAYS(0);
    }
    modulusul_t mm;
    unsigned long qq = mpz_get_ui(q);
    modul_initmod_ul(mm, qq);
    // transforming as residues
    residueul_t mat[4];
    for(int i = 0; i < 4; i++)
	modul_init(mat[i], mm);
    // be damned inefficient, but who cares?
    modul_set_int64(mat[0], A, mm);
    modul_set_int64(mat[1], B, mm);
    modul_set_int64(mat[2], C, mm);
    modul_set_int64(mat[3], D, mm);
    // Keep only one root among sigma-orbits.
    residueul_t r1, r2, r3;
    modul_init(r1, mm);
    modul_init(r2, mm);
    modul_init(r3, mm);
    residueul_t conj[ord]; // where to put conjugates
    for(int k = 0; k < ord; k++)
	modul_init(conj[k], mm);
    char used[nroots];     // used roots: non-principal conjugates
    memset(used, 0, nroots);
    for(int k = 0; k < nroots; k++){
	if(used[k]) continue;
	unsigned long rr = mpz_get_ui(roots[k]);
	modul_set_ul(r1, rr, mm);
	// build ord-1 conjugates for roots[k]
	for(int l = 0; l < ord; l++){
	    if(modul_intequal_ul(r1, qq)){
		// FIXME: sigma(oo) = A/C
		ASSERT_ALWAYS(0);
	    }
	    // denominator: C*r1+D
	    modul_mul(r2, mat[2], r1, mm);
	    modul_add(r2, r2, mat[3], mm);
	    if(modul_is0(r2, mm)){
		// FIXME: sigma(r1) = oo
		ASSERT_ALWAYS(0);
	    }
	    modul_inv(r3, r2, mm);
	    // numerator: A*r1+B
	    modul_mul(r1, mat[0], r1, mm);
	    modul_add(r1, r1, mat[1], mm);
	    modul_mul(r1, r3, r1, mm);
	    modul_set(conj[l], r1, mm);
	}
#if 0 // debug. 
	printf("orbit for %lu: %lu", qq, rr);
	for(int l = 0; l < ord-1; l++)
	    printf(" -> %lu", conj[l][0]);
	printf("\n");
#endif
	// check: sigma^ord(r1) should be rr
	ASSERT_ALWAYS(modul_intequal_ul(r1, rr));
	// look at roots
	for(int l = k+1; l < nroots; l++){
	    unsigned long ss = mpz_get_ui(roots[l]);
	    modul_set_ul(r2, ss, mm);
	    for(int i = 0; i < ord-1; i++)
		if(modul_equal(r2, conj[i], mm)){
		    ASSERT_ALWAYS(used[l] == 0);
		    // l is some conjugate, we erase it
		    used[l] = (char)1;
		    break;
		}
	}
    }
    // now, compact roots
    int kk = 0;
    for(int k = 0; k < nroots; k++)
	if(used[k] == 0){
	    if(k > kk)
		mpz_set(roots[kk], roots[k]);
	    kk++;
	}
    ASSERT_ALWAYS(kk == (nroots/ord));
    nroots = kk;
    for(int k = 0; k < ord; k++)
	modul_clear(conj[k], mm);
    for(int i = 0; i < 4; i++)
	modul_clear(mat[i], mm);
    modul_clear(r1, mm);
    modul_clear(r2, mm);
    modul_clear(r3, mm);
    modul_clearmod(mm);
    return nroots;
}

static void adwg(FILE *output, const char *comment, int *cpt,
		 relation &rel, int64_t a, int64_t b){
    if(b < 0) { a = -a; b = -b; }
    rel.a = a; rel.b = (uint64_t)b;
    rel.print(output, comment);
    *cpt += 1;
}

/* removing p^vp from the list of factors in rel. */
static void remove_galois_factors(relation &rel, int p, int vp){
    int ok = 0;

    for(int side = 0 ; side < 2 ; side++){
	for(unsigned int i = 0 ; i < rel.sides[side].size(); i++)
	    if(mpz_cmp_ui(rel.sides[side][i].p, p) == 0){
		ok = 1;
		ASSERT_ALWAYS(rel.sides[side][i].e >= vp);
		rel.sides[side][i].e -= vp;
	    }
    }
    /* indeed, p was present */
    ASSERT_ALWAYS(ok == 1);
}

/* adding p^vp to the list of factors in rel. */
static void add_galois_factors(relation &rel, int p, int vp){
    int ok[2] = {0, 0};

    for(int side = 0 ; side < 2 ; side++){
	for(unsigned int i = 0 ; i < rel.sides[side].size(); i++)
	    if(mpz_cmp_ui(rel.sides[side][i].p, p) == 0){
		ok[side] = 1;
		rel.sides[side][i].e += vp;
	    }
    }
    // FIXME: are we sure this is safe?
    for(int side = 0 ; side < 2 ; side++)
	if(ok[side] == 0)
	    /* we must add p^vp */
	    for(int i = 0; i < vp; i++)
		rel.add(side, p);
}

/* adding relations on the fly in Galois cases */
static void add_relations_with_galois(const char *galois, FILE *output, 
				      const char *comment, int *cpt,
				      relation &rel){
    int64_t a0, b0, a1, b1, a2, b2, a3, b3, a5, b5, aa, bb, a;
    uint64_t b;
    int d;

    a = rel.a; b = rel.b; // should be obsolete one day
    // (a0, b0) = sigma^0((a, b)) = (a, b)
    a0 = rel.a; b0 = (int64_t)rel.b;
    if(strcmp(galois, "autom2.1") == 0)
	// remember, 1/x is for plain autom
	// 1/y is for special Galois: x^4+1 -> DO NOT DUPLICATE RELATIONS!
	// (a-b/x) = 1/x*(-b+a*x)
	adwg(output, comment, cpt, rel, -b0, -a0);
    else if(strcmp(galois, "autom3.1") == 0){
	// modified, not checked
	a1 = -b0; b1 = a0-b0;
	adwg(output, comment, cpt, rel, a1, b1);
	a2 = -b1; b2 = a1-b1;
	adwg(output, comment, cpt, rel, a2, b2);
    }
    else if(strcmp(galois, "autom3.2") == 0){
	b1 = (int64_t)b;
	a2 = b1; b2 = -a-b1;
	adwg(output, comment, cpt, rel, a2, b2);
	a3 = b2; b3 = -a2-b2;
	adwg(output, comment, cpt, rel, a3, b3);
    }
    else if(strcmp(galois, "autom4.1") == 0){
	// FIXME: rewrite and check
	a1 = a; b1 = (int64_t)b;
	// tricky: sig^2((a, b)) = (2b, -2a) ~ (b, -a)
	aa = b1; bb = -a1;
	if(bb < 0){ aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
	// same factorization as for (a, b)
	rel.print(output, comment);
	*cpt += 1;
	// sig((a, b)) = (-(a+b), a-b)
	aa = -(a1+b1);
	bb = a1-b1;
	int am2 = a1 & 1, bm2 = b1 & 1;
	if(am2+bm2 == 1){
	    // (a, b) = (1, 0) or (0, 1) mod 2
	    // aa and bb are odd, aa/bb = 1 mod 2
	    // we must add "2,2" in front of f and g
	    add_galois_factors(rel, 2, 2);
	}
	else{
	    // (a, b) = (1, 1), aa and bb are even
	    // we must remove "2,2" in front of f and g
	    // taken from relation.cpp
	    remove_galois_factors(rel, 2, 2);
	    // remove common powers of 2
	    do {
		aa >>= 1;
		bb >>= 1;
	    } while((aa & 1) == 0 && (bb & 1) == 0);
	}
	if(bb < 0){ aa = -aa; bb = -bb; }
	rel.a = aa; rel.b = (uint64_t)bb;
	rel.print(output, comment);
	*cpt += 1;
	// sig^3((a, b)) = sig((b, -a)) = (a-b, a+b)
	aa = -aa; // FIXME: check!
	if(aa < 0){ aa = -aa; bb = -bb; }
	rel.a = bb; rel.b = (uint64_t)aa;
	rel.print(output, comment);
	*cpt += 1;
    }
    else if(strcmp(galois, "autom6.1") == 0){
	// fact do not change
	adwg(output, comment, cpt, rel, a0 + b0, -a0); // (a2, b2)
	adwg(output, comment, cpt, rel, b0, -(a0+b0)); // (a4, b4)

	// fact do change
        a1 = -(2*a0+b0); b1= a0-b0;
	d = 0;
	while(((a1 % 3) == 0) && ((b1 % 3) == 0)){
	    a1 /= 3;
	    b1 /= 3;
	    d++;
	}
	fprintf(output, "# d1=%d\n", d);
	a3 =-(2*b0+a0); b3 = 2*a0+b0;
	a5 = a0-b0;     b5 = 2*b0+a0;
	if(d == 0)
	    // we need to add 3^3
	    add_galois_factors(rel, 3, 3);
	else
	    // we need to remove 3^3
	    remove_galois_factors(rel, 3, 3);
	adwg(output, comment, cpt, rel, a1, b1); // (a1/3^d, b1/3^d)
	for(int i = 0; i < d; i++){
	    a3 /= 3;
	    b3 /= 3;
	    a5 /= 3;
	    b5 /= 3;
	}
	adwg(output, comment, cpt, rel, a3, b3); // (a3/3^d, b3/3^d)
	adwg(output, comment, cpt, rel, a5, b5); // (a5/3^d, b5/3^d)
    }
}

/* {{{ Populating the todo list */
/* See below in main() for documentation about the q-range and q-list
 * modes */
/* These functions return non-zero if the todo list is not empty */
int las_todo_feed_qrange(las_info_ptr las, param_list pl)
{
    /* If we still have entries in the stack, don't add more now */
    if (!las->todo->empty())
        return 1;

    const unsigned long push_at_least_this_many = 1;

    mpz_ptr q0 = las->todo_q0;
    mpz_ptr q1 = las->todo_q1;

    int qside = las->config_base->side;

    mpz_t * roots;
    mpz_poly_ptr f = las->cpoly->pols[qside];
    roots = (mpz_t *) malloc (f->deg * sizeof (mpz_t));
    for(int i = 0 ; i < f->deg  ; i++) {
        mpz_init(roots[i]);
    }

    if (mpz_cmp_ui(q0, 0) == 0) {
        parse_command_line_q0_q1(las->todo, q0, q1, pl, qside);
        if (las->random_sampling) {
            /* For random sampling, it's important that for all integers in
             * the range [q0, q1[, their nextprime() is within the range, and
             * that at least one such has roots mod f. Make sure that
             * this is the case.
             */
            mpz_t q, q1_orig;
            mpz_init(q);
            mpz_init_set(q1_orig, q1);
            /* we need to know the limit of the q range */
            for(unsigned long i = 1 ; ; i++) {
                mpz_sub_ui(q, q1, i);
                next_legitimate_specialq(q, q, 0);
                if (mpz_cmp(q, q1) >= 0)
                    continue;
                if (mpz_poly_roots (roots, f, q) > 0)
                    break;
                /* small optimization: avoid redoing root finding
                 * several times */
                mpz_set (q1, q);
                i = 1;
            }
            /* now q is the largest prime < q1 with f having roots mod q */
            mpz_add_ui (q1, q, 1);
            /* so now if we pick an integer in [q0, q1[, then its nextprime()
             * will be in [q0, q1_orig[, which is what we look for,
             * really.
             */
            if (mpz_cmp(q0, q1) > 0) {
                gmp_fprintf(stderr, "Error: range [%Zd,%Zd[ contains no prime with roots mod f\n", q0, q1_orig);
                exit(EXIT_FAILURE);
            }
            mpz_clear(q);
            mpz_clear(q1_orig);
        }
    }

    /* Otherwise we're going to process the next few sq's and put them
     * into the list */
    /* The loop processes all special-q in [q, q1]. On loop entry, the value
       in q is required to be a legitimate special-q, and will be added to
       the stack. */
    if (!las->random_sampling) {
        /* handy aliases */
        mpz_ptr q = q0;

        /* If nq_max is specified, then q1 has no effect, even though it
         * has been set equal to q */
        for ( ; las->todo->size() < push_at_least_this_many &&
                (las->nq_max < UINT_MAX || mpz_cmp(q, q1) < 0) &&
                las->nq_pushed < las->nq_max ; )
        {
            int nroots = mpz_poly_roots (roots, f, q);
            if (nroots == 0) {
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# polynomial has no roots for q = %Zu\n", q);
            }

            if (las->galois != NULL)
                nroots = skip_galois_roots(nroots, q, roots, las->galois);

            for(int i = 0 ; i < nroots && las->nq_pushed < las->nq_max; i++) {
                las->nq_pushed++;
                las_todo_push(las->todo, q, roots[nroots-1-i], qside);
            }

            next_legitimate_specialq(q, q, 1);
        }
    } else {
        /* we don't care much about being truly uniform here */
        mpz_t q;
        mpz_init(q);
        for ( ; las->todo->size() < push_at_least_this_many && las->nq_pushed < las->nq_max ; ) {
            mpz_sub(q, q1, q0);
            mpz_urandomm(q, las->rstate, q);
            mpz_add(q, q, q0);
            next_legitimate_specialq(q, q, 0);
            int nroots = mpz_poly_roots (roots, f, q);
            if (!nroots) continue;
            if (las->galois != NULL)
                nroots = skip_galois_roots(nroots, q, roots, las->galois);
            unsigned long i = gmp_urandomm_ui(las->rstate, nroots);
            las->nq_pushed++;
            las_todo_push(las->todo, q, roots[i], qside);
        }
        mpz_clear(q);
    }

    for(int i = 0 ; i < f->deg  ; i++) {
        mpz_clear(roots[i]);
    }
    free(roots);
    return las->todo->size();
}

/* Format of a file with a list of special-q (-todo option):
 *   - Comments are allowed (start line with #)
 *   - Blank lines are ignored
 *   - Each valid line must have the form
 *       s q r
 *     where s is the side (0 or 1) of the special q, and q and r are as usual.
 */
int las_todo_feed_qlist(las_info_ptr las, param_list pl)
{
    if (!las->todo->empty())
        return 1;

    char line[1024];
    FILE * f = las->todo_list_fd;
    /* The fgets call below is blocking, so flush las->output here just to
     * be sure. */
    fflush(las->output);
    char * x;
    for( ; ; ) {
        x = fgets(line, sizeof(line), f);
        /* Tolerate comments and blank lines */
        if (x == NULL) return 0;
        if (*x == '#') continue;
        for( ; *x && isspace(*x) ; x++) ;
        if (!*x) continue;
        break;
    }

    /* We have a new entry to parse */
    mpz_t p, r;
    int side = -1;
    mpz_init(p);
    mpz_init(r);
    int rc;
    switch(*x++) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
                   x--;
                   errno = 0;
                   side = strtoul(x, &x, 0);
                   ASSERT_ALWAYS(!errno);
                   ASSERT_ALWAYS(side < 2);
                   break;
        default:
                   fprintf(stderr, "%s: parse error at %s\n",
                           param_list_lookup_string(pl, "todo"), line);
                   /* We may as well default on the command-line switch */
                   exit(1);
    }

    int nread1 = 0;
    int nread2 = 0;

    mpz_set_ui(r, 0);
    for( ; *x && !isdigit(*x) ; x++) ;
    rc = gmp_sscanf(x, "%Zi%n %Zi%n", p, &nread1, r, &nread2);
    ASSERT_ALWAYS(rc == 1 || rc == 2); /* %n does not count */
    x += (rc==1) ? nread1 : nread2;
    ASSERT_ALWAYS(mpz_probab_prime_p(p, 2));
    {
        mpz_poly_ptr f = las->cpoly->pols[side];
        /* specifying the rational root as <0
         * means that it must be recomputed. Putting 0 does not have this
         * effect, since it is a legitimate value after all.
         */
        if (rc < 2 || (f->deg == 1 && rc == 2 && mpz_cmp_ui(r, 0) < 0)) {
            // For rational side, we can compute the root easily.
            ASSERT_ALWAYS(f->deg == 1);
            int nroots = mpz_poly_roots (&r, f, p);
            ASSERT_ALWAYS(nroots == 1);
        }
    }

    for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));
    las_todo_push(las->todo, p, r, side);
    mpz_clear(p);
    mpz_clear(r);
    return 1;
}


int las_todo_feed(las_info_ptr las, param_list pl)
{
    if (!las->todo->empty())
        return 1;
    if (las->todo_list_fd)
        return las_todo_feed_qlist(las, pl);
    else
        return las_todo_feed_qrange(las, pl);
}
/* }}} */

/* Only when tracing. This function gets called once per
 * special-q only. Here we compute the two norms corresponding to
 * the traced (a,b) pair, and start by dividing out the special-q
 * from the one where it should divide */
void init_trace_k(sieve_info_srcptr si, param_list pl)
{
    struct trace_ab_t ab;
    struct trace_ij_t ij;
    struct trace_Nx_t Nx;
    int have_trace_ab = 0, have_trace_ij = 0, have_trace_Nx = 0;

    const char *abstr = param_list_lookup_string(pl, "traceab");
    if (abstr != NULL) {
        if (sscanf(abstr, "%" SCNd64",%" SCNu64, &ab.a, &ab.b) == 2)
            have_trace_ab = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceab %s\n",
                     abstr);
            exit (EXIT_FAILURE);
        }
    }

    const char *ijstr = param_list_lookup_string(pl, "traceij");
    if (ijstr != NULL) {
        if (sscanf(ijstr, "%d,%u", &ij.i, &ij.j) == 2) {
            have_trace_ij = 1;
        } else {
            fprintf (stderr, "Invalid value for parameter: -traceij %s\n",
                     ijstr);
            exit (EXIT_FAILURE);
        }
    }

    const char *Nxstr = param_list_lookup_string(pl, "traceNx");
    if (Nxstr != NULL) {
        if (sscanf(Nxstr, "%u,%u", &Nx.N, &Nx.x) == 2)
            have_trace_Nx = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceNx %s\n",
                     Nxstr);
            exit (EXIT_FAILURE);
        }
    }

    trace_per_sq_init(si, have_trace_Nx ? &Nx : NULL,
        have_trace_ab ? &ab : NULL, 
        have_trace_ij ? &ij : NULL);
}

/* {{{ apply_buckets */
template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
void
apply_one_update (unsigned char * const S,
        const bucket_update_t<1, HINT> * const u,
        const unsigned char logp, where_am_I_ptr w)
{
  WHERE_AM_I_UPDATE(w, h, u->hint);
  WHERE_AM_I_UPDATE(w, x, u->x);
  sieve_increase(S + (u->x), logp, w);
}

template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
void
apply_one_bucket (unsigned char *S,
        const bucket_array_t<1, HINT> &BA, const int i,
        const fb_part *fb, where_am_I_ptr w)
{
  WHERE_AM_I_UPDATE(w, p, 0);

  for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
    const bucket_update_t<1, HINT> *it = BA.begin(i, i_slice);
    const bucket_update_t<1, HINT> * const it_end = BA.end(i, i_slice);
    const slice_index_t slice_index = BA.get_slice_index(i_slice);
    const unsigned char logp = fb->get_slice(slice_index)->get_logp();

    const bucket_update_t<1, HINT> *next_align;
    if (sizeof(bucket_update_t<1, HINT>) == 4) {
      next_align = (bucket_update_t<1, HINT> *) (((size_t) it + 0x3F) & ~((size_t) 0x3F));
      if (UNLIKELY(next_align > it_end)) next_align = it_end;
    } else {
      next_align = it_end;
    }

    while (it != next_align)
      apply_one_update<HINT> (S, it++, logp, w);

    while (it + 16 <= it_end) {
      uint64_t x0, x1, x2, x3, x4, x5, x6, x7;
      uint16_t x;
#ifdef HAVE_SSE2
      _mm_prefetch(((unsigned char *) it)+256, _MM_HINT_NTA);
#endif
      x0 = ((uint64_t *) it)[0];
      x1 = ((uint64_t *) it)[1];
      x2 = ((uint64_t *) it)[2];
      x3 = ((uint64_t *) it)[3];
      x4 = ((uint64_t *) it)[4];
      x5 = ((uint64_t *) it)[5];
      x6 = ((uint64_t *) it)[6];
      x7 = ((uint64_t *) it)[7];
      it += 16;
      __asm__ __volatile__ (""); /* To be sure all x? are read together in one operation */
#ifdef CADO_LITTLE_ENDIAN
#define INSERT_2_VALUES(X) do {						\
	(X) >>= 16; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32;							\
	WHERE_AM_I_UPDATE(w, x, X); sieve_increase(S + (X), logp, w);	\
      } while (0);
#else
#define INSERT_2_VALUES(X) do {						\
	x = (uint16_t) (X);						\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
      } while (0);

#endif
      INSERT_2_VALUES(x0); INSERT_2_VALUES(x1); INSERT_2_VALUES(x2); INSERT_2_VALUES(x3);
      INSERT_2_VALUES(x4); INSERT_2_VALUES(x5); INSERT_2_VALUES(x6); INSERT_2_VALUES(x7);
    }
    while (it != it_end)
      apply_one_update<HINT> (S, it++, logp, w);
  }
}

// Create the two instances, the longhint_t being specialized.
template 
void apply_one_bucket<shorthint_t> (unsigned char *S,
        const bucket_array_t<1, shorthint_t> &BA, const int i,
        const fb_part *fb, where_am_I_ptr w);

template <>
void apply_one_bucket<longhint_t> (unsigned char *S,
        const bucket_array_t<1, longhint_t> &BA, const int i,
        const fb_part *fb, where_am_I_ptr w) {
  WHERE_AM_I_UPDATE(w, p, 0);

  // There is only one slice.
  slice_index_t i_slice = 0;
  const bucket_update_t<1, longhint_t> *it = BA.begin(i, i_slice);
  const bucket_update_t<1, longhint_t> * const it_end = BA.end(i, i_slice);

  // FIXME: Computing logp for each and every entry seems really, really
  // inefficient. Could we add it to a bucket_update_t of "longhint"
  // type?
  while (it != it_end) {
    slice_index_t index = it->index;
    const unsigned char logp = fb->get_slice(index)->get_logp();
    apply_one_update<longhint_t> (S, it++, logp, w);
  }
}


/* {{{ Trial division */
typedef struct {
    uint64_t *fac;
    int n;
} factor_list_t;

#define FL_MAX_SIZE 200

void factor_list_init(factor_list_t *fl) {
    fl->fac = (uint64_t *) malloc (FL_MAX_SIZE * sizeof(uint64_t));
    ASSERT_ALWAYS(fl->fac != NULL);
    fl->n = 0;
}

void factor_list_clear(factor_list_t *fl) {
    free(fl->fac);
}

static void 
factor_list_add(factor_list_t *fl, const uint64_t p)
{
  ASSERT_ALWAYS(fl->n < FL_MAX_SIZE);
  fl->fac[fl->n++] = p;
}

bool
factor_list_contains(const factor_list_t *fl, const uint64_t p)
{
  for (int i = 0; i < fl->n; i++) {
    if (fl->fac[i] == p)
      return true;
  }
  return false;
}

// print a comma-separated list of factors.
// returns the number of factor printed (in particular, a comma is needed
// after this output only if the return value is non zero)
int factor_list_fprint(FILE *f, factor_list_t fl) {
    int i;
    for (i = 0; i < fl.n; ++i) {
        if (i) fprintf(f, ",");
        fprintf(f, "%" PRIx64, fl.fac[i]);
    }
    return i;
}


static const int bucket_prime_stats = 0;
static long nr_bucket_primes = 0;
static long nr_bucket_longhints = 0;
static long nr_div_tests = 0;
static long nr_composite_tests = 0;
static long nr_wrap_was_composite = 0;
/* The entries in BP must be sorted in order of increasing x */
static void
divide_primes_from_bucket (factor_list_t *fl, mpz_t norm, const unsigned int N, const unsigned int x,
                           bucket_primes_t *BP, const int very_verbose)
{
  while (!BP->is_end()) {
      const bucket_update_t<1, primehint_t> prime = BP->get_next_update();
      if (prime.x > x)
        {
          BP->rewind_by_1();
          break;
        }
      if (prime.x == x) {
          if (bucket_prime_stats) nr_bucket_primes++;
          const unsigned long p = prime.p;
          if (very_verbose) {
              verbose_output_vfprint(0, 1, gmp_vfprintf,
                                     "# N = %u, x = %d, dividing out prime hint p = %lu, norm = %Zd\n",
                                     N, x, p, norm);
          }
          /* If powers of a prime p get bucket-sieved and more than one such
              power hits, then the second (and later) hints will find a
              cofactor that already had all powers of p divided out by the
              loop below that removes multiplicities. Thus, if a prime does
              not divide, we check whether it was divided out before (and thus
              is in the factors list) */
          if (UNLIKELY(!mpz_divisible_ui_p (norm, p))) {
              if(!factor_list_contains(fl, p)) {
                  verbose_output_print(1, 0,
                           "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                           p, N, x);
                  abort();
              } else {
                  verbose_output_print(0, 2,
                           "# Note: p = %lu does not divide at (N,x) = (%u,%d), was divided out before\n",
                           p, N, x);
              }
          } else do {
              /* Remove powers of prime divisors */
              factor_list_add(fl, p);
              mpz_divexact_ui (norm, norm, p);
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}


/* The entries in BP must be sorted in order of increasing x */
static void
divide_hints_from_bucket (factor_list_t *fl, mpz_t norm, const unsigned int N, const unsigned int x,
                          bucket_array_complete *purged, const fb_factorbase *fb, const int very_verbose)
{
  while (!purged->is_end()) {
      const bucket_update_t<1, longhint_t> complete_hint = purged->get_next_update();
      if (complete_hint.x > x)
        {
          purged->rewind_by_1();
          break;
        }
      if (complete_hint.x == x) {
          if (bucket_prime_stats) nr_bucket_longhints++;
          const fb_slice_interface *slice = fb->get_slice(complete_hint.index);
          ASSERT_ALWAYS(slice != NULL);
          const unsigned long p = slice->get_prime(complete_hint.hint);
          if (very_verbose) {
              const unsigned char k = slice->get_k(complete_hint.hint);
              verbose_output_print(0, 1,
                                   "# N = %u, x = %d, dividing out slice hint, "
                                   "index = %lu offset = %lu ",
                                   N, x, (unsigned long) complete_hint.index,
                                   (unsigned long) complete_hint.hint);
              if (slice->is_general()) {
                  verbose_output_print(0, 1, "(general)");
              } else {
                  verbose_output_print(0, 1, "(%d roots)", slice->get_nr_roots());
              }
              verbose_output_vfprint(0, 1, gmp_vfprintf,
                                     ", q = %lu^%hhu, norm = %Zd\n",
                                     p, k, norm);
          }
          if (UNLIKELY(!mpz_divisible_ui_p (norm, p))) {
              if (!factor_list_contains(fl, p)) {
                  verbose_output_print(1, 0,
                           "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                           p, N, x);
                  abort();
              } else {
                  verbose_output_print(0, 2,
                           "# Note: p = %lu does not divide at (N,x) = (%u,%d), was divided out before\n",
                           p, N, x);
              }
          } else do {
              factor_list_add(fl, p);
              mpz_divexact_ui (norm, norm, p);
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}


NOPROFILE_STATIC void
trial_div (factor_list_t *fl, mpz_t norm, const unsigned int N, unsigned int x,
           const bool handle_2, bucket_primes_t *primes,
           bucket_array_complete *purged,
	   trialdiv_divisor_t *trialdiv_data,
           int64_t a, uint64_t b,
           const fb_factorbase *fb)
{
#ifdef TRACE_K
    const int trial_div_very_verbose = trace_on_spot_ab(a,b);
#else
    const int trial_div_very_verbose = 0;
#endif
    int nr_factors;
    fl->n = 0; /* reset factor list */

    if (trial_div_very_verbose) {
        verbose_output_start_batch();
        verbose_output_print(TRACE_CHANNEL, 0, "# trial_div() entry, N = %u, x = %d, a = %" PRId64 ", b = %" PRIu64 ", norm = ", N, x, a, b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "%Zd\n", norm);
    }

    // handle 2 separately, if it is in fb
    if (handle_2) {
        int bit = mpz_scan1(norm, 0);
        int i;
        for (i = 0; i < bit; ++i) {
            fl->fac[fl->n] = 2;
            fl->n++;
        }
        if (trial_div_very_verbose)
            verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# x = %d, dividing out 2^%d, norm = %Zd\n", x, bit, norm);
        mpz_tdiv_q_2exp(norm, norm, bit);
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket (fl, norm, N, x, primes, trial_div_very_verbose);
    divide_hints_from_bucket (fl, norm, N, x, purged, fb, trial_div_very_verbose);
    if (trial_div_very_verbose)
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# x = %d, after dividing out bucket/resieved norm = %Zd\n", x, norm);

    do {
      /* Trial divide primes with precomputed tables */
#define TRIALDIV_MAX_FACTORS 32
      int i;
      unsigned long factors[TRIALDIV_MAX_FACTORS];
      if (trial_div_very_verbose) {
          verbose_output_print(TRACE_CHANNEL, 0, "# Trial division by");
          for (i = 0; trialdiv_data[i].p != 1; i++)
              verbose_output_print(TRACE_CHANNEL, 0, " %lu", trialdiv_data[i].p);
          verbose_output_print(TRACE_CHANNEL, 0, "\n# Factors found: ");
      }

      nr_factors = trialdiv (factors, norm, trialdiv_data, TRIALDIV_MAX_FACTORS);

      for (i = 0; i < MIN(nr_factors, TRIALDIV_MAX_FACTORS); i++)
      {
          if (trial_div_very_verbose)
              verbose_output_print (TRACE_CHANNEL, 0, " %lu", factors[i]);
          factor_list_add (fl, factors[i]);
      }
      if (trial_div_very_verbose) {
          verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "\n# After trialdiv(): norm = %Zd\n", norm);
      }
    } while (nr_factors == TRIALDIV_MAX_FACTORS + 1);

    if (trial_div_very_verbose)
        verbose_output_end_batch();
}
/* }}} */

#ifdef  DLP_DESCENT
/* This returns true only if this descent node is now done, either based
 * on the new relation we have registered, or because the previous
 * relation is better anyway */
bool register_contending_relation(las_info_srcptr las, sieve_info_srcptr si, relation & rel)
{
    if (las->tree->must_avoid(rel)) {
        verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: we have already used this relation, avoiding\n");
        return true;
    }

    /* compute rho for all primes, even on the rational side */
    rel.fixup_r(true);

    descent_tree::candidate_relation contender;
    contender.rel = rel;
    double time_left = 0;

    for(int side = 0 ; side < 2 ; side++) {
        for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
            relation::pr const& v(rel.sides[side][i]);
            if (mpz_cmp(si->doing->p, v.p) == 0)
                continue;
            unsigned long p = mpz_get_ui(v.p);
            if (mpz_fits_ulong_p(v.p)) {
                unsigned long r = mpz_get_ui(v.r);
                if (las->dlog_base->is_known(side, p, r))
                    continue;
            }

            unsigned int n = mpz_sizeinbase(v.p, 2);
            int k = (n <= las->max_hint_bitsize[side] ? las->hint_lookups[side][n] : -1);
            if (k < 0) {
                /* This is not worrysome per se. We just do
                 * not have the info in the descent hint table,
                 * period.
                 */
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: cannot estimate refactoring time for relation involving %d%c (%Zd,%Zd)\n", n, sidenames[side][0], v.p, v.r);
                time_left = INFINITY;
            } else {
                if (std::isfinite(time_left))
                    time_left += las->hint_table[k]->expected_time;
            }
            contender.outstanding.push_back(std::make_pair(side, v));
        }
    }
    verbose_output_print(0, 1, "# [descent] This relation entails an additional time of %.2f for the smoothing process (%zu children)\n",
            time_left, contender.outstanding.size());

    /* when we're re-examining this special-q because of a previous
     * failure, there's absolutely no reason to hurry up on a relation */
    contender.set_time_left(time_left, si->doing->iteration ? INFINITY : general_grace_time_ratio);

    return las->tree->new_candidate_relation(contender);
}
#endif /* DLP_DESCENT */

/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */
// FIXME NOPROFILE_STATIC
int
factor_survivors (thread_data *th, int N, where_am_I_ptr w MAYBE_UNUSED)
{
    las_info_srcptr las = th->las;
    sieve_info_ptr si = th->si;
    las_report_ptr rep = th->rep;
    int cpt = 0;
    int surv = 0, copr = 0;
    mpz_t norm[NB_POLYS_MAX];
    factor_list_t factors[NB_POLYS_MAX];
    mpz_array_t *lps[NB_POLYS_MAX] = { NULL, };
    uint32_array_t *lps_m[NB_POLYS_MAX] = { NULL, }; /* corresponding multiplicities */
    bucket_primes_t primes[2] = {bucket_primes_t(BUCKET_REGION), bucket_primes_t(BUCKET_REGION)};
    bucket_array_complete purged[2] = {bucket_array_complete(BUCKET_REGION), bucket_array_complete(BUCKET_REGION)};
    uint32_t cof_bitsize[2] = {0, 0}; /* placate gcc */
    const unsigned int first_j = N << (LOG_BUCKET_REGION - si->conf->logI);
    const unsigned long nr_lines = 1U << (LOG_BUCKET_REGION - si->conf->logI);
    unsigned char * S[2] = {th->sides[0].bucket_region, th->sides[1].bucket_region};

    for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
        lps[side] = alloc_mpz_array (1);
        lps_m[side] = alloc_uint32_array (1);

        factor_list_init(&factors[side]);
        mpz_init (norm[side]);
    }

#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# When entering factor_survivors for bucket %u, "
                "S[0][%u]=%u, S[1][%u]=%u\n",
                trace_Nx.N, trace_Nx.x, S[0][trace_Nx.x], trace_Nx.x, S[1][trace_Nx.x]);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n",
                traced_norms[0], traced_norms[1]);
    }
#endif  /* }}} */

    if (las->verbose >= 2)
        th->update_checksums();

    /* This is the one which gets the merged information in the end */
    unsigned char * SS = S[0];

#ifdef TRACE_K /* {{{ */
    sieve_side_info_ptr side0 = si->sides[0];
    sieve_side_info_ptr side1 = si->sides[1];
    for (int x = 0; x < 1 << LOG_BUCKET_REGION; x++) {
        if (trace_on_spot_Nx(N, x)) {
            verbose_output_print(TRACE_CHANNEL, 0,
                    "# side0->Bound[%u]=%u, side1t->Bound[%u]=%u\n",
                    S[0][trace_Nx.x], S[0][x] <= side0->bound ? 0 : side0->bound,
                    S[1][trace_Nx.x], S[1][x] <= side0->bound ? 0 : side1->bound);
        }
    }
#endif /* }}} */

    for (unsigned int j = 0; j < nr_lines; j++)
    {
        unsigned char * const both_S[2] = {
            S[0] + (j << si->conf->logI), 
            S[1] + (j << si->conf->logI)
        };
        const unsigned char both_bounds[2] = {
            si->sides[0]->bound,
            si->sides[1]->bound,
        };
        surv += search_survivors_in_line(both_S, both_bounds,
                                         si->conf->logI, j + first_j, N,
                                         si->j_div, si->conf->unsieve_thresh,
                                         si->us);
        /* Make survivor search create a list of x-coordinates that survived
           instead of changing sieve array? More localized accesses in
           purge_bucket() that way */
    }

    /* Copy those bucket entries that belong to sieving survivors and
       store them with the complete prime */
    /* FIXME: choose a sensible size here */

    for(int side = 0 ; side < 2 ; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        // From N we can deduce the bucket_index. They are not the same
        // when there are multiple-level buckets.
        uint32_t bucket_index = N % si->nb_buckets[1];

        const bucket_array_t<1, shorthint_t> *BA =
            th->ws->cbegin_BA<1, shorthint_t>(side);
        const bucket_array_t<1, shorthint_t> * const BA_end =
            th->ws->cend_BA<1, shorthint_t>(side);
        for (; BA != BA_end; BA++)  {
            purged[side].purge(*BA, bucket_index, SS);
        }

        /* Add entries coming from downsorting, if any */
        const bucket_array_t<1, longhint_t> *BAd =
            th->ws->cbegin_BA<1, longhint_t>(side);
        const bucket_array_t<1, longhint_t> * const BAd_end =
            th->ws->cend_BA<1, longhint_t>(side);
        for (; BAd != BAd_end; BAd++)  {
            purged[side].purge(*BAd, bucket_index, SS);
        }

        /* Resieve small primes for this bucket region and store them 
           together with the primes recovered from the bucket updates */
        resieve_small_bucket_region (&primes[side], N, SS,
                th->si->sides[side]->rsd, th->sides[side].rsdpos, si, w);

        /* Sort the entries to avoid O(n^2) complexity when looking for
           primes during trial division */
        purged[side].sort();
        primes[side].sort();
    }

    /* Scan array one long word at a time. If any byte is <255, i.e. if
       the long word is != 0xFFFF...FF, examine the bytes 
       FIXME: We can use SSE to scan 16 bytes at a time, but have to make 
       sure that SS is 16-aligned first, thus currently disabled. */
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
    const size_t together = sizeof(__m128i);
    __m128i ones128 = (__m128i) {-1,-1};
    const __m128i * SS_lw = (const __m128i *)SS;
#else
    const int together = sizeof(unsigned long);
    const unsigned long * SS_lw = (const unsigned long *)SS;
#endif

    for ( ; (unsigned char *) SS_lw < SS + BUCKET_REGION; SS_lw++) {
#ifdef DLP_DESCENT
      if (las->tree->must_take_decision())
	break;
#endif
#ifdef TRACE_K
        size_t trace_offset = (const unsigned char *) SS_lw - SS;
        if ((unsigned int) N == trace_Nx.N && (unsigned int) trace_offset <= trace_Nx.x && 
            (unsigned int) trace_offset + together > trace_Nx.x) {
            verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u has value %u\n",
                    trace_Nx.x, trace_Nx.N, SS[trace_Nx.x]);
        }
#endif
#if defined(HAVE_SSE41) && defined(SSE_SURVIVOR_SEARCH)
        if (LIKELY(_mm_testc_si128(*SS_lw, ones128)))
            continue;
#else
        if (LIKELY(*SS_lw == (unsigned long)(-1L)))
            continue;
#endif
        size_t offset = (const unsigned char *) SS_lw - SS;
        for (size_t x = offset; x < offset + together; ++x) {
            if (SS[x] == 255) continue;

            th->rep->survivor_sizes[S[0][x]][S[1][x]]++;
            
            /* For factor_leftover_norm, we need to pass the information of the
             * sieve bound. If a cofactor is less than the square of the sieve
             * bound, it is necessarily prime. we implement this by keeping the
             * log to base 2 of the sieve limits on each side, and compare the
             * bitsize of the cofactor with their double.
             */
            int64_t a;
            uint64_t b;

            // Compute algebraic and rational norms.
            NxToAB (&a, &b, N, x, si);
#ifdef TRACE_K
            if (trace_on_spot_ab(a, b))
              verbose_output_print(TRACE_CHANNEL, 0, "# about to start cofactorization for (%"
                       PRId64 ",%" PRIu64 ")  %zu %u\n", a, b, x, SS[x]);
#endif
            /* since a,b both even were not sieved, either a or b should
             * be odd. However, exceptionally small norms, even without
             * sieving, may fall below the report bound (see tracker
             * issue #15437). Therefore it is safe to continue here. */
            // ASSERT((a | b) & 1);
            if (UNLIKELY(((a | b) & 1) == 0))
            {
                th->rep->both_even++;
                continue;
#if 0
                verbose_output_print (1, 0,
                        "# Error: a and b both even for N = %d, x = %d,\n"
                        "i = %d, j = %d, a = %ld, b = %lu\n",
                        N, x, ((x + N*BUCKET_REGION) & (si->I - 1))
                        - (si->I >> 1),
                        (x + N*BUCKET_REGION) >> si->conf->logI,
                        (long) a, (unsigned long) b);
                abort();
#endif
            }

            /* Since the q-lattice is exactly those (a, b) with
               a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
            /* FIXME: fast divisibility test here! */
            /* Dec 2014: on a c90, it takes 0.1 % of total sieving time*/
            if (b == 0 || (mpz_cmp_ui(si->doing->p, b) <= 0 && b % mpz_get_ui(si->doing->p) == 0))
                continue;

            copr++;

            int pass = 1;

            int i;
            unsigned int j;
#ifdef DLP_MNFSL
	    bool ismooth[NB_POLYS_MAX];
	    int nbsmooth = 0;
#endif
            for(int side = 0 ; pass && side < si->cpoly->nb_polys ; side++) {
                // Trial divide norm on side 'side'
                /* Compute the norms using the polynomials transformed to 
                   i,j-coordinates. The transformed polynomial on the 
                   special-q side is already divided by q */
                NxToIJ (&i, &j, N, x, si);
		mpz_poly_homogeneous_eval_siui (norm[side], si->sides[side]->fij, i, j);

#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0,
                            gmp_vfprintf, "# start trial division for norm=%Zd ", norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            "on %s side for (%" PRId64 ",%" PRIu64 ")\n", sidenames[side], a, b);
                }
#endif
                verbose_output_print(1, 2, "FIXME %s, line %d\n", __FILE__, __LINE__);
                const bool handle_2 = true; /* FIXME */
                const fb_factorbase *fb = th->si->sides[side]->fb;
                trial_div (&factors[side], norm[side], N, x,
                        handle_2,
                        &primes[side], &purged[side],
                        si->sides[side]->trialdiv_data,
                        a, b, fb);

                pass = check_leftover_norm (norm[side], si, side);
#ifdef TRACE_K
                if (trace_on_spot_ab(a, b)) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                            "# checked leftover norm=%Zd", norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            " on %s side for (%" PRId64 ",%" PRIu64 "): %d\n",
                            sidenames[side], a, b, pass);
                }
#endif
#ifdef DLP_MNFSL
		ismooth[side] = (pass == 0 ? false : true);
		if(pass == 1)
		    nbsmooth++;
		else if(pass == 0){
		    if(side > 0)
			pass = 1; /* humf, what a "trick" to force the loop */
		}
#endif
            }
#ifdef DLP_MNFSL
	    pass = (nbsmooth > 1 ? 1 : 0);
#endif
            if (pass == 0)
		continue;

#ifdef BATCH
	    verbose_output_start_batch ();
	    cofac_list_add ((cofac_list_t*) las->L, a, b, norm[0], norm[1],
			    si->qbasis.q);
	    verbose_output_end_batch ();
	    continue; /* we will deal with all cofactors at the end of las */
#endif
#if 0 /* activate here to create a cofac file for the 'smooth' binary */
	    gmp_printf ("LOG %ld %lu %Zd %Zd %Zd\n", a, b, norm[0], norm[1],
			si->qbasis.q);
#endif

            if (cof_stats == 1)
            {
                cof_bitsize[0] = mpz_sizeinbase (norm[0], 2);
                cof_bitsize[1] = mpz_sizeinbase (norm[1], 2);
		/* no need to use a mutex here: either we use one thread only
		   to compute the cofactorization data and if several threads
		   the order is irrelevant. The only problem that can happen
		   is when two threads increase the value at the same time,
		   and it is increased by 1 instead of 2, but this should
		   happen rarely. */
		cof_call[cof_bitsize[0]][cof_bitsize[1]] ++;
            }
	    rep->ttcof -= microseconds_thread ();
            pass = factor_both_leftover_norms(norm, lps, lps_m, si);
	    rep->ttcof += microseconds_thread ();
#ifdef TRACE_K
            if (trace_on_spot_ab(a, b) && pass == 0) {
              verbose_output_print(TRACE_CHANNEL, 0,
                      "# factor_leftover_norm failed for (%" PRId64 ",%" PRIu64 "), ", a, b);
              verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                      "remains %Zd, %Zd unfactored\n", norm[0], norm[1]);
            }
#endif
            if (pass <= 0) continue; /* a factor was > 2^lpb, or some
                                        factorization was incomplete */

            /* yippee: we found a relation! */

            if (cof_stats == 1) /* learning phase */
                cof_succ[cof_bitsize[0]][cof_bitsize[1]] ++;
	    
            // ASSERT (bin_gcd_int64_safe (a, b) == 1);

            relation rel(a, b);

            /* Note that we explicitly do not bother about storing r in
             * the relations below */
            for (int side = 0; side < 2; side++) {
                for(int i = 0 ; i < factors[side].n ; i++)
                    rel.add(side, factors[side].fac[i], 0);

                for (unsigned int i = 0; i < lps[side]->length; ++i)
                    rel.add(side, lps[side]->data[i], 0);
            }

            rel.add(si->conf->side, si->doing->p, 0);

            rel.compress();

#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                verbose_output_print(TRACE_CHANNEL, 0, "# Relation for (%"
                        PRId64 ",%" PRIu64 ") printed\n", a, b);
            }
#endif
            {
                int do_check = th->las->suppress_duplicates;
                /* note that if we have large primes which don't fit in
                 * an unsigned long, then the duplicate check will
                 * quickly return "no".
                 */
                int is_dup = do_check
                    && relation_is_duplicate(rel, las->nb_threads, si);
                const char *comment = is_dup ? "# DUPE " : "";
                FILE *output;
                if (!is_dup)
                    cpt++;
                verbose_output_start_batch();   /* lock I/O */
                if (prepend_relation_time) {
                    verbose_output_print(0, 1, "(%1.4f) ", seconds() - tt_qstart);
                }
                verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n",
                        i, j, S[0][x], S[1][x]);
                for (size_t i_output = 0;
                     (output = verbose_output_get(0, 0, i_output)) != NULL;
                     i_output++) {
                    rel.print(output, comment);
		    // once filtering is ok for all Galois cases, 
		    // this entire block would have to disappear
		    if(las->galois != NULL)
			// adding relations on the fly in Galois cases
			add_relations_with_galois(las->galois, output, comment,
						  &cpt, rel);
		}
                verbose_output_end_batch();     /* unlock I/O */
            }

            /* Build histogram of lucky S[x] values */
            th->rep->report_sizes[S[0][x]][S[1][x]]++;

#ifdef  DLP_DESCENT
            if (register_contending_relation(las, si, rel))
                break;
#endif  /* DLP_DESCENT */
        }
    }

    verbose_output_print(0, 3, "# There were %d survivors in bucket %d\n", surv, N);
    th->rep->survivors1 += surv;
    th->rep->survivors2 += copr;

    for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
        mpz_clear(norm[side]);
        factor_list_clear(&factors[side]);
        clear_uint32_array (lps_m[side]);
        clear_mpz_array (lps[side]);
    }

    return cpt;
}

/* }}} */

/****************************************************************************/

MAYBE_UNUSED static inline void subusb(unsigned char *S1, unsigned char *S2, ssize_t offset)
{
  int ex = (unsigned int) S1[offset] - (unsigned int) S2[offset];
  if (UNLIKELY(ex < 0)) S1[offset] = 0; else S1[offset] = ex;	     
}

/* S1 = S1 - S2, with "-" in saturated arithmetical,
 * and memset(S2, 0, EndS1-S1).
 */
void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2) {
#ifndef HAVE_SSE2
  ssize_t mysize = EndS1 - S1;
  unsigned char *cS2 = S2;
  while (S1 < EndS1) {
    subusb(S1,S2,0);
    subusb(S1,S2,1);
    subusb(S1,S2,2);
    subusb(S1,S2,3);
    subusb(S1,S2,4);
    subusb(S1,S2,5);
    subusb(S1,S2,6);
    subusb(S1,S2,7);
    S1 += 8; S2 += 8;
  }
  memset(cS2, 0, mysize);
#else
  __m128i *S1i = (__m128i *) S1, *EndS1i = (__m128i *) EndS1, *S2i = (__m128i *) S2,
    z = _mm_setzero_si128();
  while (S1i < EndS1i) {
    __m128i x0, x1, x2, x3;
    __asm__ __volatile__
      ("prefetcht0 0x1000(%0)\n"
       "prefetcht0 0x1000(%1)\n"
       "movdqa (%0),%2\n"
       "movdqa 0x10(%0),%3\n"
       "movdqa 0x20(%0),%4\n"
       "movdqa 0x30(%0),%5\n"
       "psubusb (%1),%2\n"
       "psubusb 0x10(%1),%3\n"
       "psubusb 0x20(%1),%4\n"
       "psubusb 0x30(%1),%5\n"
       "movdqa %6,(%1)\n"
       "movdqa %6,0x10(%1)\n"
       "movdqa %6,0x20(%1)\n"
       "movdqa %6,0x30(%1)\n"
       "movdqa %2,(%0)\n"
       "movdqa %3,0x10(%0)\n"
       "movdqa %4,0x20(%0)\n"
       "movdqa %5,0x30(%0)\n"
       "add $0x40,%0\n"
       "add $0x40,%1\n"
       : "+&r"(S1i), "+&r"(S2i), "=&x"(x0), "=&x"(x1), "=&x"(x2), "=&x"(x3) : "x"(z));
    /* I prefer use ASM than intrinsics to be sure each 4 instructions which
     * use exactly a cache line are together. I'm 99% sure it's not useful...
     * but it's more beautiful :-)
     */
    /*
    __m128i x0, x1, x2, x3;
    _mm_prefetch(S1i + 16, _MM_HINT_T0); _mm_prefetch(S2i + 16, _MM_HINT_T0);
    x0 = _mm_load_si128(S1i + 0);         x1 = _mm_load_si128(S1i + 1);
    x2 = _mm_load_si128(S1i + 2);         x3 = _mm_load_si128(S1i + 3);
    x0 = _mm_subs_epu8(S2i[0], x0);       x1 = _mm_subs_epu8(S2i[1], x1);
    x2 = _mm_subs_epu8(S2i[2], x2);       x3 = _mm_subs_epu8(S2i[3], x3);
    _mm_store_si128(S2i + 0, z);          _mm_store_si128(S1i + 1, z);
    _mm_store_si128(S2i + 2, z);          _mm_store_si128(S1i + 3, z);
    _mm_store_si128(S1i + 0, x0);         _mm_store_si128(S1i + 1, x1);
    _mm_store_si128(S1i + 2, x2);         _mm_store_si128(S1i + 3, x3);
    S1i += 4; S2i += 4;
    */
  }
#endif 
}

/* Move above ? */
/* {{{ process_bucket_region
 * th->id gives the number of the thread: it is supposed to deal with the set
 * of bucket_regions corresponding to that number, ie those that are
 * congruent to id mod nb_thread.
 *
 * The other threads are accessed by combining the thread pointer th and
 * the thread id: the i-th thread is at th - id + i
 */
void * process_bucket_region(thread_data *th)
{
    where_am_I w MAYBE_UNUSED;
    las_info_srcptr las = th->las;
    sieve_info_ptr si = th->si;
    uint32_t first_region0_index = th->first_region0_index;
    if (si->toplevel == 1) {
        first_region0_index = 0;
    }

    WHERE_AM_I_UPDATE(w, si, si);

    las_report_ptr rep = th->rep;

    unsigned char * S[2];

    unsigned int skiprows = (BUCKET_REGION >> si->conf->logI)*(las->nb_threads-1);

    /* This is local to this thread */
    for(int side = 0 ; side < 2 ; side++) {
        thread_side_data &ts = th->sides[side];
        S[side] = ts.bucket_region;
    }

    unsigned char *SS = th->SS;
    memset(SS, 0, BUCKET_REGION);

    /* loop over appropriate set of sieve regions */
    for (uint32_t ii = 0; ii < si->nb_buckets[1]; ii ++)
      {
        uint32_t i = first_region0_index + ii;
        if ((i % las->nb_threads) != (uint32_t)th->id)
            continue;
        WHERE_AM_I_UPDATE(w, N, i);

        if (recursive_descent) {
            /* For the descent mode, we bail out as early as possible. We
             * need to do so in a multithread-compatible way, though.
             * Therefore the following access is mutex-protected within
             * las->tree. */
            if (las->tree->must_take_decision())
                break;
        } else if (exit_after_rel_found) {
            if (rep->reports)
                break;
        }

        for (int side = 0; side < 2; side++)
          {
            WHERE_AM_I_UPDATE(w, side, side);

            sieve_side_info_ptr s = si->sides[side];
            thread_side_data &ts = th->sides[side];
        
            /* Init norms */
            rep->tn[side] -= seconds_thread ();
#ifdef SMART_NORM
	    init_norms_bucket_region(S[side], i, si, side, 1);
#else
	    init_norms_bucket_region(S[side], i, si, side, 0);
#endif
            // Invalidate the first row except (1,0)
            if (side == 0 && i == 0) {
                int pos10 = 1+((si->I)>>1);
                unsigned char n10 = S[side][pos10];
                memset(S[side], 255, si->I);
                S[side][pos10] = n10;
            }
            rep->tn[side] += seconds_thread ();
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              verbose_output_print(TRACE_CHANNEL, 0, "# After %s init_norms_bucket_region, N=%u S[%u]=%u\n",
                       sidenames[side], w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            /* Apply buckets */
            rep->ttbuckets_apply -= seconds_thread();
            const bucket_array_t<1, shorthint_t> *BA =
                th->ws->cbegin_BA<1, shorthint_t>(side);
            const bucket_array_t<1, shorthint_t> * const BA_end =
                th->ws->cend_BA<1, shorthint_t>(side);
            for (; BA != BA_end; BA++)  {
                apply_one_bucket(SS, *BA, ii, ts.fb->get_part(1), w);
            }

            /* Apply downsorted buckets, if necessary. */
            if (si->toplevel > 1) {
                const bucket_array_t<1, longhint_t> *BAd =
                    th->ws->cbegin_BA<1, longhint_t>(side);
                const bucket_array_t<1, longhint_t> * const BAd_end =
                    th->ws->cend_BA<1, longhint_t>(side);
                for (; BAd != BAd_end; BAd++)  {
                    // FIXME: the updates could come from part 3 as well,
                    // not only part 2.
                    ASSERT_ALWAYS(si->toplevel <= 2);
                    apply_one_bucket(SS, *BAd, ii, ts.fb->get_part(2), w);
                }
            }

	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds_thread();

            /* Sieve small primes */
            sieve_small_bucket_region(SS, i, s->ssd, ts.ssdpos, si, side, w);
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w->N))
              verbose_output_print(TRACE_CHANNEL, 0,
                      "# Final value on %s side, N=%u rat_S[%u]=%u\n",
                      sidenames[side], w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
        }

        /* Factor survivors */
        rep->ttf -= seconds_thread ();
        rep->reports += factor_survivors (th, i, w);
        rep->ttf += seconds_thread ();

        /* Reset resieving data */
        for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
            sieve_side_info_ptr s = si->sides[side];
            thread_side_data &ts = th->sides[side];
            small_sieve_skip_stride(s->ssd, ts.ssdpos, skiprows, si);
            int * b = s->fb_parts_x->rs;
            memcpy(ts.rsdpos, ts.ssdpos + b[0], (b[1]-b[0]) * sizeof(int64_t));
        }
    }

    return NULL;
}/*}}}*/

/* {{{ las_report_accumulate_threads_and_display
 * This function does three distinct things.
 *  - accumulates the timing reports for all threads into a collated report
 *  - display the per-sq timing relative to this report, and the given
 *    timing argument (in seconds).
 *  - merge the per-sq report into a global report
 */
void las_report_accumulate_threads_and_display(las_info_ptr las,
        sieve_info_ptr si, las_report_ptr report,
        thread_workspaces *ws, double qt0)
{
    /* Display results for this special q */
    las_report rep;
    las_report_init(rep);
    sieve_checksum checksum_post_sieve[2];
    
    ws->accumulate(rep, checksum_post_sieve);

    verbose_output_print(0, 2, "# ");
    /* verbose_output_print(0, 2, "%lu survivors after rational sieve,", rep->survivors0); */
    verbose_output_print(0, 2, "%lu survivors after algebraic sieve, ", rep->survivors1);
    verbose_output_print(0, 2, "coprime: %lu\n", rep->survivors2);
    verbose_output_print(0, 2, "# Checksums over sieve region: after all sieving: %u, %u\n", checksum_post_sieve[0].get_checksum(), checksum_post_sieve[1].get_checksum());
    verbose_output_vfprint(0, 1, gmp_vfprintf, "# %lu relation(s) for %s (%Zd,%Zd)\n", rep->reports, sidenames[si->conf->side], si->doing->p, si->doing->r);
    double qtts = qt0 - rep->tn[0] - rep->tn[1] - rep->ttf;
    if (rep->both_even) {
        verbose_output_print(0, 1, "# Warning: found %lu hits with i,j both even (not a bug, but should be very rare)\n", rep->both_even);
    }
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (dont_print_tally && las->nb_threads > 1) {
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [tally available only in mono-thread]\n", qt0);
    } else {
	//convert ttcof from microseconds to seconds
	rep->ttcof *= 0.000001;
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [norm %1.4f+%1.4f, sieving %1.4f"
	    " (%1.4f + %1.4f + %1.4f),"
            " factor %1.4f (%1.4f + %1.4f)]\n", qt0,
            rep->tn[0],
            rep->tn[1],
            qtts,
            rep->ttbuckets_fill,
            rep->ttbuckets_apply,
	    qtts-rep->ttbuckets_fill-rep->ttbuckets_apply,
	    rep->ttf, rep->ttf - rep->ttcof, rep->ttcof);
    }
    las_report_accumulate(report, rep);
    las_report_clear(rep);
}/*}}}*/


/*************************** main program ************************************/


static void declare_usage(param_list pl)
{
  param_list_usage_header(pl,
  "In the names and in the descriptions of the parameters, below there are often\n"
  "aliases corresponding to the convention that 0 is the rational side and 1\n"
  "is the algebraic side. If the two sides are algebraic, then the word\n"
  "'rational' just means the side number 0. Note also that for a rational\n"
  "side, the factor base is recomputed on the fly (or cached), and there is\n"
  "no need to provide a fb0 parameter.\n"
  );

  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "fb0",   "factor base file on the rational side");
  param_list_decl_usage(pl, "fb1",   "(alias fb) factor base file on the algebraic side");
  param_list_decl_usage(pl, "fb2",   "factor base file on side 2");
  param_list_decl_usage(pl, "fbc",  "factor base cache file");
  param_list_decl_usage(pl, "q0",   "left bound of special-q range");
  param_list_decl_usage(pl, "q1",   "right bound of special-q range");
  param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
  param_list_decl_usage(pl, "v",    "(switch) verbose mode, also prints sieve-area checksums");
  param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
  param_list_decl_usage(pl, "t",   "number of threads to use");
  param_list_decl_usage(pl, "ratq", "[deprecated, alias to --sqside 0] (switch) use rational special-q");
  param_list_decl_usage(pl, "sqside", "put special-q on this side");

  param_list_decl_usage(pl, "I",    "set sieving region to 2^I");
  param_list_decl_usage(pl, "skew", "(alias S) skewness");
  param_list_decl_usage(pl, "lim0", "factor base bound on side 0");
  param_list_decl_usage(pl, "lim1", "factor base bound on side 1");
  param_list_decl_usage(pl, "lim2", "factor base bound on side 2");
  param_list_decl_usage(pl, "lpb0", "set large prime bound on side 0 to 2^lpb0");
  param_list_decl_usage(pl, "lpb1", "set large prime bound on side 1 to 2^lpb1");
  param_list_decl_usage(pl, "lpb2", "set large prime bound on side 2 to 2^lpb1");
  param_list_decl_usage(pl, "mfb0", "set rational cofactor bound on side 0 to 2^mfb0");
  param_list_decl_usage(pl, "mfb1", "set rational cofactor bound on side 1 to 2^mfb1");
  param_list_decl_usage(pl, "mfb2", "set rational cofactor bound on side 2 to 2^mfb2");
  param_list_decl_usage(pl, "lambda0", "lambda value on side 0");
  param_list_decl_usage(pl, "lambda1", "lambda value on side 1");
  param_list_decl_usage(pl, "lambda2", "lambda value on side 2");
  param_list_decl_usage(pl, "powlim0", "limit on powers on side 0");
  param_list_decl_usage(pl, "powlim1", "limit on powers on side 1");
  param_list_decl_usage(pl, "powlim2", "limit on powers on side 2");
  param_list_decl_usage(pl, "ncurves0", "controls number of curves on side 0");
  param_list_decl_usage(pl, "ncurves1", "controls number of curves on side 1");
  param_list_decl_usage(pl, "ncurves2", "controls number of curves on side 2");
  param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
  param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh");
  param_list_decl_usage(pl, "bkthresh1", "2-level bucket-sieve primes p >= bkthresh1");
  param_list_decl_usage(pl, "unsievethresh", "Unsieve all p > unsievethresh where p|gcd(a,b)");

  param_list_decl_usage(pl, "allow-largesq", "(switch) allows large special-q, e.g. for a DL descent");
  param_list_decl_usage(pl, "exit-early", "once a relation has been found, go to next special-q (value==1), or exit (value==2)");
  param_list_decl_usage(pl, "stats-stderr", "(switch) print stats to stderr in addition to stdout/out file");
  param_list_decl_usage(pl, "stats-cofact", "write statistics about the cofactorization step in file xxx");
  param_list_decl_usage(pl, "file-cofact", "provide file with strategies for the cofactorization step");
  param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
  param_list_decl_usage(pl, "prepend-relation-time", "prefix all relation produced with time offset since beginning of special-q processing");
  param_list_decl_usage(pl, "ondemand-siever-config", "(switch) defer initialization of siever precomputed structures (one per special-q side) to time of first actual use");
  param_list_decl_usage(pl, "dup", "(switch) suppress duplicate relations");
  param_list_decl_usage(pl, "galois", "(switch) for reciprocal polynomials, sieve only half of the q's");
#ifdef TRACE_K
  param_list_decl_usage(pl, "traceab", "Relation to trace, in a,b format");
  param_list_decl_usage(pl, "traceij", "Relation to trace, in i,j format");
  param_list_decl_usage(pl, "traceNx", "Relation to trace, in N,x format");
  param_list_decl_usage(pl, "traceout", "Output file for trace output, default: stderr");
#endif
  param_list_decl_usage(pl, "random-sample", "Sample this number of special-q's at random, within the range [q0,q1]");
  param_list_decl_usage(pl, "seed", "Use this seed for the random sampling of special-q's (see random-sample)");
  param_list_decl_usage(pl, "nq", "Process this number of special-q's and stop");
#ifdef  DLP_DESCENT
  param_list_decl_usage(pl, "descent-hint-table", "filename with tuned data for the descent, for each special-q bitsize");
  param_list_decl_usage(pl, "recursive-descent", "descend primes recursively");
  /* given that this option is dangerous, we enable it only for
   * las_descent
   */
  param_list_decl_usage(pl, "never-discard", "Disable the discarding process for special-q's. This is dangerous. See bug #15617");
  param_list_decl_usage(pl, "grace-time-ratio", "Fraction of the estimated further descent time which should be spent processing the current special-q, to find a possibly better relation");
  las_dlog_base::declare_parameter_usage(pl);
#endif /* DLP_DESCENT */
#ifdef BATCH
  param_list_decl_usage (pl, "batch0", "side-0 batch file");
  param_list_decl_usage (pl, "batch1", "side-1 batch file");
#endif
  verbose_decl_usage(pl);
}

int main (int argc0, char *argv0[])/*{{{*/
{
    las_info las;
    double t0, tts, wct;
    unsigned long nr_sq_processed = 0;
    unsigned long nr_sq_discarded = 0;
    int allow_largesq = 0;
    int never_discard = 0;      /* only enabled for las_descent */
    double totJ = 0.0;
    int argc = argc0;
    char **argv = argv0;
    double max_full = 0.;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    param_list pl;
    param_list_init(pl);

    declare_usage(pl);

    /* Passing NULL is allowed here. Find value with
     * param_list_parse_switch later on */
    param_list_configure_switch(pl, "-v", NULL);
    param_list_configure_switch(pl, "-ratq", NULL);
    param_list_configure_switch(pl, "-ondemand-siever-config", NULL);
    param_list_configure_switch(pl, "-allow-largesq", &allow_largesq);
    param_list_configure_switch(pl, "-stats-stderr", NULL);
    param_list_configure_switch(pl, "-prepend-relation-time", &prepend_relation_time);
    param_list_configure_switch(pl, "-dup", NULL);
    //    param_list_configure_switch(pl, "-galois", NULL);
    param_list_configure_alias(pl, "skew", "S");
    // TODO: All these aliases should disappear, one day.
    // This is just legacy.
    param_list_configure_alias(pl, "fb1", "fb");
    param_list_configure_alias(pl, "lim0", "rlim");
    param_list_configure_alias(pl, "lim1", "alim");
    param_list_configure_alias(pl, "lpb0", "lpbr");
    param_list_configure_alias(pl, "lpb1", "lpba");
    param_list_configure_alias(pl, "mfb0", "mfbr");
    param_list_configure_alias(pl, "mfb1", "mfba");
    param_list_configure_alias(pl, "lambda0", "rlambda");
    param_list_configure_alias(pl, "lambda1", "alambda");
    param_list_configure_alias(pl, "powlim0", "rpowlim");
    param_list_configure_alias(pl, "powlim1", "apowlim");
#ifdef  DLP_DESCENT
    param_list_configure_switch(pl, "-recursive-descent", &recursive_descent);
    param_list_configure_switch(pl, "-never-discard", &never_discard);
#endif

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_int(pl, "exit-early", &exit_after_rel_found);
#if DLP_DESCENT
    param_list_parse_double(pl, "grace-time-ratio", &general_grace_time_ratio);
#endif
#ifdef BATCH
    FILE *batch[2] = {NULL, NULL};
    const char *batch0_file, *batch1_file;
    if ((batch0_file = param_list_lookup_string (pl, "batch0")) == NULL)
      {
        fprintf (stderr, "Error: -batch0 is missing\n");
        param_list_print_usage (pl, NULL, stderr);
        exit (EXIT_FAILURE);
      }
    if ((batch1_file = param_list_lookup_string (pl, "batch1")) == NULL)
      {
        fprintf (stderr, "Error: -batch1 is missing\n");
        param_list_print_usage (pl, NULL, stderr);
        exit (EXIT_FAILURE);
      }
#endif

    las_info_init(las, pl);    /* side effects: prints cmdline and flags */

    /* We have the following dependency chain (not sure the account below
     * is exhaustive).
     *
     * q0 -> q0d (double) -> si->sides[*]->{scale,logmax}
     * q0 -> (I, lpb, lambda) for the descent
     * 
     * scale -> logp's in factor base.
     *
     * I -> splittings of the factor base among threads.
     *
     * This is probably enough to justify having separate sieve_info's
     * for the given sizes.
     */

    if (!param_list_parse_switch(pl, "-ondemand-siever-config")) {
        /* Create a default siever instance among las->sievers if needed */
        if (las->default_config)
            get_sieve_info_from_config(las, las->default_config, pl);

        /* Create all siever configurations from the preconfigured hints */
        /* This can also be done dynamically if needed */
        if (las->hint_table) {
            for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
                siever_config_ptr sc = h->conf;
                get_sieve_info_from_config(las, sc, pl);
            }
        }
    }

    tune_las_memset();
    
    /* We allocate one more workspace per kind of bucket than there are
       threads, so that threads have some freedom in avoiding the fullest
       bucket array. With only one thread, no balancing needs to be done,
       so we use only one bucket array. */
    // FIXME: We can't do this. Some parts of the code rely on the fact
    // that nr_workspaces == las->nb_threads. For instance,
    // process_bucket_region uses las->nb_threads while it should
    // sometimes use nr_workspaces.
    const size_t nr_workspaces = las->nb_threads;
    thread_workspaces *workspaces = new thread_workspaces(nr_workspaces, 2, las);

    las_report report;
    las_report_init(report);

    t0 = seconds ();
    wct = wct_seconds();
    verbose_output_print(0, 1, "#\n");

    where_am_I w MAYBE_UNUSED;
    WHERE_AM_I_UPDATE(w, las, las);

    /* This is used only for the descent. It's harmless otherwise. */
    las->tree = new descent_tree();

    /* {{{ Doc on todo list handling
     * The function las_todo_feed behaves in different
     * ways depending on whether we're in q-range mode or in q-list mode.
     *
     * q-range mode: the special-q's to be handled are specified as a
     * range. Then, whenever the las->todo list almost runs out, it is
     * refilled if possible, up to the limit q1 (-q0 & -rho just gives a
     * special case of this).
     *
     * q-list mode: the special-q's to be handled are always read from a
     * file. Therefore each new special-q to be handled incurs a
     * _blocking_ read on the file, until EOF. This mode is also used for
     * the descent, which has the implication that the read occurs if and
     * only if the todo list is empty. }}} */

    /* pop() is achieved by sieve_info_pick_todo_item */
    for( ; las_todo_feed(las, pl) ; ) {
        if (las_todo_pop_closing_brace(las->todo)) {
            las->tree->done_node();
            if (las->tree->depth() == 0) {
                if (recursive_descent) {
                    /* BEGIN TREE / END TREE are for the python script */
                    fprintf(las->output, "# BEGIN TREE\n");
                    las->tree->display_last_tree(las->output);
                    fprintf(las->output, "# END TREE\n");
                }
                las->tree->visited.clear();
            }
            continue;
        }

        siever_config current_config;
        memcpy(current_config, las->config_base, sizeof(siever_config));
        {
            const las_todo_entry * const next_todo = las->todo->top();
            current_config->bitsize = mpz_sizeinbase(next_todo->p, 2);
            current_config->side = next_todo->side;
        }

        /* Do we have a hint table with specifically tuned parameters,
         * well suited to this problem size ? */
        if (las->hint_table) {
            for(descent_hint_ptr h = las->hint_table[0] ; h->conf->bitsize ; h++) {
                siever_config_ptr sc = h->conf;
                if (!siever_config_match(sc, current_config)) continue;
                verbose_output_print(0, 1, "# Using existing sieving parameters from hint list for q~2^%d on the %s side [%d%c]\n", sc->bitsize, sidenames[sc->side], sc->bitsize, sidenames[sc->side][0]);
                memcpy(current_config, sc, sizeof(siever_config));
            }
        }

        {
            const las_todo_entry * const next_todo = las->todo->top();
            if (next_todo->iteration) {
                verbose_output_print(0, 1, "#\n# NOTE: we are re-playing this special-q because of %d previous failed attempt(s)\n#\n", next_todo->iteration);
                /* update sieving parameters here */
                current_config->sides[0]->lpb += next_todo->iteration;
                current_config->sides[1]->lpb += next_todo->iteration;
                /* TODO: do we update the mfb or not ? */
            }
        }

        /* Maybe create a new siever ? */
        sieve_info_ptr si = get_sieve_info_from_config(las, current_config, pl);
        WHERE_AM_I_UPDATE(w, si, si);

        sieve_info_pick_todo_item(si, las->todo);

        las->tree->new_node(si->doing);
        las_todo_push_closing_brace(las->todo, si->doing->depth);

        /* Check whether q is larger than the large prime bound.
         * This can create some problems, for instance in characters.
         * By default, this is not allowed, but the parameter
         * -allow-largesq is a by-pass to this test.
         */
        if (!allow_largesq) {
            if ((int)mpz_sizeinbase(si->doing->p, 2) >
                    si->conf->sides[si->conf->side]->lpb) {
                fprintf(stderr, "ERROR: The special q (%d bits) is larger than the "
                        "large prime bound on %s side (%d bits).\n",
                        (int)mpz_sizeinbase(si->doing->p, 2),
                        sidenames[si->conf->side],
                        si->conf->sides[si->conf->side]->lpb);
                fprintf(stderr, "       You can disable this check with "
                        "the -allow-largesq argument,\n");
                fprintf(stderr, "       It is for instance useful for the "
                        "descent.\n");
                exit(EXIT_FAILURE);
            }
        }

        double qt0 = seconds();
        tt_qstart = seconds();

        if (SkewGauss (si->qbasis, si->doing->p, si->doing->r, si->cpoly->skew) != 0)
            continue;
        si->qbasis.set_q(si->doing->p);

        /* check |a0|, |a1| < 2^31 if we use fb_root_in_qlattice_31bits */
#ifndef SUPPORT_LARGE_Q
        if (si->qbasis.a0 <= INT64_C(-2147483648) || INT64_C(2147483648) <= si->qbasis.a0 ||
            si->qbasis.a1 <= INT64_C(-2147483648) || INT64_C(2147483648) <= si->qbasis.a1)
          {
            fprintf (stderr, "Error, too large special-q, define SUPPORT_LARGE_Q. Skipping this special-q.\n");
            continue;
          }
#endif

        /* FIXME: maybe we can discard some special q's if a1/a0 is too large,
         * see http://www.mersenneforum.org/showthread.php?p=130478 
         *
         * Just for correctness at the moment we're doing it in very
         * extreme cases, see bug 15617
         */
        if (sieve_info_adjust_IJ(si, las->nb_threads) == 0) {
            if (never_discard) {
                si->J = las->nb_threads << (LOG_BUCKET_REGION - si->conf->logI);
            } else {
                nr_sq_discarded++;
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# " HILIGHT_START "Discarding %s q=%Zd; rho=%Zd;" HILIGHT_END,
                                       sidenames[si->conf->side], si->doing->p, si->doing->r);
                verbose_output_print(0, 1, " a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "; raw_J=%u;\n", 
                                     si->qbasis.a0, si->qbasis.b0, si->qbasis.a1, si->qbasis.b1, si->J);
                continue;
            }
        }


        verbose_output_vfprint(0, 1, gmp_vfprintf, "# " HILIGHT_START "Sieving %s q=%Zd; rho=%Zd;" HILIGHT_END,
                               sidenames[si->conf->side], si->doing->p, si->doing->r);

        verbose_output_print(0, 1, " a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 ";",
                             si->qbasis.a0, si->qbasis.b0, si->qbasis.a1, si->qbasis.b1);
        if (si->doing->depth) {
            verbose_output_print(0, 1, " # within descent, currently at depth %d", si->doing->depth);
        }
        verbose_output_print(0, 1, "\n");
        nr_sq_processed ++;

        /* checks the value of J,
         * precompute the skewed polynomials of f(x) and g(x), and also
         * their floating-point versions */
        sieve_info_update (si, las->nb_threads, nr_workspaces);
        totJ += (double) si->J;
        verbose_output_print(0, 2, "# I=%u; J=%u\n", si->I, si->J);
        if (las->verbose >= 2) {
            verbose_output_print (0, 1, "# f_0'(x) = ");
            mpz_poly_fprintf(las->output, si->sides[0]->fij);
            verbose_output_print (0, 1, "# f_1'(x) = ");
            mpz_poly_fprintf(las->output, si->sides[1]->fij);
        }

#ifdef TRACE_K
        init_trace_k(si, pl);
#endif

        /* WARNING. We're filling the report info for thread 0 for
         * ttbuckets_fill, while in fact the cost is over all threads.
         * The problem is always the fact that we don't have a proper
         * per-thread timer. So short of a better solution, this is what
         * we do. True, it's not clean.
         * (we need this data to be caught by
         * las_report_accumulate_threads_and_display further down, hence
         * this hack).
         */

        report->ttbuckets_fill -= seconds();

        /* Allocate buckets */
        workspaces->pickup_si(si);

        thread_pool *pool = new thread_pool(las->nb_threads);

        /* Fill in buckets on both sides at top level */
        fill_in_buckets_both(*pool, *workspaces, si);

        /* Check that buckets are not more than full.
         * Due to templates and si->toplevel being not constant, need a
         * switch. (really?)
         */
        switch(si->toplevel) {
            case 1:
                max_full = std::max(max_full,
                        workspaces->buckets_max_full<1, shorthint_t>());
                break;
            case 2:
                max_full = std::max(max_full,
                        workspaces->buckets_max_full<2, shorthint_t>());
                break;
            case 3:
                max_full = std::max(max_full,
                        workspaces->buckets_max_full<3, shorthint_t>());
                break;
            default:
                ASSERT_ALWAYS(0);
        }
        ASSERT_ALWAYS(max_full <= 1.0 ||
                fprintf (stderr, "max_full=%f, see #14987\n", max_full) == 0);
        
        report->ttbuckets_fill += seconds();

        /* Prepare small sieve and re-sieve */
	/* FIXME: not sure we need to do this for all sides */
        for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
            sieve_side_info_ptr s = si->sides[side];

            small_sieve_init(s->ssd, las, s->fb_smallsieved, si, side);
            small_sieve_info("small sieve", side, s->ssd);

            small_sieve_extract_interval(s->rsd, s->ssd, s->fb_parts_x->rs);
            small_sieve_info("resieve", side, s->rsd);

            // Initialiaze small sieve data at the first region of level 0
            // TODO: multithread this? Probably useless...
            for (int i = 0; i < las->nb_threads; ++i) {
                thread_data * th = &workspaces->thrs[i];
                sieve_side_info_ptr s = si->sides[side];
                thread_side_data &ts = th->sides[side];

                uint32_t my_row0 = (BUCKET_REGION >> si->conf->logI) * th->id;
                ts.ssdpos = small_sieve_start(s->ssd, my_row0, si);
                ts.rsdpos = small_sieve_copy_start(ts.ssdpos,
                        s->fb_parts_x->rs);
            }
        }

        if (si->toplevel == 1) {
            /* Process bucket regions in parallel */
            workspaces->thread_do(&process_bucket_region);
        } else {
            // Prepare plattices at internal levels
            // TODO: this could be multi-threaded
            plattice_x_t max_area = plattice_x_t(si->J)<<si->conf->logI;
            plattice_enumerate_area<1>::value =
                MIN(max_area, plattice_x_t(BUCKET_REGION_2));
            plattice_enumerate_area<2>::value =
                MIN(max_area, plattice_x_t(BUCKET_REGION_3));
            plattice_enumerate_area<3>::value = max_area;
            precomp_plattice_t precomp_plattice;
	    /* FIXME: not needed for MNFSL (?) but for MNFSQ */
            for (int side = 0; side < 2; ++side) {
                for (int level = 1; level < si->toplevel; ++level) {
                    const fb_part * fb = si->sides[side]->fb->get_part(level);
                    const fb_slice_interface *slice;
                    for (slice_index_t slice_index = fb->get_first_slice_index();
                            (slice = fb->get_slice(slice_index)) != NULL; 
                            slice_index++) {  
                        precomp_plattice[side][level].push_back(
                                slice->make_lattice_bases(si->qbasis, si->conf->logI));
                    }
                }
            }


            // Visit the downsorting tree depth-first.
            // If toplevel = 1, then this is just processing all bucket
            // regions.
            uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
            for (uint32_t i = 0; i < si->nb_buckets[si->toplevel]; i++) {
                switch (si->toplevel) {
                    case 2:
                        downsort_tree<1>(i, i*BRS[2]/BRS[1],
                                *workspaces, *pool, si, precomp_plattice);
                        break;
                    case 3:
                        downsort_tree<2>(i, i*BRS[3]/BRS[1],
                                *workspaces, *pool, si, precomp_plattice);
                        break;
                    default:
                        ASSERT_ALWAYS(0);
                }
            }

            // Cleanup precomputed lattice bases.
	    /* FIXME: same comment for MNFSL vs. MNFSQ */
            for(int side = 0 ; side < 2 ; side++) {
                for (int level = 1; level < si->toplevel; ++level) {
                    std::vector<plattices_vector_t*> &V =
                        precomp_plattice[side][level];
                    for (std::vector<plattices_vector_t *>::iterator it =
                            V.begin();
                            it != V.end();
                            it++) {
                        delete *it;
                    }
                }
            }
        }

        delete pool;

        // Cleanup smallsieve data
        for (int i = 0; i < las->nb_threads; ++i) {
            for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
                thread_data * th = &workspaces->thrs[i];
                thread_side_data &ts = th->sides[side];
                free(ts.ssdpos);
                free(ts.rsdpos);
            }
        }


#ifdef  DLP_DESCENT
        descent_tree::candidate_relation const& winner(las->tree->current_best_candidate());
        if (winner) {
            /* Even if not going for recursion, store this as being a
             * winning relation. This is useful for preparing the hint
             * file, and also for the initialization of the descent.
             */
            las->tree->take_decision();
            verbose_output_start_batch();
            FILE * output;
            for (size_t i = 0;
                    (output = verbose_output_get(0, 0, i)) != NULL;
                    i++) {
                winner.rel.print(output, "Taken: ");
            }
            verbose_output_end_batch();
            {
                las_todo_entry const& me(*si->doing);
                unsigned int n = mpz_sizeinbase(me.p, 2);
                verbose_output_start_batch();
                verbose_output_print (0, 1, "# taking path: ");
                for(int i = 0 ; i < me.depth ; i++) {
                    verbose_output_print (0, 1, " ");
                }
                verbose_output_print (0, 1, "%d%c ->", n, sidenames[me.side][0]);
                for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                    int side = winner.outstanding[i].first;
                    relation::pr const& v(winner.outstanding[i].second);
                    unsigned int n = mpz_sizeinbase(v.p, 2);
                    verbose_output_print (0, 1, " %d%c", n, sidenames[side][0]);
                }
                if (winner.outstanding.empty()) {
                    verbose_output_print (0, 1, " done");
                }
                verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", me.side,me.p,me.r);
                verbose_output_end_batch();
            }
            if (recursive_descent) {
                /* reschedule the possibly still missing large primes in the
                 * todo list */
                for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                    int side = winner.outstanding[i].first;
                    relation::pr const& v(winner.outstanding[i].second);
                    unsigned int n = mpz_sizeinbase(v.p, 2);
                    verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] " HILIGHT_START "pushing %s (%Zd,%Zd) [%d%c]" HILIGHT_END " to todo list\n", sidenames[side], v.p, v.r, n, sidenames[side][0]);
                    las_todo_push_withdepth(las->todo, v.p, v.r, side, si->doing->depth + 1);
                }
            }
        } else {
            las_todo_entry const& me(*si->doing);
            las->tree->mark_try_again(me.iteration + 1);
            unsigned int n = mpz_sizeinbase(me.p, 2);
            verbose_output_print (0, 1, "# taking path: %d%c -> loop (#%d)", n, sidenames[me.side][0], me.iteration + 1);
            verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", me.side,me.p,me.r);
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Failed to find a relation for " HILIGHT_START "%s (%Zd,%Zd) [%d%c]" HILIGHT_END " (iteration %d). Putting back to todo list.\n", sidenames[me.side], me.p, me.r, n, sidenames[me.side][0], me.iteration);
            las_todo_push_withdepth(las->todo, me.p, me.r, me.side, me.depth + 1, me.iteration + 1);
        }
#endif  /* DLP_DESCENT */

        /* clear */
        for(int side = 0 ; side < si->cpoly->nb_polys ; side++) {
            small_sieve_clear(si->sides[side]->ssd);
            small_sieve_clear(si->sides[side]->rsd);
        }
        qt0 = seconds() - qt0;
        las_report_accumulate_threads_and_display(las, si, report, workspaces, qt0);

#ifdef TRACE_K
        trace_per_sq_clear(si);
#endif
        if (exit_after_rel_found > 1 && report->reports > 0)
            break;
      } // end of loop over special q ideals.

    if (recursive_descent) {
        verbose_output_print(0, 1, "# Now displaying again the results of all descents\n");
        las->tree->display_all_trees(las->output);
    }
    delete las->tree;

#ifdef BATCH
    unsigned long lim[2] = {las->default_config->sides[0]->lim,
			    las->default_config->sides[1]->lim};
    int lpb[2] = {las->default_config->sides[0]->lpb,
		  las->default_config->sides[1]->lpb};
    int split = 1;
    batch[0] = fopen (batch0_file, "r");
    if (batch[0] == NULL)
      {
	create_batch_file (batch0_file, lim[0], 1UL << lpb[0],
			   las->cpoly->pols[0], split);
	batch[0] = fopen (batch0_file, "r");
      }
    ASSERT_ALWAYS(batch[0] != NULL);
    batch[1] = fopen (batch1_file, "r");
    if (batch[1] == NULL)
      {
	create_batch_file (batch1_file, lim[1], 1UL << lpb[1],
			   las->cpoly->pols[1], split);
	batch[1] = fopen (batch1_file, "r");
      }
    ASSERT_ALWAYS(batch[1] != NULL);
    double tcof_batch = seconds ();
    cofac_list_realloc (las->L, las->L->size);
    verbose_output_print (2, 1, "# Total %lu pairs of cofactors\n",
			  las->L->size);
    report->reports = find_smooth (las->L, lpb, lim, batch, las->verbose);
    verbose_output_print (2, 1, "# Detected %lu smooth cofactors\n",
			  report->reports);
    factor (las->L, report->reports, las->cpoly, lpb[0], lpb[1], las->verbose);
    tcof_batch = seconds () - tcof_batch;
    report->ttcof += tcof_batch;
    /* add to ttf since the remaining time will be computed as ttf - ttcof */
    report->ttf += tcof_batch;
#endif

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;
    verbose_output_print (2, 1, "# Average J=%1.0f for %lu special-q's, max bucket fill %f\n",
            totJ / (double) nr_sq_processed, nr_sq_processed, max_full);
    verbose_output_print (2, 1, "# Discarded %lu special-q's out of %u pushed\n",
            nr_sq_discarded, las->nq_pushed);
    tts = t0;
    tts -= report->tn[0];
    tts -= report->tn[1];
    tts -= report->ttf;

    if (las->verbose)
        facul_print_stats (las->output);

    /*{{{ Display tally */
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (bucket_prime_stats) {
        verbose_output_print(2, 1, "# nr_bucket_primes = %lu, nr_div_tests = %lu, nr_composite_tests = %lu, nr_wrap_was_composite = %lu\n",
                 nr_bucket_primes, nr_div_tests, nr_composite_tests, nr_wrap_was_composite);
    }

    if (dont_print_tally && las->nb_threads > 1) 
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [tally available only in mono-thread]\n", t0);
    else
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [norm %1.2f+%1.1f, sieving %1.1f"
                " (%1.1f + %1.1f + %1.1f),"
                " factor %1.1f (%1.1f + %1.1f)]\n", t0,
                report->tn[0],
                report->tn[1],
                tts,
                report->ttbuckets_fill,
                report->ttbuckets_apply,
                tts-report->ttbuckets_fill-report->ttbuckets_apply,
		report->ttf, report->ttf - report->ttcof, report->ttcof);

    verbose_output_print (2, 1, "# Total elapsed time %1.2fs, per special-q %gs, per relation %gs\n",
                 wct, wct / (double) nr_sq_processed, wct / (double) report->reports);
    const long peakmem = PeakMemusage();
    if (peakmem > 0)
        verbose_output_print (2, 1, "# PeakMemusage (MB) = %ld \n",
                peakmem >> 10);
    verbose_output_print (2, 1, "# Total %lu reports [%1.3gs/r, %1.1fr/sq]\n",
            report->reports, t0 / (double) report->reports,
            (double) report->reports / (double) nr_sq_processed);


    print_worst_weight_errors();

    /*}}}*/

    //{{{ print the stats of the cofactorization.
    if (cof_stats == 1) {
	int mfb0 = las->default_config->sides[0]->mfb;
	int mfb1 = las->default_config->sides[1]->mfb;
	for (int i = 0; i <= mfb0; i++) {
	    for (int j = 0; j <= mfb1; j++)
		fprintf (cof_stats_file, "%u %u %u %u\n", i, j, cof_call[i][j],
			 cof_succ[i][j]);
	    free (cof_call[i]);
	    free (cof_succ[i]);
	}
	free (cof_call);
	free (cof_succ);
	fclose (cof_stats_file);
    }
    //}}}
    delete workspaces;

    las_report_clear(report);

    las_info_clear(las);

    // Clear the ecm addition chains that are stored as global variables
    // to avoid re-computations.
    free_saved_chains();

    param_list_clear(pl);

    return 0;
}/*}}}*/

