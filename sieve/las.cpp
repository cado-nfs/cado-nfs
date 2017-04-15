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
#include "threadpool.hpp"
#include "fb.hpp"
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "relation.hpp"
#include "ecm/facul.h"
#include "bucket.hpp"
#include "trialdiv.h"
#include "las-config.h"
#include "las-types.hpp"
#include "las-coordinates.hpp"
#include "las-debug.hpp"
#include "las-duplicate.hpp"
#include "las-report-stats.hpp"
#include "las-norms.hpp"
#include "las-unsieve.hpp"
#include "las-arith.hpp"
#include "las-plattice.hpp"
#include "las-qlattice.hpp"
#include "las-smallsieve.hpp"
#include "las-descent-trees.hpp"
#include "las-cofactor.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-threads.hpp"
#include "las-todo.hpp"
#include "memusage.h"
#include "tdict.hpp"
#ifdef  DLP_DESCENT
#include "las-dlog-base.hpp"
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
int allow_largesq = 0;
int adjust_strategy = 0;

double general_grace_time_ratio = DESCENT_DEFAULT_GRACE_TIME_RATIO;


double tt_qstart;

/*****************************/

/* siever_config stuff */

struct has_same_config {/*{{{*/
    siever_config const & sc;
    has_same_config(siever_config const & sc) : sc(sc) {}
    bool operator()(siever_config const& o) const { return o == sc; }
    bool operator()(sieve_info const& o) const { return (*this)(o.conf); }
};
/*}}}*/
struct has_same_config_q {/*{{{*/
    siever_config const & sc;
    has_same_config_q(siever_config const & sc) : sc(sc) {}
    bool operator()(siever_config const& o) const { return sc.has_same_config_q(o); }
    bool operator()(sieve_info const& o) const { return (*this)(o.conf); }
};
/*}}}*/
struct has_same_sieving {/*{{{*/
    siever_config const & sc;
    has_same_sieving(siever_config const & sc) : sc(sc) {}
    bool operator()(siever_config const& o) const { return sc.has_same_sieving(o); }
    bool operator()(sieve_info const& o) const { return (*this)(o.conf); }
};
/*}}}*/
struct has_same_fb_parameters {/*{{{*/
    siever_config const & sc;
    has_same_fb_parameters(siever_config const & sc) : sc(sc) {}
    bool operator()(siever_config const& o) const { return sc.has_same_fb_parameters(o); }
    bool operator()(sieve_info const& o) const { return (*this)(o.conf); }
};
/*}}}*/
struct has_same_cofactoring {/*{{{*/
    siever_config const & sc;
    has_same_cofactoring(siever_config const & sc) : sc(sc) {}
    bool operator()(siever_config const& o) const { return sc.has_same_cofactoring(o); }
    bool operator()(sieve_info const& o) const { return (*this)(o.conf); }
};
/*}}}*/
void siever_config_display(siever_config const & sc)/*{{{*/
{
    if (sc.bitsize == 0) return;

    verbose_output_print(0, 1, "# Sieving parameters for q~2^%d on side %d\n",
            sc.bitsize, sc.side);
    /* Strive to keep these output lines untouched */
    verbose_output_print(0, 1,
	    "# Sieving parameters: lim0=%lu lim1=%lu lpb0=%d lpb1=%d\n",
	    sc.sides[0].lim, sc.sides[1].lim,
            sc.sides[0].lpb, sc.sides[1].lpb);
    verbose_output_print(0, 1,
	    "#                     mfb0=%d mfb1=%d\n",
	    sc.sides[0].mfb, sc.sides[1].mfb);
    if (sc.sides[0].lambda != 0 || sc.sides[1].lambda != 0) {
        verbose_output_print(0, 1,
                "#                     lambda0=%1.1f lambda1=%1.1f\n",
            sc.sides[0].lambda, sc.sides[1].lambda);
    }
}/*}}}*/

siever_config las_info::get_config_for_q(las_todo_entry const & doing) const /*{{{*/
{
    // arrange so that we don't have the same header line as the one
    // which prints the q-lattice basis
    verbose_output_vfprint(0, 1, gmp_vfprintf,
                         "#\n"
                         "# "
                         "Now sieving side-%d q=%Zd; rho=%Zd\n",
                         doing.side,
                         (mpz_srcptr) doing.p,
                         (mpz_srcptr) doing.r);
    siever_config config = config_base;
    config.bitsize = mpz_sizeinbase(doing.p, 2);
    config.side = doing.side;

    /* Do we have a hint table with specifically tuned parameters,
     * well suited to this problem size ? */
    for(unsigned int i = 0 ; i < hint_table.size() ; i++) {
        siever_config const & sc(hint_table[i].conf);
        if (!sc.has_same_config_q(config)) continue;
        verbose_output_print(0, 1, "# Using parameters from hint list for q~2^%d on side %d [%d@%d]\n", sc.bitsize, sc.side, sc.bitsize, sc.side);
        config = sc;
    }

    /* Check whether q is larger than the large prime bound.
     * This can create some problems, for instance in characters.
     * By default, this is not allowed, but the parameter
     * -allow-largesq is a by-pass to this test.
     */
    if (!allow_largesq) {
        if ((int)mpz_sizeinbase(doing.p, 2) >
                config.sides[config.side].lpb) {
            fprintf(stderr, "ERROR: The special q (%d bits) is larger than the "
                    "large prime bound on side %d (%d bits).\n",
                    (int) mpz_sizeinbase(doing.p, 2),
                    config.side,
                    config.sides[config.side].lpb);
            fprintf(stderr, "       You can disable this check with "
                    "the -allow-largesq argument,\n");
            fprintf(stderr, "       It is for instance useful for the "
                    "descent.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (doing.iteration) {
        verbose_output_print(0, 1, "#\n# NOTE:"
                " we are re-playing this special-q because of"
                " %d previous failed attempt(s)\n", doing.iteration);
        /* update sieving parameters here */
        double ratio = double(config.sides[0].mfb) /
            double(config.sides[0].lpb);
        config.sides[0].lpb += doing.iteration;
        config.sides[0].mfb = ratio*config.sides[0].lpb;
        ratio = double(config.sides[1].mfb) /
            double(config.sides[1].lpb);
        config.sides[1].lpb += doing.iteration;
        config.sides[1].mfb = ratio*config.sides[1].lpb;
        verbose_output_print(0, 1,
                "# NOTE: current values of lpb/mfb: %d,%d %d,%d\n#\n", 
                config.sides[0].lpb,
                config.sides[0].mfb,
                config.sides[1].lpb,
                config.sides[1].mfb);
    }

    return config;
}/*}}}*/

/* sieve_info stuff */

/* This function creates a new sieve_info structure, taking advantage of
 * structures which might already exist in las.sievers
 *  - for sieving, if factor base parameters are similar, we're going to
 *    share_factor_bases()
 *  - for cofactoring, if large prime bounds and mfbs are similar, we're
 *    going to reuse the strategies.
 * 
 * The siever_config structure to be passed to this function is not
 * permitted to lack anything.
 *
 * This function differs from las_info::get_sieve_info_from_config(),
 * since the latters also registers the returned object within
 * las.sievers (while the function here only *reads* this list).
 */
sieve_info::sieve_info(las_info & las, siever_config const & sc, param_list pl)/*{{{*/
{
    cpoly = las.cpoly;
    conf = sc;

    I = 1UL << sc.logI_adjusted;

    std::list<sieve_info>::iterator psi;
    
    /*** Sieving ***/

    psi = find_if(las.sievers.begin(), las.sievers.end(), has_same_fb_parameters(sc));

    if (psi != las.sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy factor base data from previous siever\n");
        share_factor_bases(other);
    } else {
        verbose_output_print(0, 1, "# bucket_region = %" PRIu64 "\n",
                BUCKET_REGION);
        init_factor_bases(las, pl);
        for (int side = 0; side < 2; side++) {
            print_fb_statistics(side);
        }
    }

    // Now that fb have been initialized, we can set the toplevel.
    toplevel = MAX(sides[0].fb->get_toplevel(), sides[1].fb->get_toplevel());

    /* Initialize the number of buckets */

    /* If LOG_BUCKET_REGION == sc.logI, then one bucket (whose size is the
     * L1 cache size) is actually one line. This changes some assumptions
     * in sieve_small_bucket_region and resieve_small_bucket_region, where
     * we want to differentiate on the parity on j.
     */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= conf.logI_adjusted);

    /* set the maximal value of the number of buckets. This directly
     * depends on A */
    uint32_t XX[FB_MAX_PARTS] = { 0, NB_BUCKETS_2, NB_BUCKETS_3, 0};
    uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
    XX[toplevel] = 1 +
      ((UINT64_C(1) << conf.logA) - UINT64_C(1)) / BRS[toplevel];
    for (int i = toplevel+1; i < FB_MAX_PARTS; ++i)
        XX[i] = 0;
    if (toplevel > 1 && XX[toplevel] == 1) {
        XX[toplevel-1] = 1 +
            ((UINT64_C(1) << conf.logA) - UINT64_C(1))
            / BRS[toplevel-1];
        ASSERT_ALWAYS(XX[toplevel-1] != 1);
    }
    for (int i = 0; i < FB_MAX_PARTS; ++i) {
        nb_buckets_max[i] = XX[i];
        nb_buckets[i] = XX[i];
    }

    for(int side = 0 ; side < 2 ; side++) {
	init_trialdiv(side); /* Init refactoring stuff */
        init_fb_smallsieved(side);
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
    }

    /*** Cofactoring ***/

    psi = find_if(las.sievers.begin(), las.sievers.end(), has_same_cofactoring(sc));

    if (psi != las.sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy cofactoring strategies from previous siever\n");
        strategies = other.strategies;
    } else {
        init_strategies(pl);
    }
}
/*}}}*/

void sieve_info::update (size_t nr_workspaces)/*{{{*/
{
    sieve_info & si(*this);

#ifdef SMART_NORM
        /* Compute the roots of the polynomial F(i,1) and the roots of its
         * inflection points d^2(F(i,1))/d(i)^2. Used in
         * init_smart_degree_X_norms_bucket_region_internal.  */
        for(int side = 0 ; side < 2 ; side++) {
            if (si.cpoly->pols[side]->deg >= 2)
                init_norms_roots (si, side);
        }
#endif

  /* update number of buckets at toplevel */
  uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
  si.nb_buckets[si.toplevel] = 1 +
      ((((uint64_t)si.J) << si.conf.logI_adjusted) - UINT64_C(1)) / 
      BRS[si.toplevel];
  // maybe there is only 1 bucket at toplevel and less than 256 at
  // toplevel-1, due to a tiny J.
  if (si.toplevel > 1) {
      if (si.nb_buckets[si.toplevel] == 1) {
        si.nb_buckets[si.toplevel-1] = 1 +
            ((((uint64_t)si.J) << si.conf.logI_adjusted) - UINT64_C(1)) /
            BRS[si.toplevel-1]; 
        // we forbid skipping two levels.
        ASSERT_ALWAYS(si.nb_buckets[si.toplevel-1] != 1);
      } else {
        si.nb_buckets[si.toplevel-1] = BRS[si.toplevel]/BRS[si.toplevel-1];
      }
  }

  /* Update the slices of the factor base according to new log base */
  for(int side = 0 ; side < 2 ; side++) {
      /* The safety factor controls by how much a single slice should fill a
         bucket array at most, i.e., with .5, a single slice should never fill
         a bucket array more than half-way. */
      const double safety_factor = .5;
      sieve_info::side_info & sis(si.sides[side]);
      double max_weight[FB_MAX_PARTS];
      for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++) {
        max_weight[i_part] = sis.max_bucket_fill_ratio[i_part] / nr_workspaces
            * safety_factor;
      }
      sis.fb->make_slices(sis.scale * LOG_SCALE, max_weight);
  }
}/*}}}*/

/* las_info stuff */

#ifdef  DLP_DESCENT
void las_info::init_hint_table(param_list_ptr pl)/*{{{*/
{
    const char * filename = param_list_lookup_string(pl, "descent-hint-table");
    if (filename == NULL) return;
    char line[1024];
    FILE * f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        /* There's no point in proceeding, since it would really change
         * the behaviour of the program to do so */
        exit(1);
    }
    max_hint_bitsize[0] = 0;
    max_hint_bitsize[1] = 0;
    for(;;) {
        char * x = fgets(line, sizeof(line), f);
        double t;
        unsigned long z;
        /* Tolerate comments and blank lines */
        if (x == NULL) break;
        for( ; *x && isspace(*x) ; x++) ;
        if (*x == '#') continue;
        if (!*x) continue;

        descent_hint h;
        siever_config & sc(h.conf);

        /* start with the global defaults */
        sc = config_base;

        z = strtoul(x, &x, 10);
        ASSERT_ALWAYS(z > 0);
        sc.bitsize = z;
        switch(*x++) {
            case '@' :
                sc.side = strtoul(x, &x, 0);
                ASSERT_ALWAYS(sc.side < 2);
                break;
            default:
                fprintf(stderr, "%s: parse error at %s\n", filename, line);
                exit(1);
        }
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h.expected_time = t;
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h.expected_success = t;
        for( ; *x && isspace(*x) ; x++) ;
        char letter MAYBE_UNUSED = *x;
        for( ; *x && !isdigit(*x) ; x++) ;
        z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
        if (letter == 'I') {
            sc.logA = 2*z-1;
        } else if (letter == 'A') {
            sc.logA = z;
        } else {
            fprintf(stderr, "%s: parse error (want I= or A=) at %s\n", filename, line);
            exit(EXIT_FAILURE);
        }
        
        for(int s = 0 ; s < 2 ; s++) {
            for( ; *x && isspace(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc.sides[s].lim = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc.sides[s].lpb = z;
            /* recognize this as a double. If it's < 10, we'll consider
             * this means lambda */
            {
                for( ; *x && !isdigit(*x) ; x++) ;
                double t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                if (t < 10) {
                    sc.sides[s].lambda = t;
                    sc.sides[s].mfb = t * sc.sides[s].lpb;
                    /* Then no "lambda" is allowed */
                    continue;
                } else {
                    sc.sides[s].mfb = t;
                }
            }
            if (*x == ',') {
                for( ; *x && !isdigit(*x) ; x++) ;
                t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                sc.sides[s].lambda = t;
            } else {
                /* this means "automatic" */
                sc.sides[s].lambda = 0;
            }
        }
        for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));

        if (sc.bitsize > max_hint_bitsize[sc.side])
            max_hint_bitsize[sc.side] = sc.bitsize;

        hint_table.push_back(h);
    }
    if (hint_table.empty()) {
        fprintf(stderr, "%s: no data ??\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Allocate the quick lookup tables */

    for(int s = 0 ; s < 2 ; s++) {
        unsigned int n = max_hint_bitsize[s] + 1;
        hint_lookups[s] = (int *)malloc(n * sizeof(int));
        for(unsigned int i = 0 ; i < n ; i++) {
            hint_lookups[s][i] = -1;
        }
    }
    for(unsigned int i = 0 ; i < hint_table.size() ; i++) {
        descent_hint & h(hint_table[i]);
        siever_config & sc = h.conf;
        int d = hint_lookups[sc.side][sc.bitsize];
        if (d >= 0) {
            fprintf(stderr, "Error: two hints found for %d@%d\n",
                    sc.bitsize, sc.side);
            exit(1);
        }
        hint_lookups[sc.side][sc.bitsize] = i;
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
            ASSERT_ALWAYS(sc.sides[s]->lim <= cpoly->pols[s]->lim);
        }
        */
    }

    fclose(f);
}/*}}}*/

void las_info::clear_hint_table()
{
    free(hint_lookups[0]);
    free(hint_lookups[1]);
}

#endif /* DLP_DESCENT */

/* {{{ Parse default siever config (fill all possible fields). Return
 * true if the parsed siever config is complete and can be used without
 * per-special-q info. */
bool parse_default_siever_config(siever_config & sc, param_list_ptr pl)
{
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
        sc.bitsize = mpz_sizeinbase(q0, 2);
    }
    mpz_clear(q0);
    sc.side = 1; // Legacy.
    if (!param_list_parse_int(pl, "sqside", &sc.side)) {
        verbose_output_print(0, 1, "# Warning: sqside not given, "
                "assuming side 1 for backward compatibility.\n");
    }
    param_list_parse_double(pl, "lambda0", &(sc.sides[0].lambda));
    param_list_parse_double(pl, "lambda1", &(sc.sides[1].lambda));
    int complete = 1;
    complete &= param_list_parse_ulong(pl, "lim0", &(sc.sides[0].lim));
    complete &= param_list_parse_ulong(pl, "lim1", &(sc.sides[1].lim));
    if (param_list_lookup_string(pl, "A")) {
        complete &= param_list_parse_int  (pl, "A",    &(sc.logA));
        if (param_list_lookup_string(pl, "I")) {
            fprintf(stderr, "# -A and -I are incompatible\n");
            exit(EXIT_FAILURE);
        }
        // This is useful for a default config
        sc.logI_adjusted = (sc.logA+1)/2;
    } else if (param_list_lookup_string(pl, "I")) {
        int I;
        complete &= param_list_parse_int  (pl, "I", &I);
        sc.logA = 2 * I - 1;
        verbose_output_print(0, 1, "# Interpreting -I %d as meaning -A %d\n", I, sc.logA);
        // This is useful for a default config
        sc.logI_adjusted = I;
    }

    if (sc.sides[0].lim > 2147483647UL || sc.sides[1].lim > 2147483647UL)
    {
        fprintf (stderr, "Error, lim0/lim1 must be < 2^31 (bug #21094)\n");
        exit (EXIT_FAILURE);
    }

    complete &= param_list_parse_int(pl, "lpb0",  &(sc.sides[0].lpb));
    complete &= param_list_parse_int(pl, "mfb0",  &(sc.sides[0].mfb));
    complete &= param_list_parse_int(pl, "lpb1",  &(sc.sides[1].lpb));
    complete &= param_list_parse_int(pl, "mfb1",  &(sc.sides[1].mfb));
    if (!complete) {
        verbose_output_print(0, 1, "# default siever configuration is incomplete ; required parameters are I, lim[01], lpb[01], mfb[01]\n");

    }

    /* Parse optional siever configuration parameters */
    sc.td_thresh = 1024;	/* default value */
    sc.skipped = 1;	/* default value */
    param_list_parse_uint(pl, "tdthresh", &(sc.td_thresh));
    param_list_parse_uint(pl, "skipped", &(sc.skipped));

    sc.unsieve_thresh = 100;
    if (param_list_parse_uint(pl, "unsievethresh", &(sc.unsieve_thresh))) {
        verbose_output_print(0, 1, "# Un-sieving primes > %u\n",
                sc.unsieve_thresh);
    }

    // XXX note that when the sieving range size varies with the
    // special-q, we need to accept that the bucket threshold varies,
    // too.
    //
    // As a consequence, we should make it possible to specify extra
    // parameters in the hint file format (not just I,lim,lpb,mfb).
    sc.bucket_thresh = 1 << ((sc.logA+1)/2);	/* default value */
    sc.bucket_thresh1 = 0;	/* default value */
    sc.bk_multiplier = 0.0;
    /* overrides default only if parameter is given */
    param_list_parse_ulong(pl, "bkthresh", &(sc.bucket_thresh));
    param_list_parse_ulong(pl, "bkthresh1", &(sc.bucket_thresh1));
    param_list_parse_double(pl, "bkmult", &(sc.bk_multiplier));

    const char *powlim_params[2] = {"powlim0", "powlim1"};
    for (int side = 0; side < 2; side++) {
        if (!param_list_parse_ulong(pl, powlim_params[side],
                    &sc.sides[side].powlim)) {
            sc.sides[side].powlim = sc.bucket_thresh - 1;
            verbose_output_print(0, 1,
                    "# Using default value of %lu for -%s\n",
                    sc.sides[side].powlim, powlim_params[side]);
        }
    }

    const char *ncurves_params[2] = {"ncurves0", "ncurves1"};
    for (int side = 0; side < 2; side++)
        if (!param_list_parse_int(pl, ncurves_params[side],
                    &sc.sides[side].ncurves))
            sc.sides[side].ncurves = -1;

    long dupqmax[2] = {0, 0};
    param_list_parse_long_and_long(pl, "dup-qmax", dupqmax, ",");
    sc.sides[0].qmax = dupqmax[0];
    sc.sides[1].qmax = dupqmax[1];

    return complete;
}
/* }}} */

las_info::las_info(param_list_ptr pl)/*{{{*/
{
    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{
    nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "t", &nb_threads);
    if (nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    output = stdout;
    outputname = param_list_lookup_string(pl, "out");
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }
    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */

    galois = galois_automorphism_get(param_list_lookup_string(pl, "galois"));
    verbose = param_list_parse_switch(pl, "-v");
    suppress_duplicates = param_list_parse_switch(pl, "-dup");

    verbose_interpret_parameters(pl);
    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, output, verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, output, 1);
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
    param_list_print_command_line(output, pl);
    las_display_config_flags();
    /*  Parse polynomial */
    cado_poly_init(cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }
    if (!cado_poly_read(cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    // sc.skewness = cpoly->skew;
    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    param_list_parse_double(pl, "skew", &(cpoly->skew));
    if (cpoly->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	cado_poly_clear(cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    gmp_randinit_default(rstate);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);
    // }}}

    // ----- default config and adaptive configs {{{
    // perhaps the set of parameters passed is not complete.
    if (parse_default_siever_config(config_base, pl)) {
        default_config_ptr = &config_base;
        sievers.push_back(sieve_info(*this, config_base, pl));
    }
    // }}}

    // ----- stuff roughly related to the descent {{{
#ifdef DLP_DESCENT
    int has_descent_hint = param_list_lookup_string(pl, "descent-hint-table") != NULL;
#else
    const int has_descent_hint = 0;
#endif

    if (!default_config_ptr && !has_descent_hint) {
        fprintf(stderr,
                "Error: no default config set, and no hint table either\n");
        exit(EXIT_FAILURE);
    }

    memset(max_hint_bitsize, 0, sizeof(max_hint_bitsize));
    memset(hint_lookups, 0, sizeof(hint_lookups));
    descent_helper = NULL;
#ifdef  DLP_DESCENT
    init_hint_table(pl);

    dlog_base = new las_dlog_base();
    dlog_base->lookup_parameters(pl);
    dlog_base->read();
#endif
    // }}}

    // ----- todo list and such {{{
    nq_pushed = 0;
    nq_max = UINT_MAX;
    random_sampling = 0;
    if (param_list_parse_uint(pl, "random-sample", &nq_max)) {
        random_sampling = 1;
    } else if (param_list_parse_uint(pl, "nq", &nq_max)) {
        if (param_list_lookup_string(pl, "q1") || param_list_lookup_string(pl, "rho")) {
            fprintf(stderr, "Error: argument -nq is incompatible with -q1 or -rho\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        todo_list_fd = fopen(filename, "r");
        if (todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            cado_poly_clear(cpoly);
            param_list_clear(pl);
            exit(EXIT_FAILURE);
        }
    } else {
        todo_list_fd = NULL;
    }
    // }}}

    // ----- batch mode {{{
    batch = param_list_parse_switch(pl, "-batch");
    batch_print_survivors = param_list_parse_switch(pl, "-batch-print-survivors");
    cofac_list_init (L);
    // }}} 
    
    init_cof_stats(pl);
}/*}}}*/

las_info::~las_info()/*{{{*/
{
    // ----- general operational flags {{{
    if (outputname)
        fclose_maybe_compressed(output, outputname);
    gmp_randclear(rstate);
    verbose_output_clear();
    cado_poly_clear(cpoly);
    // }}}

    // ----- default config and adaptive configs: nothing

    // ----- stuff roughly related to the descent {{{
#ifdef  DLP_DESCENT
    clear_hint_table();
    delete dlog_base;
#endif
    // }}}
 
    // ----- todo list and such {{{
    if (todo_list_fd) {
        fclose(todo_list_fd);
        todo_list_fd = NULL;
    }
    // }}}
 
    // ----- batch mode: very little
    cofac_list_clear (L);

    clear_cof_stats();
}/*}}}*/

// {{{ las_info::{init,clear,print}_cof_stats
void las_info::init_cof_stats(param_list_ptr pl)
{
    const char *statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (statsfilename != NULL) { /* a file was given */
        if (default_config_ptr == NULL) {
            fprintf(stderr, "Error: option stats-cofact works only "
                    "with a default config\n");
            exit(EXIT_FAILURE);
#ifdef DLP_DESCENT
        } else if (param_list_lookup_string(pl, "descent-hint-table")) {
            verbose_output_print(0, 1, "# Warning: option stats-cofact "
                    "only applies to the default siever config\n");
#endif
        }
        siever_config const & sc(*default_config_ptr);

        cof_stats_file = fopen (statsfilename, "w");
        if (cof_stats_file == NULL)
        {
            fprintf (stderr, "Error, cannot create file %s\n",
                    statsfilename);
            exit (EXIT_FAILURE);
        }
        //allocate cof_call and cof_succ
        int mfb0 = sc.sides[0].mfb;
        int mfb1 = sc.sides[1].mfb;
        cof_call = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        cof_succ = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        for (int i = 0; i <= mfb0; i++) {
            cof_call[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            cof_succ[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            for (int j = 0; j <= mfb1; j++)
                cof_call[i][j] = cof_succ[i][j] = 0;
        }
    } else {
        cof_stats_file = NULL;
        cof_call = NULL;
        cof_succ = NULL;
    }
}
void las_info::print_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(default_config_ptr);
    siever_config const & sc0(*default_config_ptr);
    int mfb0 = sc0.sides[0].mfb;
    int mfb1 = sc0.sides[1].mfb;
    for (int i = 0; i <= mfb0; i++) {
        for (int j = 0; j <= mfb1; j++)
            fprintf (cof_stats_file, "%u %u %u %u\n", i, j,
                    cof_call[i][j],
                    cof_succ[i][j]);
    }
}

void las_info::clear_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(default_config_ptr);
    siever_config const & sc0(*default_config_ptr);
    for (int i = 0; i <= sc0.sides[0].mfb; i++) {
        free (cof_call[i]);
        free (cof_succ[i]);
    }
    free (cof_call);
    free (cof_succ);
    fclose (cof_stats_file);
    cof_stats_file = NULL;
}
//}}}


/* Look for an existing sieve_info in las.sievers with configuration matching
   that in sc; if none exists, create one. */
sieve_info & get_sieve_info_from_config(las_info & las, siever_config const & sc, param_list pl)/*{{{*/
{
    std::list<sieve_info>::iterator psi;
#if 0
#if 0
    psi = find_if(las.sievers.begin(), las.sievers.end(), has_same_config_q_A_logI(sc));
#else
    psi = find_if(las.sievers.begin(), las.sievers.end(), has_same_config_q_logI(sc));
#endif
#endif
    psi = find_if(las.sievers.begin(), las.sievers.end(), has_same_config(sc));
    if (psi != las.sievers.end()) {
        siever_config_display(sc);
        return *psi;
    }
    las.sievers.push_back(sieve_info(las, sc, pl));
    sieve_info & si(las.sievers.back());
    verbose_output_print(0, 1, "# Creating new sieve configuration for q~2^%d on side %d (logI=%d)\n",
            sc.bitsize, sc.side, si.conf.logI_adjusted);
    siever_config_display(sc);
    return las.sievers.back();
}/*}}}*/

void las_todo_push_withdepth(las_info & las, mpz_srcptr p, mpz_srcptr r, int side, int depth, int iteration = 0)/*{{{*/
{
    las.todo.push(las_todo_entry(p, r, side, depth, iteration));
}
void las_todo_push(las_info & las, mpz_srcptr p, mpz_srcptr r, int side)
{
    las_todo_push_withdepth(las, p, r, side, 0);
}
void las_todo_push_closing_brace(las_info & las, int depth)
{
    las.todo.push(las_todo_entry(-1, depth));
}
las_todo_entry las_todo_pop(las_info & las)
{
    las_todo_entry r = las.todo.top();
    las.todo.pop();
    return r;
}
int las_todo_pop_closing_brace(las_info & las)
{
    if (las.todo.top().side >= 0)
        return 0;
    las.todo.pop();
    return 1;
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
parse_command_line_q0_q1(las_info & las, mpz_ptr q0, mpz_ptr q1, param_list pl, const int qside)
{
    ASSERT_ALWAYS(param_list_parse_mpz(pl, "q0", q0));
    if (param_list_parse_mpz(pl, "q1", q1)) {
        next_legitimate_specialq(q0, q0, 0);
        return;
    }

    /* We don't have -q1. If we have -rho, we sieve only <q0, rho>. */
    mpz_t t;
    mpz_init(t);
    if (param_list_parse_mpz(pl, "rho", t)) {
        las_todo_push(las, q0, t, qside);
        /* Set empty interval [q0 + 1, q0] as special-q interval */
        mpz_set(q1, q0);
        mpz_add_ui (q0, q0, 1);
    } else {
    /* If we don't have -rho, we sieve only q0, but all roots of it.
       If -q0 does not give a legitimate special-q value, advance to the
       next legitimate one. */
        mpz_set(t, q0);
        next_legitimate_specialq(q0, q0, 0);
        mpz_set(q1, q0);
    }
    mpz_clear(t);
}

static int
skip_galois_roots(
        galois_automorphism_srcptr sigma,
        cxx_mpz * roots,
        int nroots,
        mpz_srcptr q)
{
    if (nroots % sigma->order) {
        gmp_fprintf(stderr, "Number of roots (%d) modulo q=%Zd is not divisible by order(sigma)=%d. Don't know how to interpret -galois.\n", nroots, q, sigma->order);
        ASSERT_ALWAYS(0);
    }

    std::vector<bool> met(nroots, false);

    for(int i = 0 ; i < nroots ; i++) {
        if (met[i]) continue;
        cxx_mpz r = roots[i];
        for(int k = 1 ; k < sigma->order ; k++) {
            galois_automorphism_apply_root(sigma, r, r, q);
            int found = std::find(roots, roots + nroots, r) - roots;
            if (found == i) {
                /* root in a subfield */
                break;
            }
            if (found < nroots)
                met[found] = true;
        }
    }
    /* then we discard all those which we haven't met as conjugate of
     * some other "first found" root */
    int j = 0;
    for(int i = 0 ; i < nroots ; i++) {
        if (!met[i]) {
            if (i != j)
                std::swap(roots[i], roots[j]);
            j++;
        }
    }
    return j;
}

/* adding relations on the fly in Galois cases */
static void add_relations_with_galois(
        cado_poly_srcptr cpoly,
        galois_automorphism_srcptr sigma,
        FILE *output, const char *comment, int *cpt, relation &rel)
{
    ASSERT_ALWAYS(rel.ab().size() == 2);
    int64_t a = mpz_get_int64(rel[0]);
    int64_t b = mpz_get_int64(rel[1]);

    /* Remember that the relation got printed once, already */
    for(int i = 1 ; i < sigma->order ; i++) {
        int e = galois_automorphism_apply_ab(sigma, &a, &b);
        mpz_set_uint64(rel[0], a);
        mpz_set_int64(rel[1], b);
        if (e) {
            /* in case sigma does not have a transitive action on all
             * roots, we need to take this into account. I suppose that
             * scaling will suffice.
             */
            rel.addv(0, sigma->factor, cpoly->pols[0]->deg / sigma->order * e);
            rel.addv(1, sigma->factor, cpoly->pols[1]->deg / sigma->order * e);
            rel.compress();
        }
        *cpt += 1;
        rel.print(output, comment);
    }
}

/* {{{ Populating the todo list */
/* See below in main() for documentation about the q-range and q-list
 * modes */
/* These functions return non-zero if the todo list is not empty */
int las_todo_feed_qrange(las_info & las, param_list pl)
{
    /* If we still have entries in the stack, don't add more now */
    if (!las.todo.empty())
        return 1;

    const unsigned long push_at_least_this_many = 1;

    mpz_ptr q0 = las.todo_q0;
    mpz_ptr q1 = las.todo_q1;

    int qside = las.config_base.side;

    mpz_poly_ptr f = las.cpoly->pols[qside];
    cxx_mpz roots[MAX_DEGREE];

    if (mpz_cmp_ui(q0, 0) == 0) {
        parse_command_line_q0_q1(las, q0, q1, pl, qside);
        if (las.random_sampling) {
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
                if (mpz_poly_roots ((mpz_t*) roots, f, q) > 0)
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
    if (!las.random_sampling) {
        /* handy aliases */
        mpz_ptr q = q0;

        /* If nq_max is specified, then q1 has no effect, even though it
         * has been set equal to q */
        for ( ; las.todo.size() < push_at_least_this_many &&
                (las.nq_max < UINT_MAX || mpz_cmp(q, q1) < 0) &&
                las.nq_pushed < las.nq_max ; )
        {
            int nroots = mpz_poly_roots ((mpz_t*)roots, f, q);
            if (nroots == 0) {
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# polynomial has no roots for q = %Zu\n", q);
            }

            if (las.galois != NULL)
                nroots = skip_galois_roots(las.galois, roots, nroots, q);

            for(int i = 0 ; i < nroots && las.nq_pushed < las.nq_max; i++) {
                las.nq_pushed++;
                las_todo_push(las, q, roots[nroots-1-i], qside);
            }

            next_legitimate_specialq(q, q, 1);
        }
    } else {
        /* we don't care much about being truly uniform here */
        mpz_t q;
        mpz_init(q);
        for ( ; las.todo.size() < push_at_least_this_many && las.nq_pushed < las.nq_max ; ) {
            mpz_sub(q, q1, q0);
            mpz_urandomm(q, las.rstate, q);
            mpz_add(q, q, q0);
            next_legitimate_specialq(q, q, 0);
            int nroots = mpz_poly_roots ((mpz_t*)roots, f, q);
            if (!nroots) continue;
            if (las.galois != NULL)
                nroots = skip_galois_roots(las.galois, roots, nroots, q);
            unsigned long i = gmp_urandomm_ui(las.rstate, nroots);
            las.nq_pushed++;
            las_todo_push(las, q, roots[i], qside);
        }
        mpz_clear(q);
    }

    return las.todo.size();
}

/* Format of a file with a list of special-q (-todo option):
 *   - Comments are allowed (start line with #)
 *   - Blank lines are ignored
 *   - Each valid line must have the form
 *       s q r
 *     where s is the side (0 or 1) of the special q, and q and r are as usual.
 */
int las_todo_feed_qlist(las_info & las, param_list pl)
{
    if (!las.todo.empty())
        return 1;

    char line[1024];
    FILE * f = las.todo_list_fd;
    /* The fgets call below is blocking, so flush las.output here just to
     * be sure. */
    fflush(las.output);
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
        mpz_poly_ptr f = las.cpoly->pols[side];
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
    las_todo_push(las, p, r, side);
    mpz_clear(p);
    mpz_clear(r);
    return 1;
}


int las_todo_feed(las_info & las, param_list pl)
{
    if (!las.todo.empty())
        return 1;
    if (las.todo_list_fd)
        return las_todo_feed_qlist(las, pl);
    else
        return las_todo_feed_qrange(las, pl);
}
/* }}} */

/* Only when tracing. This function gets called once per
 * special-q only. Here we compute the two norms corresponding to
 * the traced (a,b) pair, and start by dividing out the special-q
 * from the one where it should divide */
void init_trace_k(sieve_info const & si, param_list pl)
{
    struct trace_ab_t ab;
    struct trace_ij_t ij;
    struct trace_Nx_t Nx;
    int have_trace_ab = 0, have_trace_ij = 0, have_trace_Nx = 0;

    const char *abstr = param_list_lookup_string(pl, "traceab");
    if (abstr != NULL) {
        if (sscanf(abstr, "%" SCNu64",%" SCNd64, &ab.a, &ab.b) == 2)
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
        const unsigned char logp, where_am_I & w)
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
        const fb_part *fb, where_am_I & w)
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
#if defined(__ICC) || defined(__INTEL_COMPILER)
      /* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=45414 */
      _mm_prefetch(((const char *) it)+256, _MM_HINT_NTA);
#else
      _mm_prefetch(((unsigned char *) it)+256, _MM_HINT_NTA);
#endif
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
        const fb_part *fb, where_am_I & w);

template <>
void apply_one_bucket<longhint_t> (unsigned char *S,
        const bucket_array_t<1, longhint_t> &BA, const int i,
        const fb_part *fb, where_am_I & w) {
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
           uint64_t a, int64_t b,
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
bool register_contending_relation(las_info const & las, sieve_info const & si, relation & rel)
{
    if (las.tree.must_avoid(rel)) {
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
            relation::pr const & v(rel.sides[side][i]);
            if (mpz_cmp(si.doing.p, v.p) == 0)
                continue;
            unsigned long p = mpz_get_ui(v.p);
            if (mpz_fits_ulong_p(v.p)) {
                unsigned long r = mpz_get_ui(v.r);
                if (las.dlog_base->is_known(side, p, r))
                    continue;
            }

            unsigned int n = mpz_sizeinbase(v.p, 2);
            int k = (n <= las.max_hint_bitsize[side] ? las.hint_lookups[side][n] : -1);
            if (k < 0) {
                /* This is not worrysome per se. We just do
                 * not have the info in the descent hint table,
                 * period.
                 */
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: cannot estimate refactoring time for relation involving %d@%d (%Zd,%Zd)\n", n, side, (mpz_srcptr) v.p, (mpz_srcptr) v.r);
                time_left = INFINITY;
            } else {
                if (std::isfinite(time_left))
                    time_left += las.hint_table[k].expected_time;
            }
            contender.outstanding.push_back(std::make_pair(side, v));
        }
    }
    verbose_output_print(0, 1, "# [descent] This relation entails an additional time of %.2f for the smoothing process (%zu children)\n",
            time_left, contender.outstanding.size());

    /* when we're re-examining this special-q because of a previous
     * failure, there's absolutely no reason to hurry up on a relation */
    contender.set_time_left(time_left, si.doing.iteration ? INFINITY : general_grace_time_ratio);

    return las.tree.new_candidate_relation(contender);
}
#endif /* DLP_DESCENT */

struct factor_survivors_data {
    thread_data * th;
    std::vector<uint32_t> survivors;
    std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors2;
    unsigned char * SS;
    int N;
    where_am_I & w;
    int cpt;
    int copr;
    uint32_t cof_bitsize[2];

    struct side_data {
        unsigned char * S;
        bucket_primes_t primes;
        bucket_array_complete purged;
        side_data() :
            primes(bucket_primes_t(BUCKET_REGION)),
            purged(bucket_array_complete(BUCKET_REGION))
        {}
    };

    side_data sdata[2];

    factor_survivors_data(thread_data *th, int N, where_am_I & w)
        : th(th), N(N), w(w)
    {
        cpt = 0;
        copr = 0;
        for(int side = 0 ; side < 2 ; side++) {
            sdata[side].S = th->sides[side].bucket_region;
            cof_bitsize[side]=0;
        }
        /* This is the one which gets the merged information in the end */
        SS = sdata[0].S;
    }
    ~factor_survivors_data() {
    }
    void search_survivors (timetree_t & timer);
    void convert_survivors (timetree_t & timer);
    void purge_buckets (timetree_t & timer);
    void cofactoring (timetree_t & timer);
};

void factor_survivors_data::search_survivors(timetree_t & timer)
{
    CHILD_TIMER(timer, __func__);
    las_info const & las(*th->plas);
    sieve_info & si(*th->psi);
    const unsigned int first_j = N << (LOG_BUCKET_REGION - si.conf.logI_adjusted);
    const unsigned long nr_lines = 1U << (LOG_BUCKET_REGION - si.conf.logI_adjusted);

#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# When entering factor_survivors for bucket %u, "
                "S[0][%u]=%u, S[1][%u]=%u\n",
                trace_Nx.N, trace_Nx.x, sdata[0].S[trace_Nx.x], trace_Nx.x, sdata[1].S[trace_Nx.x]);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n",
                traced_norms[0], traced_norms[1]);
    }
#endif  /* }}} */

    if (las.verbose >= 2)
        th->update_checksums();


#ifdef TRACE_K /* {{{ */
    sieve_info::side_info & side0(si.sides[0]);
    sieve_info::side_info & side1(si.sides[1]);
    for (int x = 0; x < 1 << LOG_BUCKET_REGION; x++) {
        if (trace_on_spot_Nx(N, x)) {
            verbose_output_print(TRACE_CHANNEL, 0,
                    "# side0.Bound[%u]=%u, side1.Bound[%u]=%u\n",
                    sdata[0].S[trace_Nx.x],
                    sdata[0].S[x] <= side0.bound ? 0 : side0.bound,
                    sdata[1].S[trace_Nx.x],
                    sdata[1].S[x] <= side0.bound ? 0 : side1.bound);
            /* Hmmm, is that right ? side0.bound 3 times, really ? XXX */
        }
    }
#endif /* }}} */

    survivors.reserve(128);
    for (unsigned int j = 0; j < nr_lines; j++)
    {
        unsigned char * const both_S[2] = {
            sdata[0].S + (j << si.conf.logI_adjusted), 
            sdata[1].S + (j << si.conf.logI_adjusted)
        };
        const unsigned char both_bounds[2] = {
            si.sides[0].bound,
            si.sides[1].bound,
        };
        size_t old_size = survivors.size();
        search_survivors_in_line(both_S, both_bounds,
                                 si.conf.logI_adjusted, j + first_j, N,
                                 si.j_div, si.conf.unsieve_thresh,
                                 si.us, survivors);
        /* Survivors written by search_survivors_in_line() have index
           relative to their j-line. We need to convert to index within
           the bucket region by adding line offsets. */
        for (size_t i_surv = old_size; i_surv < survivors.size(); i_surv++)
            survivors[i_surv] += j << si.conf.logI_adjusted;
    }
}

void factor_survivors_data::convert_survivors(timetree_t & timer)
{
    CHILD_TIMER(timer, __func__);
    /* Convert data type of list from uint32_t to the correct br_index_t */
    survivors2.reserve(survivors.size());
    for (size_t i = 0; i < survivors.size(); i++) {
        survivors2.push_back(survivors[i]);
    }
    survivors.clear();
}

void factor_survivors_data::purge_buckets(timetree_t & timer)
{
    CHILD_TIMER(timer, __func__);

    sieve_info & si(*th->psi);

    /* Copy those bucket entries that belong to sieving survivors and
       store them with the complete prime */
    /* FIXME: choose a sensible size here */

    for(int side = 0 ; side < 2 ; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        // From N we can deduce the bucket_index. They are not the same
        // when there are multiple-level buckets.
        uint32_t bucket_index = N % si.nb_buckets[1];

        SIBLING_TIMER(timer, "purge buckets");

        const bucket_array_t<1, shorthint_t> *BA =
            th->ws->cbegin_BA<1, shorthint_t>(side);
        const bucket_array_t<1, shorthint_t> * const BA_end =
            th->ws->cend_BA<1, shorthint_t>(side);
        for (; BA != BA_end; BA++)  {
#if defined(HAVE_SSE2) && defined(SMALLSET_PURGE)
            sdata[side].purged.purge(*BA, bucket_index, SS, survivors2);
#else
            sdata[side].purged.purge(*BA, bucket_index, SS);
#endif
        }

        /* Add entries coming from downsorting, if any */
        const bucket_array_t<1, longhint_t> *BAd =
            th->ws->cbegin_BA<1, longhint_t>(side);
        const bucket_array_t<1, longhint_t> * const BAd_end =
            th->ws->cend_BA<1, longhint_t>(side);
        for (; BAd != BAd_end; BAd++)  {
            sdata[side].purged.purge(*BAd, bucket_index, SS);
        }

        SIBLING_TIMER(timer, "resieve");
        /* Resieve small primes for this bucket region and store them 
           together with the primes recovered from the bucket updates */
        resieve_small_bucket_region (&sdata[side].primes, N, SS,
                th->psi->sides[side].rsd, th->sides[side].rsdpos, si, w);

        SIBLING_TIMER(timer, "sort primes in purged buckets");
        /* Sort the entries to avoid O(n^2) complexity when looking for
           primes during trial division */
        sdata[side].purged.sort();
        sdata[side].primes.sort();
    }

#ifdef TRACE_K
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u has value %u\n",
                trace_Nx.x, trace_Nx.N, SS[trace_Nx.x]);
    }
#endif
}

void factor_survivors_data::cofactoring (timetree_t & timer)
{
    CHILD_TIMER(timer, __func__);

    las_info const & las(*th->plas);
    sieve_info & si(*th->psi);
    las_report_ptr rep = th->rep;

    mpz_t norm[2];
    factor_list_t factors[2];
    mpz_array_t *lps[2] = { NULL, };
    uint32_array_t *lps_m[2] = { NULL, }; /* corresponding multiplicities */

    for(int side = 0 ; side < 2 ; side++) {
        lps[side] = alloc_mpz_array (1);
        lps_m[side] = alloc_uint32_array (1);
        factor_list_init(&factors[side]);
        mpz_init (norm[side]);
    }

#ifdef SUPPORT_LARGE_Q
        mpz_t az, bz;
        mpz_init(az);
        mpz_init(bz);
#endif

    for (size_t i_surv = 0 ; i_surv < survivors2.size(); i_surv++) {
#ifdef DLP_DESCENT
      if (las.tree.must_take_decision())
	break;
#endif
      const size_t x = survivors2[i_surv];
      ASSERT_ALWAYS (SS[x] != 255);

        th->rep->survivor_sizes[sdata[0].S[x]][sdata[1].S[x]]++;
        
        /* For factor_leftover_norm, we need to pass the information of the
         * sieve bound. If a cofactor is less than the square of the sieve
         * bound, it is necessarily prime. we implement this by keeping the
         * log to base 2 of the sieve limits on each side, and compare the
         * bitsize of the cofactor with their double.
         */
        int64_t a;
        int64_t b;

        SIBLING_TIMER(timer, "check_coprime");

        NxToAB (&a, &b, N, x, si);
#ifdef SUPPORT_LARGE_Q
        NxToABmpz (az, bz, N, x, si);
#endif
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
#ifndef SUPPORT_LARGE_Q
        if (UNLIKELY(((a | b) & 1) == 0))
#else
        if (UNLIKELY(mpz_even_p(az) && mpz_even_p(bz)))
#endif
        {
            th->rep->both_even++;
            continue;
        }

        /* Since the q-lattice is exactly those (a, b) with
           a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
        /* FIXME: fast divisibility test here! */
        /* Dec 2014: on a c90, it takes 0.1 % of total sieving time*/
#ifndef SUPPORT_LARGE_Q
        if (b == 0 || (mpz_cmp_ui(si.doing.p, b) <= 0 && b % mpz_get_ui(si.doing.p) == 0))
#else
        if ((mpz_cmp_ui(bz, 0) == 0) || 
            (mpz_cmp(si.doing.p, bz) <= 0 &&
             mpz_divisible_p(bz, si.doing.p)))
#endif
        {
            continue;
        }

        copr++;

        int pass = 1;

        int i;
        unsigned int j;
        for(int side = 0 ; pass && side < 2 ; side++) {

            SIBLING_TIMER(timer, "recompute complete norm");

            // Trial divide norm on side 'side'
            /* Compute the norms using the polynomials transformed to 
               i,j-coordinates. The transformed polynomial on the 
               special-q side is already divided by q */
            NxToIJ (&i, &j, N, x, si);
		mpz_poly_homogeneous_eval_sisi (norm[side], si.sides[side].fij, i, j);

#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                verbose_output_vfprint(TRACE_CHANNEL, 0,
                        gmp_vfprintf, "# start trial division for norm=%Zd ", norm[side]);
                verbose_output_print(TRACE_CHANNEL, 0,
                        "on side %d for (%" PRId64 ",%" PRIu64 ")\n", side, a, b);
            }
#endif

            SIBLING_TIMER(timer, "trial division");

            verbose_output_print(1, 2, "FIXME %s, line %d\n", __FILE__, __LINE__);
            const bool handle_2 = true; /* FIXME */
            trial_div (&factors[side], norm[side], N, x,
                    handle_2,
                    &sdata[side].primes,
                    &sdata[side].purged,
                    si.sides[side].trialdiv_data.get(),
                    a, b,
                    th->psi->sides[side].fb.get());

            SIBLING_TIMER(timer, "check_leftover_norm");

            pass = check_leftover_norm (norm[side], si, side);
#ifdef TRACE_K
            if (trace_on_spot_ab(a, b)) {
                verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                        "# checked leftover norm=%Zd", norm[side]);
                verbose_output_print(TRACE_CHANNEL, 0,
                        " on side %d for (%" PRId64 ",%" PRIu64 "): %d\n",
                        side, a, b, pass);
            }
#endif
        }

        if (!pass) continue;

        if (las.batch_print_survivors) {
#ifndef SUPPORT_LARGE_Q
            gmp_printf("%" PRId64 " %" PRIu64 " %Zd %Zd\n", a, b,
                    norm[0], norm[1]);
#else
            gmp_printf("%Zd %Zd %Zd %Zd\n", az, bz, norm[0], norm[1]);
#endif
            cpt++;
            continue;
        }

        if (las.batch)
        {
            verbose_output_start_batch ();
            cofac_list_add ((cofac_list_t*) las.L, a, b, norm[0], norm[1],
                    si.doing.side, si.doing.p);
            verbose_output_end_batch ();
            cpt++;
            continue; /* we deal with all cofactors at the end of las */
        }

        if (las.cof_stats_file) {
            cof_bitsize[0] = mpz_sizeinbase (norm[0], 2);
            cof_bitsize[1] = mpz_sizeinbase (norm[1], 2);
            /* no need to use a mutex here: either we use one thread only
               to compute the cofactorization data and if several threads
               the order is irrelevant. The only problem that can happen
               is when two threads increase the value at the same time,
               and it is increased by 1 instead of 2, but this should
               happen rarely. */
            las.cof_call[cof_bitsize[0]][cof_bitsize[1]] ++;
        }

        SIBLING_TIMER(timer, "factor_both_leftover_norms");

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
        SIBLING_TIMER(timer, "print relations");

        if (las.cof_stats_file) /* learning phase */
            las.cof_succ[cof_bitsize[0]][cof_bitsize[1]] ++;
	    
        // ASSERT (bin_gcd_int64_safe (a, b) == 1);

#ifndef SUPPORT_LARGE_Q
        relation rel(a, b);
#else
        relation rel(az, bz);
#endif

        /* Note that we explicitly do not bother about storing r in
         * the relations below */
        for (int side = 0; side < 2; side++) {
            for(int i = 0 ; i < factors[side].n ; i++)
                rel.add(side, factors[side].fac[i], 0);

            for (unsigned int i = 0; i < lps[side]->length; ++i)
                rel.add(side, lps[side]->data[i], 0);
        }

        rel.add(si.conf.side, si.doing.p, 0);

        rel.compress();

#ifdef TRACE_K
        if (trace_on_spot_ab(a, b)) {
            verbose_output_print(TRACE_CHANNEL, 0, "# Relation for (%"
                    PRId64 ",%" PRIu64 ") printed\n", a, b);
        }
#endif
        {
            int do_check = th->plas->suppress_duplicates;
            /* note that if we have large primes which don't fit in
             * an unsigned long, then the duplicate check will
             * quickly return "no".
             */
            int is_dup = do_check
                && relation_is_duplicate(rel, las.nb_threads, si);
            const char *comment = is_dup ? "# DUPE " : "";
            FILE *output;
            if (!is_dup)
                cpt++;
            verbose_output_start_batch();   /* lock I/O */
            if (prepend_relation_time) {
                verbose_output_print(0, 1, "(%1.4f) ", seconds() - tt_qstart);
            }
            verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n",
                    i, j, sdata[0].S[x], sdata[1].S[x]);
            for (size_t i_output = 0;
                 (output = verbose_output_get(0, 0, i_output)) != NULL;
                 i_output++) {
                rel.print(output, comment);
		    // once filtering is ok for all Galois cases, 
		    // this entire block would have to disappear
		    if(las.galois != NULL)
			// adding relations on the fly in Galois cases
			add_relations_with_galois(si.cpoly, las.galois, output, comment,
						  &cpt, rel);
		}
            verbose_output_end_batch();     /* unlock I/O */
        }

        /* Build histogram of lucky S[x] values */
        th->rep->report_sizes[sdata[0].S[x]][sdata[1].S[x]]++;

#ifdef  DLP_DESCENT
        if (register_contending_relation(las, si, rel))
            break;
#endif  /* DLP_DESCENT */
    }

#ifdef SUPPORT_LARGE_Q
        mpz_clear(az);
        mpz_clear(bz);
#endif

    for(int side = 0 ; side < 2 ; side++) {
        mpz_clear(norm[side]);
        factor_list_clear(&factors[side]);
        clear_uint32_array (lps_m[side]);
        clear_mpz_array (lps[side]);
    }
}

/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */
// FIXME NOPROFILE_STATIC
int
factor_survivors (timetree_t & timer, thread_data *th, int N, where_am_I & w MAYBE_UNUSED)
{
    CHILD_TIMER(timer, __func__);
    factor_survivors_data F(th, N, w);

    F.search_survivors(timer);
    F.convert_survivors(timer);

    F.purge_buckets(timer);

    F.cofactoring(timer);

    verbose_output_print(0, 3, "# There were %d survivors in bucket %d\n", F.copr, N);
    th->rep->survivors1 += F.copr;
    th->rep->survivors2 += F.copr;
    return F.cpt;
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
void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2) {/*{{{*/
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
}/*}}}*/

/* Move above ? */
/* {{{ process_bucket_region
 * th->id gives the number of the thread: it is supposed to deal with the set
 * of bucket_regions corresponding to that number, ie those that are
 * congruent to id mod nb_thread.
 *
 * The other threads are accessed by combining the thread pointer th and
 * the thread id: the i-th thread is at th - id + i
 */
void * process_bucket_region(timetree_t & timer, thread_data *th)
{
    ACTIVATE_TIMER(timer);
    CHILD_TIMER(timer, __func__);

    where_am_I w MAYBE_UNUSED;
    las_info const & las(*th->plas);
    sieve_info & si(*th->psi);
    uint32_t first_region0_index = th->first_region0_index;
    if (si.toplevel == 1) {
        first_region0_index = 0;
    }

    WHERE_AM_I_UPDATE(w, psi, &si);

    las_report_ptr rep = th->rep;

    unsigned char * S[2];

    unsigned int skiprows = (BUCKET_REGION >> si.conf.logI_adjusted)*(las.nb_threads-1);

    /* This is local to this thread */
    for(int side = 0 ; side < 2 ; side++) {
        thread_side_data &ts = th->sides[side];
        S[side] = ts.bucket_region;
    }

    unsigned char *SS = th->SS;
    memset(SS, 0, BUCKET_REGION);

    /* loop over appropriate set of sieve regions */
    for (uint32_t ii = 0; ii < si.nb_buckets[1]; ii ++)
      {
        uint32_t i = first_region0_index + ii;
        if ((i % las.nb_threads) != (uint32_t)th->id)
            continue;
        WHERE_AM_I_UPDATE(w, N, i);

        if (recursive_descent) {
            /* For the descent mode, we bail out as early as possible. We
             * need to do so in a multithread-compatible way, though.
             * Therefore the following access is mutex-protected within
             * las.tree. */
            if (las.tree.must_take_decision())
                break;
        } else if (exit_after_rel_found) {
            if (rep->reports)
                break;
        }

        for (int side = 0; side < 2; side++)
          {
            WHERE_AM_I_UPDATE(w, side, side);

            sieve_info::side_info & s(si.sides[side]);
            thread_side_data & ts = th->sides[side];
        
            SIBLING_TIMER(timer, "init norms");
            /* Init norms */
            rep->tn[side] -= seconds_thread ();
#ifdef SMART_NORM
	    init_norms_bucket_region(S[side], i, si, side, 1);
#else
	    init_norms_bucket_region(S[side], i, si, side, 0);
#endif
            // Invalidate the first row except (1,0)
            if (side == 0 && i == 0) {
                int pos10 = 1+((si.I)>>1);
                unsigned char n10 = S[side][pos10];
                memset(S[side], 255, si.I);
                S[side][pos10] = n10;
            }
            rep->tn[side] += seconds_thread ();
#if defined(TRACE_K) 
            if (trace_on_spot_N(w.N))
              verbose_output_print(TRACE_CHANNEL, 0, "# After side %d init_norms_bucket_region, N=%u S[%u]=%u\n",
                       side, w.N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            SIBLING_TIMER(timer, "apply buckets");
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
            if (si.toplevel > 1) {
                SIBLING_TIMER(timer, "apply downsorted buckets");
                const bucket_array_t<1, longhint_t> *BAd =
                    th->ws->cbegin_BA<1, longhint_t>(side);
                const bucket_array_t<1, longhint_t> * const BAd_end =
                    th->ws->cend_BA<1, longhint_t>(side);
                for (; BAd != BAd_end; BAd++)  {
                    // FIXME: the updates could come from part 3 as well,
                    // not only part 2.
                    ASSERT_ALWAYS(si.toplevel <= 2);
                    apply_one_bucket(SS, *BAd, ii, ts.fb->get_part(2), w);
                }
            }

            SIBLING_TIMER(timer, "S minus S (1)");

	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
            rep->ttbuckets_apply += seconds_thread();

            SIBLING_TIMER(timer, "small sieve");
            /* Sieve small primes */
            sieve_small_bucket_region(SS, i, s.ssd, ts.ssdpos, si, side, w);

            SIBLING_TIMER(timer, "S minus S (2)");
	    SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
            if (trace_on_spot_N(w.N))
              verbose_output_print(TRACE_CHANNEL, 0,
                      "# Final value on side %d, N=%u rat_S[%u]=%u\n",
                      side, w.N, trace_Nx.x, S[side][trace_Nx.x]);
#endif

            BOOKKEEPING_TIMER(timer);
        }

        /* Factor survivors */
        rep->ttf -= seconds_thread ();
        rep->reports += factor_survivors (timer, th, i, w);
        rep->ttf += seconds_thread ();

        SIBLING_TIMER(timer, "reposition small (re)sieve data");
        /* Reset resieving data */
        for(int side = 0 ; side < 2 ; side++) {
            sieve_info::side_info & s(si.sides[side]);
            thread_side_data & ts = th->sides[side];
            small_sieve_skip_stride(s.ssd, ts.ssdpos, skiprows, si);
            int * b = s.fb_parts_x->rs;
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
void las_report_accumulate_threads_and_display(las_info & las,
        sieve_info & si, las_report_ptr report,
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
    verbose_output_vfprint(0, 1, gmp_vfprintf, "# %lu %s for side-%d (%Zd,%Zd)\n",
              rep->reports,
              las.batch ? "survivor(s) saved" : "relation(s)",
              si.conf.side,
              (mpz_srcptr) si.doing.p,
              (mpz_srcptr) si.doing.r);
    double qtts = qt0 - rep->tn[0] - rep->tn[1] - rep->ttf;
    if (rep->both_even) {
        verbose_output_print(0, 1, "# Warning: found %lu hits with i,j both even (not a bug, but should be very rare)\n", rep->both_even);
    }
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (dont_print_tally && las.nb_threads > 1) {
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


static void declare_usage(param_list pl)/*{{{*/
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
  param_list_decl_usage(pl, "fbc",  "factor base cache file (not yet functional)");
  param_list_decl_usage(pl, "q0",   "left bound of special-q range");
  param_list_decl_usage(pl, "q1",   "right bound of special-q range");
  param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
  param_list_decl_usage(pl, "v",    "(switch) verbose mode, also prints sieve-area checksums");
  param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
  param_list_decl_usage(pl, "t",   "number of threads to use");
  param_list_decl_usage(pl, "sqside", "put special-q on this side");

  param_list_decl_usage(pl, "I",    "set sieving region to 2^I times J");
  param_list_decl_usage(pl, "A",    "set sieving region to 2^A");
  param_list_decl_usage(pl, "skew", "(alias S) skewness");
  param_list_decl_usage(pl, "lim0", "factor base bound on side 0");
  param_list_decl_usage(pl, "lim1", "factor base bound on side 1");
  param_list_decl_usage(pl, "lpb0", "set large prime bound on side 0 to 2^lpb0");
  param_list_decl_usage(pl, "lpb1", "set large prime bound on side 1 to 2^lpb1");
  param_list_decl_usage(pl, "mfb0", "set rational cofactor bound on side 0 2^mfb0");
  param_list_decl_usage(pl, "mfb1", "set rational cofactor bound on side 1 2^mfb1");
  param_list_decl_usage(pl, "lambda0", "lambda value on side 0");
  param_list_decl_usage(pl, "lambda1", "lambda value on side 1");
  param_list_decl_usage(pl, "powlim0", "limit on powers on side 0");
  param_list_decl_usage(pl, "powlim1", "limit on powers on side 1");
  param_list_decl_usage(pl, "ncurves0", "controls number of curves on side 0");
  param_list_decl_usage(pl, "ncurves1", "controls number of curves on side 1");
  param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
  param_list_decl_usage(pl, "skipped", "primes below this bound are not sieved at all");
  param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh");
  param_list_decl_usage(pl, "bkthresh1", "2-level bucket-sieve primes p >= bkthresh1");
  param_list_decl_usage(pl, "bkmult", "multiplier to use for taking margin in the bucket allocation\n");
  param_list_decl_usage(pl, "unsievethresh", "Unsieve all p > unsievethresh where p|gcd(a,b)");

  param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that bounday is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=2");
  param_list_decl_usage(pl, "allow-largesq", "(switch) allows large special-q, e.g. for a DL descent");
  param_list_decl_usage(pl, "exit-early", "once a relation has been found, go to next special-q (value==1), or exit (value==2)");
  param_list_decl_usage(pl, "stats-stderr", "(switch) print stats to stderr in addition to stdout/out file");
  param_list_decl_usage(pl, "stats-cofact", "write statistics about the cofactorization step in file xxx");
  param_list_decl_usage(pl, "file-cofact", "provide file with strategies for the cofactorization step");
  param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
  param_list_decl_usage(pl, "prepend-relation-time", "prefix all relation produced with time offset since beginning of special-q processing");
  param_list_decl_usage(pl, "ondemand-siever-config", "(switch) defer initialization of siever precomputed structures (one per special-q side) to time of first actual use");
  param_list_decl_usage(pl, "dup", "(switch) suppress duplicate relations");
  param_list_decl_usage(pl, "dup-qmax", "limits of q-sieving for 2-sided duplicate removal");
  param_list_decl_usage(pl, "batch", "(switch) use batch cofactorization");
  param_list_decl_usage(pl, "batch0", "side-0 batch file");
  param_list_decl_usage(pl, "batch1", "side-1 batch file");
  param_list_decl_usage(pl, "batchlpb0", "large prime bound on side 0 to be considered by batch cofactorization");
  param_list_decl_usage(pl, "batchlpb1", "large prime bound on side 1 to be considered by batch cofactorization");
  param_list_decl_usage(pl, "batchmfb0", "cofactor bound on side 0 to be considered after batch cofactorization");
  param_list_decl_usage(pl, "batchmfb1", "cofactor bound on side 1 to be considered after batch cofactorization");
  param_list_decl_usage(pl, "batch-print-survivors", "(switch) just print survivros for an external cofactorization");
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
  verbose_decl_usage(pl);
  tdict_decl_usage(pl);
}/*}}}*/

int main (int argc0, char *argv0[])/*{{{*/
{
    double t0, tts, wct;
    unsigned long nr_sq_processed = 0;
    unsigned long nr_sq_discarded = 0;
    int never_discard = 0;      /* only enabled for las_descent */
    double totJ = 0.0;
    double totlogI = 0.0;
    int argc = argc0;
    char **argv = argv0;
    double max_full = 0.;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    cxx_param_list pl;

    declare_usage(pl);

    /* Passing NULL is allowed here. Find value with
     * param_list_parse_switch later on */
    param_list_configure_switch(pl, "-v", NULL);
    param_list_configure_switch(pl, "-ondemand-siever-config", NULL);
    param_list_configure_switch(pl, "-allow-largesq", &allow_largesq);
    param_list_configure_switch(pl, "-stats-stderr", NULL);
    param_list_configure_switch(pl, "-prepend-relation-time", &prepend_relation_time);
    param_list_configure_switch(pl, "-dup", NULL);
    param_list_configure_switch(pl, "-batch", NULL);
    param_list_configure_switch(pl, "-batch-print-survivors", NULL);
    //    param_list_configure_switch(pl, "-galois", NULL);
    param_list_configure_alias(pl, "skew", "S");
    // TODO: All these aliases should disappear, one day.
    // This is just legacy.
    param_list_configure_alias(pl, "fb1", "fb");
#ifdef  DLP_DESCENT
    param_list_configure_switch(pl, "-recursive-descent", &recursive_descent);
    param_list_configure_switch(pl, "-never-discard", &never_discard);
#endif
    tdict_configure_switch(pl);

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

    param_list_parse_int(pl, "adjust-strategy", &adjust_strategy);
    param_list_parse_int(pl, "exit-early", &exit_after_rel_found);
#if DLP_DESCENT
    param_list_parse_double(pl, "grace-time-ratio", &general_grace_time_ratio);
#endif

    las_info las(pl);    /* side effects: prints cmdline and flags */

    /* We have the following dependency chain (not sure the account below
     * is exhaustive).
     *
     * q0 -> q0d (double) -> si.sides[*]->{scale,logmax}
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
        /* Create a default siever instance among las.sievers if needed */
        if (las.default_config_ptr)
            get_sieve_info_from_config(las, *las.default_config_ptr, pl);

        /* Create all siever configurations from the preconfigured hints */
        /* This can also be done dynamically if needed */
        for(unsigned int i = 0 ; i < las.hint_table.size() ; i++) {
            descent_hint & h(las.hint_table[i]);
            get_sieve_info_from_config(las, h.conf, pl);
        }
        verbose_output_print(0, 1, "# Done creating cached siever configurations\n");
    }

    tune_las_memset();
    
    /* We allocate one more workspace per kind of bucket than there are
       threads, so that threads have some freedom in avoiding the fullest
       bucket array. With only one thread, no balancing needs to be done,
       so we use only one bucket array. */
    const size_t nr_workspaces = las.nb_threads + ((las.nb_threads > 1)?1:0);
    thread_workspaces *workspaces = new thread_workspaces(nr_workspaces, 2, las);

    if (las.batch) {
        ASSERT_ALWAYS(las.default_config_ptr);
        siever_config const & sc0(*las.default_config_ptr);
        int lpb[2] = { sc0.sides[0].lpb, sc0.sides[1].lpb, };
        param_list_parse_int(pl, "batchlpb0", &(lpb[0]));
        param_list_parse_int(pl, "batchlpb1", &(lpb[1]));
        for(int side = 0 ; side < 2 ; side++) {
            // the product of primes up to B takes \log2(B)-\log\log 2 /
            // \log 2 bits. The added constant is 0.5287.
            if (lpb[side] + 0.5287 >= 31 + log2(GMP_LIMB_BITS)) {
                fprintf(stderr, "Gnu MP cannot deal with primes product that large (max 37 bits, asked for batchlpb%d=%d)\n", side, lpb[side]);
                abort();
            } else if (lpb[side] + 0.5287 >= 34) {
                fprintf(stderr, "Gnu MP's mpz_inp_raw and mpz_out_raw functions are limited to integers of at most 34 bits (asked for batchlpb%d=%d)\n",side,lpb[side]);
                abort();
            }
        }
    }

    las_report report;
    las_report_init(report);

    t0 = seconds ();
    wct = wct_seconds();

    where_am_I w MAYBE_UNUSED;
    WHERE_AM_I_UPDATE(w, plas, &las);

    /* This timer is not active for now. While most of the accounting is
     * done through timer_special_q, it will be used mostly for
     * collecting the info on the time spent.
     */
    timetree_t global_timer;

    /* {{{ Doc on todo list handling
     * The function las_todo_feed behaves in different
     * ways depending on whether we're in q-range mode or in q-list mode.
     *
     * q-range mode: the special-q's to be handled are specified as a
     * range. Then, whenever the las.todo list almost runs out, it is
     * refilled if possible, up to the limit q1 (-q0 & -rho just gives a
     * special case of this).
     *
     * q-list mode: the special-q's to be handled are always read from a
     * file. Therefore each new special-q to be handled incurs a
     * _blocking_ read on the file, until EOF. This mode is also used for
     * the descent, which has the implication that the read occurs if and
     * only if the todo list is empty. }}} */

    thread_pool *pool = new thread_pool(las.nb_threads);

    for( ; las_todo_feed(las, pl) ; ) {
        if (las_todo_pop_closing_brace(las)) {
            las.tree.done_node();
            if (las.tree.depth() == 0) {
                if (recursive_descent) {
                    /* BEGIN TREE / END TREE are for the python script */
                    fprintf(las.output, "# BEGIN TREE\n");
                    las.tree.display_last_tree(las.output);
                    fprintf(las.output, "# END TREE\n");
                }
                las.tree.visited.clear();
            }
            continue;
        }

        timetree_t timer_special_q;
        double qt0 = seconds();
        tt_qstart = seconds();

        ACTIVATE_TIMER(timer_special_q);

        /* pick a new entry from the stack, and do a few sanity checks */
        las_todo_entry doing = las_todo_pop(las);

        ASSERT_ALWAYS(mpz_poly_is_root(las.cpoly->pols[doing.side], doing.r, doing.p));

        SIBLING_TIMER(timer_special_q, "skew Gauss");

        sieve_range_adjust Adj(doing, las);

        if (!Adj.SkewGauss())
            continue;

#ifndef SUPPORT_LARGE_Q
        if (!Adj.Q.fits_31bits()) { // for fb_root_in_qlattice_31bits
            fprintf (stderr,
                    "Warning, special-q basis is too skewed,"
                    " skipping this special-q."
                    " Define SUPPORT_LARGE_Q to proceed anyway.\n");
            continue;
        }
#endif

        /* Try strategies for adopting the sieving range */

        int should_discard = !Adj.sieve_info_adjust_IJ();

        if (should_discard) {
            if (never_discard) {
                Adj.set_minimum_J_anyway();
            } else {
                verbose_output_vfprint(0, 1, gmp_vfprintf,
                        "# "
                        HILIGHT_START
                        "Discarding side-%d q=%Zd; rho=%Zd;"
                        HILIGHT_END,
                        doing.side,
                        (mpz_srcptr) doing.p,
                        (mpz_srcptr) doing.r);
                verbose_output_print(0, 1,
                         " a0=%" PRId64
                        "; b0=%" PRId64
                        "; a1=%" PRId64
                        "; b1=%" PRId64
                        "; raw_J=%u;\n", 
                        Adj.Q.a0, Adj.Q.b0, Adj.Q.a1, Adj.Q.b1, Adj.J);
                nr_sq_discarded++;
                continue;
            }
        }

        /* With adjust_strategy == 2, we want to display the other
         * values, too. Also, strategy 0 wants strategy 1 to run first.
         */
        if (adjust_strategy != 1)
            Adj.sieve_info_update_norm_data_Jmax();

        if (adjust_strategy >= 2)
            Adj.adjust_with_estimated_yield();

        if (adjust_strategy >= 3) {
            /* Let's change that again. We tell the code to keep logI as
             * it is currently. */
            Adj.sieve_info_update_norm_data_Jmax(true);
        }

        siever_config conf = Adj.config();
        conf.logI_adjusted = Adj.logI;

        /* done with skew gauss ! */

        BOOKKEEPING_TIMER(timer_special_q);

        /* Maybe create a new siever ? */
        sieve_info & si(get_sieve_info_from_config(las, conf, pl));

        si.recover_per_sq_values(Adj);

        si.init_j_div();
        si.init_unsieve_data();

        /* checks the value of J,
         * precompute the skewed polynomials of f(x) and g(x), and also
         * their floating-point versions */

        totJ += si.J;
        totlogI += si.conf.logI_adjusted;


        WHERE_AM_I_UPDATE(w, psi, &si);

        las.tree.new_node(si.doing);
        las_todo_push_closing_brace(las, si.doing.depth);

        nr_sq_processed ++;
        verbose_output_vfprint(0, 1, gmp_vfprintf,
                             "# "
                             HILIGHT_START
                             "Sieving side-%d q=%Zd; rho=%Zd;"
                             HILIGHT_END,
                             si.conf.side,
                             (mpz_srcptr) si.doing.p,
                             (mpz_srcptr) si.doing.r);

        verbose_output_print(0, 1, " a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "; J=%u;",
                             si.qbasis.a0, si.qbasis.b0,
                             si.qbasis.a1, si.qbasis.b1,
                             si.J);
        if (si.doing.depth) {
            verbose_output_print(0, 1, " # within descent, currently at depth %d", si.doing.depth);
        }
        verbose_output_print(0, 1, "\n");
        /* TODO: maybe print that later */
        if (!mpz_probab_prime_p(doing.p, 1)) {
            verbose_output_vfprint(1, 0, gmp_vfprintf,
                    "# Warning, q=%Zd is not prime\n",
                    (mpz_srcptr) doing.p);
        }

        verbose_output_print(0, 2, "# I=%u; J=%u\n", si.I, si.J);
        if (las.verbose >= 2) {
            verbose_output_print (0, 1, "# f_0'(x) = ");
            mpz_poly_fprintf(las.output, si.sides[0].fij);
            verbose_output_print (0, 1, "# f_1'(x) = ");
            mpz_poly_fprintf(las.output, si.sides[1].fij);
        }

        /* essentially update the fij polynomials and the max log bounds */
        si.update_norm_data();

        /* Now we're ready to sieve. We have to refresh some fields
         * in the sieve_info structure, otherwise we'll be polluted by
         * the leftovers from earlier runs.
         */
        si.update(nr_workspaces);


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
    
        workspaces->thrs[0].rep->ttbuckets_fill -= seconds();

        /* Allocate buckets */
        workspaces->pickup_si(si);

        fill_in_buckets_both(timer_special_q, *pool, *workspaces, si);

        pool->accumulate(*timer_special_q.current);

        /* Check that buckets are not more than full.
         * Due to templates and si.toplevel being not constant, need a
         * switch. (really?)
         */
        switch(si.toplevel) {
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
        
        workspaces->thrs[0].rep->ttbuckets_fill += seconds();

        SIBLING_TIMER(timer_special_q, "prepare small sieve");

        /* Prepare small sieve and re-sieve */
        for(int side = 0 ; side < 2 ; side++) {
            sieve_info::side_info & s(si.sides[side]);

            small_sieve_init(s.ssd, las, s.fb_smallsieved.get(), si, side);
            small_sieve_info("small sieve", side, s.ssd);

            small_sieve_extract_interval(s.rsd, s.ssd, s.fb_parts_x->rs);
            small_sieve_info("resieve", side, s.rsd);

            // Initialize small sieve data at the first region of level 0
            // TODO: multithread this? Probably useless...
            for (int i = 0; i < las.nb_threads; ++i) {
                thread_data * th = &workspaces->thrs[i];
                sieve_info::side_info & s(si.sides[side]);
                thread_side_data & ts = th->sides[side];

                uint32_t my_row0 = (BUCKET_REGION >> si.conf.logI_adjusted) * th->id;
                ts.ssdpos = small_sieve_start(s.ssd, my_row0, si);
                ts.rsdpos = small_sieve_copy_start(ts.ssdpos,
                        s.fb_parts_x->rs);
            }
        }
        BOOKKEEPING_TIMER(timer_special_q);

        if (si.toplevel == 1) {
            SIBLING_TIMER(timer_special_q, "process_bucket_region outer container");
            /* Process bucket regions in parallel */
            workspaces->thread_do_using_pool(*pool, &process_bucket_region);
            pool->accumulate(*timer_special_q.current);
        } else {
            SIBLING_TIMER(timer_special_q, "process_bucket_region outer container (non-MT)");
            // Prepare plattices at internal levels
            // TODO: this could be multi-threaded
            plattice_x_t max_area = plattice_x_t(si.J)<<si.conf.logI_adjusted;
            plattice_enumerate_area<1>::value =
                MIN(max_area, plattice_x_t(BUCKET_REGION_2));
            plattice_enumerate_area<2>::value =
                MIN(max_area, plattice_x_t(BUCKET_REGION_3));
            plattice_enumerate_area<3>::value = max_area;
            precomp_plattice_t precomp_plattice;
            for (int side = 0; side < 2; ++side) {
                for (int level = 1; level < si.toplevel; ++level) {
                    const fb_part * fb = si.sides[side].fb->get_part(level);
                    const fb_slice_interface *slice;
                    for (slice_index_t slice_index = fb->get_first_slice_index();
                            (slice = fb->get_slice(slice_index)) != NULL; 
                            slice_index++) {  
                        precomp_plattice[side][level].push_back(
                                slice->make_lattice_bases(si.qbasis, si.conf.logI_adjusted));
                    }
                }
            }

            SIBLING_TIMER(timer_special_q, "process_bucket_region outer container (MT)");
            // Prepare plattices at internal levels

            // Visit the downsorting tree depth-first.
            // If toplevel = 1, then this is just processing all bucket
            // regions.
            uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
            for (uint32_t i = 0; i < si.nb_buckets[si.toplevel]; i++) {
                switch (si.toplevel) {
                    case 2:
                        downsort_tree<1>(timer_special_q, i, i*BRS[2]/BRS[1],
                                *workspaces, *pool, si, precomp_plattice);
                        break;
                    case 3:
                        downsort_tree<2>(timer_special_q, i, i*BRS[3]/BRS[1],
                                *workspaces, *pool, si, precomp_plattice);
                        break;
                    default:
                        ASSERT_ALWAYS(0);
                }
            }
            pool->accumulate(*timer_special_q.current);

            BOOKKEEPING_TIMER(timer_special_q);

            // Cleanup precomputed lattice bases.
            for(int side = 0 ; side < 2 ; side++) {
                for (int level = 1; level < si.toplevel; ++level) {
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

        BOOKKEEPING_TIMER(timer_special_q);

        // Cleanup smallsieve data
        for (int i = 0; i < las.nb_threads; ++i) {
            for(int side = 0 ; side < 2 ; side++) {
                thread_data * th = &workspaces->thrs[i];
                thread_side_data &ts = th->sides[side];
                free(ts.ssdpos);
                free(ts.rsdpos);
            }
        }


#ifdef  DLP_DESCENT
        SIBLING_TIMER(timer_special_q, "descent");
        descent_tree::candidate_relation const & winner(las.tree.current_best_candidate());
        if (winner) {
            /* Even if not going for recursion, store this as being a
             * winning relation. This is useful for preparing the hint
             * file, and also for the initialization of the descent.
             */
            las.tree.take_decision();
            verbose_output_start_batch();
            FILE * output;
            for (size_t i = 0;
                    (output = verbose_output_get(0, 0, i)) != NULL;
                    i++) {
                winner.rel.print(output, "Taken: ");
            }
            verbose_output_end_batch();
            {
                las_todo_entry const & me(si.doing);
                unsigned int n = mpz_sizeinbase(me.p, 2);
                verbose_output_start_batch();
                verbose_output_print (0, 1, "# taking path: ");
                for(int i = 0 ; i < me.depth ; i++) {
                    verbose_output_print (0, 1, " ");
                }
                verbose_output_print (0, 1, "%d@%d ->", n, me.side);
                for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                    int side = winner.outstanding[i].first;
                    relation::pr const & v(winner.outstanding[i].second);
                    unsigned int n = mpz_sizeinbase(v.p, 2);
                    verbose_output_print (0, 1, " %d@%d", n, side);
                }
                if (winner.outstanding.empty()) {
                    verbose_output_print (0, 1, " done");
                }
                verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", me.side,
                        (mpz_srcptr) me.p,
                        (mpz_srcptr) me.r);
                verbose_output_end_batch();
            }
            if (recursive_descent) {
                /* reschedule the possibly still missing large primes in the
                 * todo list */
                for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                    int side = winner.outstanding[i].first;
                    relation::pr const & v(winner.outstanding[i].second);
                    unsigned int n = mpz_sizeinbase(v.p, 2);
                    verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] " HILIGHT_START "pushing side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " to todo list\n", side, (mpz_srcptr) v.p, (mpz_srcptr) v.r, n, side);
                    las_todo_push_withdepth(las, v.p, v.r, side, si.doing.depth + 1);
                }
            }
        } else {
            las_todo_entry const & me(si.doing);
            las.tree.mark_try_again(me.iteration + 1);
            unsigned int n = mpz_sizeinbase(me.p, 2);
            verbose_output_print (0, 1, "# taking path: %d@%d -> loop (#%d)", n, me.side, me.iteration + 1);
            verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", me.side,
                    (mpz_srcptr) me.p,
                    (mpz_srcptr) me.r);
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Failed to find a relation for " HILIGHT_START "side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " (iteration %d). Putting back to todo list.\n", me.side,
                    (mpz_srcptr) me.p,
                    (mpz_srcptr) me.r, n, me.side, me.iteration);
            las_todo_push_withdepth(las, me.p, me.r, me.side, me.depth + 1, me.iteration + 1);
        }
#endif  /* DLP_DESCENT */

        BOOKKEEPING_TIMER(timer_special_q);

        /* clear */
        for(int side = 0 ; side < 2 ; side++) {
            small_sieve_clear(si.sides[side].ssd);
            small_sieve_clear(si.sides[side].rsd);
        }
        qt0 = seconds() - qt0;

        pool->accumulate(timer_special_q);
        timer_special_q.stop();

        if (tdict::global_enable >= 2)
            verbose_output_print (0, 1, "%s", timer_special_q.display().c_str());

        global_timer += timer_special_q;

        las_report_accumulate_threads_and_display(las, si, report, workspaces, qt0);

#ifdef TRACE_K
        trace_per_sq_clear(si);
#endif
        if (exit_after_rel_found > 1 && report->reports > 0)
            break;
      } // end of loop over special q ideals.

    delete pool;

    if (recursive_descent) {
        verbose_output_print(0, 1, "# Now displaying again the results of all descents\n");
        las.tree.display_all_trees(las.output);
    }

    global_timer.start();

    if (las.batch)
      {
        ASSERT_ALWAYS(las.default_config_ptr);
        siever_config const & sc0(*las.default_config_ptr);
        SIBLING_TIMER(global_timer, "batch cofactorization (time is wrong because of openmp)");
	const char *batch0_file, *batch1_file;
	batch0_file = param_list_lookup_string (pl, "batch0");
	batch1_file = param_list_lookup_string (pl, "batch1");
	unsigned long lim[2] = { sc0.sides[0].lim, sc0.sides[1].lim };
	int lpb[2] = { sc0.sides[0].lpb, sc0.sides[1].lpb };
	int batchlpb[2] = { lpb[0], lpb[1]};
	int batchmfb[2] = { sc0.sides[0].lpb, sc0.sides[1].lpb};
        param_list_parse_int(pl, "batchlpb0", &(batchlpb[0]));
        param_list_parse_int(pl, "batchlpb1", &(batchlpb[1]));
        param_list_parse_int(pl, "batchmfb0", &(batchmfb[0]));
        param_list_parse_int(pl, "batchmfb1", &(batchmfb[1]));
	mpz_t batchP[2];
	mpz_init (batchP[0]);
	mpz_init (batchP[1]);
	create_batch_file (batch0_file, batchP[0], lim[0], 1UL << batchlpb[0],
			   las.cpoly->pols[0], las.output, las.nb_threads);
	create_batch_file (batch1_file, batchP[1], lim[1], 1UL << batchlpb[1],
			   las.cpoly->pols[1], las.output, las.nb_threads);
	double tcof_batch = seconds ();
	cofac_list_realloc (las.L, las.L->size);

        mpz_t B[2], L[2], M[2];

        for(int side = 0 ; side < 2 ; side++) {
            mpz_init(B[side]);
            mpz_init(L[side]);
            mpz_init(M[side]);
            mpz_ui_pow_ui(B[side], 2, batchlpb[side]);
            mpz_ui_pow_ui(L[side], 2, lpb[side]);
            mpz_ui_pow_ui(M[side], 2, batchmfb[side]);
        }

	report->reports = find_smooth (las.L, batchP, B, L, M, las.output,
				       las.nb_threads);

        for(int side = 0 ; side < 2 ; side++) {
            mpz_clear (batchP[side]);
            mpz_clear (B[side]);
            mpz_clear (L[side]);
            mpz_clear (M[side]);
        }
	factor (las.L, report->reports, las.cpoly, lpb,
		las.output, las.nb_threads);
	tcof_batch = seconds () - tcof_batch;
	report->ttcof += tcof_batch;
	/* add to ttf since the remaining time will be computed as ttf-ttcof */
	report->ttf += tcof_batch;
      }

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;
    if (adjust_strategy < 2) {
        verbose_output_print (2, 1, "# Average J=%1.0f for %lu special-q's, max bucket fill %f\n",
                totJ / (double) nr_sq_processed, nr_sq_processed, max_full);
    } else {
        verbose_output_print (2, 1, "# Average logI=%1.1f for %lu special-q's, max bucket fill %f\n",
                totlogI / (double) nr_sq_processed, nr_sq_processed, max_full);
    }
    verbose_output_print (2, 1, "# Discarded %lu special-q's out of %u pushed\n",
            nr_sq_discarded, las.nq_pushed);
    tts = t0;
    tts -= report->tn[0];
    tts -= report->tn[1];
    tts -= report->ttf;

    global_timer.stop();
    if (tdict::global_enable >= 1)
        verbose_output_print (0, 1, "%s", global_timer.display().c_str());

    if (las.verbose)
        facul_print_stats (las.output);

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

    if (dont_print_tally && las.nb_threads > 1) 
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

    las.print_cof_stats();

    //
    delete workspaces;

    las_report_clear(report);

    return 0;
}/*}}}*/

