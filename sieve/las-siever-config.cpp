#include "cado.h"
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include <cctype>
#include <cerrno>
#include "las-siever-config.hpp"
#include "verbose.h"

/* siever_config stuff */

void siever_config::declare_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "sublat", "modulus for sublattice sieving");
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
    param_list_decl_usage(pl, "dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
    param_list_decl_usage(pl, "dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
}

void siever_config::display() const /*{{{*/
{
    if (bitsize == 0) return;

    verbose_output_print(0, 1, "# Sieving parameters for q~2^%d on side %d\n",
            bitsize, side);
    /* Strive to keep these output lines untouched */
    verbose_output_print(0, 1,
	    "# Sieving parameters: lim0=%lu lim1=%lu lpb0=%d lpb1=%d\n",
	    sides[0].lim, sides[1].lim,
            sides[0].lpb, sides[1].lpb);
    verbose_output_print(0, 1,
	    "#                     mfb0=%d mfb1=%d\n",
	    sides[0].mfb, sides[1].mfb);
    if (sides[0].lambda != 0 || sides[1].lambda != 0) {
        verbose_output_print(0, 1,
                "#                     lambda0=%1.1f lambda1=%1.1f\n",
            sides[0].lambda, sides[1].lambda);
    }
}/*}}}*/

/* {{{ Parse default siever config (fill all possible fields). Return
 * true if the parsed siever config is complete and can be used without
 * per-special-q info. */
bool siever_config::parse_default(siever_config & sc, param_list_ptr pl)
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

    // Sublattices?
    sc.sublat.m = 0; // no sublattices by default.
    param_list_parse_uint(pl, "sublat", &(sc.sublat.m));

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

    long dupqmin[2] = {0, 0};
    param_list_parse_long_and_long(pl, "dup-qmin", dupqmin, ",");
    sc.sides[0].qmin = dupqmin[0];
    sc.sides[1].qmin = dupqmin[1];

    long dupqmax[2] = {LONG_MAX, LONG_MAX};
    param_list_parse_long_and_long(pl, "dup-qmax", dupqmax, ",");
    sc.sides[0].qmax = dupqmax[0];
    sc.sides[1].qmax = dupqmax[1];

    /* Change 0 (not initialized) into LONG_MAX */
    for (int side = 0; side < 2; side ++)
        if (sc.sides[side].qmin == 0)
            sc.sides[side].qmin = LONG_MAX;

    /* if qmin is not given, use lim on the sqside by default */
    if (sc.sides[sc.side].qmin == LONG_MAX)
        sc.sides[sc.side].qmin = sc.sides[sc.side].lim;

    return complete;
}
/* }}} */

siever_config siever_config_pool::get_config_for_q(las_todo_entry const & doing) const /*{{{*/
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
    siever_config config = base;
    config.bitsize = mpz_sizeinbase(doing.p, 2);
    config.side = doing.side;

    /* Do we have a hint table with specifically tuned parameters,
     * well suited to this problem size ? */
    siever_config const * adapted = get_hint(config.side, config.bitsize);
    if (adapted) {
        config = *adapted;
        verbose_output_print(0, 1, "# Using parameters from hint list for q~2^%d on side %d [%d@%d]\n", config.bitsize, config.side, config.bitsize, config.side);
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

siever_config_pool::siever_config_pool(cxx_param_list & pl)/*{{{*/
{
    default_config_ptr = NULL;
    if (siever_config::parse_default(base, pl))
        default_config_ptr = &base;

    /* support both, since we've got to realize it's not that much
     * attached to sieving. */
    const char * filename = param_list_lookup_string(pl, "hint-table");
#if DLP_DESCENT
    const char * filename2 = param_list_lookup_string(pl, "descent-hint-table");
    if (!filename) {
        filename = filename2;
    }
    if (!filename) {
        if (!default_config_ptr) {
            fprintf(stderr,
                    "Error: no default config set, and no hint table either\n");
            exit(EXIT_FAILURE);
        }
        return;
    }
#endif
    if (!filename) return;

    char line[1024];
    FILE * f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        /* There's no point in proceeding, since it would really change
         * the behaviour of the program to do so */
        exit(1);
    }
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
        siever_config & sc(h);

        /* start with the global defaults */
        sc = base;

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

        key_type K(sc.side, sc.bitsize);

        if (hints.find(K) != hints.end()) {
            fprintf(stderr, "Error: two hints found for %d@%d\n",
                    sc.bitsize, sc.side);
            exit(EXIT_FAILURE);
        }

        hints[K] = h;
    }
    fclose(f);
}/*}}}*/
