#include "cado.h" // IWYU pragma: keep
#include <cctype>              // for isspace, isdigit
#include <cerrno>              // for errno
#include <climits>            // for ULONG_MAX
#include <cstdio>             // for fprintf, stderr, fclose, fgets, fopen
#include <cstdlib>            // for exit, strtoul, strtod, EXIT_FAILURE
#include <gmp.h>               // for mpz_sizeinbase
#include "cxx_mpz.hpp"   // for cxx_mpz
#include "fb-types.h"          // for fbprime_t
#include "las-multiobj-globals.hpp"     // for dlp_descent
#include "las-siever-config.hpp"
#include "las-todo-entry.hpp"  // for las_todo_entry
#include "macros.h"            // for ASSERT_ALWAYS, MAYBE_UNUSED
#include "params.h"     // param_list_parse_*
#include "verbose.h"    // verbose_output_print

/* siever_config stuff */

void siever_config::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "I",    "set sieving region to 2^I times J, with J <= 2^(I-1) ; -I x is equivalent to -A (2*x-1)");
    param_list_decl_usage(pl, "A",    "set sieving region to (at most) 2^A");

    siever_side_config::declare_usage(pl);

    param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
    param_list_decl_usage(pl, "skipped", "primes below this bound are not sieved at all");
    param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh (default 2^I)");
    param_list_decl_usage(pl, "bkthresh1", "2-level bucket-sieve primes in [bkthresh1,lim] (default=lim, meaning inactive)");
    param_list_decl_usage(pl, "bkmult", "multiplier to use for taking margin in the bucket allocation\n");
    param_list_decl_usage(pl, "unsievethresh", "Unsieve all p > unsievethresh where p|gcd(a,b)");
}

void siever_config::display(int side, unsigned int bitsize) const /*{{{*/
{
    if (bitsize == 0) return;

    verbose_output_print(0, 2, "# Sieving parameters for q~2^%d on side %d\n",
            bitsize, side);
    /* Strive to keep these output lines untouched */
    verbose_output_print(0, 2,
	    "# Sieving parameters: lim0=%lu lim1=%lu lpb0=%d lpb1=%d\n",
	    sides[0].lim, sides[1].lim,
            sides[0].lpb, sides[1].lpb);
    verbose_output_print(0, 2,
	    "#                     mfb0=%d mfb1=%d\n",
	    sides[0].mfb, sides[1].mfb);
    if (sides[0].lambda != 0 || sides[1].lambda != 0) {
        verbose_output_print(0, 2,
                "#                     lambda0=%1.1f lambda1=%1.1f\n",
            sides[0].lambda, sides[1].lambda);
    }
}/*}}}*/

/* {{{ Parse default siever config (fill all possible fields). Return
 * true if the parsed siever config is complete and can be used without
 * per-special-q info. */
bool siever_config::parse_default(siever_config & sc, cxx_param_list & pl, int n)
{
    /* The default config is not necessarily a complete bit of
     * information.
     */

    auto found = siever_side_config::parse(pl, sc.sides, n);

    bool complete = true;

    for (auto const & s : { "lim", "lpb", "mfb" }) {
        if (found[s] < 2)
            complete = false;
    }

    /*
     * Note that lim0 lim1 powlim0 powlim1 are also parsed from
     * fb_factorbase::fb_factorbase, using only the command line (and no
     * hint file, in particular). This is intended to be the "max" pair
     * of limits.
     */

    /*
     * Note: the stuff about the config being complete or not is mostly
     * rubbish now...
     */
    if (param_list_lookup_string(pl, "A")) {
        complete &= param_list_parse_int  (pl, "A",    &(sc.logA));
        if (param_list_lookup_string(pl, "I")) {
            fprintf(stderr, "# -A and -I are incompatible\n");
            exit(EXIT_FAILURE);
        }
    } else if (param_list_lookup_string(pl, "I")) {
        int I;
        complete &= param_list_parse_int  (pl, "I", &I);
        sc.logA = 2 * I - 1;
        verbose_output_print(0, 1, "# Interpreting -I %d as meaning -A %d\n", I, sc.logA);
    } else {
        complete = false;
    }

#if ULONG_BITS > 32
    for(auto const & s : sc.sides) {
        if (s.lim > 4294967295UL)
        {
            fprintf (stderr, "Error, lim0/lim1 must be < 2^32\n");
            exit (EXIT_FAILURE);
        }
    }
#endif

    if (!complete) {
        verbose_output_print(0, 1, "# default siever configuration is incomplete ; required parameters are I, lim[01], lpb[01], mfb[01]\n");

    }

    // Sublattices?
    sc.sublat_bound = 0; // no sublattices by default.
    param_list_parse_uint(pl, "sublat", &sc.sublat_bound);

    /* Parse optional siever configuration parameters */
    param_list_parse_uint(pl, "tdthresh", &(sc.td_thresh));
    param_list_parse_uint(pl, "skipped", &(sc.skipped));

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
    //
    // The logic that sets bucket_thresh to 2^I by default is done late,
    // namely when we are about to create a new slicing for the factor
    // base.
    /* overrides default only if parameter is given */
    param_list_parse_ulong(pl, "bkthresh", &(sc.bucket_thresh));
    param_list_parse_ulong(pl, "bkthresh1", &(sc.bucket_thresh1));

    if (sc.bucket_thresh) {
        for (auto & s : sc.sides) {
            if (s.powlim == ULONG_MAX) {
                s.powlim = sc.bucket_thresh - 1;
            }
            /* This message is also printed by
             * fb_factorbase::fb_factorbase
             */
            /*
               verbose_output_print(0, 2,
               "# Using default value of %lu for -%s\n",
               sc.sides[side].powlim, powlim_params[side]);
               */
        }
    }

    return complete;
}
/* }}} */

/* returns a set of thresholds that is compatible with the command-line
 * defaults that we see here, as well as the logI that we've just made
 * our mind on using.
 *
 * XXX NOTE XXX : the scale field is set to 0 by this function.
 *
 * XXX NOTE XXX : the nr_workspaces field is set to 0 by this function,
 * because the caller is expected to set that instead.
 *
 *
 */
fb_factorbase::key_type siever_config::instantiate_thresholds(int side) const
{
    fbprime_t fbb = sides[side].lim;
    fbprime_t bucket_thresh = this->bucket_thresh;
    fbprime_t bucket_thresh1 = this->bucket_thresh1;

    if (bucket_thresh == 0)
        bucket_thresh = 1UL << logI;
    if (bucket_thresh < (1UL << logI)) {
        verbose_output_print(0, 1, "# Warning: with logI = %d,"
                " we can't have %lu as the bucket threshold. Using %lu\n",
                logI,
                (unsigned long) bucket_thresh,
                1UL << logI);
        bucket_thresh = 1UL << logI;
    }

    if (bucket_thresh > fbb) bucket_thresh = fbb;
    if (bucket_thresh1 == 0 || bucket_thresh1 > fbb) bucket_thresh1 = fbb;
    if (bucket_thresh > bucket_thresh1) bucket_thresh1 = bucket_thresh;

    return fb_factorbase::key_type {
        {{bucket_thresh, bucket_thresh1, fbb, fbb}},
            td_thresh,
            skipped,
            0,
            0
    };
}
siever_config siever_config_pool::get_config_for_q(las_todo_entry const & doing) const /*{{{*/
{
    siever_config config = base;
    unsigned int bitsize = mpz_sizeinbase(doing.p, 2);
    int side = doing.side;

    /* Do we have a hint table with specifically tuned parameters,
     * well suited to this problem size ? */
    siever_config const * adapted = get_hint(side, bitsize);
    if (adapted) {
        config = *adapted;
        verbose_output_print(0, 1, "# Using parameters from hint list for q~2^%d on side %d [%d@%d]\n", bitsize, side, bitsize, side);
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

void 
siever_config_pool::declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "hint-table", "filename with per-special q sieving data");
    if (dlp_descent)
        param_list_decl_usage(pl, "descent-hint-table", "Alias to hint-table");
}
/*}}}*/

void siever_config_pool::parse_hints_file(const char * filename)/*{{{*/
{
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

        int side;
        unsigned int bitsize;

        z = strtoul(x, &x, 10);
        ASSERT_ALWAYS(z > 0);
        bitsize = z;
        switch(*x++) {
            case '@' :
                side = strtoul(x, &x, 0);
                ASSERT_ALWAYS(side < 2);
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

        key_type K(side, bitsize);

        if (hints.find(K) != hints.end()) {
            fprintf(stderr, "Error: two hints found for %d@%d\n",
                    bitsize, side);
            exit(EXIT_FAILURE);
        }

        hints[K] = h;
    }
    fclose(f);
}
/*}}}*/
siever_config_pool::siever_config_pool(cxx_param_list & pl, int nb_polys)/*{{{*/
{
    default_config_ptr = NULL;
    if (siever_config::parse_default(base, pl, nb_polys))
        default_config_ptr = &base;

    /* support both, since we've got to realize it's not that much
     * attached to sieving. */
    const char * filename = param_list_lookup_string(pl, "hint-table");
    if (dlp_descent) {
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
    }

    parse_hints_file(filename);

    /* Do checks for #30092 */

    if (default_config_ptr) {
        /* no need to do this check if we don't have a default siever
         * config!
         */
        if (base.logA < LOG_BUCKET_REGION) {
            fprintf(stderr, "Error: I=%d (or A=%d) is incompatible with LOG_BUCKET_REGION=%d. Try -B %d\n",
                    (base.logA + 1) / 2, base.logA, LOG_BUCKET_REGION,
                    base.logA);
            exit(EXIT_FAILURE);
        }
    }

    for(auto const & kh : hints) {
        siever_config const & sc(kh.second);
        if (sc.logA < LOG_BUCKET_REGION) {
            fprintf(stderr, "Error: I=%d (or A=%d) is incompatible with LOG_BUCKET_REGION=%d. Try -B %d\n",
                    (sc.logA + 1) / 2, sc.logA, LOG_BUCKET_REGION,
                    sc.logA);
            exit(EXIT_FAILURE);
        }
    }
}/*}}}*/

