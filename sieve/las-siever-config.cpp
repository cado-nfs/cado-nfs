#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cmath>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include <gmp.h>
#include "fmt/format.h"

#include "fb.hpp"
#include "fb-types.hpp"
#include "las-config.hpp"
#include "las-multiobj-globals.hpp"
#include "las-siever-config.hpp"
#include "las-side-config.hpp"
#include "las-special-q-task.hpp"
#include "las-special-q-task-tree.hpp"
#include "macros.h"
#include "params.h"
#include "verbose.h"
#include "utils_cxx.hpp"

/* siever_config stuff */

void siever_config::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "I",    "set sieving region to 2^I times J, with J <= 2^(I-1) ; -I x is equivalent to -A (2*x-1)");
    param_list_decl_usage(pl, "A",    "set sieving region to (at most) 2^A");

    siever_side_config::declare_usage(pl);

    param_list_decl_usage(pl, "tdthresh", "trial-divide primes p/r <= ththresh (r=number of roots)");
    param_list_decl_usage(pl, "skipped", "primes below this bound are not sieved at all");
    param_list_decl_usage(pl, "bkthresh", "bucket-sieve primes p >= bkthresh (default 2^I)");
#if MAX_TOPLEVEL >= 2
    param_list_decl_usage(pl, "bkthresh1", "2-level bucket-sieve primes in [bkthresh1,lim] (default=lim, meaning inactive)");
#endif
#if MAX_TOPLEVEL >= 3
    param_list_decl_usage(pl, "bkthresh2", "3-level bucket-sieve primes in [bkthresh2,lim] (default=lim, meaning inactive)");
#endif
    static_assert(MAX_TOPLEVEL == 3);
    param_list_decl_usage(pl, "bkmult", "multiplier to use for taking margin in the bucket allocation\n");
    param_list_decl_usage(pl, "unsievethresh", "Unsieve all p > unsievethresh where p|gcd(a,b)");
    param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that boundary is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=0");
}

/* {{{ Parse default siever config (fill all possible fields). Return
 * true if the parsed siever config is complete and can be used without
 * per-special-q info. */
bool siever_config::parse_default(siever_config & sc, cxx_param_list & pl, int nb_polys)
{
    /* The default config is not necessarily a complete bit of
     * information.
     */

    auto found = siever_side_config::parse(pl, sc.sides, nb_polys);

    bool complete = true;

    for (auto const & s : { "lim", "lpb", "mfb" }) {
        if (found[s] < nb_polys)
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
        verbose_fmt_print(0, 1, "# Interpreting -I {} as meaning -A {}\n",
                I, sc.logA);
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
        verbose_fmt_print(0, 1, "# default siever configuration is incomplete ; required parameters are I, lim[01], lpb[01], mfb[01]\n");

    }

    // Sublattices?
    sc.sublat_bound = 0; // no sublattices by default.
    param_list_parse_uint(pl, "sublat", &sc.sublat_bound);

    /* Parse optional siever configuration parameters */
    param_list_parse_uint(pl, "tdthresh", &(sc.td_thresh));
    param_list_parse_uint(pl, "skipped", &(sc.skipped));

    if (param_list_parse_uint(pl, "unsievethresh", &(sc.unsieve_thresh))) {
        verbose_fmt_print(0, 1, "# Un-sieving primes > {}\n",
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
#if MAX_TOPLEVEL >= 2
    param_list_parse_ulong(pl, "bkthresh1", &(sc.bucket_thresh1));
#endif
#if MAX_TOPLEVEL >= 3
    param_list_parse_ulong(pl, "bkthresh2", &(sc.bucket_thresh2));
#endif
    static_assert(MAX_TOPLEVEL == 3);

    if (sc.bucket_thresh) {
        for (auto & s : sc.sides) {
            if (s.powlim == ULONG_MAX) {
                s.powlim = sc.bucket_thresh - 1;
            }
            /* This message is also printed by
             * fb_factorbase::fb_factorbase
             */
            /*
               verbose_fmt_print(0, 2,
               "# Using default value of {} for -{}\n",
               sc.sides[side].powlim, powlim_params[side]);
               */
        }
    }
    param_list_parse_int(pl, "adjust-strategy", &sc.adjust_strategy);

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
    const fbprime_t fbb = sides[side].lim;
    decltype(fb_factorbase::key_type::thresholds) ret;

    ret.fill(fbb);

    {
        fbprime_t bucket_thresh = this->bucket_thresh;

        if (bucket_thresh == 0)
            bucket_thresh = 1UL << logI;
        if (bucket_thresh < (1UL << logI)) {
            verbose_fmt_print(0, 1, "# Warning: with logI = {},"
                    " we can't have {} as the bucket threshold. Using {}\n",
                    logI,
                    bucket_thresh,
                    1UL << logI);
            bucket_thresh = 1UL << logI;
        }
        ret[0] = std::min(bucket_thresh, fbb);
    }

#if MAX_TOPLEVEL >= 2
    {
        fbprime_t bucket_thresh1 = this->bucket_thresh1;
        if (bucket_thresh1 == 0 || bucket_thresh1 > fbb) bucket_thresh1 = fbb;
        ret[1] = std::max(bucket_thresh1, ret[0]);
    }
#endif

#if MAX_TOPLEVEL >= 3
    {
        fbprime_t bucket_thresh2 = this->bucket_thresh2;
        if (bucket_thresh2 == 0 || bucket_thresh2 > fbb) bucket_thresh2 = fbb;
        ret[2] = std::max(bucket_thresh2, ret[1]);
    }
#endif
    static_assert(MAX_TOPLEVEL == 3);

    return {
         .thresholds = ret,
         .td_thresh = td_thresh,
         .skipped = skipped
    };
}

siever_config siever_config_pool::get_config_for_q(special_q_task const & doing) const /*{{{*/
{
    siever_config config = base;
    unsigned int const bitsize = mpz_sizeinbase(doing.p, 2);
    int const side = doing.side;

    /* Do we have a hint table with specifically tuned parameters,
     * well suited to this problem size ? */
    siever_config const * adapted = get_hint(side, bitsize);
    if (adapted) {
        config = *adapted;
        verbose_fmt_print(0, 1, "# Using parameters from hint list for q~2^{} on side {} [{}]\n", bitsize, side, doing.shortname());
    }

    auto const * tt = dynamic_cast<special_q_task_tree const *>(&doing);

    if (tt != nullptr && tt->try_again) {
        verbose_fmt_print(0, 1, "#\n# NOTE:"
                " we are re-playing this special-q because of"
                " {} previous failed attempt(s)\n", tt->try_again);

        std::string parameters_info;

        int c = tt->try_again;

        if (base.adjust_strategy != 2) {
            config.adjust_strategy = 2;
            parameters_info += fmt::format(" adjust_strategy={}", config.adjust_strategy);
            c--;
        }

        const int dA = std::min(c, max_increase_logA);
        if (dA) {
            config.logA += dA;
            parameters_info += fmt::format(" A={}", config.logA);
            c -= dA;
        }

        const int dlpb = std::min(c, max_increase_lpb);
        if (dlpb) {
            for(size_t side = 0 ; side < config.sides.size() ; side++) {
                auto & s = config.sides[side];
                const double lambda = double_ratio(s.mfb, s.lpb);
                s.lpb += dlpb;
                s.mfb = lround(lambda * s.lpb);
                parameters_info += fmt::format(" lpb{}={} mfb{}={}",
                        side, s.lpb, side, s.mfb);
            }
            c -= dlpb;

            /* Now if at this point c is still >0, we want to abort. How
             * do we do this and let the task collection know that we
             * want to abort?
             */
        }

        verbose_fmt_print(0, 1,
                "# NOTE: modified parameters are: {}\n#\n", 
                parameters_info);
    }

    return config;
}
/* }}} */

void 
siever_config_pool::declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "hint-table", "filename with per-special q sieving data");
    if (dlp_descent)
        param_list_decl_usage(pl, "descent-hint-table", "Alias to hint-table");
    param_list_decl_usage(pl, "fuzzy-descent", "Allow descent steps to start over with fuzzier and fuzzier parameters, with no limit");
    param_list_decl_usage(pl, "descent-max-increase-A", "When retrying steps in the descent, limit the increase of A to this value"); //  (default %d)", max_increase_logA_default);
    param_list_decl_usage(pl, "descent-max-increase-lpb", "When retrying steps in the descent, limit the increase of lpb to this value"); //  (default %d)", max_increase_lpb_default);
}
/*}}}*/

void siever_config_pool::parse_hints_file(const char * filename)/*{{{*/
{
    if (!filename) return;

    char line[1024];
    FILE * f;
    f = fopen(filename, "r");
    if (f == nullptr) {
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
        if (x == nullptr) break;
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
        char const letter MAYBE_UNUSED = *x;
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
        
        for(auto & S : sc.sides) {
            for( ; *x && isspace(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            S.lim = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            S.lpb = z;
            /* recognize this as a double. If it's < 10, we'll consider
             * this means lambda */
            {
                for( ; *x && !isdigit(*x) ; x++) ;
                double const t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                if (t < 10) {
                    S.lambda = t;
                    S.mfb = t * S.lpb;
                    /* Then no "lambda" is allowed */
                    continue;
                } else {
                    S.mfb = t;
                }
            }
            if (*x == ',') {
                for( ; *x && !isdigit(*x) ; x++) ;
                t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                S.lambda = t;
            } else {
                /* this means "automatic" */
                S.lambda = 0;
            }
        }
        for( ; *x ; x++) {
            if (*x == '#')
                break;
            if (!isspace(*x)) {
                fprintf(stderr, "Error: found leftover data in hint file while reading line %d@%d (note that the polynomial file has %zu sides). Leftover text is %s\n", bitsize, side, sc.sides.size(), x);
                exit(EXIT_FAILURE);
            }
        }


        key_type const K(side, bitsize);

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
    param_list_parse(pl, "descent-max-increase-A", max_increase_logA);
    param_list_parse(pl, "descent-max-increase-lpb", max_increase_lpb);

    if (param_list_lookup_string(pl, "fuzzy-descent")) {
        /* well, we do set limits, still */
        max_increase_logA = 16;
        max_increase_lpb = 64;
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

