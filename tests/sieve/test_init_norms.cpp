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
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "las-config.h"
#include "las-norms.hpp"
#include "memusage.h"

int adjust_strategy = 0;

/*{{{ stuff copied from las.cpp */
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


void ensure_qrange_has_prime_ideals(cxx_mpz const & q0, cxx_mpz & q1, mpz_poly_srcptr f)
{
    /* For random sampling, it's important that for all integers in
     * the range [q0, q1[, their nextprime() is within the range, and
     * that at least one such has roots mod f. Make sure that
     * this is the case.
     */
    cxx_mpz q, q1_orig = q1;
    /* we need to know the limit of the q range */
    for(unsigned long i = 1 ; ; i++) {
        mpz_sub_ui(q, q1, i);
        next_legitimate_specialq(q, q, 0);
        if (mpz_cmp(q, q1) >= 0)
            continue;
        if (mpz_poly_roots(NULL, f, q) > 0)
            break;
        /* small optimization: avoid redoing root finding
         * several times (for all i such that nextprime(q1-i) is
         * the q we've just tested.  */
        q1 = q;
        i = 1;
    }
    /* now q is the largest prime < q1 with f having roots mod q */
    mpz_add_ui (q1, q, 1);
    /* so now if we pick an integer in [q0, q1[, then its nextprime()
     * will be in [q0, q1_orig[, which is what we look for,
     * really.
     */
    if (mpz_cmp(q0, q1) > 0) {
        gmp_fprintf(stderr, "Error: range [%Zd,%Zd[ contains no prime with roots mod f\n", (mpz_srcptr) q0, (mpz_srcptr) q1_orig);
        exit(EXIT_FAILURE);
    }
}

/*}}}*/

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
  param_list_decl_usage(pl, "skew", "(alias S) skewness");

  param_list_decl_usage(pl, "v",    "(switch) verbose mode, also prints sieve-area checksums");

  param_list_decl_usage(pl, "q0",   "left bound of special-q range");
  param_list_decl_usage(pl, "q1",   "right bound of special-q range");
  param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
  param_list_decl_usage(pl, "check-bucket",  "force checking on that particular bucket region");
  param_list_decl_usage(pl, "sqside", "put special-q on this side");
  param_list_decl_usage(pl, "random-sample", "Sample this number of special-q's at random, within the range [q0,q1]");
  param_list_decl_usage(pl, "random-seed", "Use this seed for the random sampling of special-q's (see random-sample)");
  param_list_decl_usage(pl, "nq", "Process this number of special-q's and stop");
  param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");

  param_list_decl_usage(pl, "I",    "set sieving region to 2^I times J");
  param_list_decl_usage(pl, "A",    "set sieving region to 2^A");

  siever_config::declare_usage(pl);

  param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that boundary is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=0");

  param_list_decl_usage(pl, "nfills-speed-test",    "number of bucket region norm fills to simulate per special q");
  param_list_decl_usage(pl, "norm-sides",    "on which sides we should check norms\n");
  param_list_decl_usage(pl, "norm-impls",    "which norm implementations we should check\n");
  param_list_decl_usage(pl, "hush-max-jitter",    "as its name says, do not bother reporting when jitter is below this threshold");
  param_list_decl_usage(pl, "abort-on-jitter",    "exit with failure if jitter exceeds thresholds (one per side)");

  verbose_decl_usage(pl);
}/*}}}*/

int main (int argc0, char *argv0[])/*{{{*/
{
    int argc = argc0;
    char **argv = argv0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    cxx_param_list pl;

    declare_usage(pl);

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    cxx_cado_poly cpoly;

    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }
    if (!cado_poly_read(cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	exit(EXIT_FAILURE);
    }
    cxx_mpz q0, q1, rho;
    int sqside;
    int nq_max = 1;
    int nfills_speed_test = 32;
    int hush_max_jitter = 0;
    int abort_on_jitter[2] = {INT_MAX, INT_MAX};
    int check_bucket = -1;      /* defaults to random pick, can be forced */
    unsigned long seed = 0;

    bool ok = true;

    ok = ok && param_list_parse_mpz(pl, "q0", q0);
    ok = ok && param_list_parse_int(pl, "sqside", &sqside);
    bool okrange = ok && param_list_parse_mpz(pl, "q1", q1);
    bool ok_qrho = ok && param_list_parse_mpz(pl, "rho", rho);
    param_list_parse_int(pl, "check-bucket", &check_bucket);
    param_list_parse_int(pl, "nfills-speed-test", &nfills_speed_test);
    param_list_parse_int(pl, "random-sample", &nq_max);
    param_list_parse_ulong(pl, "random-seed", &seed);
    param_list_parse_int(pl, "hush-max-jitter", &hush_max_jitter);
    param_list_parse_int_and_int(pl, "abort-on-jitter", abort_on_jitter, ",");

    if (okrange == ok_qrho) {
        fprintf(stderr, "Must provide sqside, q0, and either q1 or rho\n");
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }
    if (ok_qrho && param_list_lookup_string(pl, "random-seed")) {
        fprintf(stderr, "-rho and -random-sample are incompatible\n");
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }

    /* These two are mandatory for siever_config::parse_default ; however
     * for our application here, they're totally useless, as we're not
     * initializing a factor base */
    if (!param_list_lookup_string(pl, "lim0"))
        param_list_add_key(pl, "lim0", "0", PARAMETER_FROM_FILE);
    if (!param_list_lookup_string(pl, "lim1"))
        param_list_add_key(pl, "lim1", "0", PARAMETER_FROM_FILE);

    siever_config config_base;
    if (!siever_config::parse_default(config_base, pl)) {
        fprintf(stderr, "Error: please provide a full set of {lim,mfb,lpb}{0,1} parameters\n");
        param_list_print_usage(pl, NULL, stderr);
        exit(EXIT_FAILURE);
    }

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);

    std::vector<int> sides;
    sides.push_back(0);
    sides.push_back(1);

    int * opt_sides;
    unsigned int nopt_sides;
    if (param_list_parse_int_list_size(pl, "norm-sides", &opt_sides, &nopt_sides)) {
        sides.clear();
        for(unsigned int i = 0 ; i < nopt_sides ; i++) {
            sides.push_back(opt_sides[i]);
        }
        free(opt_sides);
    }

    std::vector<std::string> impls;
    impls.push_back("reference");
    impls.push_back("smart");

    char ** opt_impls;
    int nopt_impls;
    if (param_list_parse_string_list_alloc(pl, "norm-impls", &opt_impls, &nopt_impls, ",")) {
        impls.clear();
        for(int i = 0 ; i < nopt_impls ; i++) {
            impls.push_back(opt_impls[i]);
            free(opt_impls[i]);
        }
        free(opt_impls);
    }

    /* That's a maximum only. Currently we have only two lognorm
     * implementations defined. 
     */
#define NCODES  3
    
    ASSERT_ALWAYS(impls.size() <= NCODES);

    unsigned char * S[NCODES];
    double tt[NCODES][2];
    double ttmin[NCODES][2];
    double ttmax[NCODES][2];
    double tt2[NCODES][2];
    for(int c = 0 ; c < NCODES ; c++) {
        for(int s = 0 ; s < 2 ; s++) {
            tt[c][s] = tt2[c][s] = 0;
            ttmin[c][s] = DBL_MAX;
            ttmax[c][s] = DBL_MIN;
        }
    }
    int dd[NCODES][2];
    int ddmin[NCODES][2];
    int ddmax[NCODES][2];
    int dd2[NCODES][2];
    for(int c = 0 ; c < NCODES ; c++) {
        for(int s = 0 ; s < 2 ; s++) {
            dd[c][s] = dd2[c][s] = 0;
            ddmin[c][s] = INT_MAX;
            ddmax[c][s] = INT_MIN;
        }
    }

    if (okrange)
        ensure_qrange_has_prime_ideals(q0, q1, cpoly->pols[sqside]);

    int must_abort = 0;

    for(int qnum = 0 ; qnum < nq_max ; qnum++) {
        /* we don't care much about being truly uniform here */
        cxx_mpz q;
        if (okrange) {
            for(;;) {
                mpz_sub(q, q1, q0);
                mpz_urandomm(q, rstate, q);
                mpz_add(q, q, q0);
                next_legitimate_specialq(q, q, 0);
                cxx_mpz roots[MAX_DEGREE];
                int nroots = mpz_poly_roots ((mpz_t*)roots, cpoly->pols[sqside], q);
                if (nroots) {
                    unsigned long i = gmp_urandomm_ui(rstate, nroots);
                    rho = roots[i];
                    break;
                }
            }
        } else {
            q = q0;
        }
        las_todo_entry doing(q, rho, sqside);

        sieve_range_adjust Adj(doing, cpoly, config_base, 1);

        if (!Adj.SkewGauss())
            continue;

        /* Try strategies for adopting the sieving range */
        int should_discard = !Adj.sieve_info_adjust_IJ();

        if (should_discard) {
                verbose_output_vfprint(0, 1, gmp_vfprintf,
                        "# "
                        "Discarding side-%d q=%Zd; rho=%Zd;",
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
                continue;
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

        /* It's a bit of a hack, yes. If we tinker with I, then we are
         * varying the notion of bucket-sieved primes. So the "default"
         * setting varies, and if there's a user-supplied value, it
         * should by no means fall below the minimum admissible value.
         */
        conf.bucket_thresh = 1UL << conf.logI_adjusted;
        param_list_parse_ulong(pl, "bkthresh", &(conf.bucket_thresh));
        if (conf.bucket_thresh < (1UL << conf.logI_adjusted)) {
            verbose_output_print(0, 1, "# Warning: with logI = %d,"
                    " we can't have %lu as the bucket threshold. Using %lu\n",
                    conf.logI_adjusted,
                    conf.bucket_thresh,
                    1UL << conf.logI_adjusted);
            conf.bucket_thresh = 1UL << conf.logI_adjusted;
        }
        /* done with skew gauss ! */

        verbose_output_vfprint(0, 1, gmp_vfprintf,
                             "# "
                             "Sieving side-%d q=%Zd; rho=%Zd;",
                             conf.side,
                             (mpz_srcptr) doing.p,
                             (mpz_srcptr) doing.r);

        verbose_output_print(0, 1, " a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "; J=%u;",
                             Adj.Q.a0, Adj.Q.b0,
                             Adj.Q.a1, Adj.Q.b1,
                             Adj.J);
        verbose_output_print(0, 1, "\n");
        /* TODO: maybe print that later */
        if (!mpz_probab_prime_p(doing.p, 1)) {
            verbose_output_vfprint(1, 0, gmp_vfprintf,
                    "# Warning, q=%Zd is not prime\n",
                    (mpz_srcptr) doing.p);
        }
        verbose_output_print(0, 2, "# I=%u; J=%u\n", 1U << conf.logI_adjusted, Adj.J);

        std::shared_ptr<lognorm_base> lognorms[NCODES][2];

        for(int side : sides) {
            for(size_t c = 0 ; c < impls.size() ; c++) {
                std::string const & s(impls[c]);
                if (s == "reference") {
                    lognorms[c][side] = std::make_shared<lognorm_reference>(conf, cpoly, side, Adj.Q, Adj.J);
                } else if (s == "smart") {
                    /* For the moment we keep the "smart" code... */
                    lognorms[c][side] = std::make_shared<lognorm_smart>(conf, cpoly, side, Adj.Q, Adj.J);
                } else {
                    fprintf(stderr, "no such implementation: %s\n", s.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }

        int logI = conf.logI_adjusted;
        size_t I = 1UL << logI;
        size_t J = Adj.J;
        int B = 1 << LOG_BUCKET_REGION;
        for(size_t c = 0 ; c < impls.size() ; c++) {
            S[c] = new unsigned char[B + MEMSET_MIN];
            memset(S[c], 0, B);
        }

        /* do a correctness check */
        for(int side : sides) {
            int N = (check_bucket >= 0) ? check_bucket : gmp_urandomm_ui(rstate, iceildiv(I*J, B));     
            for(size_t c = 0 ; c < impls.size() ; c++) {
                lognorms[c][side]->fill(S[c], N);
                if (c == 0) continue;
                int dmin=INT_MAX;
                int dmax=INT_MIN;
                int idmin=-1;
                int idmax=-1;
                double d1=0;
                double d2=0;
                for(int i = 0 ; i < B ; i++) {
                    int d = (int) S[c][i] - (int) S[0][i];
                    if (d < dmin) { dmin = d; idmin = i; }
                    if (d > dmax) { dmax = d; idmax = i; }
                    d1 += d;
                    d2 += d*d;
                }
                if (dmin < ddmin[c][side]) ddmin[c][side] = dmin;
                if (dmax > ddmax[c][side]) ddmax[c][side] = dmax;
                dd[c][side] += d1;
                dd2[c][side] += d2;
                if (dmin < -hush_max_jitter || dmax > hush_max_jitter) {
                    d1 /= B;
                    d2 /= B;
                    fprintf(stderr, "Norm computation disagree for side %d (region %d, %s vs %s); min %d (@%d) max %d (@%d) avg %.1f sdev %.1f\n",
                            side, N,
                            impls[c].c_str(), impls[0].c_str(),
                            dmin, idmin, dmax, idmax, d1, sqrt(d2 - d1*d1));
                }
                if (MAX(-dmin, dmax) > abort_on_jitter[side]) {
                    must_abort = 1;
                    gmp_fprintf(stderr,
                            "###### The jitter reported above will"
                            " cause a program failure\n"
                            "###### Reproduce with:\n"
                            "###### -sqside %d -q0 %Zd -rho %Zd -check-bucket %d\n",
                            sqside, (mpz_srcptr) q, (mpz_srcptr) rho, N);
                    abort();
                }
            }
        }

        /* do a speed test. Since B is essentially fixed, there's
         * no real need to make that adaptative.
         */
        for(int side : sides) {
            gmp_randstate_t rstate2;

            for(size_t c = 0 ; c < impls.size() ; c++) {
                gmp_randinit_set(rstate2, rstate);
                double t = -(double) wct_seconds();
                for(int i = 0 ; i < nfills_speed_test ; i++) {
                    lognorms[c][side]->fill(S[c], gmp_urandomm_ui(rstate2, iceildiv(I*J, B)));
                }
                printf("# Side %d, lognorm %s code: %.3f microseconds per bucket region\n", 
                        side,
                        impls[c].c_str(),
                        1e6 * (t += wct_seconds()) / nfills_speed_test);
                gmp_randclear(rstate2);
                tt[c][side] += t;
                if (t < ttmin[c][side]) ttmin[c][side] = t;
                if (t > ttmax[c][side]) ttmax[c][side] = t;
                tt2[c][side] += t * t;
            }
        }

        for(size_t c = 0 ; c < impls.size() ; c++) {
            delete[] S[c];
        }
    }

    {
        size_t B = 1 << LOG_BUCKET_REGION;
        size_t n = B * nq_max;
        printf("\n# difference values versus %s code over %zu cells\n",
                impls[0].c_str(),
                n);
        for(int side : sides) {
            for(size_t c = 1 ; c < impls.size() ; c++) {
                double a = (double) dd[c][side] / n;
                int amin = ddmin[c][side];
                int amax = ddmax[c][side];
                double a2 = (double) dd2[c][side] / n - a*a;
                printf("# Side %d, %s: %.3f [%d - %d, sd %.3f]\n",
                        side,
                        impls[c].c_str(),
                        a, amin, amax, sqrt(a2));

            }
        }
    }

    if (nfills_speed_test) {
        size_t n = nfills_speed_test * nq_max;
        printf("\n# microseconds per bucket region [average over %zu fills, min-max over %d fills]\n", n, nfills_speed_test);
        for(int side : sides) {
            for(size_t c = 0 ; c < impls.size() ; c++) {
                double a = tt[c][side] / nq_max;
                double amin = ttmin[c][side] / nfills_speed_test;
                double amax = ttmax[c][side] / nfills_speed_test;
                double a2 = tt2[c][side] / nq_max - a*a;
                a /= nfills_speed_test;
                a2 = sqrt(a2) / nfills_speed_test;
                printf("# Side %d, %s : %.3f [%.3f - %.3f, sd %.3f]\n",
                        side,
                        impls[c].c_str(),
                        1e6 * a, 1e6 * amin, 1e6 * amax, 1e6 * a2);
            }
        }
    }

    gmp_randclear(rstate);

    return must_abort ? EXIT_FAILURE : EXIT_SUCCESS;
}/*}}}*/

