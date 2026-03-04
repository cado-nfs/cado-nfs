#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <climits>

#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#ifdef HAVE_MINGW
#include <fcntl.h>
#endif

#include <gmp.h>

#include "gmp_aux.h"
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "las-config.hpp"
#include "las-coordinates.hpp"
#include "las-norms.hpp"
#include "las-siever-config.hpp"
#include "special-q.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "params.h"
#include "rootfinder.h"
#include "timing.h"
#include "verbose.h"

static int adjust_strategy = 0;

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


static void ensure_qrange_has_prime_ideals(cxx_mpz const & q0, cxx_mpz & q1, mpz_poly_srcptr f)
{
    /* For random sampling, it's important that for all integers in
     * the range [q0, q1[, their nextprime() is within the range, and
     * that at least one such has roots mod f. Make sure that
     * this is the case.
     */
    cxx_mpz q, q1_orig = q1;
    cxx_gmp_randstate rstate;
    /* we need to know the limit of the q range */
    for(unsigned long i = 1 ; ; i++) {
        mpz_sub_ui(q, q1, i);
        next_legitimate_specialq(q, q, 0);
        if (mpz_cmp(q, q1) >= 0)
            continue;
        if (mpz_poly_roots(nullptr, f, q, rstate) > 0)
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

static void declare_usage(cxx_param_list & pl)/*{{{*/
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
  param_list_decl_usage(pl, "skew", "skewness");

  param_list_decl_usage(pl, "v",    "verbose mode, also prints sieve-area checksums");

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

  siever_config::declare_usage<NFS>(pl);

  param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that boundary is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=0");

  param_list_decl_usage(pl, "nfills-speed-test",    "number of bucket region norm fills to simulate per special q");
  param_list_decl_usage(pl, "norm-sides",    "on which sides we should check norms\n");
  param_list_decl_usage(pl, "norm-impls",    "which norm implementations we should check\n");
  param_list_decl_usage(pl, "hush-max-jitter",    "as its name says, do not bother reporting when jitter is below this threshold");
  param_list_decl_usage(pl, "abort-on-jitter",    "exit with failure if jitter exceeds thresholds (one per side)");

  verbose_decl_usage(pl);
}/*}}}*/

// coverity[root_function]
int main(int argc0, char const * argv0[])
    /*{{{*/
{
    int argc = argc0;
    const char **argv = argv0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif
    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    cxx_param_list pl;

    declare_usage(pl);
    param_list_decl_usage(pl, "log-bucket-region", "set bucket region to 2^x");
    param_list_configure_alias(pl, "log-bucket-region", "B");

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    cxx_cado_poly cpoly(pl);

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
    bool const okrange = ok && param_list_parse_mpz(pl, "q1", q1);
    bool const ok_qrho = ok && param_list_parse_mpz(pl, "rho", rho);
    param_list_parse_int(pl, "check-bucket", &check_bucket);
    param_list_parse_int(pl, "nfills-speed-test", &nfills_speed_test);
    param_list_parse_int(pl, "random-sample", &nq_max);
    param_list_parse_ulong(pl, "random-seed", &seed);
    param_list_parse_int(pl, "hush-max-jitter", &hush_max_jitter);
    param_list_parse_int_and_int(pl, "abort-on-jitter", abort_on_jitter, ",");
    param_list_parse_int(pl, "log-bucket-region", &LOG_BUCKET_REGION);
    set_LOG_BUCKET_REGION();

    if (okrange == ok_qrho) {
        fprintf(stderr, "Must provide sqside, q0, and either q1 or rho\n");
        param_list_print_usage(pl, nullptr, stderr);
        exit(EXIT_FAILURE);
    }
    if (ok_qrho && param_list_lookup_string(pl, "random-seed")) {
        fprintf(stderr, "-rho and -random-sample are incompatible\n");
        param_list_print_usage(pl, nullptr, stderr);
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
    if (!siever_config::parse_default<NFS>(config_base, pl, cpoly->nb_polys)) {
        fprintf(stderr, "Error: please provide a full set of {lim,mfb,lpb}{0,1} parameters\n");
        param_list_print_usage(pl, nullptr, stderr);
        exit(EXIT_FAILURE);
    }

    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, seed);

    std::vector<int> sides;
    if (!param_list_parse(pl, "norm-sides", sides)) {
        for(int side = 0 ; side < cpoly->nb_polys ; side++)
            sides.push_back(side);
    }

    std::vector<std::string> impls;
    if (!param_list_parse(pl, "norm-impls", impls)) {
        impls.emplace_back("reference");
        impls.emplace_back("smart");
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
    int ddmin[NCODES][2];
    int ddmax[NCODES][2];
    double dd[NCODES][2];
    double dd2[NCODES][2];
    for(int c = 0 ; c < NCODES ; c++) {
        for(int s = 0 ; s < 2 ; s++) {
            dd[c][s] = dd2[c][s] = 0;
            ddmin[c][s] = INT_MAX;
            ddmax[c][s] = INT_MIN;
        }
    }

    /* only the smart layer defines a per-object stats info (number of
     * logapprox endpoints).
     */
    size_t impl_stats[NCODES][2] = {{0,}};

    if (okrange)
        ensure_qrange_has_prime_ideals(q0, q1, cpoly->pols[sqside]);

    for(int qnum = 0 ; qnum < nq_max ; qnum++) {
        /* we don't care much about being truly uniform here */
        cxx_mpz q;
        if (okrange) {
            for(;;) {
                mpz_sub(q, q1, q0);
                mpz_urandomm(q, rstate, q);
                mpz_add(q, q, q0);
                next_legitimate_specialq(q, q, 0);
                auto roots = mpz_poly_roots(cpoly->pols[sqside], q, rstate);
                if (!roots.empty()) {
                    auto const i = gmp_urandomm_ui(rstate, roots.size());
                    rho = roots[i];
                    break;
                }
            }
        } else {
            q = q0;
        }
        special_q doing(q, rho, sqside);

        sieve_range_adjust Adj(doing, cpoly, config_base);

        /* Try strategies for adopting the sieving range */
        int const should_discard = !Adj.sieve_info_adjust_IJ();

        if (should_discard) {
                verbose_fmt_print(0, 1,
                        "# Discarding {}; raw_J={};\n",
                        Adj.Q, Adj.J);
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
        conf.logI = Adj.logI;

        /* done with skew gauss ! */

        verbose_fmt_print(0, 1, "# Sieving {}; J={};\n", Adj.Q, Adj.J);
        /* TODO: maybe print that later */
        if (!mpz_probab_prime_p(doing.p, 1)) {
            verbose_fmt_print(1, 0,
                    "# Warning, q={} is not prime\n",
                    doing.p);
        }
        verbose_fmt_print(0, 2, "# I={}; J={}", 1U << conf.logI, Adj.J);

        std::unique_ptr<lognorm_base> lognorms[NCODES][2];

        for(int const side : sides) {
            for(size_t c = 0 ; c < impls.size() ; c++) {
                std::string const & s(impls[c]);
                if (s == "reference") {
                    lognorms[c][side] = std::make_unique<lognorm_reference>(conf, cpoly, side, Adj.Q, Adj.logI, Adj.J);
                } else if (s == "smart") {
                    lognorms[c][side] = std::make_unique<lognorm_smart>(conf, cpoly, side, Adj.Q, Adj.logI, Adj.J);
                    impl_stats[c][side] += dynamic_cast<lognorm_smart*>(lognorms[c][side].get())->G.endpoints.size();
                } else {
                    fprintf(stderr, "no such implementation: %s\n", s.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }

        int const logI = conf.logI;
        size_t const I = 1UL << logI;
        size_t const J = Adj.J;
        int const B = 1 << LOG_BUCKET_REGION;
        for(size_t c = 0 ; c < impls.size() ; c++) {
            S[c] = new unsigned char[B + MEMSET_MIN];
            memset(S[c], 0, B);
        }

        /* do a correctness check */
        for(int const side : sides) {
            unsigned int const N = (check_bucket >= 0) ? check_bucket : (unsigned int) gmp_urandomm_ui(rstate, iceildiv(((uint64_t) I)*J, B));
            for(size_t c = 0 ; c < impls.size() ; c++) {
                lognorms[c][side]->fill(S[c], N);
                if (c == 0) continue;
                int dmin=INT_MAX;
                int dmax=INT_MIN;
                unsigned int xdmin = UINT_MAX;
                unsigned int xdmax = UINT_MAX;
                double d1=0;
                double d2=0;
                for(int x = 0 ; x < B ; x++) {
                    int const d = (int) S[c][x] - (int) S[0][x];
                    if (d < dmin) { dmin = d; xdmin = x; }
                    if (d > dmax) { dmax = d; xdmax = x; }
                    d1 += d;
                    d2 += d*d;
                }
                ddmin[c][side] = std::min(ddmin[c][side], dmin);
                ddmax[c][side] = std::max(ddmax[c][side], dmax);
                dd[c][side] += d1;
                dd2[c][side] += d2;
                if (dmin < -hush_max_jitter || dmax > hush_max_jitter) {
                    d1 /= B;
                    d2 /= B;
                    int imin; unsigned int jmin; double zmin;
                    convert_Nx_to_ij(imin, jmin, N, xdmin, logI);
                    zmin = (double) imin / jmin;
                    int imax; unsigned int jmax; double zmax;
                    convert_Nx_to_ij(imax, jmax, N, xdmax, logI);
                    zmax = (double) imax / jmax;

                    fprintf(stderr, "Norm computation disagree for side %d"
                            " (region %d, %s vs %s);\n",
                            side, N,
                            impls[c].c_str(), impls[0].c_str());
                        fprintf(stderr, " min %d (@%d == %d,%u ~ %.2f)\n",
                            dmin, xdmin, imin, jmin, zmin);
                        fprintf(stderr, " max %d (@%d == %d,%u ~ %.2f)\n",
                            dmax, xdmax, imax, jmax, zmax);
                        fprintf(stderr, " avg %.1f sdev %.1f\n",
                            d1, sqrt(d2 - d1*d1));
                }
                if (MAX(-dmin, dmax) > abort_on_jitter[side]) {
                    fmt::print(stderr,
                            "###### The jitter reported above will"
                            " cause a program failure\n"
                            "###### Reproduce with:\n"
                            "###### -sqside {} -q0 {} -rho {} -check-bucket {}",
                            sqside, q, rho, N);
                    abort();
                }
            }
        }

        /* do a speed test. Since B is essentially fixed, there's
         * no real need to make that adaptative.
         */
        for(int const side : sides) {
            cxx_gmp_randstate rstate2;

            for(size_t c = 0 ; c < impls.size() ; c++) {
                rstate2 = rstate;
                double t = -wct_seconds();
                for(int i = 0 ; i < nfills_speed_test ; i++) {
                    lognorms[c][side]->fill(S[c], gmp_urandomm_ui(rstate2, iceildiv(((uint64_t) I)*J, B)));
                }
                printf("# Side %d, lognorm %s code: %.3f microseconds per bucket region\n", 
                        side,
                        impls[c].c_str(),
                        1e6 * (t += wct_seconds()) / nfills_speed_test);
                tt[c][side] += t;
                ttmin[c][side] = std::min(ttmin[c][side], t);
                ttmax[c][side] = std::max(ttmax[c][side], t);
                tt2[c][side] += t * t;
            }
        }

        for(size_t c = 0 ; c < impls.size() ; c++) {
            delete[] S[c];
        }
    }

    {
        size_t const B = 1 << LOG_BUCKET_REGION;
        size_t const n = B * nq_max;
        printf("\n# difference values versus %s code over %zu cells\n",
                impls[0].c_str(),
                n);
        for(int const side : sides) {
            for(size_t c = 1 ; c < impls.size() ; c++) {
                double const a = dd[c][side] / (double) n;
                int const amin = ddmin[c][side];
                int const amax = ddmax[c][side];
                double const a2 = dd2[c][side] / (double) n - a*a;
                printf("# Side %d, %s: %.3f [%d - %d, sd %.3f]\n",
                        side,
                        impls[c].c_str(),
                        a, amin, amax, sqrt(a2));

            }
        }
    }

    if (nfills_speed_test) {
        size_t const n = nfills_speed_test * nq_max;
        for(int const side : sides) {
            for(size_t c = 0 ; c < impls.size() ; c++) {
                if (impl_stats[c][side]) {
                    fmt::print("# side {}, {}: approximation using {} lines"
                            " (avg over {} q's)\n",
                            side, impls[c],
                            double(impl_stats[c][side]) / nq_max,
                            nq_max);
                }
            }
        }
        printf("\n# microseconds per bucket region [average over %zu fills, min-max over %d fills]\n", n, nfills_speed_test);
        for(int const side : sides) {
            for(size_t c = 0 ; c < impls.size() ; c++) {
                double a = tt[c][side] / nq_max;
                double const amin = ttmin[c][side] / nfills_speed_test;
                double const amax = ttmax[c][side] / nfills_speed_test;
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

    return EXIT_SUCCESS;
}/*}}}*/

