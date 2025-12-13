#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <stdexcept>
#include <ostream>

#include <ctime>

#include <gmp.h>
#include "fmt/base.h"

#include "factor.hpp"
#include "getprime.h"
#include "params.h"     // param_list
#include "verbose.h"
#include "cxx_mpz.hpp"
#include "ecm.h"
#include "macros.h"

/* fmt/ranges.h is needed to print std::set<std::pair<cxx_mpz, int> > aka
 * fully_factor::primes, within fully_factor::print_progress
 */
#include "fmt/ranges.h" // IWYU pragma: keep

fully_factor::fully_factor(cxx_mpz n)
  : N(std::move(n))
{
    mpz_abs(cofactor, N); /* easier to work with positive number */
    if (cofactor > 1U) {
        composites.emplace_back(cofactor);
    }

    ecm_init(params);

    /* Remove power of 2 and 3 as ECM does not like them */
    if (mpz_even_p(cofactor)) {
      split_composite(composites.begin(), 2U);
    }
    if (mpz_fdiv_ui(cofactor, 3U) == 0) {
      split_composite(composites.begin(), 3U);
    }
    /* We just did trial division up to B=4 */
    B = 4U;
    BB = 16U;

    print_progress("After init");
}

/* Copy from https://members.loria.fr/PZimmermann/records/ecm/params.html */
const std::vector<std::pair<unsigned int, unsigned int>> fully_factor::ECM_DATA = {
        { 1358, 2 },  /* 30 bits */
        { 1270, 5 }, { 1629, 10 }, { 4537, 10 }, { 12322, 9 }, { 12820, 18 },
        { 21905, 21 }, { 24433, 41 }, { 32918, 66 }, { 64703, 71 },
        { 76620, 119 }, { 155247, 123 }, { 183849, 219 }, { 245335, 321 },
        { 445657, 339 }, /* 100 bits */
        { 643986, 468 }, { 1305195, 439 }, { 1305195, 818 }, { 3071166, 649 },
        { 3784867, 949 }, { 4572523, 1507 }, { 7982718, 1497 },
        { 9267681, 2399 }, { 22025673, 1826 }, { 22025673, 3159 },
        { 26345943, 4532 }, { 35158748, 6076 }, { 46919468, 8177 },
        { 47862548, 14038 }, { 153319098, 7166 }, { 153319098, 12017 },
        { 188949210, 16238 }, { 410593604, 13174 }, { 496041799, 17798 },
        { 491130495, 29584 }, /* 200 bits */
};

void
fully_factor::get_B1_ncurves(unsigned int & B1, unsigned int & ncurves,
                             cxx_mpz const & n)
{
    const unsigned int min_pbits = 30U;
    /* If n has m bits, we look for primes of m/2 bits */
    unsigned int pbits = (mpz_sizeinbase(n, 2) + 1U)/2U;
    pbits = std::max(pbits, min_pbits);
    const unsigned int index = std::min((pbits-min_pbits+4U)/5U,
            (unsigned int) ECM_DATA.size()-1U);
    B1 = ECM_DATA[index].first;
    ncurves = ECM_DATA[index].second;
}

void
fully_factor::print_progress(const char * header) const
{
    verbose_fmt_print(0, 2, "# {}:\n#   primes={}\n#   cofactor={} ({} bits)\n"
                            "#   composites={}\n", header, primes, cofactor,
                            mpz_sizeinbase(cofactor, 2), composites);
}

/* Assumes n has no prime factor < B. So if n is smaller than B^2, it is
 * necessarily a prime.
 */
bool
fully_factor::isprime(cxx_mpz const & n) const
{
    return (mpz_cmpabs_ui(n, 1U) > 0 && mpz_cmpabs(n, BB) < 0)
           || mpz_probab_prime_p(n, isprime_niter);
}

void
fully_factor::write_as_power(cxx_mpz & r, int & e, cxx_mpz const & n)
{
    /* write n as r^e */
    if (mpz_perfect_power_p(n)) {
        /* Unfortunately, GMP does not tell us which power it is... Worst
         * case is n is a power of 2.
         */
        for (e = mpz_sizeinbase(n, 2) - 1; e > 1; --e) {
            if (mpz_root(r, n, e)) {
                break;
            }
        }
        ASSERT_ALWAYS(e > 1);
    } else {
        r = n;
        e = 1;
    }
}

void
fully_factor::push_prime(cxx_mpz const & p)
{
    int exp = mpz_remove(cofactor, cofactor, p);
    ASSERT_ALWAYS(exp > 0);
    verbose_fmt_print(0, 2, "# adding {}^{} to the list of prime "
                            "factors\n", p, exp);
    primes.emplace(p, exp);
}

std::list<cxx_mpz>::iterator
fully_factor::split_composite(std::list<cxx_mpz>::iterator it,
                              cxx_mpz const & f)
{
    ASSERT_ALWAYS(mpz_divisible_p(*it, f));
    if (f > 1U) {
        cxx_mpz cf, g, r;
        int e;

        write_as_power(r, e, f); /* f = r^e */
        bool is_r_prime = isprime(r);

        if (!is_r_prime && f == *it) {
            /* if f is not a prime power and not a proper factor, do nothing */
            return ++it;
        }

        /* composites are coprime, no need to remove p from other composites. */
        mpz_remove(cf, *it, r);
        verbose_fmt_print(0, 2, "# updating coprime from {} to {}\n", *it, cf);

        if (is_r_prime) { /* f is a prime or a prime power */
            push_prime(r);
        } else {
            verbose_fmt_print(0, 2, "# adding {} to the list of composites\n",
                                    r);
            composites.push_front(r);
            /* composites should be coprime, we need to deal with the case where
             * gcd(r, cf) is not 1.
             */
            mpz_gcd (g, r, cf);
            if (g != 1U) {
                cxx_mpz g2, cf_bak = cf;
                mpz_divexact(cf, cf_bak, g);
                mpz_gcd(g2, g, cf);
                ASSERT_ALWAYS(g2 == 1U); //TODO not implemented yet
                verbose_fmt_print(0, 2, "# updating coprime from {} to {}\n",
                                        cf_bak, cf);
                split_composite(composites.begin(), g);
            }
        }

        /* Deal with cf: 1, prime, prime power or composites ? */
        if (cf == 1U) {
            return composites.erase(it);
        } else {
            write_as_power(r, e, cf);
            if (isprime(r)) { /* cf is a prime or a prime power */
                push_prime(r);
                /* *it has been completly factored, can be erased */
                return composites.erase(it);
            } else {
                if (e > 1) {
                    verbose_fmt_print(0, 2, "# updating coprime from {0} to {1}"
                                            " using {0}={1}^{2}\n", cf, r, e);
                }
                *it = r;
                return ++it;
            }
        }
    }
    return ++it;
}

void
fully_factor::trial_division_up_to(unsigned long bound)
{
    if (bound <= B) {
        /* nothing to do: trial division was already done with a larger bound */
        return;
    }

    prime_info pi;
    prime_info_init(pi);

    for (unsigned long p = 2; p < bound; p = getprime_mt(pi)) {
        for (auto it = composites.begin(); it != composites.end(); ) {
            if (mpz_divisible_ui_p(*it, p)) {
                it = split_composite(it, p);
            } else {
                ++it;
            }
        }
        if (is_complete()) {
            break;
        }
    }

    prime_info_clear(pi);

    B = bound;
    BB = B;
    mpz_mul(BB, BB, BB);

    print_progress("After trial division");
}

void fully_factor::use_hints(std::vector<cxx_mpz> const & hints)
{
    cxx_mpz g, r;
    for (auto const & f: hints) {
        mpz_gcd(g, N, f);
        if (g == 1U)
        {
            fmt::print("Warning, hint {} is coprime with N, skipping it\n", f);
            continue;
        }
        for (auto it = composites.begin(); it != composites.end(); ) {
            mpz_gcd(r, *it, g);
            if (r != 1U && r != *it) {
                verbose_fmt_print(0, 2, "# Using hint {} to factor {}\n",
                                        f, *it);
                it = split_composite(it, r);
            } else {
              ++it;
            }
        }
    }
    print_progress("After hints");
}

std::list<cxx_mpz>::iterator
fully_factor::ecmlib_wrapper(std::list<cxx_mpz>::iterator it, unsigned int B1,
                             int method, gmp_randstate_t randgen)
{
    cxx_mpz f;

    /* copy code from ecm_reset here, as it is not available in older versions
     * of gmp-ecm and ECM_VERSION is only available as a string, so it is not
     * easy to do some compile-time comparisons.
     */
    mpz_set_ui (params->sigma, 0U);
    params->B1done = ECM_DEFAULT_B1_DONE;
    mpz_set_ui (params->x, 0U);

    params->method = method;
    params->verbose = std::max(0, ecmlib_verbose);

    if (randgen != nullptr) {
        const unsigned long r = gmp_urandomb_ui(randgen, 32U);
        mpz_set_ui(method == ECM_ECM ? params->sigma : params->x, r);
    }

    char const * s = method == ECM_ECM ? "ECM"
                                       : (method == ECM_PM1 ? "PM1" : "PP1");
    verbose_fmt_print(0, 2, "{}: B1={}; target={}\n", s, B1, *it);

    const int res = ecm_factor(f, *it, B1, params);

    if (res < 0) { /* error => abort */
        throw std::runtime_error("GMP-ECM failed\n");
    } else if (res > 0 && f != *it) { /* factor found */
        verbose_fmt_print(0, 2, "# factor found: {}\n", f);
        return split_composite(it, f);
    } else {
      return ++it;
    }
}

void
fully_factor::apply_method_to_all_composites(int method,
                                          std::vector<unsigned int> const & B1s,
                                          char const * method_str,
                                          gmp_randstate_t randgen)
{
    for (unsigned int B1: B1s) {
        verbose_fmt_print(0, 2, "{}: B1={}\n", method_str, B1);
        for (auto it = composites.begin(); it != composites.end(); ) {
            it = ecmlib_wrapper(it, B1, method, randgen);
        }
        if (is_complete())
            break;
    }
    print_progress(fmt::format("After {}", method_str).c_str());
}

void
fully_factor::do_ECM_based_on_composites_size(unsigned int niter,
                                              gmp_randstate_t randgen)
{
    cxx_mpz bak;
    for (unsigned int iter = 0; iter < niter; ++iter) {
        bool made_progress = true;
        while (!is_complete() && made_progress) {
            made_progress = false;
            for (auto it = composites.begin(); it != composites.end(); ) {
                bak = *it;
                unsigned int B1, ncurves;
                get_B1_ncurves(B1, ncurves, *it);
                for (unsigned int c = 0; c < ncurves ; ++c) {
                    it = ecmlib_wrapper(it, B1, ECM_ECM, randgen);
                    if (it == composites.end() || *it != bak) {
                        /* we make some progress on the composite */
                        made_progress = true;
                        break;
                    }
                }
            }
            if (!made_progress) {
                print_progress("After more ECM");
                break;
            }
        }
    }
}

void
fully_factor::factor_using_default_strategy(unsigned long trialdiv_bound,
                                            std::vector<cxx_mpz> const & hints,
                                            unsigned long PM1_B1,
                                            unsigned long PP1_B1,
                                            unsigned int nECM_rounds,
                                            gmp_randstate_t randgen)
{
    trial_division_up_to(trialdiv_bound);
    if (!is_complete() && !hints.empty()) {
        use_hints(hints);
    }
    if (!is_complete()) {
        do_PM1(PM1_B1, randgen);
    }
    if (!is_complete()) {
        do_PP1(PP1_B1, randgen);
    }
    /* First round of ECM to try to separate medium factors */
    if (!is_complete()) {
        std::vector<unsigned int> B1s;
        for (double B1d = 150.; B1d <= ECM_DATA[0].first; B1d += sqrt(B1d)) {
            B1s.push_back((unsigned int) B1d);
        }
        apply_method_to_all_composites(ECM_ECM, B1s, "ECM 1st round", randgen);
    }
    /* Second round of ECM */
    if (!is_complete()) {
        do_ECM_based_on_composites_size(nECM_rounds, randgen);
    }
}

std::ostream &
operator<<(std::ostream & o, fully_factor const & fact)
{
    cxx_mpz r = 1U;
    bool first = true;
    for (auto const & p: fact.primes) {
        o << (first ? "" : " * ") << p.first << "^" << p.second;
        first = false;
        for (int i = 0; i < p.second; ++i) {
            mpz_mul(r, r, p.first);
        }
    }
    for (auto const & c: fact.composites) {
        o << (first ? "" : " * ") << c;
        mpz_mul(r, r, c);
        first = false;
    }
    mpz_divexact(r, fact.N, r);
    if (r != 1U && r != -1) {
        o << " * " << r;
    }
    return o;
}


static void
usage (cxx_param_list & pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


struct command_line
{
    cxx_mpz N;
    unsigned long trialdiv_bound = 1UL << 24U;
    int isprime_niter = 15; /* number of iterations for primality testing */
    unsigned long PM1_B1 = 5000;
    unsigned long PP1_B1 = 3000;
    unsigned int effort = 2U;
    std::vector<cxx_mpz> hints;
    int verbosity_level = 1; /* each -v on command line increases it by 1 */
    unsigned int seed = 0;

    static void declare_usage(cxx_param_list & pl) {
        param_list_usage_header(pl, "Try as much as possible to fully factor "
                                    "the input integer. The sign of N is not "
                                    "considered.\n");
        param_list_decl_usage(pl, "N", "number to factor");
        param_list_decl_usage(pl, "trialdiv_bound", "bound for trial division");
        param_list_decl_usage(pl, "isprime-niter",
                              "number of iterations for primality testing");
        param_list_decl_usage(pl, "seed", "seed for random generator");
        param_list_decl_usage(pl, "pm1-B1", "B1 used for P-1 algo");
        param_list_decl_usage(pl, "pp1-B1", "B1 used for P+1 algo");
        param_list_decl_usage(pl, "effort", "number of ECM round");
        param_list_decl_usage(pl, "v", "enable verbose output");
        param_list_decl_usage(pl, "hints", "list of possible factors");
        verbose_decl_usage(pl);
    }

    void configure_switches(cxx_param_list & pl)
    {
        param_list_configure_switch(pl, "-v", &verbosity_level);
    }

    void lookup_parameters(cxx_param_list & pl)
    {
        param_list_parse(pl, "N", N);
        param_list_parse(pl, "trialdiv_bound", trialdiv_bound);
        param_list_parse(pl, "isprime_niter", isprime_niter);
        param_list_parse(pl, "pm1-B1", PM1_B1);
        param_list_parse(pl, "pp1-B1", PP1_B1);
        param_list_parse(pl, "effort", effort);
        param_list_parse(pl, "hints", hints);
        if (!param_list_parse(pl, "seed", seed)) {
            seed = time(nullptr);
        }
    }

    void check_inconsistencies(const char * argv0, cxx_param_list & pl) const
    {
        if (!N) {
            fmt::print(stderr, "Error, -N should not be 0\n");
            usage(pl, argv0);
        }

        if (trialdiv_bound < 2) {
            fmt::print(stderr, "Error, -trialdiv_bound should be >= 2\n");
            usage(pl, argv0);
        } else if (trialdiv_bound < (1U << 20)) {
            fmt::print(stderr, "Warning, we recommend that -trialdiv_bound be "
                               ">= 2^20\n");
        }

        if (isprime_niter < 1) {
            fmt::print(stderr, "Error, -isprime-niter should be >= 1\n");
            usage(pl, argv0);
        }
    }
};

int main(int argc, char const * argv[])
{
    const char *argv0 = argv[0];

    cxx_param_list pl;
    command_line cmdline;

    command_line::declare_usage(pl);
    cmdline.configure_switches(pl);

    argv++, argc--;
    if (argc == 0)
      usage(pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        usage(pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    cmdline.lookup_parameters(pl);
    if (param_list_warn_unused(pl)) {
        usage(pl, argv0);
    }
    cmdline.check_inconsistencies(argv0, pl);

    verbose_output_init(2);
    verbose_output_add(0, stdout, cmdline.verbosity_level);
    verbose_output_add(1, stderr, 1);

    verbose_fmt_print(0, 1, "N={} ({} bits)\nseed={}\n",
                       cmdline.N, mpz_sizeinbase(cmdline.N, 2), cmdline.seed);

    gmp_randstate_t randgen;
    gmp_randinit_default(randgen);
    gmp_randseed_ui(randgen, cmdline.seed);

    fully_factor fact(cmdline.N);
    fact.isprime_niter = cmdline.isprime_niter;
    fact.ecmlib_verbose = cmdline.verbosity_level - 2;

    fact.factor_using_default_strategy(cmdline.trialdiv_bound,cmdline.hints,
                                       cmdline.PM1_B1, cmdline.PP1_B1,
                                       cmdline.effort, randgen);

    fmt::print("factorization: {}\n", fact);

    gmp_randclear(randgen);
    verbose_output_clear();

    return fact.is_complete() ? EXIT_SUCCESS : EXIT_FAILURE;
}
