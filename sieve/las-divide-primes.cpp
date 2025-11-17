#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <ostream>
#include <sstream>
#include <vector>

#include <gmp.h>

// bucket_primes_t::get_next_update
#include "bucket-push-update.hpp" // IWYU pragma: keep

#include "bucket.hpp"
#include "cxx_mpz.hpp"
#include "fb.hpp"
#include "las-divide-primes.hpp"
#include "las-output.hpp"
#include "las-where-am-i-proxy.hpp"
#include "macros.h"
#include "trialdiv.hpp"
#include "verbose.h"

/*  Trial division */

static void factor_list_add(factor_list_t & fl, uint64_t const p)
{
    fl.push_back(p);
}

static bool factor_list_contains(factor_list_t const & fl, uint64_t const p)
{
    return std::find(fl.begin(), fl.end(), p) != fl.end();
}

// print a comma-separated list of factors.
// returns the number of factor printed (in particular, a comma is needed
// after this output only if the return value is non zero)
int factor_list_fprint(FILE * f, factor_list_t const & fl)
{
    for (size_t i = 0; i < fl.size(); ++i) {
        if (i)
            fprintf(f, ",");
        fprintf(f, "%" PRIx64, fl[i]);
    }
    return fl.size();
}

static int const bucket_prime_stats = 0;

static long nr_bucket_primes = 0;
static long nr_bucket_longhints = 0;
static long nr_div_tests = 0;
static long nr_composite_tests = 0;
static long nr_wrap_was_composite = 0;

void display_bucket_prime_stats()
{
    if (bucket_prime_stats) {
        verbose_output_print(
            2, 1,
            "# nr_bucket_primes = %lu, nr_div_tests = %lu, nr_composite_tests "
            "= %lu, nr_wrap_was_composite = %lu\n",
            nr_bucket_primes, nr_div_tests, nr_composite_tests,
            nr_wrap_was_composite);
    }
}

/* The entries in BP must be sorted in order of increasing x */
static void divide_primes_from_bucket(factor_list_t & fl, mpz_t norm,
                                      unsigned int const N,
                                      unsigned int const x,
                                      bucket_primes_t * BP,
                                      int const very_verbose)
{
    while (!BP->is_end()) {
        bucket_update_t<1, primehint_t> const prime = BP->get_next_update();
        if (prime.x > x) {
            BP->rewind_by_1();
            break;
        }
        if (prime.x == x) {
            if (bucket_prime_stats)
                nr_bucket_primes++;
            unsigned long const p = prime.p;
            if (very_verbose) {
                verbose_output_vfprint(0, 1, gmp_vfprintf,
                                       "# N = %u, x = %d, dividing out prime "
                                       "hint p = %lu, norm = %Zd\n",
                                       N, x, p, norm);
            }
            /* If powers of a prime p get bucket-sieved and more than one such
                power hits, then the second (and later) hints will find a
                cofactor that already had all powers of p divided out by the
                loop below that removes multiplicities. Thus, if a prime does
                not divide, we check whether it was divided out before (and thus
                is in the factors list) */
            if (UNLIKELY(!mpz_divisible_ui_p(norm, p))) {
                if (!factor_list_contains(fl, p)) {
                    verbose_output_print(
                        1, 0,
                        "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                        p, N, x);
                    abort();
                } else {
                    verbose_output_print(
                        0, 2,
                        "# Note (harmless): p = %lu does not divide at (N,x) = "
                        "(%u,%d), was divided out before\n",
                        p, N, x);
                }
            } else
                do {
                    /* Remove powers of prime divisors */
                    factor_list_add(fl, p);
                    mpz_divexact_ui(norm, norm, p);
                    /* Lacking bucket-sieving for powers, we have to check for
                     * divisibility once again */
                } while (mpz_divisible_ui_p(norm, p));
        }
    }
}

/* The entries in BP must be sorted in order of increasing x */
static void divide_hints_from_bucket(factor_list_t & fl, mpz_t norm,
                                     unsigned int const N, unsigned int const x,
                                     bucket_array_complete * purged,
                                     fb_factorbase::slicing const & fbs,
                                     int const very_verbose)
{
    while (!purged->is_end()) {
        bucket_update_t<1, longhint_t> const complete_hint =
            purged->get_next_update();
        if (complete_hint.x > x) {
            purged->rewind_by_1();
            break;
        }
        if (complete_hint.x == x) {
            if (bucket_prime_stats)
                nr_bucket_longhints++;
            fb_slice_interface const & fb_slice =
                fbs[complete_hint.slice_index];
            unsigned long const p = fb_slice.get_prime(complete_hint.hint);
            if (very_verbose) {
                unsigned char const k = fb_slice.get_k(complete_hint.hint);
                verbose_output_print(
                    0, 1,
                    "# N = %u, x = %d, dividing out fb_slice hint, "
                    "index = %lu offset = %lu ",
                    N, x, (unsigned long)complete_hint.slice_index,
                    (unsigned long)complete_hint.hint);
                if (fb_slice.is_general()) {
                    verbose_output_print(0, 1, "(general)");
                } else {
                    verbose_output_print(0, 1, "(%d roots)",
                                         fb_slice.get_nr_roots());
                }
                verbose_output_vfprint(0, 1, gmp_vfprintf,
                                       ", q = %lu^%hhu, norm = %Zd\n", p, k,
                                       norm);
            }
            if (UNLIKELY(!mpz_divisible_ui_p(norm, p))) {
                if (!factor_list_contains(fl, p)) {
                    verbose_output_print(
                        1, 0,
                        "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                        p, N, x);
                    abort();
                } else {
                    verbose_output_print(
                        0, 2,
                        "# Note (harmless): p = %lu (from hint) does not "
                        "divide at (N,x) = (%u,%d), was divided out before\n",
                        p, N, x);
                }
            } else
                do {
                    factor_list_add(fl, p);
                    mpz_divexact_ui(norm, norm, p);
                } while (mpz_divisible_ui_p(norm, p));
        }
    }
}

/* Extract all known primes (from bucket and small sieves) from the norm.
 * It also removes all the tiny factors that were not resieved and are
 * therefore trial-divided. (see -ththresh parameter)
 *
 * Note: there is another function trialdiv() without underscore that
 * does just the second step.
 */
void divide_known_primes(std::vector<uint64_t> & fl, cxx_mpz & norm,
                         unsigned int const N, unsigned int x,
                         bool const handle_2, bucket_primes_t * primes,
                         bucket_array_complete * purged,
                         trialdiv_data const & td, int64_t a, uint64_t b,
                         fb_factorbase::slicing const & fbs)
{
    int const trial_div_very_verbose = extern_trace_on_spot_ab(a, b);

    if (trial_div_very_verbose) {
        verbose_output_start_batch();
        verbose_output_print(
            TRACE_CHANNEL, 0,
            "# divide_known_primes() entry, N = %u, x = %d, a = %" PRId64
            ", b = %" PRIu64 ", norm = ",
            N, x, a, b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "%Zd\n",
                               (mpz_srcptr)norm);
    }

    // handle 2 separately, if it is in fb
    if (handle_2) {
        int const bit = mpz_scan1(norm, 0);
        for (int i = 0; i < bit; ++i)
            fl.push_back(2);
        if (trial_div_very_verbose)
            verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                                   "# x = %d, dividing out 2^%d, norm = %Zd\n",
                                   x, bit, (mpz_srcptr)norm);
        mpz_tdiv_q_2exp(norm, norm, bit);
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket(fl, norm, N, x, primes, trial_div_very_verbose);
    size_t const nf_divide_primes = fl.size();

    // now remove prime hints in "purged". If we had no factor base, then
    // we really should have an empty list here.
    divide_hints_from_bucket(fl, norm, N, x, purged, fbs,
                             trial_div_very_verbose);
    size_t const nf_divide_hints = fl.size();

    if (trial_div_very_verbose)
        verbose_output_vfprint(
            TRACE_CHANNEL, 0, gmp_vfprintf,
            "# x = %d, after dividing out bucket/resieved norm = %Zd\n", x,
            (mpz_srcptr)norm);

    /* Trial divide primes with precomputed tables */

    if (trial_div_very_verbose) {
        std::ostringstream os;
        for (auto p: td)
            os << " " << p.p;
        verbose_output_print(TRACE_CHANNEL, 0, "# Trial division by%s\n",
                             os.str().c_str());
    }

    td.trial_divide(fl, norm);
    size_t const nf_td = fl.size();

    if (trial_div_very_verbose) {
        std::ostringstream os;
        size_t i = 0;
        if (i < nf_divide_primes) {
            os << " [resieved:";
            for (; i < nf_divide_primes; ++i)
                os << " " << fl[i];
            os << "]";
        }
        if (i < nf_divide_hints) {
            os << " [hints:";
            for (; i < nf_divide_hints; ++i)
                os << " " << fl[i];
            os << "]";
        }
        if (i < nf_td) {
            os << " [trialdiv:";
            for (; i < nf_td; ++i)
                os << " " << fl[i];
            os << "]";
        }
        verbose_output_print(TRACE_CHANNEL, 0, "# %zu factors found:%s\n",
                             fl.size(), os.str().c_str());
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                               "# After trialdiv(): norm = %Zd\n",
                               (mpz_srcptr)norm);
    }

    if (trial_div_very_verbose)
        verbose_output_end_batch();
}
/*  */
