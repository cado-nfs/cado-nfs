#include "cado.h" // IWYU pragma: keep

// NOLINTBEGIN(bugprone-reserved-identifier,cert-dcl37-c,cert-dcl51-cpp)
#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
// NOLINTEND(bugprone-reserved-identifier,cert-dcl37-c,cert-dcl51-cpp)

#include <cmath>
#include <cstddef>

#include <vector>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>

#include <gmp.h>

#include "cado_expression_parser.hpp"
#include "cxx_mpz.hpp"
#include "getprime.h"
#include "misc.h"
#include "runtime_numeric_cast.hpp"
#include "random_distributions.hpp"

double nprimes_interval(double p0, double p1)
{
    double s0;
    if (p0 <= 1) {
        s0 = 0;
    } else {
        const double l0 = log(p0);
#ifdef HAVE_STDCPP_MATH_SPEC_FUNCS
        s0 = std::expint(l0);
#else
        s0 = p0*(1/l0+1/pow(l0,2)+2/pow(l0,3)+6/pow(l0,4));
#endif
    }

    double s1;
    if (p1 <= 1) {
        s1 = 0;
    } else {
        const double l1 = log(p1);
#ifdef HAVE_STDCPP_MATH_SPEC_FUNCS
        s1 = std::expint(l1);
#else
        s1 = p1*(1/l1+1/pow(l1,2)+2/pow(l1,3)+6/pow(l1,4));
#endif
    }

    return s1 - s0;
}

/* returns the number of primes <= 2^n ; result is exact up to some
 * bound. Note that exact counts are known for larger values as well (see
 * test_prime_count.cpp). Anyway the output that we get is accurate to a
 * relative precision of 2^-20.
 */
double prime_pi_2exp(unsigned int n)
{
    static const double A7053[] = {
        0, 1, 2, 4, 6, 11, 18, 31,
        54, 97, 172, 309, 564, 1028, 1900, 3512,
        6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163,
        1077871, 2063689, 3957809, 7603553,
        14630843, 28192750, 54400028, 105097565,
        203280221, 393615806, 762939111, 1480206279,
        2874398515, 5586502348, 10866266172, 21151907950,
        41203088796, 80316571436, 156661034233, 305761713237,
        /*
        597116381732, 1166746786182, 2280998753949, 4461632979717,
        8731188863470, 17094432576778, 33483379603407, 65612899915304,
        128625503610475,
        */
        };
    /* there's more at https://oeis.org/A007053/b007053.txt */
    if (n < sizeof(A7053) / sizeof(A7053[0]))
        return A7053[n];
    else
        return nprimes_interval(2, ldexp(1, runtime_numeric_cast<int>(n)));
}

/* generate a random i-bit integer that roughly follows the distribution
 * of i-bit prime numbers
 */
double random_along_prime_distribution(unsigned int bits, gmp_randstate_t rstate)
{
    const double n0 = prime_pi_2exp(bits - 1);
    const double n1 = prime_pi_2exp(bits);
    double a, b, n;

    /* we assume the n-th prime is in a*n*log(n)+b, thus we want:
     *    a*n0*log(n0) + b = 2^(bits-1)
     *    a*n1*log(n1) + b = 2^bits
     */
    a = ldexp(1.0, runtime_numeric_cast<int>(bits) - 1) / (n1 * log(n1) - n0 * log(n0));
    b = ldexp(1.0, runtime_numeric_cast<int>(bits)) - a * n1 * log(n1);
    n = n0 + (n1 - n0) * random_uniform(rstate);
    return a * n * log(n) + b;
}

std::vector<unsigned long> subdivide_primes_interval(unsigned long p0, unsigned long p1, size_t n)
{
    std::vector<unsigned long> ret;
    ret.push_back(p0);
    unsigned long const previous = p0;
    double const total_count = nprimes_interval(
            (double) (p0),
            (double) (p1));
    /* by proceeding like this, we're wasting time, since p1 always
     * serves as an endpoint, so that we have a complexity which is
     * roughly n * log(p1-p0). We could have something like log(p1-p0) +
     * 2*log((p1-p0)/2) + 4 * log((p1-p0)/4) + ... which would save a
     * lower-order additive term. No big deal, really.
     */
    for(size_t i = 1 ; i < n ; i++) {
        /* find smallest p such that nprimes_interval(previous, p) >= i *
         * total_count / n ; do simple dichotomy.
         */
        double const target = i * total_count / n;
        unsigned long q0 = previous;
        unsigned long q1 = p1;
        unsigned long q = previous + (p1 - previous) / (n - i);
        for( ; q > q0 ; ) {
            double const r = nprimes_interval(
                    (double) (p0),
                    (double) (q));
            if (r < target)
                q0 = q;
            else
                q1 = q;
            q = (q0 + q1) / 2;
        }
        ret.push_back(q);
    }
    ret.push_back(p1);
    return ret;
}

/* This is meant to be used by C code */
void subdivide_primes_interval_proxy(unsigned long * r, unsigned long p0, unsigned long p1, size_t n)
{
    auto v = subdivide_primes_interval(p0, p1, n);
    std::copy(v.begin(), v.end(), r);
}

struct mpz_parser_traits {
    static constexpr const int accept_literals = 0;
    typedef cxx_mpz type;
    typedef cxx_mpz number_type;
    static void add(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_add(c, a, b);
    }
    static void sub(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_sub(c, a, b);
    }
    static void neg(cxx_mpz & c, cxx_mpz const & a) {
        mpz_neg(c, a);
    }
    static void mul(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_mul(c, a, b);
    }
    static void pow_ui(cxx_mpz & c, cxx_mpz const & a, unsigned long e) {
        mpz_pow_ui(c, a, e);
    }
    static void swap(cxx_mpz & a, cxx_mpz & b) {
        mpz_swap(a, b);
    }
    static void set(cxx_mpz & a, cxx_mpz const & z) {
        mpz_set(a, z);
    }
    static void set_literal_power(cxx_mpz &, std::string const&, unsigned long) {
        // never called. we could do some gymnastics to statically elide
        // this call, but that does not seem to be worth it.
    }
};

typedef cado_expression_parser<mpz_parser_traits> integer_parser;

cxx_mpz mpz_from_expression(const char * value)
{
    std::istringstream is(value);
    integer_parser P;
    P.tokenize(is);
    return P.parse();
}

int mpz_set_from_expression(mpz_ptr f, const char * value)
{
    try {
        mpz_set(f, mpz_from_expression(value));
    } catch (cado_expression_parser_details::token_error const & p) {
        return 0;
    } catch (cado_expression_parser_details::parse_error const & p) {
        return 0;
    }
    return 1;
}


std::vector<std::pair<cxx_mpz, int> > trial_division(cxx_mpz const& n0, unsigned long B, cxx_mpz & cofactor)/*{{{*/
{
    std::vector<std::pair<cxx_mpz, int> > res;
    prime_info pinf;

    prime_info_init (pinf);
    cxx_mpz n = n0;

    /* if n takes k bits it means that n < 2^k. If p^2>=2^k, then we
     * can break. A sufficient condition for this is
     * p^2 >= 2^(2*(ceiling(k/2)))
     * p   >= 2^(  (floor((k+1)/2)))
     * because regardless of the parity of k, this implies our
     * condition
     */
    unsigned long bound_shift = (mpz_sizeinbase(n, 2) + 1) / 2;
    for (unsigned long p = 2; p < B; p = getprime_mt (pinf)) {
        if (bound_shift < ULONG_BITS && p >> bound_shift) {
            if (mpz_cmpabs_ui(n, 1U) > 0 && mpz_cmpabs_ui(n, B) < 0) {
                res.emplace_back(n, 1);
                n = 1U;
                if (mpz_sgn(n0) < 0) {
                    mpz_neg(res.back().first, res.back().first);
                    mpz_neg(n, n);
                }
            }
            break;
        }
        if (!mpz_divisible_ui_p(n, p)) continue;
        int k = 0;
        for( ; mpz_divisible_ui_p(n, p) ; mpz_fdiv_q_ui(n, n, p), k++);
        bound_shift = (mpz_sizeinbase(n, 2) + 1) / 2;
        res.emplace_back(p, k);
    }
    cofactor = n;

    prime_info_clear (pinf); /* free the tables */
    return res;
}
/*}}}*/

/*
 * derived_filename as in misc.h has C linkage, so we can't have a C++
 * overload. OTOH, only replay.c uses the C version.
std::string derived_filename(std::string const & s, const char * what, const char * ext)
{
    char * p = derived_filename(s.c_str(), what, ext);
    std::string r(p);
    p.clear();
}
*/
