#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
#include "cado.h" // IWYU pragma: keep
#include <cmath>
#include <vector>
#include <sstream>
#include "cxx_mpz.hpp"
#include "cado_expression_parser.hpp"

#include "misc.h"
#include "getprime.h"

double nprimes_interval(double p0, double p1)
{
#ifdef HAVE_STDCPP_MATH_SPEC_FUNCS
    return std::expint(log(p1)) - std::expint(log(p0));
#else
    /* that can't be sooo wrong... */
    double l0 = log(p0);
    double l1 = log(p1);
    double s1 = p1*(1/l1+1/pow(l1,2)+2/pow(l1,3)+6/pow(l1,4));
    double s0 = p0*(1/l0+1/pow(l0,2)+2/pow(l0,3)+6/pow(l0,4));
    return s1 - s0;
#endif
}

std::vector<unsigned long> subdivide_primes_interval(unsigned long p0, unsigned long p1, size_t n)
{
    std::vector<unsigned long> ret;
    ret.push_back(p0);
    unsigned long previous = p0;
    double total_count = nprimes_interval(p0, p1);
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
        double target = i * total_count / n;
        unsigned long q0 = previous;
        unsigned long q1 = p1;
        unsigned long q = previous + (p1 - previous) / (n - i);
        for( ; q > q0 ; ) {
            double r = nprimes_interval(p0, q);
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
    void add(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_add(c, a, b);
    }
    void sub(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_sub(c, a, b);
    }
    void mul(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_mul(c, a, b);
    }
    void pow_ui(cxx_mpz & c, cxx_mpz const & a, unsigned long e) {
        mpz_pow_ui(c, a, e);
    }
    void swap(cxx_mpz & a, cxx_mpz & b) {
        mpz_swap(a, b);
    }
    void set_mpz(cxx_mpz & a, cxx_mpz const & z) {
        mpz_set(a, z);
    }
    void set_literal_power(cxx_mpz &, char, unsigned long) {
        // never called. we could do some gymnastics to statically elide
        // this call, but that does not seem to be worth it.
    }
};

typedef cado_expression_parser<mpz_parser_traits> integer_parser;

int mpz_set_from_expression(mpz_ptr f, const char * value)
{
    std::istringstream is(value);

    integer_parser P;
    try {
        P.tokenize(is);
        cxx_mpz tmp = P.parse();
        mpz_set(f, tmp);
    } catch (integer_parser::token_error const & p) {
        return 0;
    } catch (integer_parser::parse_error const & p) {
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
        if (bound_shift < ULONG_BITS && p >> bound_shift)
            break;
        if (!mpz_divisible_ui_p(n, p)) continue;
        int k = 0;
        for( ; mpz_divisible_ui_p(n, p) ; mpz_fdiv_q_ui(n, n, p), k++);
        bound_shift = (mpz_sizeinbase(n, 2) + 1) / 2;
        res.push_back(std::make_pair(cxx_mpz(p), k));
    }
    // cout << "remaining discriminant " << n << "\n";
    cofactor = n;
    prime_info_clear (pinf); /* free the tables */
    return res;
}
/*}}}*/

