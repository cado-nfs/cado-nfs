#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
#include "cado.h" // IWYU pragma: keep
#include <cmath>
#include <vector>
#include "misc.h"

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

