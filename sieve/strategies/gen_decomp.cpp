#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include "misc.h"

#include "decomp.hpp"
#include "gen_decomp.hpp"
#include "tab_decomp.hpp"

#define PREC 1000000

/* generate a random i-bit integer, distributed with density in 1/log(x). */
static double gen_random(unsigned int i)
{
    double n0 = prime_pi_2exp(i - 1);
    double n1 = prime_pi_2exp(i);
    double a, b, n;

    /* we assume the n-th prime is in a*n*log(n)+b, thus we want:
       a*n0*log(n0) + b = 2^(i-1)
       a*n1*log(n1) + b = 2^i */
    a = ldexp(1.0, i - 1) / (n1 * log(n1) - n0 * log(n0));
    b = ldexp(1.0, i) - a * n1 * log(n1);
    n = n0 + (n1 - n0) * ((double)rand() / (RAND_MAX + 1UL));
    return a * n * log(n) + b;
}

double T[256] = {
    0,
};

/* estimates number of products of primes:
   l[0]-bit * l[1]-bit * ... * l[k-1]-bit
   with product of mfb bits, each prime should be >= lim */
static double psi(unsigned int * l, unsigned int k, unsigned long mfb,
                  unsigned long lim)
{
    unsigned long ok, tot;
    double S, r, rmin, rmax;

    ok = tot = 0;
    S = 1;
    for (unsigned int i = 0; i < k; i++)
        S *= T[l[i]];

    if (S == 0.0)
        return 0;
    rmin = ldexp(1.0, (int)mfb - 1);
    rmax = ldexp(1.0, (int)mfb);
    while (tot < PREC) {
        r = 1;
        for (unsigned int i = 0; i < k; i++) {
            double p = gen_random(l[i]);
            if (p >= (double)lim)
                r *= p;
            else
                r = 0.0;
        }
        tot++;
        if (rmin <= r && r <= rmax)
            ok++;
    }
    // printf ("%lu", mfb);
    for (unsigned int i = 0; i < k; i++) {
        unsigned int j;
        for (j = 1; i >= j && l[i - j] == l[i]; j++)
            ;
        S /= (double)j;
    }
    /* for (i = 0; i < k; i++) */
    /*   printf (" %d", l[i]); */
    /* printf (" %e # %d\n", (double) ok * S / (double) tot, ++count); */
    /* fflush (stdout); */
    return (double)ok * S / (double)tot;
}

tabular_decomp generate_all_decomp(unsigned int mfb, unsigned long lim)
{
    tabular_decomp res;
    unsigned int l[256];
    unsigned int imin = ceil(log2((double)(lim + 1)));
    for (unsigned int i = 1; (i < 256) && i <= mfb; i++) {
        double p0, p1;

        p0 = ldexp(1.0, (int)i - 1);
        p0 = (p0 < (double)lim) ? (double)lim / log((double)lim)
                                : prime_pi_2exp(i - 1);
        p1 = ldexp(1.0, (int)i);
        p1 = (p1 <= (double)lim) ? (double)lim / log((double)lim) : prime_pi_2exp(i);
        T[i] = p1 - p0;
        // if (T[i] != 0) printf ("T[%lu]=%.0f\n", i, T[i]);
    }
    unsigned int kmax = floor(log(ldexp(1.0, (int)mfb)) / log((double)lim));
    for (unsigned int k = 2; k <= kmax; k++) {
        /* a product l[0]-bit * l[1]-bit * ... * l[k-1]-bit can have
           l[0] + l[1] + ... + l[k-1] - (k-1) bits */
        for (unsigned int j = 0; j < k; j++) {
            for (unsigned int i = 0; i < k - 1; i++)
                l[i] = imin;
            l[k - 1] = mfb + j - (k - 1) * imin;
            if (l[k - 1] < imin)
                continue;
            while (1) {
                res.emplace_back(psi(l, k, mfb, lim), l, l + k);

                /* next partition of mfb + j */
                unsigned int t;
                for (t = k - 1; t > 0; t--)
                    /* try to increase l[t-1] and decrease l[t] */
                    if (l[t - 1] + 1 <= l[t]) {
                        l[t - 1]++;
                        unsigned int s = 1;
                        for (unsigned int u = t; u < k; u++) {
                            s -= l[u];
                            l[u] = l[u - 1];
                            s += l[u];
                        }
                        /* the total has increased by s */
                        l[k - 1] -= s;
                        if (l[k - 2] <= l[k - 1])
                            break;
                    }
                if (t == 0)
                    break;
            }
        }
    }
    return res;
}
