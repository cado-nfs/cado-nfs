#include "cado.h"
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <gmp.h>
#include "utils.h"
#include "gf2x.h"
#include "gf2x-fft.h"
#include "math.h"

#include "lingen_mat_types.hpp"

void usage()
{
    fprintf(stderr, "Usage: ./bench_polmatmul [--nrep <k>] <N> <m> <n>\n");
    exit(1);
}

/* Counted in unsigned longs -- must exceed all caches by large, so as to
 * avoid cheating on benches */
#define DATA_POOL_SIZE  (1 << 23)

/* handy globals */
unsigned int nrep = 10;

#if 0
void bench_polmul(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    printf("gf2x_mul (%u x %u):", d1, d2);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        gf2x_mul(h, f, nw1, g, nw2);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);
    printf("\n");

    delete[] f;
    delete[] g;
    delete[] h;
}
#endif

#define NREPS_MAX       1000
#define REPEAT_TIME_MAX 1
#define do_and_noprinttime(v, line) do {		        	\
    clock_t tt = clock();						\
    unsigned int r;							\
    clock_t clocklim = tt + REPEAT_TIME_MAX * CLOCKS_PER_SEC;		\
    for(r = 0 ; r < NREPS_MAX && clock() < clocklim ; r++) {		\
        line;								\
    }									\
    v=(double) (clock()-tt)/r/CLOCKS_PER_SEC;           		\
} while (0)

#define do_and_printtime(v, line, fmt, scale) do {			\
    do_and_noprinttime(v, line);                            \
    printf(fmt, v * scale);                                                 \
} while (0)

#if 0
void bench_c128(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt0;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    c128_info_t o;

    c128_init(o, d1, d2);

    c128_t * tf = c128_alloc(o, 1);
    c128_t * tg = c128_alloc(o, 1);
    c128_t * th = c128_alloc(o, 1);

    printf("c128 (%u x %u):", d1, d2);

    tt0 = cputime();
    double v;

    do_and_printtime(v, c128_dft(o, tf, f, d1));
    do_and_printtime(v, c128_dft(o, tg, g, d2));
    do_and_printtime(v, c128_compose(o, th, tf, tg));
    do_and_printtime(v, c128_ift(o, h, d1 + d2, th));

    c128_clear(o);
    printf(" tot %f", (double) (cputime()-tt0)/nrep);
    printf("\n");
    delete[] f;
    delete[] g;
    delete[] h;
}
#endif

/* Since strassen is a recursive algorithm that falls back on matrices
 * whose sizes are halved, and since for our applications, the matrices
 * do have a large power of 2 in their dimension, we consider different
 * possible dimensions of the form (d1*2^i,d2*2^i) times (d2*2^i,d3*2^i),
 * for all values of i from 0 to some bound, and for all value of d from
 * 1 to some power of 2 (of course 0 for d makes no sense).
 */

#define BITS_IN_DIM_D   3
#define BITS_IN_DIM_I   4

/* This is in case 1 << BITS_IN_DIM_I is unacceptably large. */
#define MAX_I_FOR_TUNING        8    /* < (1 << BITS_IN_DIM_I) */

// dimension must be <= d * 2^i
// d <= (1 << BITS_IN_DIM_D)
// i < (1 << BITS_IN_DIM_I)
//
// one tuning table entry contains
// 3*BITS_IN_DIM_D+BITS_IN_DIM_I bits
//
// representing a max dimension (1 << BITS_IN_DIM_D) * (1 << ((1 <<
// BITS_IN_DIM_I) - 1)) ; however not all dimensions can be represented
// up to this size, for obvious reasons.
//
// 1,3 means at most 256        ;  6 bits per entry index
// 2,3 means at most 512        ;  9 bits per entry index
// 1,4 means at most 65536      ;  7 bits per entry index
// 2,4 means at most 131072     ;  10 bits per entry index
// 3,4 means at most 262144     ;  13 bits per entry index

// contains the value n0 such that for nbits >= n0, we'd better use strassen

#define SELECTOR_NB_INDICES (1 << (3*BITS_IN_DIM_D+BITS_IN_DIM_I))

struct my_strassen_selector {

    unsigned int strassen_threshold[SELECTOR_NB_INDICES];

    unsigned int threshold_index(unsigned int m, unsigned int n, unsigned int p) const {
        unsigned int combined = m|n|p;
        unsigned int b = ctzl(combined);
        if (b >= (1UL << BITS_IN_DIM_I)) {
            // arrange so that we 
            b = (1UL << BITS_IN_DIM_I) - 1;
        }
        m >>= b; m-=1; ASSERT(m < (1 << BITS_IN_DIM_D));
        n >>= b; n-=1; ASSERT(n < (1 << BITS_IN_DIM_D));
        p >>= b; p-=1; ASSERT(p < (1 << BITS_IN_DIM_D));

        unsigned int index = 0;
        index = m;
        index = (index << BITS_IN_DIM_D) | n;
        index = (index << BITS_IN_DIM_D) | p;
        index = (index << BITS_IN_DIM_I) | b;
        return index;
    }

    unsigned int const & threshold(unsigned int m, unsigned int n, unsigned int p) const {
        return strassen_threshold[threshold_index(m,n,p)];
    }

    unsigned int& threshold(unsigned int m, unsigned int n, unsigned int p) {
        return strassen_threshold[threshold_index(m,n,p)];
    }

    int operator()(unsigned int m, unsigned int n, unsigned int p, unsigned int nbits) const {
        int answer = nbits >= threshold(m,n,p);
        // printf("([%u,%u,%u,%u]=>%d)",m,n,p,nbits,answer);
        return answer;
    }
    void dump(const char * name) const {
        printf("#ifdef  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS\n");
        printf("#define %s_STRASSEN_THRESHOLDS_D %u\n", name, BITS_IN_DIM_D);
        printf("#define %s_STRASSEN_THRESHOLDS_I %u\n", name, BITS_IN_DIM_I);
        printf("#define %s_STRASSEN_THRESHOLDS {\t\\\n\t", name);
        unsigned int disp = 0;
        for(unsigned int i = 0 ; i < SELECTOR_NB_INDICES ; i++) {
            if (strassen_threshold[i]) {
                printf(" {%u, %d},", i, (int) strassen_threshold[i]);
                if (++disp % 4 == 0) {
                    printf("\t\\\n\t");
                }
            }
        }
        printf("}\n");
        printf("#else   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */\n");
        char * s = strdup(name);
        for(unsigned int i = 0 ; i < strlen(s) ; i++) {
            if (s[i] >= 'A' && s[i] <= 'Z')
                s[i] |= 0x20;
        }
        printf("template<> unsigned int foo<%s>::default_selector_data[] = {\n",
                s);
        free(s);
        printf("\t/* D = %u */\n", BITS_IN_DIM_D);
        printf("\t/* I = %u */\n", BITS_IN_DIM_I);
        printf("\t");
        disp = 0;
        for(unsigned int i = 0 ; i < SELECTOR_NB_INDICES ; i++) {
            if (strassen_threshold[i]) {
                printf(" [%u] = %d,", i, (int) strassen_threshold[i]);
                if (++disp % 4 == 0) {
                    printf("\t\\\n\t");
                }
            }
        }
        printf("}\n");
        printf("#endif   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */\n");
    }
};


struct logbook {
    double t1;
    double t2;
    double c;
    double i;
};

template<typename T>
struct foo {
    static logbook l;
    static my_strassen_selector s;
    static unsigned int default_selector_data[];
};

template<> logbook foo<fake_fft>::l = logbook();
template<> logbook foo<c128_fft>::l = logbook();
template<> logbook foo<gf2x_tfft>::l = logbook();
template<> my_strassen_selector foo<fake_fft>::s = my_strassen_selector();
template<> my_strassen_selector foo<c128_fft>::s = my_strassen_selector();
template<> my_strassen_selector foo<gf2x_tfft>::s = my_strassen_selector();

#define STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS
#include "strassen-thresholds.h"
/*
template<> unsigned int foo<fake_fft>::default_selector_data[] = {};
template<> unsigned int foo<c128_fft>::default_selector_data[] = {};
template<> unsigned int foo<gf2x_tfft>::default_selector_data[] = {};
*/


template<typename fft_type>
int depth(fft_type const& o, size_t d1, size_t d2, size_t d3)
{
    int k;
    my_strassen_selector& s(foo<fft_type>::s);
    size_t nbits = o.size() * sizeof(typename fft_type::t) * CHAR_BIT;
    for( k = 0 ; (((d1|d2|d3) >> k) & 1UL) == 0 ; k++) {

        if (!s(d1 >> k, d2 >> k, d3 >> k, nbits))
            break;
    }
    return k;
}

template<typename fft_type>
unsigned long nmults(fft_type const& o, size_t d1, size_t d2, size_t d3)
{
    int d = depth(o, d1, d2, d3);
    unsigned long nmuls = d1 * d2 * d3;
    for(int k = 0 ; k < d ; k++) {
        nmuls /= 8;
        nmuls *= 7;
    }
    return nmuls;
}


/* One step:
 *
 * pi matrix (m+n) * (m+n), length ceil(d/2)*m/(m+n) ; two transforms.
 * E matrix m * (m+n), length ceil(d/2)*(1+m/(m+n)) = d-ceil(d/2)*n/(m+n) ; one transform
 *
 * length of E*pi d
 * length of pi*pi d*m/(m+n)
 *
 */

static inline size_t W(size_t x) { return (x + ULONG_BITS - 1) / ULONG_BITS; }
static inline size_t I(size_t x) { return x / ULONG_BITS; }
static inline size_t R(size_t x) { return x % ULONG_BITS; }
static inline unsigned long MASK(size_t x) { return (1UL << R(x)) - 1UL; }

unsigned long * tidy_data(unsigned long * data, size_t n1)
{
    unsigned long * p = data + (random() % (DATA_POOL_SIZE/2));
    if (R(n1)) p[I(n1)]&=MASK(n1);
    return p;
}

    template<typename T>
void fft_times(double& dft1, double& dft2, double& compose, double& ift,
        T& o, unsigned long n1, unsigned long n2, unsigned long * data)
{
    typename T::t * f = o.alloc(1);
    typename T::t * g = o.alloc(1);
    /* Make sure we're timing for at least a second. */
    unsigned long n3 = n1 + n2 - 1;

    do_and_printtime(dft1, o.dft(f, tidy_data(data, n1), n1), " %.2f", 1.0e6);
    do_and_printtime(dft2, o.dft(g, tidy_data(data, n2), n2), " %.2f", 1.0e6);
    do_and_printtime(compose, o.compose(g, g, f), " %.2f", 1.0e6);
    do_and_printtime(ift, o.ift(tidy_data(data, n3), n3, g), " %.2f", 1.0e6);

    o.free(f, 1);
    o.free(g, 1);
}


void randomize(polmat & t)
{
    for(unsigned int i = 0 ; i < t.nrows ; i++) {
        for(unsigned int j = 0 ; j < t.ncols ; j++) {
            mpn_random((mp_limb_t *) t.poly(i,j), W(t.ncoef));
        }
    }
}

    template<typename fft_type>
void randomize(tpolmat<fft_type> & t)
{
    for(unsigned int i = 0 ; i < t.nrows ; i++) {
        for(unsigned int j = 0 ; j < t.ncols ; j++) {
            mpn_random((mp_limb_t *) t.poly(i,j),
                    (t.po->size() * sizeof(typename fft_type::t)) / sizeof(mp_limb_t));
        }
    }
}


#if 0
/* One step:
 *
 * pi matrix (m+n) * (m+n), length d/2*m/(m+n) ; two transforms.
 * E matrix m * (m+n), length d/2*(1+m/(m+n)) = d-d/2*n/(m+n) ; one transform
 *
 * length of E*pi d
 * length of pi*pi d*m/(m+n)
 *
 */

    template<typename T>
void bench_polmatmul(const char * s,
        unsigned int n1, unsigned int n2,
        unsigned int d1, unsigned int d2)
{
    unsigned int nc1 = d1 + 1;
    unsigned int nc2 = d2 + 1;
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);

    printf("polmatmul<%s> (dim %u x %u deg %u x %u)", s, n1, n2, d1, d2);

    polmat f(n1, n2, nc1);
    polmat g(n2, n2, nc2);
    polmat h(n1, n2, nc1 + nc2 - 1);

    randomize(f);
    randomize(g);

    T o(d1, d2);

    tpolmat<T> tf(n1,n2,o);
    tpolmat<T> tg(n2,n2,o);
    tpolmat<T> th(n1,n2,o);

    do_and_printtime(foo<T>::l.t1, transform(tf, f, o, d1));
    do_and_printtime(foo<T>::l.t2, transform(tg, g, o, d2));
    do_and_printtime(foo<T>::l.c, compose(th, tf, tg, o));
    do_and_printtime(foo<T>::l.i, itransform(h, th, o, d1 + d2));

    double t = foo<T>::l.t1 + foo<T>::l.t2 + foo<T>::l.c + foo<T>::l.i;

    printf(" tot %f", t);
    printf("\n");
}
#endif

    template<typename fft_type>
void tune_strassen1(fft_type const& base,
        unsigned int d1, unsigned int d2, unsigned int d3, size_t maxlen)
{
    my_strassen_selector& s(foo<fft_type>::s);

    unsigned int bp = ctzl(d1 | d2 | d3);
    unsigned int dd1 = d1 >> bp;
    unsigned int dd2 = d2 >> bp;
    unsigned int dd3 = d3 >> bp;
    printf(": -- tuning strassen for %ux%u * %ux%u ; max %u levels above %ux%u * %ux%u --\n",d1,d2,d2,d3,bp,dd1,dd2,dd2,dd3);
    ASSERT((dd1-1) < (1 << BITS_IN_DIM_D));
    ASSERT((dd2-1) < (1 << BITS_IN_DIM_D));
    ASSERT((dd3-1) < (1 << BITS_IN_DIM_D));
    // for the very small matrices, we just _cannot_ use strassen because
    // there is an odd dimension hanging around.
    s.threshold(dd1,dd2,dd3) = UINT_MAX;
    printf("Threshold for %ux%u * %ux%u is %u bits [index %u]\n",
            dd1, dd2, dd2, dd3, UINT_MAX,
            s.threshold_index(dd1,dd2,dd3));
    dd1 <<= 1;
    dd2 <<= 1;
    dd3 <<= 1;
    // If we know that for smaller matrices, strassen was paying off from
    // nbits on, then of course it is going to pay off at worst at this
    // size.
    // old_wt stands for old threshold in words.
    unsigned int max_wt = iceildiv(maxlen, ULONG_BITS);
    unsigned int min_wt = 0;
    unsigned int old_wt = max_wt;
    size_t earliest_good_strassen = UINT_MAX;
    for( ; dd1 <= d1 ; dd1 <<= 1, dd2 <<= 1, dd3 <<= 1) {
        printf("Tuning for %ux%u * %ux%u\n",dd1, dd2, dd2, dd3);
        if (old_wt == min_wt) {
            printf("[keeping low threshold]\n");
            s.threshold(dd1,dd2,dd3) = 0;
            continue;
        }
        // in order to avoid accumulating errors, we introduce some
        // looseness here. Setting wt1 to twice the old wt value will
        // force re-checking at the previous cutoff value. If ever this
        // cutoff value was wrong, we'll settle for one which sits above
        // in the interval.
        unsigned int wt1 = 2 * old_wt;
        unsigned int wt0 = min_wt;
        for( ; wt1 - wt0 > 1  && wt0 < max_wt; ) {
            // assuming cubic is better for 64 * t0, ans strassen better
            // for 64 * t1
            unsigned int wt = (wt1 + wt0) / 2;
            // unsigned int t1 = wt1 * ULONG_BITS;
            // unsigned int t0 = wt0 * ULONG_BITS;
            unsigned int t = wt * ULONG_BITS;
            if (t == 0) { t = 2; }
            fft_type o(t / 2, t / 2, base);
            typedef typename fft_type::t ft;
            unsigned int tr = (o.size() * sizeof(ft)) * CHAR_BIT;
            tpolmat<fft_type> tf(dd1, dd2, o);
            tpolmat<fft_type> tg(dd2, dd3, o);
            tpolmat<fft_type> th(dd1, dd3, o);

            randomize(tf);
            randomize(tg);
            printf("%u [%u]", t, tr);
            double t_cubic;
            double t_strassen;
            // force cubic
            s.threshold(dd1,dd2,dd3) = UINT_MAX;
            do_and_printtime(t_cubic, compose_inner(th, tf, tg, o, s)," %.2g",1);
            // force strassen
            s.threshold(dd1,dd2,dd3) = 0;
            do_and_printtime(t_strassen, compose_inner(th, tf, tg, o, s)," %.2g",1);
            if (t_cubic <= t_strassen) {
                printf(" cubic [%.2f]\n", t_cubic / t_strassen);
                wt0 = wt;
            } else {
                printf(" strassen [%.2f]\n", t_cubic / t_strassen);
                wt1 = wt;
                earliest_good_strassen = tr;
            }
        }
        if (wt1 > max_wt)
            wt1 = max_wt;
        old_wt = wt1;
        if (wt1 == max_wt) wt1 = UINT_MAX;
        if (wt1 == min_wt) wt1 = 0;
        if (wt1 == 0) earliest_good_strassen = 0;
        if (wt1 == UINT_MAX) earliest_good_strassen = UINT_MAX;
        s.threshold(dd1,dd2,dd3) = earliest_good_strassen;
        printf("Threshold for %ux%u * %ux%u is %u bits [index %u]\n",
                dd1, dd2, dd2, dd3, s.threshold(dd1,dd2,dd3),
                s.threshold_index(dd1,dd2,dd3));
    }
    printf("\n");
}

    template<typename fft_type>
void tune_strassen(fft_type const& base, size_t maxlen)
{
    my_strassen_selector& s(foo<fft_type>::s);
    for(unsigned int i = 1 ; i <= (1 << BITS_IN_DIM_D) ; i++) {
        for(unsigned int j = 1 ; j <= (1 << BITS_IN_DIM_D) ; j++) {
            for(unsigned int k = 1 ; k <= (1 << BITS_IN_DIM_D) ; k++) {
                tune_strassen1(base,
                        i << MAX_I_FOR_TUNING,
                        j << MAX_I_FOR_TUNING,
                        k << MAX_I_FOR_TUNING, maxlen);
            }
        }
    }
    for(int i = 0 ; i < SELECTOR_NB_INDICES ; i++) {
        unsigned int t = s.strassen_threshold[i];
        if (t) printf("%u %u\n", i, t);
    }
}

    template<typename fft_type>
void plot_compose(const char * name,
        unsigned int n1, unsigned int n2, unsigned int n3, unsigned long wt)
{
    my_strassen_selector& s(foo<fft_type>::s);

    unsigned int nx1 = n1;
    unsigned int nx2 = n2;
    unsigned int nx3 = n3;
    unsigned long nmul = 1;
    for( ; !((nx1|nx2|nx3)&1UL) ; ) {
        if (!s(nx1,nx2,nx3,wt*ULONG_BITS))
            break;
        nx1>>=1;
        nx2>>=1;
        nx3>>=1;
        nmul *= 7;
    }
    nmul *= nx1 * nx2 * nx3;
    printf("For %ux%u * %ux%u, strassen begins at %ux%u * %ux%u ; roughly %lu mults\n", n1,n2,n2,n3,nx1,nx2,nx2,nx3,nmul);

    for( ; !((nx1|nx2|nx3)&1UL) ; nx1>>=1,nx2>>=1,nx3>>=1) ;

    /* Start by a measurement of the unit time */
    fft_type o(wt/2,wt/2);
    double single_compose;
    {
        typename fft_type::t * pf = o.alloc(1);
        typename fft_type::t * pg = o.alloc(1);
        typename fft_type::t * ph = o.alloc(1);
        do_and_noprinttime(single_compose, o.compose(ph, pf, pg));
        o.clear(pf, 1);
        o.clear(pg, 1);
        o.clear(ph, 1);
    }

    /* Now bench all possible matrix timings. */
    nmul = nx1 * nx2 * nx3;
    for( ; nx1 < n1 ; ) {
        tpolmat<fft_type> tf(nx1, nx2, o);
        tpolmat<fft_type> tg(nx2, nx3, o);
        tpolmat<fft_type> th(nx1, nx3, o);

        randomize(tf);
        randomize(tg);
        double v;
        do_and_noprinttime(v, compose_inner(th, tf, tg, o, s));
        printf("%ux%u * %ux%u : %.2e = %.2f%% of %lu * %.2e\n",
                nx1,nx2,nx2,nx3,v,
                100.0 * v/nmul/single_compose,nmul,single_compose);
        nx1 <<= 1;
        nx2 <<= 1;
        nx3 <<= 1;
        nmul *= s(nx1,nx2,nx3,wt*ULONG_BITS) ? 7 : 8;
    }
}

    template<typename fft_type>
void bench_one_polmm_projected_sub(fft_type& o, unsigned long d1, unsigned long d2, unsigned long d3, unsigned long n1, unsigned long n2, unsigned long * data)
{
    double dft1, dft2, compose, ift;
    fft_times(dft1, dft2, compose, ift, o, n1, n2, data);
    double total_compose_time = compose * nmults(o, d1, d2, d3);
    int d = depth(o, d1, d2,  d3);

    printf(" --> dft: %.2f ; mul(d%d): %.2f ; ift: %.2f\n",
            (dft1 * d1 * d2 + dft2 * d2 * d3),
            d,
            total_compose_time,
            ift * d1 * d3);
}

void bench_one_polmm_projected(unsigned long d1, unsigned long d2, unsigned long d3, unsigned long n1, unsigned long n2, unsigned long * data)
{
    printf("Timings %lux%lu (%lu-bit entries)"
            " times %lux%lu (%lu-bit entries) [projected timings]\n",
            d1, d2, n1, 
            d2, d3, n2);

    c128_fft oc(n1, n2); printf("c128:");
    bench_one_polmm_projected_sub(oc, d1, d2, d3, n1, n2, data);

    fake_fft of(n1, n2); printf("fake:");
    bench_one_polmm_projected_sub(of, d1, d2, d3, n1, n2, data);

    gf2x_tfft os(n1, n2, 81); printf("tfft(%u):", 81);
    bench_one_polmm_projected_sub(os, d1, d2, d3, n1, n2, data);
}

    template<typename fft_type>
void bench_one_polmm_complete_sub(fft_type& o, unsigned long d1, unsigned long d2, unsigned long d3, unsigned long n1, unsigned long n2)
{
    my_strassen_selector& s(foo<fft_type>::s);

    polmat f(d1, d2, n1);
    polmat g(d2, d2, n2);
    polmat h(d1, d3, n1 + n2 - 1);
    tpolmat<fft_type> tf(d1, d2, o);
    tpolmat<fft_type> tg(d2, d3, o);
    tpolmat<fft_type> th(d1, d3, o);

    randomize(f);
    randomize(g);
    randomize(tf);
    randomize(tg);

    logbook& l(foo<fft_type>::l);

    do_and_printtime(l.t1, transform(tf, f, o, n1), " dft1: %.2f", 1);
    do_and_printtime(l.t2, transform(tg, g, o, n2), " dft2: %.2f", 1);

    do_and_printtime(l.c, compose_inner(th, tf, tg, o, s), " compose: %.2f", 1);

#if 0
    printf(" [");
    /* Do composition by checking how deep strassen might pay off. */
    for(unsigned int k = 0 ; (((d1 | d2 | d3) >> k) & 1UL) == 0 ; k++) {
        double naive, strassen;
        // Force naive
        s.threshold(d1>>k,d2>>k,d3>>k) = UINT_MAX;
        do_and_printtime(naive, compose_inner(th, tf, tg, o, s));
        // Force strassen
        s.threshold(d1>>k,d2>>k,d3>>k) = 0;
        do_and_printtime(strassen, compose_inner(th, tf, tg, o, s));
        if (naive < strassen) {
            l.c = naive;
            break;
        }
        l.c = strassen;
    }
    printf("]");
#endif

    do_and_printtime(l.i, itransform(h, th, o, n1 + n2 - 1), " ift: %.2f", 1);
    printf("\n");
}

void bench_one_polmm_complete(unsigned long d1, unsigned long d2, unsigned long d3, unsigned long n1, unsigned long n2)
{
    printf("Timings %lux%lu (%lu-bit entries)"
            " times %lux%lu (%lu-bit entries) [complete timings]\n",
            d1, d2, n1, 
            d2, d3, n2);

    c128_fft oc(n1, n2); printf("c128:");
    bench_one_polmm_complete_sub(oc, d1, d2, d3, n1, n2);

    fake_fft of(n1, n2); printf("fake:");
    bench_one_polmm_complete_sub(of, d1, d2, d3, n1, n2);

    gf2x_tfft os(n1, n2, 81); printf("tfft(%u):", 81);
    bench_one_polmm_complete_sub(os, d1, d2, d3, n1, n2);
}


int main(int argc, char * argv[])
{
    /* This is a bench program. We want the output quick.  */
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    unsigned long m = 0, n = 0, N = 0;
    argc--,argv++;

    for( ; argc ; ) {
        if (strcmp(argv[0],"--nrep") == 0) {
            argc--,argv++;
            if (!argc) usage();
            nrep = atol(argv[0]);
            argc--,argv++;
            continue;
        }
        if (N == 0) {
            N = atol(argv[0]);
            argc--,argv++;
            continue;
        }
        if (m == 0) {
            m = atol(argv[0]);
            argc--,argv++;
            continue;
        }
        if (n == 0) {
            n = atol(argv[0]);
            argc--,argv++;
            continue;
        }
        usage();
    }
    // if (N == 0 || n == 0) usage();

    unsigned long b = m + n;

#if 0
    /* tune strassen for all sizes... Takes a looong time */
#define MAX_CONSIDERED_THRESHOLD        1000000
    {
        size_t d = (N * b / m / n);
        size_t n1 = d;
        size_t n2 = d * m / b;
        c128_fft oc(n1, n2); printf("c128");
        tune_strassen(oc, MAX_CONSIDERED_THRESHOLD);
        fake_fft of(n1, n2); printf("fake");
        tune_strassen(of, MAX_CONSIDERED_THRESHOLD);
        gf2x_tfft os(n1, n2, 81); printf("tfft(%u)", 81);
        tune_strassen(of, MAX_CONSIDERED_THRESHOLD);
    }
#endif
#if 1
    /* Tune for E * pi */
    {
        size_t d = (N * b / m / n);

        unsigned long dl = d-d/2;
        unsigned long pi_l_len = dl*m/b;   // always <= dl
        unsigned long chop = dl - pi_l_len;

        size_t n1 = d - chop;
        size_t n2 = pi_l_len;
        printf("Top-level multiplications E*pi: len%zu, %lux%lu * len%zu, %lux%lu\n",
                n1, m, b, n2, b, b);

        c128_fft oc(n1, n2); printf("c128");
        tune_strassen1(oc, m, b, b, n1 + n2 - 1);
        fake_fft of(n1, n2); printf("fake");
        tune_strassen1(of, m, b, b, n1 + n2 - 1);
        gf2x_tfft os(n1, n2, 81); printf("tfft(%u)", 81);
        tune_strassen1(os, m, b, b, n1 + n2 - 1);

        printf("Options for composition at top level E*pi\n");
        printf("c128: %u levels of strassen, %lu pol.muls\n",
                depth(oc,m,b,b), nmults(oc,m,b,b));
        printf("fake: %u levels of strassen, %lu pol.muls\n",
                depth(of,m,b,b), nmults(of,m,b,b));
        printf("gf2x: %u levels of strassen, %lu pol.muls\n",
                depth(os,m,b,b), nmults(os,m,b,b));
    }
    /* Tune for pi * pi */
    {
        size_t d = (N * b / m / n);

        unsigned long dl = d-d/2;
        unsigned long pi_l_len = dl*m/b;   // always <= dl
        size_t n1 = pi_l_len;
        size_t n2 = pi_l_len;

        c128_fft oc(n1, n2); printf("c128");
        tune_strassen1(oc, b, b, b, n1 + n2 - 1);
        fake_fft of(n1, n2); printf("fake");
        tune_strassen1(of, b, b, b, n1 + n2 - 1);
        gf2x_tfft os(n1, n2, 81); printf("tfft(%u)", 81);
        tune_strassen1(os, b, b, b, n1 + n2 - 1);
        printf("Options for composition at top level pi*pi\n");
        printf("c128: %u levels of strassen, %lu pol.muls\n",
                depth(oc,b,b,b), nmults(oc,b,b,b));
        printf("fake: %u levels of strassen, %lu pol.muls\n",
                depth(of,b,b,b), nmults(of,b,b,b));
        printf("gf2x: %u levels of strassen, %lu pol.muls\n",
                depth(os,b,b,b), nmults(os,b,b,b));
    }
    foo<c128_fft>::s.dump("C128");
    foo<fake_fft>::s.dump("FAKE");
    foo<gf2x_tfft>::s.dump("GF2X_TFFT");
#endif

    /* Seed the random state. Ugly at will. */
    extern char             __gmp_rands_initialized;
    extern gmp_randstate_t  __gmp_rands;
    __gmp_rands_initialized = 1;
    gmp_randinit_default (__gmp_rands);
    gmp_randseed_ui(__gmp_rands, time(NULL));

    /* Allocate an area which is very significantly larger than any of
     * the cache sizes */
    unsigned long * data ;
    data = (unsigned long *) malloc(DATA_POOL_SIZE * sizeof(unsigned long));
    for(unsigned int i = 0 ; i < DATA_POOL_SIZE ; i++) data[i] = random();

    printf("Some timings are in microseconds\n");
    printf("Matrix size %lu\n", N);
    // unsigned long c = m - n;

    for(unsigned long level = 0 ; ; level++) {
        unsigned long d = (N * b / m / n) >> level;
        if (d <= 64) {
            printf("[%lu] is below threshold (length %lu)\n", level, d);
            break;
        }

        if (DATA_POOL_SIZE * ULONG_BITS / 2 < d) {
            fprintf(stderr, "Please increa DATA_POOL_SIZE to at least %zu\n",
                    W(2*d));
            exit(1);
        }

        printf("[%lu] length of E is %lu (%lu MB)\n", level, d, d*m*b>>23);
        unsigned long dl = d-d/2;
        unsigned long pi_l_len = dl*m/b;   // always <= dl
        printf("[%lu] Top-level length of pi_left is %lu (%lu MB)\n", level, pi_l_len, pi_l_len*b*b>>23);
        unsigned long chop = dl - pi_l_len;
        printf("[%lu] Number of chopped of bits at top level is %lu\n", level, chop);
        printf("[%lu] Top-level degree of truncated E is %lu\n", level, d - chop);
        printf("[%lu] Degree of product E'*pi_left is %lu\n", level, d-chop + pi_l_len);
        printf("[%lu] Degree of product pi_left*pi_right is %lu\n", level, d*m/b);

        bench_one_polmm_projected(m, b, b, d-chop, pi_l_len, data);
        bench_one_polmm_projected(b, b, b, pi_l_len, pi_l_len, data);
        bench_one_polmm_complete(m, b, b, d-chop, pi_l_len);
        bench_one_polmm_complete(b, b, b, pi_l_len, pi_l_len);

#if 0
        /* Start by benching degree d polynomial multiplication */
        bench_polmul(d, d);

        /* Then time for cantor, separating the different steps */
        bench_c128(d, d);

        /* Now the matrix versions */
        bench_polmatmul<c128_fft>("c128_fft", n, n, d, d);
        bench_polmatmul<fake_fft>("fake_fft", n, n, d, d);
        bench_polmatmul<gf2x_tfft>("gf2x_tfft", n, n, d, d);

        /* Now do these again, but for unbalanced computations */

        bench_polmul(d, d/4);
        bench_c128(d, d/4);
        bench_polmatmul<c128_fft>("c128_fft", n/2, n, d, d/4);
        bench_polmatmul<fake_fft>("fake_fft", n/2, n, d, d/4);
        bench_polmatmul<gf2x_tfft>("gf2x_tfft", n/2, n, d, d/4);
#endif
    }

    free(data);
}


