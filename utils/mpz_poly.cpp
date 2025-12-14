/*
   These files (mpz_poly.*) implement arithmetics of polynomials whose
   coefficients are in multiprecision integers (using mpz_t from GNU MP).

   We use them in sqrt/algsqrt.c to represent rings of integers.

   Please see the file README_POLYNOMIALS for more details and
   comparisons to other poly formats.
*/

#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#ifdef HAVE_OPENMP
#include <type_traits>
#endif
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "cado_expression_parser.hpp"
#include "cxx_mpz.hpp"
#include "double_poly.h"
#include "gmp_aux.h"
#include "lll.h"
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_polymodF.h"
#include "mpz_poly_parallel.hpp"
#include "portability.h"
#include "rootfinder.h"
#include "runtime_numeric_cast.hpp"
/* and just because we expose a proxy to usp.c's root finding... */
#include "usp.h"
#include "utils_cxx.hpp"

// scan-headers: stop here
#ifdef MPZ_POLY_TIMINGS
#include "timing.h"
#include <ctime>
#else
#ifndef CADO_MPZ_POLY_H
#error "please include mpz_poly.h first"
#endif
#endif

/* A note on the parallel interface.
 *
 * The interface in mpz_poly.h has some of its functions duplicated in
 * mpz_poly_parallel.hpp. Those functions can be used either with or
 * without openmp.
 *
 * A function such as ::mpz_poly_blah is implemented here as a simple
 * trampoline to mpz_poly_notparallel_info::mpz_poly_blah.
 *
 * Then both mpz_poly_notparallel_info::mpz_poly_blah and
 * mpz_poly_parallel_info::mpz_poly_blah are created as instances of a
 * template member functions, so that the "parallel" criterion can
 * resolve statically. (static resolution is the main reason why I
 * didn't go for virtual functions)
 *
 * The code to implement mpz_poly_blah is therefore:
 *
 * void mpz_poly_blah(...)
 * {
 *      mpz_poly_notparallel_info().mpz_poly_blah(...);
 * }
 * template<typename inf>
 * void mpz_poly_parallel_interface<inf>::mpz_poly_blah(...)
 * {
 *      then std::is_same<inf, mpz_poly_notparallel_info>::value is a
 *      compile-time constant
 * }
 *
 */

// compute timings for Q^a mod(f,p): square, multiplication, reduction mod f
// PZ piece of code:
// in mpz_poly.h: add #define MPZ_POLY_TIMINGS
// beware: these timers are not thread-safe.
#ifdef MPZ_POLY_TIMINGS
static double timer[3] = {0.0, 0.0, 0.0};
static unsigned long calls[3] = {0, 0, 0};
#endif

#define TIMER_MUL 0
#define TIMER_SQR 1
#define TIMER_RED 2

#ifdef MPZ_POLY_TIMINGS
#define START_TIMER double t = seconds_thread()
#define RESTART_TIMER t = seconds_thread()
#define END_TIMER(x) add_timer(x, seconds_thread() - t)

/* flag=0: computing roots of f mod p
        1: checking irreducibility of f
        2: LLL reduction
        3: computing MurphyE */
static void add_timer(int flag, double t)
{
    timer[flag] += t;
    calls[flag] += 1;
}
#else
#define START_TIMER
#define RESTART_TIMER
#define END_TIMER(x)
#endif

#ifdef MPZ_POLY_TIMINGS
void print_timings_pow_mod_f_mod_p()
{
    printf(" (mul %.2fs/%lu, square %.2fs/%lu, red mod (f,p) %.2fs/%lu)",
           timer[TIMER_MUL], calls[TIMER_MUL], timer[TIMER_SQR],
           calls[TIMER_SQR], timer[TIMER_RED], calls[TIMER_RED]);
}
#endif

/* This sets a polynomial coefficient, even if that means reallocation.
 *
 * A note about a few design choices:
 *
 * - should we update the degree when we return a read-write accessor to
 *   a coefficient that is above the current degree?
 *   No, because we don't know yet what coefficient we're going to put
 *   anyway, and if we do that we may end up breaking consistency of the
 *   object.
 * - should we return zero on all coefficients above the current degree?
 *   It's not totally easy, but there are more use cases that hinge on
 *   the idea that we don't: this makes it possible to set all
 *   coefficients, possibly depending on each other, and _then_, outside
 *   the loop, set the degree with cleandeg.
 */
mpz_ptr mpz_poly_coeff(mpz_poly_ptr f, int i)
{
    mpz_poly_realloc(f, i + 1);
    return f->_coeff[i];
}

mpz_srcptr mpz_poly_coeff_const(mpz_poly_srcptr f, int i)
{
    static cxx_mpz zero = 0;
    if ((unsigned int)i < f->alloc)
        return f->_coeff[i];
    else
        return zero;
}

/* --------------------------------------------------------------------------
   Static functions
   -------------------------------------------------------------------------- */

/* Return f=g*h, where g has degree r, and h has degree s. */
static int mpz_poly_mul_basecase(mpz_t * f, mpz_t * g, int r, mpz_t * h, int s)
{
    int i, j;
    ASSERT(f != g && f != h);
    for (i = 0; i <= r + s; i++)
        mpz_set_ui(f[i], 0);
    for (i = 0; i <= r; ++i)
        for (j = 0; j <= s; ++j)
            mpz_addmul(f[i + j], g[i], h[j]);
    return r + s;
}

/* f <- g^2 where g has degree r */
static int mpz_poly_sqr_basecase(mpz_t * f, mpz_t * g, int r)
{
    int i, j;
    ASSERT(f != g);
    for (i = 0; i <= 2 * r; i++)
        mpz_set_ui(f[i], 0);
    for (i = 0; i <= r; ++i)
        for (j = 0; j < i; ++j)
            mpz_addmul(f[i + j], g[i], g[j]);
    for (i = 0; i < 2 * r; i++)
        mpz_mul_2exp(f[i], f[i], 1);
    for (i = 0; i < r; i++) {
        mpz_mul(f[2 * r], g[i], g[i]);
        mpz_add(f[2 * i], f[2 * i], f[2 * r]);
    }
    mpz_mul(f[2 * r], g[r], g[r]);
    return 2 * r;
}

/* To interpolate a polynomial of degree <= MAX_TC_DEGREE, we use the following
   MAX_TC_DEGREE evaluation points, plus the point at infinity.
   The first point should always be 0. */
static long const tc_points[MAX_TC_DEGREE] = {0,  1, -1, 2, -2, 3, -3, 4, -4, 5,
                                              -5, 6, -6, 7, -7, 8, -8, 9, -9};

#if GNUC_VERSION_ATLEAST(7, 0, 0)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#endif
/* Given f[0]...f[t] that contain respectively f(0), ..., f(t),
   put in f[0]...f[t] the coefficients of f. Assumes t <= MAX_TC_DEGREE.

   In the square root, with an algebraic polynomial of degree d,
   we have to multiply polynomials of degree d-1, thus we have t=2(d-1).
*/
static void mpz_poly_mul_tc_interpolate(mpz_t * f, int t)
{
    /* The alignment constraint below is bogus. We've seen this be
     * mandatory with g++ 12.2.1 20220924 on alpine linux, which seems to
     * very aggressively use simd rewriting on this part of the code, and
     * believe that M is always aligned. It is not, and forcing it to be
     * aligned has the effect of avoiding a segfault. But it definitely
     * looks like a compiler bug.
     */
    int64_t M[MAX_TC_DEGREE + 1][MAX_TC_DEGREE + 1] ATTR_ALIGNED(32), g, h;
    int i, j, k, l;
    /* G[] is the list of gcd's that appear in the forward Gauss loop, in the
       order they appear (they don't depend on t, since we start with the low
       triangular submatrix for t-1) */
    static int64_t const G[] = {1,
                                1,
                                2,
                                6,
                                24,
                                120,
                                720,
                                5040,
                                40320,
                                362880,
                                3628800,
                                39916800,
                                479001600,
                                6227020800,
                                87178291200,
                                1307674368000,
                                20922789888000,
                                355687428096000};

    ASSERT(t <= MAX_TC_DEGREE); /* Ensures that all M[i][j] fit in uint64_t,
                                   and similarly for all intermediate
                                   computations on M[i][j]. This avoids the
                                   use of mpz_t to store the M[i][j]. */
    ASSERT_ALWAYS(t >= 0);

    /* initialize M[i][j] = tc_points[i]^j */
    for (i = 0; i <= t; i++)
        for (j = 0; j <= t; j++) {
            if (i < t)
                M[i][j] = (j == 0) ? 1 : tc_points[i] * M[i][j - 1];
            else /* evaluation point is infinity */
                M[i][j] = (j == t);
        }

    /* Forward Gauss: zero the under-diagonal coefficients while going down.
       Since the last point of evaluation is infinity, the last row is already
       reduced, thus we go up to i=t-1. */
    for (i = 1; i < t; i++) {
        for (j = 0, l = 0; j < i; j++) {
            g = G[l++]; /* same as gcd_int64 (M[i][j], M[j][j]) */
            h = M[i][j] / g;
            g = M[j][j] / g;
            /* f[i] <- g*f[i] - h*f[j] */
            ASSERT(g > 0);
            if (g != 1)
                mpz_mul_uint64(f[i], f[i], g);
            mpz_submul_int64(f[i], f[j], h);
            for (k = j; k <= t; k++)
                M[i][k] = g * M[i][k] - h * M[j][k];
        }
    }

    /* now zero upper-diagonal coefficients while going up */
    for (i = t - 1; i >= 0; i--) {
        for (j = i + 1; j <= t; j++)
            /* f[i] = f[i] - M[i][j] * f[j] */
            mpz_submul_int64(f[i], f[j], M[i][j]);
        ASSERT(mpz_divisible_uint64_p(f[i], M[i][i]));
        mpz_divexact_uint64(f[i], f[i], M[i][i]);
    }
}
#if GNUC_VERSION_ATLEAST(7, 0, 0)
#pragma GCC diagnostic pop
#endif

/* v <- g(i) */
static void mpz_poly_mul_eval_si(mpz_t v, mpz_t * g, int r, long i)
{
    mpz_set(v, g[r]);
    for (int j = r - 1; j >= 0; j--) {
        mpz_mul_si(v, v, i);
        mpz_add(v, v, g[j]);
    }
}

/* Generic Toom-Cook implementation: stores in f[0..r+s] the coefficients
   of g*h, where g has degree r and h has degree s, and their coefficients
   are in g[0..r] and h[0..s].
   Assumes f differs from g and h, and f[0..r+s] are already allocated.
   Returns the degree of f.
*/
template <typename inf>
static int mpz_poly_mul_tc(inf & inf_arg, mpz_t * f, mpz_t * g, int r,
                           mpz_t * h, int s)
{
    int t;

    if ((r == -1) || (s == -1)) /* g or h is 0 */
        return -1;

    if (r < s)
        return mpz_poly_mul_tc(inf_arg, f, h, s, g, r);

    /* now r >= s */

    if (s == 0) {
        for (int i = 0; i <= r; i++)
            mpz_mul(f[i], g[i], h[0]);
        return r;
    }

    t = r + s; /* degree of f, which has t+1 coefficients */

    if (t == 2) /* necessary r = s = 1, use vanilla Karatsuba with 3 MUL,
                   this is always optimal */
    {
        /* TODO: does it make sense to do the two last MULs below in
         * parallel if entries are large enough ? */
        mpz_add(f[0], g[0], g[1]);
        mpz_add(f[2], h[0], h[1]);

        mpz_mul(f[1], f[0], f[2]);

        mpz_mul(f[0], g[0], h[0]);
        mpz_mul(f[2], g[1], h[1]);

        mpz_sub(f[1], f[1], f[0]);

        mpz_sub(f[1], f[1], f[2]);
        return 2;
    }

    if (t == 3) /* r = 2, s = 1, the following code in 4 MUL is optimal */
    {
        /* TODO: do part of the stuff (if relevant) in parallel if entries
         * are large enough ? */
        mpz_add(f[1], g[2], g[0]);
        mpz_add(f[2], f[1], g[1]); /* g(1) */
        mpz_sub(f[1], f[1], g[1]); /* g(-1) */
        mpz_add(f[0], h[0], h[1]); /* h(1) */
        mpz_mul(f[2], f[2], f[0]); /* f(1) = g(1)*h(1) */
        mpz_sub(f[0], h[0], h[1]); /* h(-1) */
        mpz_mul(f[0], f[1], f[0]); /* f(-1) = g(-1)*h(-1) */
        mpz_sub(f[1], f[2], f[0]); /* f(1) - f(-1) = 2*g2*h1+2*g1*h0+2*g0*h1 */
        mpz_add(f[2], f[2], f[0]); /* f(1) + f(-1) = 2*g2*h0+2*g1*h1+2*g0*h0 */
        mpz_div_2exp(f[1], f[1], 1); /* g2*h1 + g1*h0 + g0*h1 */
        mpz_div_2exp(f[2], f[2], 1); /* g2*h0 + g1*h1 + g0*h0 */
        mpz_mul(f[0], g[0], h[0]);
        mpz_sub(f[2], f[2], f[0]); /* g2*h0 + g1*h1 */
        mpz_mul(f[3], g[2], h[1]);
        mpz_sub(f[1], f[1], f[3]); /* g1*h0 + g0*h1 */
        return 3;
    }

    if (r == 2) /* Necessarily s = 2 since r >= s and the cases s = 0 and
                   (r,s) = (2,1) have already been treated. */
    {
        /* TODO: do part of the stuff (if relevant) in parallel if entries
         * are large enough ? */
        /* we use here the code from Appendix Z from "Towards Optimal Toom-Cook
           Multiplication for Univariate and Multivariate Polynomials in
           Characteristic 2 and 0" from Marco Bodrato, WAIFI 2007, LNCS 4547,
           with the change: W3 -> f[1], W2 -> f[3], W1 -> f[2]. */
        mpz_add(f[0], g[2], g[0]);
        mpz_add(f[4], h[2], h[0]);
        mpz_sub(f[3], f[0], g[1]);
        mpz_sub(f[2], f[4], h[1]);
        mpz_add(f[0], f[0], g[1]);
        mpz_add(f[4], f[4], h[1]);
        mpz_mul(f[1], f[3], f[2]);
        mpz_mul(f[2], f[0], f[4]);
        mpz_add(f[0], f[0], g[2]);
        mpz_mul_2exp(f[0], f[0], 1);
        mpz_sub(f[0], f[0], g[0]);
        mpz_add(f[4], f[4], h[2]);
        mpz_mul_2exp(f[4], f[4], 1);
        mpz_sub(f[4], f[4], h[0]);
        mpz_mul(f[3], f[0], f[4]);
        mpz_mul(f[0], g[0], h[0]);
        mpz_mul(f[4], g[2], h[2]);
        /* interpolation */
        mpz_sub(f[3], f[3], f[1]);
        ASSERT(mpz_divisible_ui_p(f[3], 3));
        mpz_divexact_ui(f[3], f[3], 3);
        mpz_sub(f[1], f[2], f[1]);
        ASSERT(mpz_divisible_ui_p(f[1], 2));
        mpz_div_2exp(f[1], f[1], 1); /* exact */
        mpz_sub(f[2], f[2], f[0]);
        mpz_sub(f[3], f[3], f[2]);
        ASSERT(mpz_divisible_ui_p(f[3], 2));
        mpz_div_2exp(f[3], f[3], 1); /* exact */
        mpz_submul_ui(f[3], f[4], 2);
        mpz_sub(f[2], f[2], f[1]);
        mpz_sub(f[2], f[2], f[4]);
        mpz_sub(f[1], f[1], f[3]);
        return 4;
    }

    if (t > MAX_TC_DEGREE) {
        /* naive product */
        /* currently we have to resort to this for larger degree, because
         * the generic Toom implementation is bounded in degree.
         */
        return mpz_poly_mul_basecase(f, g, r, h, s);
    }

    /* now t <= MAX_TC_DEGREE */

    /* store g(tc_points[i])*h(tc_points[i]) in f[i] for 0 <= i <= t */

    ASSERT(tc_points[0] == 0);

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= t; i++) {
        if (i == 0) /* evaluate at 0 */
            mpz_mul(f[0], g[0], h[0]);
        else if (i == t) /* evaluate at infinity */
            mpz_mul(f[t], g[r], h[s]);
        else {
            ASSERT_FOR_STATIC_ANALYZER(i < t);
            mpz_t tmp;
            mpz_init(tmp);
            /* f[i] <- g(i) */
            mpz_poly_mul_eval_si(f[i], g, r, tc_points[i]);
            /* f[t] <- h(i) */
            mpz_poly_mul_eval_si(tmp, h, s, tc_points[i]);
            /* f[i] <- g(i)*h(i) */
            mpz_mul(f[i], f[i], tmp);
            mpz_clear(tmp);
        }
    }

    mpz_poly_mul_tc_interpolate(f, t);

    return t;
}

/* f2*x^2+f1*x+f0 = (g1*x+g0)^2 using Karatsuba and 3 SQR */
static int mpz_poly_sqr_tc2(mpz_t * f, mpz_t * g)
{
    mpz_mul(f[0], g[0], g[0]);
    mpz_mul(f[2], g[1], g[1]);
    mpz_add(f[1], g[0], g[1]);
    mpz_mul(f[1], f[1], f[1]);
    mpz_sub(f[1], f[1], f[0]);
    mpz_sub(f[1], f[1], f[2]);
    return 2;
}

#if 1
/* Algorithm SQR3 from Asymmetric Squaring Formulae by Chung and Hasan,
   with 4 SQR and 1 MUL. */
static int mpz_poly_sqr_tc3(mpz_t * c, mpz_t * a)
{
    mpz_mul(c[0], a[0], a[0]); /* c0 = S0 = a0^2 */
    mpz_add(c[1], a[2], a[0]); /* c1 = a2 + a0 */
    mpz_sub(c[2], c[1], a[1]); /* c2 = a2 - a1 + a0 */
    mpz_add(c[1], c[1], a[1]); /* c1 = a2 + a1 + a0 */
    mpz_mul(c[1], c[1], c[1]); /* S1 = (a2 + a1 + a0)^2 */
    mpz_mul(c[2], c[2], c[2]); /* S2 = (a2 - a1 + a0)^2 */
    mpz_mul(c[3], a[2], a[1]);
    mpz_mul_2exp(c[3], c[3], 1); /* S3 = 2*a1*a2 */
    mpz_add(c[2], c[1], c[2]);
    mpz_tdiv_q_2exp(c[2], c[2], 1); /* T1 = (S1 + S2) / 2 */
    mpz_sub(c[1], c[1], c[2]);      /* S1 - T1 */
    mpz_sub(c[1], c[1], c[3]);      /* c1 = S1 - T1 - S3 */
    mpz_mul(c[4], a[2], a[2]);      /* c4 = S4 = a2^2 */
    mpz_sub(c[2], c[2], c[4]);      /* T1 - S4 */
    mpz_sub(c[2], c[2], c[0]);      /* c2 = T1 - S4 - S0 */
    return 4;
}
#else
/* (g2*x^2+g1*x+g0)^2 = g2^2*x^4 + (2*g2*g1)*x^3 +
   (g1^2+2*g2*g0)*x^2 + (2*g1*g0)*x + g0^2 in 3 SQR + 2 MUL.
   Experimentally this is faster for less than 4096 bits than the
   algorithm in 5 SQR from mpz_poly_mul_tc_interpolate. */
static int mpz_poly_sqr_tc3(mpz_t * f, mpz_t * g)
{
    mpz_mul(f[4], g[2], g[2]); /* g2^2 */
    mpz_mul(f[3], g[2], g[1]);
    mpz_mul_2exp(f[3], f[3], 1); /* 2*g2*g1 */
    mpz_mul(f[1], g[1], g[0]);
    mpz_mul_2exp(f[1], f[1], 1); /* 2*g1*g0 */
    mpz_mul(f[0], g[0], g[0]);   /* g0^2 */
    mpz_add(f[2], g[2], g[1]);
    mpz_add(f[2], f[2], g[0]);
    mpz_mul(f[2], f[2], f[2]); /* (g2+g1+g0)^2 */
    mpz_sub(f[2], f[2], f[4]);
    mpz_sub(f[2], f[2], f[0]); /* g1^2 + 2*g2*g0 + 2*g2*g1 + 2*g1*g0 */
    mpz_sub(f[2], f[2], f[3]); /* g1^2 + 2*g2*g0 + 2*g1*g0 */
    mpz_sub(f[2], f[2], f[1]); /* g1^2 + 2*g2*g0 */
    return 4;
}
#endif

#if 1
/* 2 levels of Karatsuba: 9 SQR */
static int mpz_poly_sqr_tc4(mpz_t * f, mpz_t * g)
{
    mpz_t t;

    mpz_init(t);
    /* f4*x^2+f3*x+f2 <- [(g3+g1)*x + (g2+g0)]^2 */
    mpz_add(f[0], g[0], g[2]);
    mpz_add(f[1], g[1], g[3]);
    mpz_poly_sqr_tc2(f + 2, f);
    /* save f2 into t */
    mpz_swap(t, f[2]);
    /* f2*x^2+f1*x+f0 <- (g1*x+g0)^2 */
    mpz_poly_sqr_tc2(f, g);
    /* subtract f2*x^2+f1*x+f0 from f4*x^2+f3*x+t */
    mpz_sub(f[4], f[4], f[2]);
    mpz_sub(f[3], f[3], f[1]);
    mpz_sub(t, t, f[0]);
    /* re-add t into f2 */
    mpz_add(f[2], f[2], t);
    /* save f4 into t */
    mpz_swap(t, f[4]);
    /* f6*x^2+f5*x+f4 <- (g3*x+g2)^2 */
    mpz_poly_sqr_tc2(f + 4, g + 2);
    /* subtract f6*x^2+f5*x+f4 from t*x^2+f3*x+f2 */
    mpz_sub(t, t, f[6]);
    mpz_sub(f[3], f[3], f[5]);
    mpz_sub(f[2], f[2], f[4]);
    /* re-add t into f4 */
    mpz_add(f[4], f[4], t);
    mpz_clear(t);
    return 6;
}
#else
static int mpz_poly_sqr_tc4(mpz_t * f, mpz_t * g)
{
    /* (g3*x^3+g2*x^2+g1*x+g0)^2 = g3^2*x^6 + (2*g3*g2)*x^5
       + (g2^2+2*g3*g1)*x^4 + (2*g3*g0+2*g2*g1)*x^3
       + (g1^2+2*g2*g0)*x^2 + (2*g1*g0)*x + g0^2 in 4 SQR + 6 MUL.
       Experimentally this is faster for less than 4096 bits than the
       algorithm in 7 SQR from mpz_poly_mul_tc_interpolate. */
    mpz_mul(f[6], g[3], g[3]); /* g3^2 */
    mpz_mul(f[5], g[3], g[2]);
    mpz_mul_2exp(f[5], f[5], 1); /* 2*g3*g2 */
    mpz_mul(f[1], g[1], g[0]);
    mpz_mul_2exp(f[1], f[1], 1); /* 2*g1*g0 */
    mpz_mul(f[2], g[1], g[1]);   /* g1^2 */
    mpz_mul(f[0], g[2], g[0]);
    mpz_addmul_ui(f[2], f[0], 2); /* g1^2+2*g2*g0 */
    mpz_mul(f[4], g[2], g[2]);    /* g2^2 */
    mpz_mul(f[0], g[3], g[1]);
    mpz_addmul_ui(f[4], f[0], 2); /* g2^2+2*g3*g1 */
    mpz_mul(f[3], g[3], g[0]);
    mpz_mul(f[0], g[2], g[1]);
    mpz_add(f[3], f[3], f[0]);
    mpz_mul_2exp(f[3], f[3], 1); /* 2*g3*g0+2*g2*g1 */
    mpz_mul(f[0], g[0], g[0]);   /* g0^2 */
    return 6;
}
#endif

/* Same as mpz_poly_mul_tc for the squaring: store in f[0..2r] the coefficients
   of g^2, where g has degree r, and its coefficients are in g[0..r].
   Assumes f differs from g, and f[0..2r] are already allocated.
   Returns the degree of f.
*/
template <typename inf>
static int mpz_poly_sqr_tc(inf & /* unused */, mpz_t * f, mpz_t * g, int r)
{
    int t;
    size_t nbits;

    if (r == -1) /* g is 0 */
        return -1;

    if (r == 0) {
        mpz_mul(f[0], g[0], g[0]);
        return 0;
    }

    if (r == 1) /* 3 SQR: always optimal */
        return mpz_poly_sqr_tc2(f, g);

    /* we assume the number of bits of f[0] is representative of the size of
       the coefficients of f */
    nbits = mpz_sizeinbase(f[0], 2);

    if (r == 2 && nbits < 4096)
        return mpz_poly_sqr_tc3(f, g);

    if (r == 3 && nbits < 4096)
        return mpz_poly_sqr_tc4(f, g);

    t = 2 * r; /* product has degree t thus t+1 coefficients */
    if (t > MAX_TC_DEGREE) {
        /* naive product */
        /* currently we have to resort to this for larger degree, because
         * the generic toom implementation is bounded in degree.
         */
        return mpz_poly_sqr_basecase(f, g, r);
    }

    ASSERT(t <= MAX_TC_DEGREE);

    /* store g(tc_points[i])^2 in f[i] for 0 <= i <= t */

    ASSERT(tc_points[0] == 0);

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= t; i++) {
        if (i == 0) /* evaluate at 0 */
            mpz_mul(f[0], g[0], g[0]);
        else if (i == t) /* evaluate at infinity */
            mpz_mul(f[t], g[r], g[r]);
        else {
            ASSERT_FOR_STATIC_ANALYZER(i < t);
            /* f[i] <- g(i) */
            mpz_poly_mul_eval_si(f[i], g, r, tc_points[i]);
            mpz_mul(f[i], f[i], f[i]);
        }
    }

    mpz_poly_mul_tc_interpolate(f, t);

    return t;
}

/* --------------------------------------------------------------------------
   Public functions
   -------------------------------------------------------------------------- */

/* Management of the structure, set and print coefficients. */

/* Allocate a polynomial that can contain 'd+1' coefficients and set to zero.
   We allow d < 0, which is equivalent to d = -1.

   XXX mpz_poly_init takes a degree ; mpz_poly_realloc takes a number of coeffs
 */
void mpz_poly_init(mpz_poly_ptr f, int d)
{
    f->deg = -1;
    if (d < 0) {
        f->alloc = 0;
        f->_coeff = (mpz_t *)nullptr;
    } else {
        int i;
        f->alloc = d + 1;
        f->_coeff = (mpz_t *)malloc((d + 1) * sizeof(mpz_t));
        FATAL_ERROR_CHECK(f->_coeff == nullptr, "not enough memory");
        for (i = 0; i <= d; ++i)
            mpz_init_set_ui(f->_coeff[i], 0);
    }
}

/* realloc f to (at least) nc coefficients
   XXX mpz_poly_init takes a degree ; mpz_poly_realloc takes a number of coeffs
   */
void mpz_poly_realloc(mpz_poly_ptr f, unsigned int nc)
{
    ASSERT_ALWAYS(nc <= (unsigned int)INT_MAX);
    if (f->alloc < nc) {
        checked_realloc(f->_coeff, nc);
        for (unsigned int i = f->alloc; i < nc; i++)
            mpz_init_set_ui(f->_coeff[i], 0);
        f->alloc = nc;
    }
}

/* Copy f to g, where g must be initialized (but there is no need it has
   enough allocated coefficients, since mpz_poly_setcoeff reallocates if
   needed). */
void mpz_poly_set(mpz_poly_ptr g, mpz_poly_srcptr f)
{
    int i;

    if (g == f)
        return; /* nothing to do */

    g->deg = f->deg;
    for (i = f->deg; i >= 0; --i)
        mpz_poly_setcoeff(g, i, f->_coeff[i]);
}

void mpz_poly_set_double_poly(mpz_poly_ptr g, double_poly_srcptr f)
{
    mpz_poly_realloc(g, f->deg + 1);
    g->deg = f->deg;
    for (int i = f->deg; i >= 0; --i)
        mpz_set_d(g->_coeff[i], f->coeff[i]);
}

/* Init polynomial rel and set it to a - b*x */
void mpz_poly_init_set_ab(mpz_poly_ptr rel, int64_t a, uint64_t b)
{
    mpz_poly_init(rel, 1);
    mpz_poly_set_ab(rel, a, b);
}

void mpz_poly_set_ab(mpz_poly_ptr rel, int64_t a, uint64_t b)
{
    mpz_poly_set_zero(rel);
    mpz_poly_setcoeff_int64(rel, 0, a);
    mpz_poly_setcoeff_int64(rel, 1, -runtime_numeric_cast<int64_t>(b));
}

/* rel <- a-b*x */
void mpz_poly_init_set_mpz_ab(mpz_poly_ptr rel, mpz_srcptr a, mpz_srcptr b)
{
    mpz_poly_init(rel, 1);
    mpz_poly_set_mpz_ab(rel, a, b);
}

void mpz_poly_set_mpz_ab(mpz_poly_ptr rel, mpz_srcptr a, mpz_srcptr b)
{
    mpz_poly_set_zero(rel);
    mpz_poly_setcoeff(rel, 0, a);
    mpz_t mb;
    mpz_init(mb);
    mpz_neg(mb, b);
    mpz_poly_setcoeff(rel, 1, mb);
    mpz_clear(mb);
}

/* swap f and g */
void mpz_poly_swap(mpz_poly_ptr f, mpz_poly_ptr g)
{
    std::swap(f->alloc, g->alloc);
    std::swap(f->deg, g->deg);
    std::swap(f->_coeff, g->_coeff);
}

/* Free polynomial f in mpz_poly. */
void mpz_poly_clear(mpz_poly_ptr f)
{
    for (unsigned i = 0; i < f->alloc; ++i)
        mpz_clear(f->_coeff[i]);
    if (f->_coeff != nullptr)
        free(f->_coeff);
    f->_coeff = nullptr; /* to avoid a double-free */
    memset(f, 0, sizeof(mpz_poly));
    f->deg = -1;
    f->alloc = 0; /* to avoid a double-free */
}

namespace
{
struct urandomm {
    using argtype = mpz_srcptr;
    argtype k;
    explicit urandomm(argtype k)
        : k(k)
    {
        ASSERT_ALWAYS(mpz_sgn(k) > 0);
    }
    void fetch_half(cxx_mpz & h) const { mpz_div_2exp(h, k, 1); }
    void operator()(mpz_ptr x, gmp_randstate_t state) const
    {
        mpz_urandomm(x, state, k);
    }
};
struct urandomm_ui {
    using argtype = unsigned long;
    argtype k;
    explicit urandomm_ui(argtype k)
        : k(k)
    {
        ASSERT_ALWAYS(k > 0);
    }
    void fetch_half(cxx_mpz & h) const { h = k / 2; }
    void operator()(mpz_ptr x, gmp_randstate_t state) const
    {
        mpz_set_ui(x, gmp_urandomm_ui(state, k));
    }
};

struct urandomb {
    using argtype = int;
    argtype k;
    explicit urandomb(argtype k)
        : k(k)
    {
        ASSERT_ALWAYS(k > 0);
    }
    void fetch_half(cxx_mpz & h) const
    {
        h = 1;
        mpz_mul_2exp(h, h, k - 1);
    }
    void operator()(mpz_ptr x, gmp_randstate_t state) const
    {
        mpz_urandomb(x, state, k);
    }
};
struct rrandomb {
    /* note that mpz_rrandomb(..., ..., k) only gives k bits of entropy
     * (return value between 2^(k-1) and 2^k-1).
     */
    using argtype = int;
    argtype k;
    explicit rrandomb(argtype k)
        : k(k)
    {
        ASSERT_ALWAYS(k > 0);
    }
    void fetch_half(cxx_mpz & h) const
    {
        h = 1;
        mpz_mul_2exp(h, h, k - 1);
    }
    void operator()(mpz_ptr x, gmp_randstate_t state) const
    {
        /* get 2^k <= x < 2^(k+1) */
        mpz_rrandomb(x, state, k + 1);
        /* subtract 2^k, i.e.  clear the k-th bit */
        mpz_clrbit(x, k);
        /* we now have 0 <= x < 2^k */
    }
};
/* Put random coefficients of k bits in a polynomial
 * Ensure the coefficient of degree d is not zero.
 */
template <typename R>
void mpz_poly_set_random_internal(mpz_poly_ptr f, int d,
                                  gmp_randstate_ptr state, R const & r,
                                  int flags)
{
    /* Note: this code used to implicitly assume exact as being set to
     * true */
    cxx_mpz u;

    bool const exact = !(flags & mpz_poly_random_flags::MPZ_POLY_DEGREE_UPPER_BOUND);
    bool const is_signed = flags & mpz_poly_random_flags::MPZ_POLY_SIGNED_COEFFICIENTS;
    bool const monic = flags & mpz_poly_random_flags::MPZ_POLY_MONIC;

    if (is_signed)
        r.fetch_half(u);

    for (int i = 0; i <= d; i++) {
        mpz_ptr fi = mpz_poly_coeff(f, i);
        do {
            r(fi, state);
            if (is_signed)
                mpz_sub(fi, fi, u);
        } while (i == d && exact && mpz_cmp_ui(fi, 0) == 0);
    }
    mpz_poly_cleandeg(f, d);
    if (monic) {
        if (f->deg == -1) {
            mpz_poly_set_xi(f, 0);
        } else {
            mpz_set_ui(f->_coeff[f->deg], 1);
        }
    }
}
} // namespace

void mpz_poly_set_randomb(mpz_poly_ptr f, int d,
                                  gmp_randstate_ptr state, int k, int flags)
{
    if (flags & mpz_poly_random_flags::MPZ_POLY_RRANDOM)
        mpz_poly_set_random_internal(f, d, state, rrandomb(k), flags);
    else
        mpz_poly_set_random_internal(f, d, state, urandomb(k), flags);
}

void mpz_poly_set_randomm(mpz_poly_ptr f, int d,
                          gmp_randstate_ptr state, mpz_srcptr N, int flags)
{
    ASSERT_ALWAYS(!(flags & mpz_poly_random_flags::MPZ_POLY_RRANDOM));
    mpz_poly_set_random_internal(f, d, state, urandomm(N), flags);
}

void mpz_poly_set_randomm_ui(mpz_poly_ptr f, int d, gmp_randstate_ptr state,
                              unsigned long m, int flags)
{
    ASSERT_ALWAYS(!(flags & mpz_poly_random_flags::MPZ_POLY_RRANDOM));
    mpz_poly_set_random_internal(f, d, state, urandomm_ui(m), flags);
}

/* removed mpz_poly_set_deg, as for all purposes there is no reason to
 * not use the more robust mpz_poly_cleandeg */

/* Find polynomial degree. */
void mpz_poly_cleandeg(mpz_poly_ptr f, int deg)
{
    ASSERT(deg >= -1);
    if ((unsigned int)(deg + 1) >= f->alloc)
        deg = (int)f->alloc - 1;
    while ((deg >= 0) && mpz_sgn(f->_coeff[deg]) == 0)
        deg--;
    f->deg = deg;
}

/* Sets f to the polynomial of degree d, of coefficients
   given by coeffs. */
void mpz_poly_setcoeffs(mpz_poly_ptr f, mpz_t * coeffs, int d)
{
    for (int i = d; i >= 0; --i)
        mpz_poly_setcoeff(f, i, coeffs[i]);
    mpz_poly_cleandeg(f, d);
}

/* Set signed long int coefficient for the i-th term. */
void mpz_poly_setcoeffs_si(mpz_poly_ptr f, long int const * h, int d)
{
    for (int i = d; i >= 0; --i)
        mpz_poly_setcoeff_si(f, i, h[i]);
    mpz_poly_cleandeg(f, d);
}

/* Set unsigned long int coefficient for the i-th term. */
void mpz_poly_setcoeffs_ui(mpz_poly_ptr f, unsigned long int const * h, int d)
{
    for (int i = d; i >= 0; --i)
        mpz_poly_setcoeff_ui(f, i, h[i]);
    mpz_poly_cleandeg(f, d);
}

/* Set a zero polynomial. */
void mpz_poly_set_zero(mpz_poly_ptr f)
{
    f->deg = -1;
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff(mpz_poly_ptr f, int i, mpz_srcptr z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* Set signed int coefficient for the i-th term. */
void mpz_poly_setcoeff_si(mpz_poly_ptr f, int i, long z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set_si(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* Set unsigned int coefficient for the i-th term. */
void mpz_poly_setcoeff_ui(mpz_poly_ptr f, int i, unsigned long z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set_ui(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* Set int64 coefficient for the i-th term. */
void mpz_poly_setcoeff_int64(mpz_poly_ptr f, int i, int64_t z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set_int64(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* Set uint64 coefficient for the i-th term. */
void mpz_poly_setcoeff_uint64(mpz_poly_ptr f, int i, uint64_t z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set_uint64(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* f[i] <- z */
void mpz_poly_setcoeff_double(mpz_poly_ptr f, int i, double z)
{
    mpz_poly_realloc(f, i + 1);
    /* ensure consistency of inserted coefficients */
    for (int j = f->deg + 1; j < i; j++)
        mpz_set_ui(f->_coeff[j], 0);
    mpz_set_d(f->_coeff[i], z);
    if (i >= f->deg)
        mpz_poly_cleandeg(f, i);
}

/* x^i is often useful */
void mpz_poly_set_xi(mpz_poly_ptr f, int i)
{
    mpz_poly_realloc(f, i + 1);
    for (int j = 0; j <= i; j++) {
        mpz_set_ui(f->_coeff[j], j == i);
    }
    f->deg = i;
}

void mpz_poly_set_mpz(mpz_poly_ptr f, mpz_srcptr z)
{
    mpz_poly_realloc(f, 1);
    mpz_set(f->_coeff[0], z);
    mpz_poly_cleandeg(f, 0);
}

void mpz_poly_set_ui(mpz_poly_ptr f, unsigned long z)
{
    mpz_poly_realloc(f, 1);
    mpz_set_ui(f->_coeff[0], z);
    mpz_poly_cleandeg(f, 0);
}

void mpz_poly_set_si(mpz_poly_ptr f, long z)
{
    mpz_poly_realloc(f, 1);
    mpz_set_si(f->_coeff[0], z);
    mpz_poly_cleandeg(f, 0);
}

/* g <- quo (f, x^i) */
void mpz_poly_div_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i)
{
    if (f->deg < i) {
        mpz_poly_set_zero(g);
        return;
    }
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }
    if (g == f) {
        /* rotate the coefficients, don't do any freeing of the mpz's: we
         * assume that we might have a use for them later.
         * (new mpz_t[] does obviously not call mpz_init)
         */
        auto const temp = std::unique_ptr<mpz_t[]>(new mpz_t[i]);
        memcpy(temp.get(), g->_coeff, i * sizeof(mpz_t));
        memmove(g->_coeff, g->_coeff + i, (g->deg + 1 - i) * sizeof(mpz_t));
        memcpy(g->_coeff + (g->deg + 1 - i), temp.get(), i * sizeof(mpz_t));
        g->deg -= i;
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg - i);
    for (int j = i; j <= f->deg; j++) {
        mpz_set(g->_coeff[j - i], f->_coeff[j]);
    }
    g->deg = f->deg - i;
}

/* g <- f * x^i */
void mpz_poly_mul_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i)
{
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg + i);

    if (g == f) {
        for (int j = g->deg; j >= 0; j--) {
            mpz_swap(g->_coeff[j + i], g->_coeff[j]);
        }
    } else {
        for (int j = g->deg; j >= 0; j--) {
            mpz_set(g->_coeff[j + i], g->_coeff[j]);
        }
    }
    for (int j = 0; j < i; j++) {
        mpz_set_ui(g->_coeff[j], 0);
    }
    g->deg = f->deg + i;
}

/* g <- f * (x + a) */
void mpz_poly_mul_xplusa(mpz_poly_ptr g, mpz_poly_srcptr f, mpz_srcptr a)
{
    mpz_t aux;
    mpz_init_set_ui(aux, 0);
    mpz_poly_realloc(g, f->deg + 2);
    for (int i = 0; i <= f->deg; i++) {
        /* aux is is coeff of degree [i-1]. We want
         * (coeff_i, aux) <- (coeff_{i-1} + a * coeff_i, coeff_i)
         *                   (aux + a * coeff_i, coeff_i)
         */
        if (a)
            mpz_addmul(aux, f->_coeff[i], a);
        mpz_swap(g->_coeff[i], aux);
        if (f != g) {
            mpz_set(aux, f->_coeff[i]);
        }
    }
    /* last coefficient */
    mpz_swap(g->_coeff[f->deg + 1], aux);
    /* This is just as valid as it was for f */
    g->deg = f->deg + 1;
    mpz_clear(aux);
}

/* return the valuation of f */
int mpz_poly_valuation(mpz_poly_srcptr f)
{
    int n = 0;
    ASSERT(f->deg >= 0);
    for (; n < f->deg && mpz_cmp_ui(f->_coeff[n], 0) == 0; n++)
        ;
    return n;
}

/* return true if f is zero or a multiple of x^deg(f) */
int mpz_poly_is_monomial_multiple(mpz_poly_srcptr f)
{
    return mpz_poly_degree(f) == -1 ||
           mpz_poly_valuation(f) == mpz_poly_degree(f);
}

/* return true if f == x^deg(f) */
int mpz_poly_is_monomial(mpz_poly_srcptr f)
{
    return mpz_poly_degree(f) >= 0 &&
           mpz_poly_valuation(f) == mpz_poly_degree(f) &&
           mpz_cmp_ui(mpz_poly_lc(f), 1) == 0;
}

int mpz_poly_asprintf(char ** res, mpz_poly_srcptr f)
{
    std::string const s = cxx_mpz_poly(f).print_poly("x");
    *res = strdup(s.c_str());
    return (int)s.size();
}

/* Print coefficients of f. */
void mpz_poly_fprintf(FILE * fp, mpz_poly_srcptr f)
{
    fmt::print(fp, "{}\n", cxx_mpz_poly(f));
}

/* Print f of degree d with the following format
    f0<sep>f1<sep>...<sep>fd\n
   Print only '\n' if f = 0 (ie deg(f) = -1)
*/
void mpz_poly_fprintf_coeffs(FILE * fp, mpz_poly_srcptr f, char const * sep)
{
    fmt::print(fp, "{}\n", mpz_poly_coeff_list(cxx_mpz_poly(f), sep));
}

/* Read a polynomial printed using mpz_poly_fprintf_coeffs, with the same
   separator. */
void mpz_poly_fscanf_coeffs(FILE * fp, mpz_poly_ptr f, char const * sep)
{
    int c, deg = -1;
    mpz_t z;

    mpz_init(z);
    while ((c = getc(fp)) != '\n') {
        ungetc(c, fp);
        int const ret = gmp_fscanf(fp, "%Zd", z);
        ASSERT_ALWAYS(ret == 1);
        deg++;
        mpz_poly_setcoeff(f, deg, z);
        c = getc(fp);
        if (c == '\n')
            break;
        for (char const * s = sep; *s; s++)
            ASSERT_ALWAYS(c == *s);
    }
    mpz_clear(z);
}

/* Print f of degree d with the following format
    <pre><letter>0: f0\n
    <pre><letter>1: f1\n
    ...
    <pre><letter>d: fd\n
   Print nothing if f = 0 (ie deg(f) = -1)
*/
void mpz_poly_fprintf_cado_format(FILE * fp, mpz_poly_srcptr f,
                                  char const letter, char const * prefix)
{
    for (int i = 0; i <= f->deg; i++) {
        if (prefix)
            fputs(prefix, fp);
        gmp_fprintf(fp, "%c%d: %Zd\n", letter, i, f->_coeff[i]);
    }
}

void mpz_poly_asprintf_cado_format(char ** pstr, mpz_poly_srcptr f,
                                   char const letter, char const * prefix)
{
    size_t size = 10;
    char * str = (char *)malloc(size);
    size_t p = 0;
    for (int i = 0; i <= f->deg; i++) {
        for (size_t n = SIZE_MAX;;) {
            n = gmp_snprintf(str + p, size - p, "%s%c%d: %Zd\n",
                             prefix ? prefix : "", letter, i, f->_coeff[i]);
            if (n < size - p) {
                p += n;
                break;
            }
            size *= 2;
            checked_realloc(str, size);
        }
    }
    *pstr = str;
}

void mpz_poly_print_raw(mpz_poly_srcptr f)
{
    cxx_mpz_poly F;
    mpz_poly_set(F, f);
    std::string const s = F.print_poly("x");
    printf("%s\n", s.c_str());
}

/* -------------------------------------------------------------------------- */

/* Tests and comparison functions */

/* return 0 if a and b are equal,
 * -1 if a is "smaller" and 1 if b is "bigger", for some arbitrary
 * ordering, defined as the limit of the evaluation in x as x tends to
 * +infinity
 * Assumes a and b are normalized.
 */
int mpz_poly_cmp(mpz_poly_srcptr a, mpz_poly_srcptr b)
{
    for (int i = std::max(a->deg, b->deg); i >= 0; i--) {
        int r;
        if (i > b->deg)
            r = mpz_sgn(a->_coeff[i]);
        else if (i > a->deg)
            r = -mpz_sgn(b->_coeff[i]);
        else
            r = mpz_cmp(a->_coeff[i], b->_coeff[i]);
        if (r)
            return r;
    }
    return 0;
}

int mpz_poly_cmp_mpz(mpz_poly_srcptr a, mpz_srcptr b)
{
    for (int i = a->deg ; i >= 1; i--) {
        int const r = mpz_sgn(a->_coeff[i]);
        if (r)
            return r;
    }
    if (a->deg < 0) {
        return -mpz_sgn(b);
    } else {
        return mpz_cmp(a->_coeff[0], b);
    }
}

/* return 1 if f is normalized, i.e. f[deg] != 0, or the null polynomial.  */
int mpz_poly_normalized_p(mpz_poly_srcptr f)
{
    return (f->deg == -1) || mpz_sgn(f->_coeff[f->deg]) != 0;
}

/* return 1 if f is monic, i.e. f[deg] == 1, return 0 otherwise (null
 * polynomial is considered monic).
 */
int mpz_poly_is_monic(mpz_poly_srcptr f)
{
    if (f->deg == -1)
        return 1;
    else if (mpz_cmp_ui(f->_coeff[f->deg], 1) == 0)
        return 1;
    else
        return 0;
}

/* -------------------------------------------------------------------------- */

/* Set f=-g.
   Note: f can be the same as g. */
void mpz_poly_neg(mpz_poly_ptr f, mpz_poly_srcptr g)
{
    mpz_poly_realloc(f, g->deg + 1);
    for (int i = 0; i <= g->deg; i++) {
        mpz_neg(f->_coeff[i], g->_coeff[i]);
    }
    f->deg = g->deg;
}

/* Set f=g+h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_add(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h)
{
    int i, maxdeg;
    mpz_t z;
    mpz_init(z);
    maxdeg = std::max(g->deg, h->deg);
    mpz_poly_realloc(f, maxdeg + 1);
    for (i = 0; i <= maxdeg; i++) {
        if (i <= g->deg)
            mpz_set(z, g->_coeff[i]);
        else
            mpz_set_ui(z, 0);
        if (i <= h->deg)
            mpz_add(z, z, h->_coeff[i]);
        mpz_set(f->_coeff[i], z);
    }
    f->deg = maxdeg;
    mpz_clear(z);
    mpz_poly_cleandeg(f, maxdeg);
}

/* Set f=g-h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_sub(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h)
{
    mpz_poly_notparallel_info().mpz_poly_sub(f, g, h);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_sub(mpz_poly_ptr f,
                                                    mpz_poly_srcptr g,
                                                    mpz_poly_srcptr h) const
{
    int const maxdeg = std::max(g->deg, h->deg);
    mpz_poly_realloc(f, maxdeg + 1);
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= maxdeg; i++) {
        if (i <= g->deg && i <= h->deg)
            mpz_sub(f->_coeff[i], g->_coeff[i], h->_coeff[i]);
        else if (i <= g->deg)
            mpz_set(f->_coeff[i], g->_coeff[i]);
        else
            mpz_neg(f->_coeff[i], h->_coeff[i]);
    }
    f->deg = maxdeg;
    mpz_poly_cleandeg(f, maxdeg);
}

/* Set g=f+a where a is unsigned long. */
void mpz_poly_add_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_add_ui(g->_coeff[0], f->_coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->_coeff[0], a);
        g->deg = 0;
    }
}

/* Set g=f-a where a is unsigned long. */
void mpz_poly_sub_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_sub_ui(g->_coeff[0], f->_coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->_coeff[0], a);
        mpz_neg(g->_coeff[0], g->_coeff[0]);
        g->deg = 0;
    }
}

void mpz_poly_add_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a)
{
    mpz_poly_set(f, g);
    if (f->deg == -1) {
        mpz_poly_set_mpz(f, a);
    } else {
        mpz_add(f->_coeff[0], f->_coeff[0], a);
        if (f->deg == 0)
            mpz_poly_cleandeg(f, 0);
    }
}

void mpz_poly_sub_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a)
{
    mpz_poly_set(f, g);
    if (f->deg == -1) {
        mpz_poly_set_mpz(f, a);
        mpz_neg(f->_coeff[0], f->_coeff[0]);
    } else {
        mpz_sub(f->_coeff[0], f->_coeff[0], a);
        if (f->deg == 0)
            mpz_poly_cleandeg(f, 0);
    }
}

/* Set f=g-h (mod m). */
void mpz_poly_sub_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h,
                          mpz_srcptr m)
{
    mpz_poly_notparallel_info().mpz_poly_sub_mod_mpz(f, g, h, m);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_sub_mod_mpz(mpz_poly_ptr f,
                                                            mpz_poly_srcptr g,
                                                            mpz_poly_srcptr h,
                                                            mpz_srcptr m) const
{
    mpz_poly_sub(f, g, h);
    mpz_poly_mod_mpz(f, f, m, nullptr);
}

/* Set f=g*h. Note: f might equal g or h.
   Assumes the g[g->deg] and h[h->deg[] are not zero. */
void mpz_poly_mul(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h)
{
    mpz_poly_notparallel_info().mpz_poly_mul(f, g, h);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_mul(mpz_poly_ptr f,
                                                    mpz_poly_srcptr g,
                                                    mpz_poly_srcptr h) const
{
    if (f == h || f == g) {
        mpz_poly aux;
        mpz_poly_init(aux, -1);
        mpz_poly_mul(aux, g, h);
        mpz_poly_set(f, aux);
        mpz_poly_clear(aux);
        return;
    }

    /* now f differs from g and h */

    if ((g->deg == -1) || (h->deg == -1)) {
        f->deg = -1;
        return;
    }

    mpz_poly_realloc(f, g->deg + h->deg + 1);

    if (g == h) /* this is a square */
    {
        f->deg = mpz_poly_sqr_tc((inf &)*this, f->_coeff, g->_coeff, g->deg);
        return;
    }

    if (g->deg == 0) {
        mpz_poly_mul_mpz(f, h, g->_coeff[0]);
        return;
    }

    if (h->deg == 0) {
        mpz_poly_mul_mpz(f, g, h->_coeff[0]);
        return;
    }

    ASSERT(mpz_cmp_ui(g->_coeff[g->deg], 0) != 0);
    ASSERT(mpz_cmp_ui(h->_coeff[h->deg], 0) != 0);
    ASSERT(f != g);
    ASSERT(f != h);

#if 1
    f->deg = mpz_poly_mul_tc((inf &)*this, f->_coeff, g->_coeff, g->deg,
                             h->_coeff, h->deg);
#else /* segmentation, this code has problem with huge runs, for example       \
         degree 5 with lifting to 631516975 bits */
    {
        mpz_t G, H;
        size_t sg, sh, s;
        int i;

        mpz_init(G);
        mpz_init(H);
        sg = mpz_poly_sizeinbase(g, 2);
        sh = mpz_poly_sizeinbase(h, 2);
        /* the +1 accounts for a possible sign */
        for (s = sg + sh + 1, i = (g->deg >= h->deg) ? h->deg + 1 : g->deg + 1;
             i > 1; i = (i + 1) / 2, s++)
            ;
        mpz_set(G, g->_coeff[g->deg]);
        for (i = g->deg - 1; i >= 0; i--) {
            mpz_mul_2exp(G, G, s);
            mpz_add(G, G, g->_coeff[i]);
        }
        /* sanity check: G should have sizeinbase(lc(g))+d*s bits (or -1) */
        size_t size_g = mpz_sizeinbase(g->_coeff[g->deg], 2) + g->deg * s;
        ASSERT(mpz_sizeinbase(G, 2) == size_g ||
               mpz_sizeinbase(G, 2) == size_g - 1);
        mpz_set(H, h->_coeff[h->deg]);
        for (i = h->deg - 1; i >= 0; i--) {
            mpz_mul_2exp(H, H, s);
            mpz_add(H, H, h->_coeff[i]);
        }
        /* sanity check: H should have sizeinbase(lc(h))+d*s bits (or -1) */
        size_t size_h = mpz_sizeinbase(h->_coeff[h->deg], 2) + h->deg * s;
        ASSERT(mpz_sizeinbase(H, 2) == size_h ||
               mpz_sizeinbase(H, 2) == size_h - 1);
        size_g = mpz_sizeinbase(G, 2);
        size_h = mpz_sizeinbase(H, 2);
        /* sanity check: we verify that the product agrees both mod B and B-1 */
        mp_limb_t g0 = mpz_getlimbn(G, 0);
        mpz_mul(G, G, H);
        ASSERT(mpz_getlimbn(G, 0) == g0 * mpz_getlimbn(H, 0));
        ASSERT(mpz_sizeinbase(G, 2) == size_g + size_h ||
               mpz_sizeinbase(G, 2) == size_g + size_h - 1);
        for (i = 0; i < g->deg + h->deg; i++) {
            mpz_fdiv_r_2exp(f->_coeff[i], G, s);
            if (mpz_sizeinbase(f->_coeff[i], 2) == s) {
                mpz_cdiv_r_2exp(f->_coeff[i], G, s);
                mpz_cdiv_q_2exp(G, G, s);
            } else
                mpz_fdiv_q_2exp(G, G, s);
            ASSERT(mpz_sizeinbase(f->_coeff[i], 2) < s);
        }
        mpz_set(f->_coeff[i], G);
        mpz_clear(G);
        mpz_clear(H);
        f->deg = g->deg + h->deg;
    }
#endif
    /* there is no need to run mpz_poly_cleandeg since g[g->deg] <> 0
       and h[h->deg] <> 0 */
    ASSERT(mpz_cmp_ui(f->_coeff[f->deg], 0) != 0);
}

/* Set Q=a*P, where a is an mpz_t */
void mpz_poly_mul_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a)
{
    mpz_poly_notparallel_info().mpz_poly_mul_mpz(f, g, a);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_mul_mpz(mpz_poly_ptr Q,
                                                        mpz_poly_srcptr P,
                                                        mpz_srcptr a) const
{
    mpz_poly_realloc(Q, P->deg + 1);
#ifdef HAVE_OPENMP
#pragma omp parallel if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    {
        mpz_t aux;
        mpz_init(aux);
#ifdef HAVE_OPENMP
#pragma omp for
#endif
        for (int i = 0; i <= P->deg; i++) {
            mpz_mul(aux, P->_coeff[i], a);
            mpz_set(Q->_coeff[i], aux);
        }
        mpz_clear(aux);
    }
    mpz_poly_cleandeg(Q, P->deg);
}

/* Set Q=P/a, where a is an mpz_t. Assume a divides the content of P (the
   division is done with mpz_divexact). Otherwise the result is not correct. */
void mpz_poly_divexact_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a)
{
    mpz_poly_notparallel_info().mpz_poly_divexact_mpz(f, g, a);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_divexact_mpz(mpz_poly_ptr Q,
                                                             mpz_poly_srcptr P,
                                                             mpz_srcptr a) const
{
    mpz_poly_realloc(Q, P->deg + 1);
#ifdef HAVE_OPENMP
#pragma omp parallel if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    {
        mpz_t aux;
        mpz_init(aux);
#ifdef HAVE_OPENMP
#pragma omp for
#endif
        for (int i = 0; i <= P->deg; i++) {
            mpz_divexact(aux, P->_coeff[i], a);
            mpz_set(Q->_coeff[i], aux);
        }
        mpz_clear(aux);
    }
    mpz_poly_cleandeg(Q, P->deg);
}

/* Test whether a divides P, where a is an mpz_t. */
int mpz_poly_divisible_mpz(mpz_poly_srcptr P, mpz_srcptr a)
{
    for (int i = 0; i <= P->deg; ++i)
        if (!mpz_divisible_p(P->_coeff[i], a))
            return 0;
    return 1;
}

/* Set ft(x) = f(x+k) where k is an mpz_t, ft and f can be the same poly. */
void mpz_poly_translation(mpz_poly_ptr ft, mpz_poly_srcptr f, mpz_srcptr k)
{
    int i, j;
    int const d = f->deg;

    mpz_poly_set(ft, f);
    for (i = d - 1; i >= 0; i--)
        for (j = i; j < d; j++)
            mpz_addmul(ft->_coeff[j], ft->_coeff[j + 1], k);
}

/* Set fr = f + k * x^t * g, with t >= 0
 * fr and f can be the same poly */
void mpz_poly_rotation(mpz_poly_ptr fr, mpz_poly_srcptr f, mpz_poly_srcptr g,
                       mpz_srcptr k, int t)
{
    mpz_poly_set(fr, f);
    ASSERT_ALWAYS(t >= 0);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_addmul(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}



void mpz_poly_reverse_rotation(mpz_poly_ptr fr, mpz_poly_srcptr f,
                               mpz_poly_srcptr g, mpz_srcptr k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_submul(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

void mpz_poly_rotation_si(mpz_poly_ptr fr, mpz_poly_srcptr f, mpz_poly_srcptr g,
                          long int k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_addmul_si(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

void mpz_poly_reverse_rotation_si(mpz_poly_ptr fr, mpz_poly_srcptr f,
                                  mpz_poly_srcptr g, long int k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_submul_si(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

void mpz_poly_rotation_ui(mpz_poly_ptr fr, mpz_poly_srcptr f, mpz_poly_srcptr g,
                          unsigned long int k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_addmul_ui(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

void mpz_poly_reverse_rotation_ui(mpz_poly_ptr fr, mpz_poly_srcptr f,
                                  mpz_poly_srcptr g, unsigned long int k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_submul_ui(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

/* Set h = fr + k * x^t * g with t >= 0 */
void mpz_poly_rotation_int64(mpz_poly_ptr fr, mpz_poly_srcptr f,
                             mpz_poly_srcptr g, int64_t const k, int t)
{
    mpz_poly_set(fr, f);
    mpz_poly_realloc(fr, t + g->deg + 1);
    /* If f is not yet such that t+deg(g) <= deg(f), it is important that
     * all newly added coefficients are set to zero. realloc does not do
     * this (well, maybe it's a bug, really) */
    for (int i = fr->deg + 1; i <= t + g->deg; i++) {
        mpz_set_ui(fr->_coeff[i], 0);
    }
    for (int i = 0; i <= g->deg; i++)
        mpz_addmul_int64(fr->_coeff[i + t], g->_coeff[i], k);
    mpz_poly_cleandeg(fr, std::max(fr->deg, t + g->deg));
}

/* Set f = f + k * g */
void mpz_poly_addmul_si(mpz_poly_ptr f, mpz_poly_srcptr g, long k)
{
    mpz_poly_realloc(f, g->deg + 1);
    for (int i = f->deg + 1; i <= g->deg; i++)
        mpz_set_ui(f->_coeff[i], 0);
    for (int i = 0; i <= g->deg; i++)
        mpz_addmul_si(f->_coeff[i], g->_coeff[i], k);
    mpz_poly_cleandeg(f, std::max(f->deg, g->deg));
}

/* Set f = k * g such that deg(g) <= deg(f) (this assumption is not
   checked). */
void mpz_poly_mul_si(mpz_poly_ptr f, mpz_poly_srcptr g, long k)
{
    mpz_poly_realloc(f, g->deg + 1);
    for (int i = 0; i <= g->deg; i++)
        mpz_mul_si(f->_coeff[i], g->_coeff[i], k);
    mpz_poly_cleandeg(f, g->deg);
}

/* Set f = g / k such that k divides g
   (this assumption is not checked). */
void mpz_poly_divexact_ui(mpz_poly_ptr f, mpz_poly_srcptr g, unsigned long k)
{
    mpz_poly_realloc(f, g->deg + 1);
    for (int i = 0; i <= g->deg; i++)
        mpz_divexact_ui(f->_coeff[i], g->_coeff[i], k);
    mpz_poly_cleandeg(f, g->deg);
}

/* h=rem(h, f) mod N, f not necessarily monic, N not necessarily prime */
/* Coefficients of f must be reduced mod N
 * Coefficients of h need not be reduced mod N on input, but are reduced
 * on output.
 *
 * If division fails, stores a non-trivial factor of N in factor. This is
 * not done if factor==NULL.
 */
static int mpz_poly_pseudodiv_r(mpz_poly_ptr h, mpz_poly_srcptr f, mpz_srcptr N,
                                mpz_ptr factor)
{
    int const d = f->deg;
    int dh = h->deg;
    mpz_t tmp, inv;

    mpz_init_set_ui(inv, 1);
    mpz_init_set_ui(tmp, 1);

    mpz_set(tmp, f->_coeff[d]);
    if (mpz_cmp_ui(tmp, 1) != 0) {
        /* inv is 1/f[d] mod N */
        if (!mpz_invert(inv, tmp, N)) {
            if (factor != nullptr)
                mpz_gcd(factor, tmp, N);
            mpz_clear(inv);
            mpz_clear(tmp);
            return 0;
        }
    }

    while (dh >= d) {
        /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
        if (mpz_cmp_ui(inv, 1) != 0) {
            mpz_mul(h->_coeff[dh], h->_coeff[dh], inv);
            mpz_ndiv_r(h->_coeff[dh], h->_coeff[dh], N);
        }

        for (int i = 0; i < d; i++) {
            mpz_mul(tmp, h->_coeff[dh], f->_coeff[i]);
            mpz_mod(tmp, tmp, N);
            mpz_sub(h->_coeff[dh - d + i], h->_coeff[dh - d + i], tmp);
            mpz_ndiv_r(h->_coeff[dh - d + i], h->_coeff[dh - d + i], N);
        }

        do {
            dh--;
        } while (dh >= 0 && mpz_divisible_p(h->_coeff[dh], N));

        h->deg = dh;
    }

    mpz_clear(inv);
    mpz_clear(tmp);
    return 1;
}

/* h=rem(h, f) mod p, f not necessarily monic. */
/* returns 0 if an inverse of the leading coeff could not be found. */
int mpz_poly_div_r_mod_mpz_clobber(mpz_poly_ptr h, mpz_poly_srcptr f,
                                   mpz_srcptr p)
{
    return mpz_poly_pseudodiv_r(h, f, p, nullptr);
}

/*
   computes q, r such that f = q*g + r mod p, with deg(r) < deg(g)
   and p in mpz_t
   q and r must be allocated!
   f == r or f == q are allowed.
*/
int mpz_poly_div_qr_mod_mpz(mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f,
                            mpz_poly_srcptr g, mpz_srcptr p)
{
    int k, j, df = f->deg, dg = g->deg, dq = df - dg;
    mpz_t lg, invlg;

    if (df < dg) /* f is already reduced mod g */
    {
        mpz_poly_set (r, f);
        mpz_poly_set_zero (q);
        return 1;
    }

    /* now df >= dg */
    mpz_poly_realloc(q, dq + 1);

    mpz_init(lg);
    mpz_init_set_ui(invlg, 1);

    mpz_poly_set(r, f);
    q->deg = dq;

    mpz_set(lg, g->_coeff[dg]);
    mpz_mod(lg, lg, p);
    /* invlg = 1/g[dg] mod p */
    if (mpz_cmp_ui(lg, 1) != 0)
        if (!mpz_invert(invlg, lg, p)) {
            mpz_clear(lg);
            mpz_clear(invlg);
            return 0;
        }

    for (k = df - dg; k >= 0; k--) {
        mpz_mul(q->_coeff[k], r->_coeff[k + dg], invlg);
        mpz_mod(q->_coeff[k], q->_coeff[k], p);
        for (j = dg + k; j >= k; j--) {
            mpz_submul(r->_coeff[j], q->_coeff[k], g->_coeff[j - k]);
            mpz_mod(r->_coeff[j], r->_coeff[j], p);
        }
    }
    mpz_poly_cleandeg(r, r->deg);

    mpz_clear(invlg);
    mpz_clear(lg);
    return 1;
}

/* This also computes q and r such that f = q * g + r, but over Z, not
 * modulo a prime. Also, we do not assume that g is monic. Of course, if
 * it is not, then most often the result will be undefined (over Z). We
 * return something well-defined if q and r happen to be integer
 * polynomials.
 * We return 0 if this is not the case (in which case q and r are
 * undefined).
 * r==f is allowed.
 *
 * TODO: do the parallel version
 */
int mpz_poly_div_qr(mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f,
                    mpz_poly_srcptr g)
{
    int k, j, df = f->deg, dg = g->deg, dq = df - dg;

    if (df < dg) /* f is already reduced mod g */
    {
        mpz_poly_set_zero(q);
        mpz_poly_set(r, f);
        return 1;
    }

    /* now df >= dg */
    /* realloc takes a number of coefficients, not a degree */
    mpz_poly_realloc(q, dq + 1);

    mpz_poly_set(r, f);
    q->deg = dq;

    mpz_srcptr lg = g->_coeff[dg];

    for (k = df - dg; k >= 0; k--) {
        if (!mpz_divisible_p(r->_coeff[k + dg], lg))
            return 0;
        mpz_divexact(q->_coeff[k], r->_coeff[k + dg], lg);
        for (j = dg + k; j >= k; j--) {
            mpz_submul(r->_coeff[j], q->_coeff[k], g->_coeff[j - k]);
        }
    }
    mpz_poly_cleandeg(r, r->deg);
    return 1;
}

int mpz_poly_div_r(mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g)
{
    cxx_mpz_poly quo;
    return mpz_poly_div_qr(quo, r, f, g);
}

int mpz_poly_mod(mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g)
{
    return mpz_poly_div_r(r, f, g);
}

/* q=divexact(h, f) mod p, f not necessarily monic.
   Assumes lc(h) <> 0 mod p.
   Clobbers h. */
/* Coefficients of f must be reduced mod p
 * Coefficients of h need not be reduced mod p
 * Coefficients of q are reduced mod p
 */
static int mpz_poly_divexact_clobber(mpz_poly_ptr q, mpz_poly_ptr h,
                                     mpz_poly_srcptr f, mpz_srcptr p)
{
    int i, d = f->deg, dh = h->deg;
    mpz_t t, aux;

    mpz_init(t);
    mpz_init(aux);
    ASSERT(d >= 0);
    ASSERT(dh >= 0);
    ASSERT(dh >= d);
    ASSERT(mpz_divisible_p(h->_coeff[dh], p) == 0);

    mpz_poly_realloc(q, dh + 1 - d);
    q->deg = dh - d;
    /* t is 1/f[d] mod p */
    mpz_set(aux, f->_coeff[d]);
    mpz_mod(aux, aux, p);
    if (!mpz_invert(t, aux, p)) {
        mpz_clear(t);
        mpz_clear(aux);
        return 0;
    }

    while (dh >= d) {

        /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
        if (mpz_cmp_ui(t, 1) != 0) {
            mpz_mul(h->_coeff[dh], h->_coeff[dh], t);
            mpz_mod(h->_coeff[dh], h->_coeff[dh], p);
        }
        mpz_set(q->_coeff[dh - d], h->_coeff[dh]);
        mpz_mod(q->_coeff[dh - d], q->_coeff[dh - d], p);

        /* we only need to update the coefficients of degree >= d of h,
           i.e., we want i >= 2d - dh */
        for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
            mpz_mul(aux, h->_coeff[dh], f->_coeff[i]);
            mpz_sub(h->_coeff[dh - d + i], h->_coeff[dh - d + i], aux);
            mpz_mod(h->_coeff[dh - d + i], h->_coeff[dh - d + i], p);
        }
        dh--;
    }
    /* since lc(h) <> 0 mod p, q is normalized */

    mpz_clear(t);
    mpz_clear(aux);
    return 1;
}

/* q <- divexact(h, f) mod p, f not necessarily monic. */
/* Coefficients of f must be reduced mod p
 * Coefficients of h need not be reduced mod p
 * Coefficients of q are reduced mod p
 */
int mpz_poly_divexact(mpz_poly_ptr q, mpz_poly_srcptr h, mpz_poly_srcptr f,
                      mpz_srcptr p)
{
    cxx_mpz_poly hh = h;
    return mpz_poly_divexact_clobber(q, hh, f, p);
}

/* Set f=g/2 (mod m), where f might equal g.
   Assumes m is odd. */
/* If coefficients of g are reduced mod m, then coefficients of f are
 * reduced.
 */
void mpz_poly_div_2_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a)
{
    mpz_poly_notparallel_info().mpz_poly_div_2_mod_mpz(f, g, a);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_div_2_mod_mpz(
    mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr m) const
{
    ASSERT_ALWAYS(mpz_scan1(m, 0) == 0);

    mpz_poly_realloc(f, g->deg + 1);

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = g->deg; i >= 0; --i) {
        if (mpz_scan1(g->_coeff[i], 0) == 0) /* g[i] is odd */
        {
            mpz_add(f->_coeff[i], g->_coeff[i], m);
            mpz_div_2exp(f->_coeff[i], f->_coeff[i], 1);
        } else
            mpz_div_2exp(f->_coeff[i], g->_coeff[i], 1);
    }
}

/* Set res=f(x). Assumes res and x are different variables. */
void mpz_poly_eval(mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x)
{
    int i, d;
    ASSERT_ALWAYS(res != x);
    d = f->deg;
    if (d == -1) {
        mpz_set_ui(res, 0);
        return;
    }
    mpz_set(res, f->_coeff[d]);
    for (i = d - 1; i >= 0; --i) {
        mpz_mul(res, res, x);
        mpz_add(res, res, f->_coeff[i]);
    }
}

/* Set res=f(x). Assumes res and x are different variables. */
void mpz_poly_eval_poly(mpz_poly_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr x)
{
    int i, d;
    ASSERT_ALWAYS(res != x);
    ASSERT_ALWAYS(res != f);
    d = f->deg;
    if (d == -1) {
        mpz_poly_set_zero(res);
        return;
    }
    mpz_poly_set_mpz(res, f->_coeff[d]);
    for (i = d - 1; i >= 0; --i) {
        mpz_poly_mul(res, res, x);
        mpz_poly_add_mpz(res, res, f->_coeff[i]);
    }
}

/* Set res=f(x) where x is an unsigned long. */
void mpz_poly_eval_ui(mpz_ptr res, mpz_poly_srcptr f, unsigned long x)
{
    int const d = f->deg;

    mpz_set(res, f->_coeff[d]);
    for (int i = d - 1; i >= 0; i--) {
        mpz_mul_ui(res, res, x);
        mpz_add(res, res, f->_coeff[i]);
    }
}

/* Set res=f'(x), where x is an unsigned long */
void mpz_poly_eval_diff_ui(mpz_ptr res, mpz_poly_srcptr f, unsigned long x)
{
    int const d = f->deg;

    mpz_mul_ui(res, f->_coeff[d], d);
    for (int i = d - 1; i >= 1; i--) {
        mpz_mul_ui(res, res, x);
        mpz_addmul_ui(res, f->_coeff[i], i); /* res <- res + i*f[i] */
    }
}

void mpz_poly_eval_diff(mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x)
{
    int const d = f->deg;
    ASSERT_ALWAYS(res != x);

    mpz_mul_ui(res, f->_coeff[d], d);
    for (int i = d - 1; i >= 1; i--) {
        mpz_mul(res, res, x);
        mpz_addmul_ui(res, f->_coeff[i], i); /* res <- res + i*f[i] */
    }
}

/* Set res=f'(x), where x is a polynomial */
void mpz_poly_eval_diff_poly(mpz_poly_ptr res, mpz_poly_srcptr f,
                             mpz_poly_srcptr x)
{
    ASSERT_ALWAYS(res != x);
    ASSERT_ALWAYS(res != f);
    int const d = f->deg;
    mpz_poly_realloc(res, f->deg * x->deg + 1);
    mpz_t t;
    mpz_init(t);
    mpz_mul_ui(t, f->_coeff[d], d);
    mpz_poly_add_mpz(res, res, t);
    for (int i = d - 1; i >= 1; i--) {
        mpz_poly_mul(res, res, x);
        mpz_mul_ui(t, f->_coeff[i], i); /* res <- res + i*f[i] */
        mpz_poly_add_mpz(res, res, t);
    }
    mpz_clear(t);
}

/* Return 1 if poly(root) % modulus == 0, return 0 otherwise */
/* Coefficients of f(x) need not be reduced mod m */
int mpz_poly_is_root(mpz_poly_srcptr poly, mpz_srcptr root, mpz_srcptr modulus)
{
    cxx_mpz x;
    mpz_poly_eval_mod_mpz(x, poly, root, modulus);
    return mpz_cmp_ui(x, 0) == 0;
}

/* Set res=f(x) (mod m).  Assume res and x are different variables. */
/* Coefficients of f(x) need not be reduced mod m.
 * The computed value res is reduced mod m
 */
void mpz_poly_eval_mod_mpz(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x,
                           mpz_srcptr m)
{
    int i, d;

    d = f->deg;
    if (d == -1) {
        mpz_set_ui(res, 0);
        return;
    }
    mpz_mod(res, f->_coeff[d], m);
    for (i = d - 1; i >= 0; --i) {
        mpz_mul(res, res, x);
        mpz_add(res, res, f->_coeff[i]);
        mpz_mod(res, res, m);
    }
}

/* This evaluates several polynomials at the same point w. It is possible
 * to do fewer modular reductions in this case.
 *
 * When k==1, we use mpz_poly_eval_mod_mpz instead, since it's faster
 * then.
 */
/* Coefficients of f[i](x) need not be reduced mod m.
 * The computed values r[i] are reduced mod m
 */
void mpz_poly_eval_several_mod_mpz(mpz_ptr * r, mpz_poly_srcptr * f, int k,
                                   mpz_srcptr x, mpz_srcptr m)
{
    int i;

    if (k == 1) {
        mpz_poly_eval_mod_mpz(r[0], f[0], x, m);
        return;
    }

    mpz_t w;
    mpz_init(w);
    int maxdeg = -1;
    for (int j = 0; j < k; j++) {
        if (f[j]->deg >= 0)
            mpz_set(r[j], f[j]->_coeff[0]);
        else
            mpz_set_ui(r[j], 0);
        maxdeg = std::max(maxdeg, f[j]->deg);
    }

    mpz_set(w, x);
    for (int j = 0; j < k; j++) {
        if (f[j]->deg >= 1)
            mpz_addmul(r[j], w, f[j]->_coeff[1]);
    }
    for (i = 2; i <= maxdeg; i++) {
        mpz_mul(w, w, x);
        mpz_mod(w, w, m);
        for (int j = 0; j < k; j++)
            if (f[j]->deg >= i)
                mpz_addmul(r[j], w, f[j]->_coeff[i]);
    }
    for (int j = 0; j < k; j++) {
        mpz_mod(r[j], r[j], m);
    }
    mpz_clear(w);
}

/* Set Q = P/lc(P) (mod m). Q and P might be identical. */
/* Coefficients of P need not be reduced mod m
 * Coefficients of Q are reduced mod m
 */
void mpz_poly_makemonic_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr m)
{
    int i;
    mpz_t aux;
    mpz_init(aux);
    for (i = P->deg; i >= 0; i--) {
        mpz_mod(aux, P->_coeff[i], m);
        if (mpz_cmp_ui(aux, 0) != 0)
            break;
    }
    /* i is the degree of the leading monomial */
    Q->deg = i;
    if (i < 0) {
        /* if i == -1, then Q is the zero polynomial, there's nothing to do */
        mpz_clear(aux);
        return;
    }

    mpz_t aux2;
    mpz_init(aux2);
    mpz_invert(aux2, aux, m);
    for (i = 0; i < Q->deg; ++i) {
        mpz_mul(aux, aux2, P->_coeff[i]);
        mpz_mod(aux, aux, m);
        mpz_poly_setcoeff(Q, i, aux);
    }
    mpz_clear(aux2);

    /* we can directly set the leading coefficient to 1 */
    mpz_poly_setcoeff_si(Q, Q->deg, 1);
    mpz_clear(aux);
}

/* Algorithm 2.5 from "Modern Computer Arithmetic" */
void barrett_precompute_inverse(mpz_ptr invm, mpz_srcptr m)
{
    size_t const n = mpz_sizeinbase(m, 2);
    ASSERT_ALWAYS(mpz_cmp_ui(m, 0) > 0);
    /* with B = 2^n, we have B/2 <= m < B */
    mpz_set_ui(invm, 0);
    mpz_setbit(invm, 2 * n);   /* invm = B^2 */
    mpz_tdiv_q(invm, invm, m); /* floor(B^2/m) */
}

/* r <- a mod m */
static void mpz_mod_barrett(mpz_ptr r, mpz_srcptr a, mpz_srcptr m,
                            mpz_srcptr invm)
{
    size_t const n = mpz_sizeinbase(m, 2);
    mpz_srcptr r_or_a = a;

    while (mpz_sizeinbase(r_or_a, 2) > n + 1) {
        mpz_t a1;
        size_t const sr = mpz_sizeinbase(r_or_a, 2);
        mpz_init(a1);
        /* if sr <= 2n we consider the sr-n most significant bits of r,
           otherwise we take the n most significant bits */
        mpz_tdiv_q_2exp(a1, r_or_a, (sr <= 2 * n) ? n : sr - n);
        mpz_mul(a1, a1, invm);
        mpz_tdiv_q_2exp(a1, a1, n);
        mpz_mul(a1, a1, m);
        /* if sr > 2*n we have to multiply by 2^(sr-2n) */
        if (sr >= 2 * n)
            mpz_mul_2exp(a1, a1, sr - 2 * n);
        mpz_sub(r, r_or_a, a1);
        r_or_a = r;
        mpz_clear(a1);
    }
    /* now r_or_a has at most n bits */
    while (mpz_cmpabs(r_or_a, m) >= 0) {
        if (mpz_cmp_ui(r_or_a, 0) > 0)
            mpz_sub(r, r_or_a, m);
        else
            mpz_add(r, r_or_a, m);
        r_or_a = r;
    }
    /* now |r_or_a| < m */
    if (mpz_cmp_ui(r_or_a, 0) < 0) {
        mpz_add(r, r_or_a, m);
        r_or_a = r;
    }
    /* now 0 <= r_or_a < m */
    if (r_or_a == a && r != a)
        mpz_set(r, r_or_a);
}

/* Coefficients of A need not be reduced mod m
 * Coefficients of R are reduced mod m
 * If invm = NULL, use mpz_mod.
 * Otherwise use mpz_mod_barrett, where invm should have been precomputed
 * using barrett_precompute_inverse (invm, m).
 */
int mpz_poly_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m,
                     mpz_srcptr invm)
{
    return mpz_poly_notparallel_info().mpz_poly_mod_mpz(R, A, m, invm);
}
template <typename inf>
int mpz_poly_parallel_interface<inf>::mpz_poly_mod_mpz(mpz_poly_ptr R,
                                                       mpz_poly_srcptr A,
                                                       mpz_srcptr m,
                                                       mpz_srcptr invm) const
{
    /* reduce lower coefficients */
    mpz_poly_realloc(R, A->deg + 1);
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= A->deg; ++i)
        if (invm == nullptr)
            mpz_mod(R->_coeff[i], A->_coeff[i], m);
        else
            mpz_mod_barrett(R->_coeff[i], A->_coeff[i], m, invm);

    mpz_poly_cleandeg(R, A->deg);
    return R->deg;
}

/* reduce non-negative coefficients in [0, p-1], negative ones in [1-p, -1] */
int mpz_poly_mod_mpz_lazy(mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m)
{
    mpz_poly_realloc(R, A->deg + 1);
    for (int i = 0; i <= A->deg; ++i) {
        if (mpz_sgn(A->_coeff[i]) >= 0)
            mpz_mod(R->_coeff[i], A->_coeff[i], m);
        else {
            mpz_neg(R->_coeff[i], A->_coeff[i]);
            mpz_mod(R->_coeff[i], R->_coeff[i], m);
            mpz_neg(R->_coeff[i], R->_coeff[i]);
        }
    }

    mpz_poly_cleandeg(R, A->deg);
    return R->deg;
}

/* Reduce R[d]*x^d + ... + R[0] mod f[df]*x^df + ... + f[0] modulo m.
   Return the degree of the remainder.
   Coefficients of f must be reduced mod m on input.
   Coefficients of R need not be reduced mod m on input, but are reduced
   on output.
   If invf is not NULL, it should be 1/m mod lc(f). */
int mpz_poly_mod_f_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m,
                           mpz_srcptr invf, mpz_srcptr invm)
{
    return mpz_poly_notparallel_info().mpz_poly_mod_f_mod_mpz(R, f, m, invf,
                                                              invm);
}
template <typename inf>
int mpz_poly_parallel_interface<inf>::mpz_poly_mod_f_mod_mpz(
    mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf,
    mpz_srcptr invm) const
{
    cxx_mpz c;
    std::unique_ptr<cxx_mpz> aux;
    mpz_srcptr invf_local = nullptr;
    size_t size_f, size_R;

    if (f == nullptr) {
        /* This is really an obscure feature. We should get rid of it.
         * It's used, though.
         */
        goto reduce_R;
    }

    /* The invf parameter is ignored if f is monic */
    if (!mpz_poly_is_monic(f)) {
        if (invf) {
            invf_local = invf;
        } else {
            aux = std::make_unique<cxx_mpz>();
            /* aux = 1/m mod lc(f) */
            mpz_invert(*aux, m, f->_coeff[f->deg]);
            invf_local = *aux;
        }
    }

    size_f = mpz_poly_size(f);
    size_R = mpz_poly_size(R);

    // FIXME: write a subquadratic variant
    while (R->deg >= f->deg) {
        /* Here m is large (thousand to million bits) and lc(f) is small
         * (typically one word). We first subtract lambda * m * x^(dR-df)
         * ---which is zero mod m--- to R such that the new coefficient
         * of degree dR is divisible by lc(f), i.e., lambda = lc(R)/m
         * mod lc(f). Then if c = (lc(R) - lambda * m) / lc(f), we
         * subtract c * x^(dR-df) * f.
         * Of course, if f is monic, we don't have to do that.
         */
        if (invf_local) {
            mpz_mod(c, R->_coeff[R->deg], f->_coeff[f->deg]); /* lc(R) mod lc(f) */
            mpz_mul(c, c, invf_local);
            mpz_mod(c, c, f->_coeff[f->deg]); /* lc(R)/m mod lc(f) */
            mpz_submul(R->_coeff[R->deg], m, c); /* lc(R) - m * (lc(R) / m mod lc(f)) */
        }
        ASSERT(mpz_divisible_p(R->_coeff[R->deg], f->_coeff[f->deg]));
        mpz_divexact(c, R->_coeff[R->deg], f->_coeff[f->deg]);
        /* If R[deg] has initially size 2n, and f[deg] = O(1), then c has size
           2n here. However, in the equal-degree factorization, even if f[deg]
           = O(1), the lower coefficients of f might have n bits. Thus we decide
           to reduce whenever the total size exceeds 2n. */
        /* FIXME: commit a85444984 changed the line below but not the
         * comment above. */
        size_t const size_c = mpz_size(c);
        if (size_c + size_f > (3 * size_R) / 2)
            mpz_mod(c, c, m);
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
        for (int i = R->deg - 1; i >= R->deg - f->deg; --i)
            mpz_submul(R->_coeff[i], c, f->_coeff[f->deg - R->deg + i]);
        R->deg--;
    }

reduce_R:
    mpz_poly_mod_mpz(R, R, m, invm);

    return R->deg;
}

/* Stores in g the polynomial linked to f by:
 *  - alpha is a root of f if and only if lc(f)*alpha is a root of g
 *  - g is monic
 *
 * This has nothing to do with makemonic_mpz! The operation here
 * works in Z[x]
 */
void mpz_poly_to_monic(mpz_poly_ptr g, mpz_poly_srcptr f)
{
    mpz_t fd, temp;

    if (f->deg == -1) {
        mpz_poly_set_zero(g);
        return;
    }

    mpz_init(fd);
    mpz_init(temp);

    mpz_poly_set(g, f);
    for (int k = 0; k < g->deg; k++) {
        mpz_set(temp, mpz_poly_coeff_const(g, k));
        for (int j = 1; j <= g->deg - 1 - k; j++) {
            mpz_mul(temp, temp, mpz_poly_lc(f));
        }
        mpz_poly_setcoeff(g, k, temp);
    }
    mpz_poly_setcoeff_ui(g, g->deg, 1);

    mpz_clear(temp);
    mpz_clear(fd);
}

/*  Reduce frac (= num / denom) mod F mod m ,
    i.e. compute num * denom^-1 mod F mod m .
    The return value is in num, denom is set to constant polynomial 1
    */
/* TODO: We must state the input / output requirements with respect to
 * reduction mod m */
void mpz_poly_reduce_frac_mod_f_mod_mpz(mpz_poly_ptr num, mpz_poly_ptr denom,
                                        mpz_poly_srcptr F, mpz_srcptr m)
{
    mpz_poly_notparallel_info().mpz_poly_reduce_frac_mod_f_mod_mpz(num, denom,
                                                                   F, m);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_reduce_frac_mod_f_mod_mpz(
    mpz_poly_ptr num, mpz_poly_ptr denom, mpz_poly_srcptr F, mpz_srcptr m) const
{
    if (denom->deg == 0) {
        mpz_t inv;
        mpz_init(inv);
        mpz_set(inv, mpz_poly_coeff_const(denom, 0)); /* inv <- denom[0] */
        mpz_invert(inv, inv, m);                      /* inv <- denom[0]^-1 */
        mpz_poly_mul_mpz(num, num, inv);              /* num <- num * inv */
        mpz_poly_mod_mpz(num, num, m, nullptr); /* num <- num * inv mod m */
        mpz_clear(inv);
    } else {
        mpz_poly g, U, V;
        mpz_poly_init(g, 0);
        mpz_poly_init(U, 0);
        mpz_poly_init(V, 0);
        mpz_poly_xgcd_mpz(g, F, denom, U, V, m);
        mpz_poly_mul(num, num, V);
        mpz_poly_mod_f_mod_mpz(num, F, m, nullptr, nullptr);
        mpz_poly_clear(g);
        mpz_poly_clear(U);
        mpz_poly_clear(V);
    }
    mpz_poly_set_zero(denom);
    mpz_poly_setcoeff_si(denom, 0, 1);
}

/* Q = P1*P2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Coefficients of P1 and P2 need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If invf is not NULL, it is 1/m mod lc(f). */
void mpz_poly_mul_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P1,
                                mpz_poly_srcptr P2, mpz_poly_srcptr f,
                                mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm)
{
    mpz_poly_notparallel_info().mpz_poly_mul_mod_f_mod_mpz(Q, P1, P2, f, m,
                                                           invf, invm);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_mul_mod_f_mod_mpz(
    mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f,
    mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm) const
{
    int const d1 = P1->deg;
    int const d2 = P2->deg;
    int d = d1 + d2;
    mpz_poly R;

    mpz_poly_init(R, d);
#ifdef MPZ_POLY_TIMINGS
    if (mpz_fits_sint_p(f->_coeff[0])) {
        START_TIMER;
        d = mpz_poly_mul_tc((inf &)*this, R->_coeff, P1->_coeff, d1, P2->_coeff,
                            d2);
        mpz_poly_cleandeg(R, d);
        END_TIMER(TIMER_MUL);
        // reduce mod f
        RESTART_TIMER;
        mpz_poly_mod_f_mod_mpz(R, f, m, invf, invm);
        END_TIMER(TIMER_RED);
    } else
#endif
    {
        d = mpz_poly_mul_tc((inf &)*this, R->_coeff, P1->_coeff, d1, P2->_coeff,
                            d2);
        mpz_poly_cleandeg(R, d);
        // reduce mod f
        mpz_poly_mod_f_mod_mpz(R, f, m, invf, invm);
    }
    mpz_poly_set(Q, R);
    mpz_poly_clear(R);
}

/* Q = P1*P2 mod f, assuming f is monic */
void mpz_poly_mul_mod_f(mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2,
                        mpz_poly_srcptr f)
{
    mpz_poly_notparallel_info().mpz_poly_mul_mod_f(Q, P1, P2, f);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_mul_mod_f(
    mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2,
    mpz_poly_srcptr f) const
{
    mpz_poly_mul(Q, P1, P2);
    mpz_poly_div_r(Q, Q, f);
}

/* Q = P^2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Coefficients of P need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If not NULL, invf = 1/m mod lc(f). */
void mpz_poly_sqr_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P,
                                mpz_poly_srcptr f, mpz_srcptr m,
                                mpz_srcptr invf, mpz_srcptr invm)
{
    mpz_poly_notparallel_info().mpz_poly_sqr_mod_f_mod_mpz(Q, P, f, m, invf,
                                                           invm);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_sqr_mod_f_mod_mpz(
    mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr m,
    mpz_srcptr invf, mpz_srcptr invm) const
{
    int const d1 = P->deg;
    int d = d1 + d1;
    mpz_poly R;

    mpz_poly_init(R, d);

    /* Fast squaring in 2d1+1 squares, i.e., 2d-1 squares.
       For d=5, this gives 9 squares. */
    // compute timing only if f has short coefficients
#ifdef MPZ_POLY_TIMINGS
    if (mpz_fits_sint_p(f->_coeff[0])) {
        START_TIMER;
        d = mpz_poly_sqr_tc((inf &)*this, R->_coeff, P->_coeff, d1);
        mpz_poly_cleandeg(R, d);
        END_TIMER(TIMER_SQR);
        // reduce mod f
        RESTART_TIMER;
        mpz_poly_mod_f_mod_mpz(R, f, m, invf, invm);
        END_TIMER(TIMER_RED);
    } else
#endif
    {
        d = mpz_poly_sqr_tc((inf &)*this, R->_coeff, P->_coeff, d1);
        mpz_poly_cleandeg(R, d);
        // reduce mod f
        mpz_poly_mod_f_mod_mpz(R, f, m, invf, invm);
    }

    mpz_poly_set(Q, R);
    mpz_poly_clear(R);
}

/* Affects the derivative of f to df. */
void mpz_poly_derivative(mpz_poly_ptr df, mpz_poly_srcptr f)
{
    if (f->deg <= 0) {
        df->deg = -1;
        return;
    }

    /* This is a no-op if df == f, since f->deg < f->deg + 1 */
    mpz_poly_realloc(df, f->deg - 1 + 1);

    for (int n = 0; n <= f->deg - 1; n++)
        mpz_mul_si(df->_coeff[n], f->_coeff[n + 1], n + 1);

    df->deg = f->deg - 1;
}

/* B = A^n */
void mpz_poly_pow_ui(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n) /*{{{*/
{
    if (n == 0) {
        mpz_poly_set_xi(B, 0);
        return;
    }
    if (A->deg < 0) {
        mpz_poly_set_zero(B);
        return;
    }
    if (B == A) {
        mpz_poly C;
        mpz_poly_init(C, n * A->deg);
        mpz_poly_pow_ui(C, A, n);
        mpz_poly_swap(B, C);
        mpz_poly_clear(C);
        return;
    }
    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    mpz_poly_set(B, A);
    for (; k >>= 1;) {
        mpz_poly_mul(B, B, B);
        if (n & k)
            mpz_poly_mul(B, B, A);
    }
} /*}}}*/

/* B = A^n mod f, assuming f is monic */
void mpz_poly_pow_ui_mod_f(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n,
                           mpz_poly_srcptr f) /*{{{*/
{
    if (n == 0) {
        mpz_poly_set_xi(B, 0);
        return;
    }
    if (A->deg < 0) {
        mpz_poly_set_zero(B);
        return;
    }
    if (B == A || B == f) {
        mpz_poly C;
        mpz_poly_init(C, f->deg - 1);
        mpz_poly_pow_ui_mod_f(C, A, n, f);
        mpz_poly_swap(B, C);
        mpz_poly_clear(C);
        return;
    }
    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    mpz_poly_realloc(B, n * A->deg + 1);
    mpz_poly_set(B, A);
    mpz_poly Q;
    mpz_poly_init(Q, 3 * f->deg - 3);
    for (; k >>= 1;) {
        mpz_poly_mul(Q, B, B);
        mpz_poly_swap(Q, B);
        if (n & k) {
            mpz_poly_mul(Q, B, A);
            mpz_poly_swap(Q, B);
        }
        mpz_poly_div_r(B, B, f);
    }
    mpz_poly_clear(Q);
} /*}}}*/

/* Q = P^a mod f, mod p (f is the algebraic polynomial, non monic)
 * f may be NULL, in case there is only reduction mod p.
 * Coefficients of f must be reduced mod m on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void mpz_poly_pow_mod_f_mod_ui(mpz_poly_ptr Q, mpz_poly_srcptr P,
                               mpz_poly_srcptr f, mpz_srcptr a, unsigned long p)
{
    mpz_poly_notparallel_info().mpz_poly_pow_mod_f_mod_ui(Q, P, f, a, p);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_pow_mod_f_mod_ui(
    mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a,
    unsigned long p) const
{
    mpz_t m;
    mpz_init_set_ui(m, p);
    mpz_poly_pow_mod_f_mod_mpz(Q, P, f, a, m);
    mpz_clear(m);
}

/* Coefficients of f must be reduced mod m on input
 * f may be NULL, in case there is only reduction mod p.
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void mpz_poly_pow_ui_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P,
                                   mpz_poly_srcptr f, unsigned long a,
                                   mpz_srcptr p)
{
    mpz_poly_notparallel_info().mpz_poly_pow_ui_mod_f_mod_mpz(Q, P, f, a, p);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_pow_ui_mod_f_mod_mpz(
    mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, unsigned long a,
    mpz_srcptr p) const
{
    mpz_t az;
    mpz_init_set_ui(az, a);
    mpz_poly_pow_mod_f_mod_mpz(Q, P, f, az, p);
    mpz_clear(az);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t
 * f may be NULL, in case there is only reduction mod p.
 * Coefficients of f must be reduced mod p on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p
 */
void mpz_poly_pow_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P,
                                mpz_poly_srcptr f, mpz_srcptr a, mpz_srcptr p)
{
    mpz_poly_notparallel_info().mpz_poly_pow_mod_f_mod_mpz(Q, P, f, a, p);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_pow_mod_f_mod_mpz(
    mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a,
    mpz_srcptr p) const
{
    int k = mpz_sizeinbase(a, 2), l, L = 0, j;
    mpz_poly R, *T = nullptr;
    mpz_t invf;

    if (mpz_cmp_ui(a, 0) == 0) {
        mpz_poly_set_xi(Q, 0);
        return;
    }

    mpz_poly_init(R, f ? 2 * f->deg : -1);

    // Initialize R to P
    mpz_poly_set(R, P);

    /* compute invf = 1/p mod lc(f) */
    if (f != nullptr) {
        mpz_init(invf);
        mpz_invert(invf, p, f->_coeff[f->deg]);
    }

    /* We use base-2^l exponentiation with sliding window,
       thus we need to precompute P^2, P^3, ..., P^(2^l-1).
       The expected average cost k squarings plus k/(l+1) + 2^(l-1) multiplies.
    */

    for (l = 1; k / (l + 1) + (1 << (l - 1)) > k / (l + 2) + (1 << l); l++)
        ;
    /* this gives (roughly) l=1 for k < 8, l=2 for k < 27, l=3 for k < 84,
       l=4 for k < 245 */
    if (l > 1) {
        L = 1 << (l - 1);
        T = (mpz_poly *)malloc(L * sizeof(mpz_poly));
        /* we store P^2 in T[0], P^3 in T[1], ..., P^(2^l-1) in T[L-1] */
        for (j = 0; j < L; j++)
            mpz_poly_init(T[j], f ? 2 * f->deg : -1);
        mpz_poly_sqr_mod_f_mod_mpz(T[0], R, f, p, invf, nullptr); /* P^2 */
        mpz_poly_mul_mod_f_mod_mpz(T[1], T[0], R, f, p, invf,
                                   nullptr); /* P^3 */
        for (j = 2; j < L; j++)
            mpz_poly_mul_mod_f_mod_mpz(T[j], T[j - 1], T[0], f, p, invf,
                                       nullptr);
    }

    // Horner
    for (k -= 2; k >= 0;) {
        while (k >= 0 && mpz_tstbit(a, k) == 0) {
            mpz_poly_sqr_mod_f_mod_mpz(R, R, f, p, invf, nullptr);
            k--;
        }
        if (k < 0)
            break;
        j = mpz_scan1(a, (k >= l) ? k - (l - 1) : 0);
        /* if l is 1 then j==k*/
        ASSERT_ALWAYS(l > 1 || j == k);
        /* new window starts at bit k, and ends at bit j <= k */
        int e = 0;
        while (k >= j) {
            mpz_poly_sqr_mod_f_mod_mpz(R, R, f, p, invf, nullptr);
            e = 2 * e + mpz_tstbit(a, k);
            k--;
        }
        /* if l is 1 then e == 1 */
        ASSERT_ALWAYS(l > 1 || e == 1);
        mpz_poly_mul_mod_f_mod_mpz(R, R, (e == 1) ? P : T[e / 2], f, p, invf,
                                   nullptr);
    }

    mpz_poly_swap(Q, R);
    mpz_poly_clear(R);
    if (f != nullptr)
        mpz_clear(invf);
    if (l > 1) {
        for (k = 0; k < L; k++)
            mpz_poly_clear(T[k]);
        free(T);
    }
}

/* Return a list of polynomials P[0], P[1], ..., P[l] such that
   P0 = Q[l-1] + p^K[1]*P[l]
   Q[l-1] = Q[l-2] + p^K[2]*P[l-1]
   ...
   Q[2] = Q[1] + p^K[l-1]*P[2]
   Q[1] = Q[0] + p*P[1]         where K[l] = 1
   Q[0] = P[0]
   ...
   With all coefficients of P[i] smaller than p^(K[l-i]-K[l-(i-1)]).
   Assume K[l]=1, K[l-1]=2, ..., K[i] = 2*K[i+1] or 2*K[i+1]-1.
   P0 = P[0] + p*P[1] + p^2*P[2] + ... + p^K[1]*P[l] < p^K[0]
   The end of the list is P[l+1]=0.
   Assume l > 0.
*/
mpz_poly * mpz_poly_base_modp_init(mpz_poly_srcptr P0, unsigned long p,
                                   unsigned long * K, int l)
{
    return mpz_poly_notparallel_info().mpz_poly_base_modp_init(P0, p, K, l);
}
template <typename inf>
mpz_poly * mpz_poly_parallel_interface<inf>::mpz_poly_base_modp_init(
    mpz_poly_srcptr P0, unsigned long p, unsigned long * K, int l) const
{
    mpz_poly * P;
    mpz_t * pk;

    ASSERT_ALWAYS(l > 0);
    ASSERT_ALWAYS(K[l] == 1);

    /* initialize pk[i] = p^K[l-i] for 0 <= i < l */
    pk = (mpz_t *)malloc(l * sizeof(mpz_t));
    FATAL_ERROR_CHECK(pk == nullptr, "not enough memory");
    mpz_init_set_ui(pk[0], p);
    /* this loop cannot be parallelized, since pk[i] depends on pk[i-1] */
    for (int i = 1; i < l; i++) {
        mpz_init(pk[i]);
        mpz_mul(pk[i], pk[i - 1], pk[i - 1]);
        if (K[l - i] & 1) {
            ASSERT_ALWAYS(K[l - i] == 2 * K[l - i + 1] - 1);
            mpz_div_ui(pk[i], pk[i], p);
        } else
            ASSERT_ALWAYS(K[l - i] == 2 * K[l - i + 1]);
    }

    /* now decompose P0: we need P[0], P[1] for factor p, P[2] for p^2,
       ..., P[l] for p^K[1], and one for the end of list,
       thus l+2 polynomials */
    P = (mpz_poly *)malloc((l + 2) * sizeof(mpz_poly));
    FATAL_ERROR_CHECK(P == nullptr, "not enough memory");
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i < l + 2; i++)
        mpz_poly_init(P[i], P0->deg);
    /* P[l+1] is initialized to 0 by mpz_poly_init */

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= P0->deg; i++)
        mpz_tdiv_qr(P[l]->_coeff[i], P[l - 1]->_coeff[i], P0->_coeff[i],
                    pk[l - 1]);
    mpz_poly_cleandeg(P[l], P0->deg);

    /* now go down */
    for (int j = l - 1; j >= 1; j--) {
        /* reduce P[j] into P[j] and P[j-1] */
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
        for (int i = 0; i <= P0->deg; i++)
            mpz_tdiv_qr(P[j]->_coeff[i], P[j - 1]->_coeff[i], P[j]->_coeff[i],
                        pk[j - 1]);
        mpz_poly_cleandeg(P[j], P0->deg);
    }
    mpz_poly_cleandeg(P[0], P0->deg);

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i < l; i++)
        mpz_clear(pk[i]);
    free(pk);

    return P;
}

/* a <- a + pk*P[k] */
void mpz_poly_base_modp_lift(mpz_poly_ptr a, mpz_poly * P, int k, mpz_srcptr pk)
{
    mpz_poly_notparallel_info().mpz_poly_base_modp_lift(a, P, k, pk);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_base_modp_lift(
    mpz_poly_ptr a, mpz_poly * P, int k, mpz_srcptr pk) const
{
    /* first check P[k] exists and is not zero */
    if (P[k]->deg == -1)
        return;

    mpz_poly_realloc(a, P[k]->deg + 1);

    int const imax = (P[k]->deg < a->deg) ? P[k]->deg : a->deg;
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = 0; i <= imax; i++)
        mpz_addmul(a->_coeff[i], P[k]->_coeff[i], pk);

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
    for (int i = imax + 1; i <= P[k]->deg; i++)
        mpz_mul(a->_coeff[i], P[k]->_coeff[i], pk);

    mpz_poly_cleandeg(a, (a->deg >= P[k]->deg) ? a->deg : P[k]->deg);
}

void mpz_poly_base_modp_clear(mpz_poly * P, int l)
{
    for (int i = 0; i < l + 2; i++)
        mpz_poly_clear(P[i]);
    free(P);
}

/* return the maximal size of the coefficients of f in base b */
size_t mpz_poly_sizeinbase(mpz_poly_srcptr f, int b)
{
    size_t S = 0, s;
    int i;
    int const d = f->deg;

    for (i = 0; i <= d; i++) {
        s = mpz_sizeinbase(f->_coeff[i], b);
        if (s > S)
            S = s;
    }
    return S;
}

/* return the maximal limb-size of the coefficients of f */
size_t mpz_poly_size(mpz_poly_srcptr f)
{
    size_t S = 0, s;
    int i;
    int const d = f->deg;

    for (i = 0; i <= d; i++) {
        s = mpz_size(f->_coeff[i]);
        if (s > S)
            S = s;
    }
    return S;
}

void mpz_poly_infinity_norm(mpz_ptr in, mpz_poly_srcptr f)
{
    if (f->deg == -1) {
        mpz_set_ui(in, 0);
    } else {
        mpz_abs(in, f->_coeff[0]);
        for (int i = 1; i <= f->deg; i++) {
            if (mpz_cmpabs(f->_coeff[i], in) > 0)
                mpz_abs(in, f->_coeff[i]);
        }
    }
}

/* return the total size (in bytes) to store the polynomial f */
size_t mpz_poly_totalsize(mpz_poly_srcptr f)
{
    int i;
    size_t s = 0;

    for (i = 0; i <= f->deg; i++)
        s += mpz_size(f->_coeff[i]);
    return s * sizeof(mp_limb_t);
}

/* f=gcd(f, g) mod p, with p in mpz_t */
/* clobbers g */
/* Coefficients of f and g need not be reduced mod p on input.
 * Coefficients of f are reduced mod p on output */
static void mpz_poly_gcd_mpz_clobber(mpz_poly_ptr f, mpz_poly_ptr g,
                                     mpz_srcptr p)
{
    /* First reduce mod p */
    mpz_poly_mod_mpz(f, f, p, nullptr);
    mpz_poly_mod_mpz(g, g, p, nullptr);
    while (g->deg >= 0) {
        mpz_poly_div_r_mod_mpz_clobber(f, g, p);
        /* now deg(f) < deg(g): swap f and g */
        mpz_poly_swap(f, g);
    }
}

/* f <- gcd(a, b) mod p. */
/* Coefficients of a and b need not be reduced mod p
 * Coefficients of f are reduced mod p */
void mpz_poly_gcd_mpz(mpz_poly_ptr f, mpz_poly_srcptr a, mpz_poly_srcptr b,
                      mpz_srcptr p)
{
    mpz_poly hh;
    if (f == b) {
        mpz_poly_init(hh, a->deg);
        mpz_poly_set(hh, a);
    } else {
        mpz_poly_init(hh, b->deg);
        mpz_poly_set(hh, b);
        mpz_poly_set(f, a); /* will do nothing if f = a */
    }
    mpz_poly_gcd_mpz_clobber(f, hh, p);
    mpz_poly_clear(hh);
}

/* Attempt to compute the f=gcd(f, g) mod N, where N is not necessarily a
 * prime. If at some point a division fails, this gives a proper factor
 * of N that is put in the corresponding argument.
 * The return value tells whether the process was successful (1 means
 * that no inversion failed, 0 means that a factor was found).
 *
 * If a factor is found, and if the parameter "factor" is not NULL, then
 * the encountered factored is stored there.
 *
 * WARNING: this function destroys its input.
 */
/* Coefficients of f and g need not be reduced mod p on input.
 * Coefficients of f are reduced mod p on output */
int mpz_poly_pseudogcd_mpz(mpz_poly_ptr f, mpz_poly_ptr g, mpz_srcptr N,
                           mpz_ptr factor)
{
    mpz_poly_mod_mpz(f, f, N, nullptr);
    mpz_poly_mod_mpz(g, g, N, nullptr);
    while (g->deg >= 0) {
        int const ret = mpz_poly_pseudodiv_r(f, g, N, factor);
        if (!ret)
            return ret;
        /* now deg(f) < deg(g): swap f and g */
        mpz_poly_swap(f, g);
    }
    // success: all inversions mod N worked.
    return 1;
}

/* computes d = gcd(f, g) = u*f + v*g mod p, with p in mpz_t */
/* Coefficients of f and g need not be reduced mod p.
 * Coefficients of d, u, v are reduced mod p
 *
 * Note that this is likely to fail quite miserably if p is not a prime
 * number.
 */
void mpz_poly_xgcd_mpz(mpz_poly_ptr d, mpz_poly_srcptr f, mpz_poly_srcptr g,
                       mpz_poly_ptr u, mpz_poly_ptr v, mpz_srcptr p)
{
    if (f->deg < g->deg) {
        mpz_poly_xgcd_mpz(d, g, f, v, u, p);
        return;
    }

    cxx_mpz_poly u0, v0, r0;
    cxx_mpz_poly u1, v1, r1;

    mpz_poly_set_ui(u0, 1);
    mpz_poly_set_ui(v0, 0);
    mpz_poly_set(r0, f);
    mpz_poly_set_ui(u1, 0);
    mpz_poly_set_ui(v1, 1);
    mpz_poly_set(r1, g);

    mpz_poly_mod_mpz(r0, r0, p, nullptr);
    mpz_poly_mod_mpz(r1, r1, p, nullptr);

    while (r1->deg >= 0) {
        cxx_mpz_poly q, tmp;

        /* q, r0 := r0 div r1 mod p
         * yes, replacing the dividend by the remainder works */
        int const ok = mpz_poly_div_qr_mod_mpz(q, r0, r0, r1, p);

        /* if this fails, then we need a pseudo_xgcd_mpz */
        ASSERT_ALWAYS(ok);
        mpz_poly_swap(r0, r1);

        /* u0 := u0 - q * u1 mod p */
        mpz_poly_mul(tmp, q, u1);
        mpz_poly_sub_mod_mpz(u0, u0, tmp, p);
        mpz_poly_swap(u0, u1);

        /* v0 := v0 - q * v1 mod p */
        mpz_poly_mul(tmp, q, v1);
        mpz_poly_sub_mod_mpz(v0, v0, tmp, p);
        mpz_poly_swap(v0, v1);
    }

    /* make monic */
    if (mpz_cmp_ui(r0->_coeff[r0->deg], 1) != 0) {
        cxx_mpz inv;
        mpz_invert(inv, r0->_coeff[r0->deg], p);
        mpz_poly_mul_mpz(r0, r0, inv);
        mpz_poly_mod_mpz(r0, r0, p, nullptr);
        mpz_poly_mul_mpz(u0, u0, inv);
        mpz_poly_mod_mpz(u0, u0, p, nullptr);
        mpz_poly_mul_mpz(v0, v0, inv);
        mpz_poly_mod_mpz(v0, v0, p, nullptr);
    }

    mpz_poly_swap(u, u0);
    mpz_poly_swap(v, v0);
    mpz_poly_swap(d, r0);
}

/* Homographic transform on polynomials */
/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
*/
cxx_mpz_poly cxx_mpz_poly::homography(std::array<int64_t, 4> const & H) const
{
    cxx_mpz_poly const & F(*this);
    cxx_mpz_poly Fij = F;
    cxx_mpz f0;
    int const d = F->deg;

    mpz_t * fij = Fij->_coeff;

    /* g holds the coefficients of (b0*i+b1)^l */
    std::vector<cxx_mpz> g(d + 1);

    /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
       F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
       We use that formula recursively. */

    mpz_set_ui(g[0], 1); /* g = 1 */

    for (int k = d - 1; k >= 0; k--) {
        /* invariant: we have already translated coefficients of degree > k,
           in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
           g[0..d - (k+1)]:
           f[k] <- a1*f[k+1]
           ...
           f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
           ...
           f[d] <- a0*f[d] */
        mpz_swap(f0, fij[k]); /* save the new constant coefficient */
        mpz_mul_si(fij[k], fij[k + 1], H[2]);
        for (int l = k + 1; l < d; l++) {
            mpz_mul_si(fij[l], fij[l], H[0]);
            mpz_addmul_si(fij[l], fij[l + 1], H[2]);
        }
        mpz_mul_si(fij[d], fij[d], H[0]);

        /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
           g[d-k] = b0*g[d-k-1]
           ...
           g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
           ...
           g[0] = b1*g[0]
           */
        mpz_mul_si(g[d - k], g[d - k - 1], H[1]);
        for (int l = d - k - 1; l > 0; l--) {
            mpz_mul_si(g[l], g[l], H[3]);
            mpz_addmul_si(g[l], g[l - 1], H[1]);
        }
        mpz_mul_si(g[0], g[0], H[3]);

        /* now g has degree d-k, and we add f0*g */
        for (int l = 0; l <= d - k; l++)
            mpz_addmul(fij[l + k], g[l], f0);
    }

    mpz_poly_cleandeg(Fij, Fij->deg);

    return Fij;
}

/* Linear transform on polynomials: return the polynomial f(u*x+v) */
cxx_mpz_poly cxx_mpz_poly::linear_transform(cxx_mpz const & u,
                                            cxx_mpz const & v) const
{
    cxx_mpz_poly L{v, u}; /* L = (ux + v) */
    cxx_mpz_poly Fij;
    int const d = degree();

    Fij = coeff(d);

    for (int k = d-1; k >= 0; --k) {
        Fij *= L;
        Fij += coeff(k);
    }
    return Fij;
}

/* v <- f(i,j), where f is homogeneous of degree d */
void mpz_poly_homogeneous_eval_siui(mpz_ptr v, mpz_poly_srcptr f,
                                    int64_t const i, uint64_t const j)
{
    unsigned int k = f->deg;
    ASSERT(k > 0);
    mpz_set(v, f->_coeff[f->deg]);
    mpz_mul_si(v, f->_coeff[k], i);
    cxx_mpz jpow;
    mpz_set_uint64(jpow, j);
    mpz_addmul(v, f->_coeff[--k], jpow); /* v = i*f[d] + j*f[d-1] */
    for (; k-- > 0;) {
        /* this test will be resolved at compile time by most compilers */
        if ((uint64_t)ULONG_MAX >=
            UINT64_MAX) { /* hardcode since this function is critical in las */
            mpz_mul_si(v, v, i);
            mpz_mul_ui(jpow, jpow, j);
        } else {
            mpz_mul_int64(v, v, i);
            mpz_mul_uint64(jpow, jpow, j);
        }
        mpz_addmul(v, f->_coeff[k], jpow);
    }
}

/* put in c the content of f */
void mpz_poly_content(mpz_ptr c, mpz_poly_srcptr F)
{
    int i;
    mpz_t * f = F->_coeff;
    int const d = F->deg;

    if (d == -1) {
        mpz_set_ui(c, 0);
        return;
    }

    mpz_set(c, f[0]);
    for (i = 1; i <= d; i++)
        mpz_gcd(c, c, f[i]);
    mpz_abs(c, c);
}

int mpz_poly_has_trivial_content(mpz_poly_srcptr F)
{
    int i;
    mpz_t * f = F->_coeff;
    int const d = F->deg;
    mpz_t c;
    mpz_init_set(c, f[0]);
    mpz_abs(c, c);
    for (i = 1; i <= d && mpz_cmp_ui(c, 1) > 0; i++) {
        mpz_gcd(c, c, f[i]);
        mpz_abs(c, c);
    }
    int const res = mpz_cmp_ui(c, 1) == 0;
    mpz_clear(c);
    return res;
}

/* return non-zero if the polynomial has non-trivial content */
int mpz_poly_divide_by_content(mpz_poly_ptr F)
{
    cxx_mpz c;
    mpz_poly_content(c, F);
    if (mpz_cmp_ui(c, 1) == 0)
        return 0;
    for (int i = 0; i <= F->deg; i++)
        mpz_fdiv_q(F->_coeff[i], F->_coeff[i], c);
    return 1;
}

/*
 * Compute the pseudo division of a and b such that
 *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * b must be non-zero
 * set q=NULL if the quotient is not needed
 * a==r is supported
 */
void mpz_poly_pseudo_division(mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr a,
                              mpz_poly_srcptr b)
{
    ASSERT_ALWAYS(b->deg != -1);

    if (q)
        mpz_poly_set_zero(q);
    mpz_poly_set(r, a);

    if (a->deg < b->deg)
        return;

    int const m = a->deg;
    int const n = b->deg;
    int e = m - n + 1;

    cxx_mpz d = mpz_poly_lc(b);

    while (r->deg >= n) {
        cxx_mpz_poly s;
        mpz_poly_setcoeff(s, r->deg - n, mpz_poly_lc(r));
        if (q) {
            mpz_poly_mul_mpz(q, q, d);
            mpz_poly_add(q, q, s);
        }

        mpz_poly_mul_mpz(r, r, d);
        mpz_poly_mul(s, b, s);
        mpz_poly_sub(r, r, s);
        e--;
    }

    ASSERT(e >= 0);

    mpz_pow_ui(d, d, (unsigned long int)e);

    if (q)
        mpz_poly_mul_mpz(q, q, d);
    mpz_poly_mul_mpz(r, r, d);
}

/*
 * Like mpz_poly_pseudo_division, but give only the remainder.
 *
 * a==r is supported
 */
void mpz_poly_pseudo_remainder(mpz_poly_ptr r, mpz_poly_srcptr a,
                               mpz_poly_srcptr b)
{
    mpz_poly_pseudo_division(nullptr, r, a, b);
}

/*
 * Compute the resultant of p and q and put the result in res.
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that the polynomials are normalized.
 */
void mpz_poly_resultant(mpz_ptr res, mpz_poly_srcptr p, mpz_poly_srcptr q)
{
    if (p->deg == -1 || q->deg == -1) {
        mpz_set_ui(res, 0);
        return;
    }

    ASSERT(mpz_cmp_ui(p->_coeff[p->deg], 0) != 0);
    ASSERT(mpz_cmp_ui(q->_coeff[q->deg], 0) != 0);

    long int s = 1;
    mpz_t g;
    mpz_t h;
    mpz_t t;
    mpz_t tmp;
    int d;
    mpz_poly r;
    mpz_poly a;
    mpz_poly b;

    mpz_init(g);
    mpz_init(h);
    mpz_init(t);
    mpz_init(tmp);
    mpz_poly_init(r, -1);
    mpz_poly_init(a, p->deg);
    mpz_poly_init(b, q->deg);

    mpz_poly_set(a, p);
    mpz_poly_set(b, q);

    mpz_poly_content(g, a);
    mpz_poly_content(h, b);

    mpz_poly_divexact_mpz(a, a, g);
    mpz_poly_divexact_mpz(b, b, h);

#ifndef NDEBUG
    mpz_poly_content(tmp, a);
    ASSERT(mpz_cmp_ui(tmp, 1) == 0);
    mpz_poly_content(tmp, b);
    ASSERT(mpz_cmp_ui(tmp, 1) == 0);
#endif // NDEBUG

    mpz_pow_ui(t, g, (unsigned long int)b->deg);
    mpz_pow_ui(tmp, h, (unsigned long int)a->deg);
    mpz_mul(t, t, tmp);

    mpz_set_ui(g, 1);
    mpz_set_ui(h, 1);

    if (a->deg < b->deg) {
        mpz_poly_swap(a, b);

        if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
            s = -1;
        }
    }

    while (b->deg > 0) {
        d = a->deg - b->deg;

        if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
            s = -s;
        }

        mpz_poly_pseudo_remainder(r, a, b);
        mpz_poly_set(a, b);

        ASSERT(d >= 0);

        mpz_pow_ui(tmp, h, (unsigned long int)d);
        mpz_mul(tmp, g, tmp);

        mpz_poly_divexact_mpz(b, r, tmp);

        mpz_set(g, mpz_poly_lc(a));

#ifdef NDEBUG
        if (d == 0) {
            ASSERT(mpz_cmp_ui(h, 1) == 0);
        }
#endif // NDEBUG
        mpz_pow_ui(h, h, (unsigned long int)(d - 1));
        mpz_pow_ui(tmp, g, (unsigned long int)d);
        mpz_divexact(h, tmp, h);
    }

    // Prevent an error if b = 0.
    if (b->deg == -1) {
        mpz_set_ui(res, 0);
    } else {
        ASSERT(a->deg > 0);
        ASSERT(b->deg == 0);

        mpz_pow_ui(h, h, (unsigned long int)(a->deg - 1));

        ASSERT(a->deg >= 0);

        mpz_pow_ui(tmp, b->_coeff[0], (unsigned long int)a->deg);
        mpz_divexact(h, tmp, h);

        mpz_mul_si(t, t, s);
        mpz_mul(h, h, t);
        mpz_set(res, h);
    }

    mpz_clear(g);
    mpz_clear(h);
    mpz_clear(t);
    mpz_clear(tmp);
    mpz_poly_clear(a);
    mpz_poly_clear(b);
    mpz_poly_clear(r);
}

void mpz_poly_discriminant(mpz_ptr res, mpz_poly_srcptr f)
{
    mpz_poly df;
    mpz_poly_init(df, f->deg);
    mpz_poly_derivative(df, f);
    mpz_poly_resultant(res, f, df);
    ASSERT(mpz_divisible_p(res, mpz_poly_lc(f)));
    mpz_divexact(res, res, mpz_poly_lc(f));
    mpz_poly_clear(df);
}

/* returns non-zero iff f is square-free in Z[x] */
int mpz_poly_squarefree_p(mpz_poly_srcptr f)
{
    mpz_poly df;
    mpz_t res;
    int ret;

    mpz_poly_init(df, f->deg);
    mpz_poly_derivative(df, f);
    mpz_init(res);
    mpz_poly_resultant(res, f, df);
    ret = mpz_cmp_ui(res, 0);
    mpz_clear(res);
    mpz_poly_clear(df);
    return ret;
}

/**
 * Quick-and-dirty test if the polynomial f is irreducible over Z.
 * This function is extracted from dlpolyselect.c written by PZ
 *
 * Let p be a prime such that f has a root r modulo p.
 * Search by LLL a small linear combination between 1, r, ..., r^(d-1).
 * If f is not irreducible, then p will be root of a factor of degree <= d-1,
 * which will yield a linear dependency of the same order as the coefficients
 * of f, otherwise if f is irreducible the linear dependency will be of order
 * p^(1/d).
 *
 * \return 0 by default, 1 if the poly is irreducible.
 *
 * this function is without any warranty and at your own risk
 */
int mpz_poly_is_irreducible_z(mpz_poly_srcptr f)
{
    mpz_t p;
    int const d = f->deg;
    ASSERT_ALWAYS(d >= -1 && d < INT_MAX);
    int const dg = d - 1;
    ASSERT_ALWAYS(dg + 2 > 0);
    int i, j, nr;
    mpz_t a, b, det, r, *roots;
    size_t normf;
    int ret = 0; // init
    mat_Z g;
    gmp_randstate_t rstate;

    mpz_init(p);
    gmp_randinit_default(rstate);

    mpz_poly_infinity_norm(p, f);
    normf = mpz_sizeinbase(p, 2);
    /* The following table might be useful to optimize the value of MARGIN:
       degree d    infinity_norm(f)    optimal MARGIN
          3               8                 5
          4               4                 5
    */
#define MARGIN 16
    /* add some margin bits */
    mpz_mul_2exp(p, p, MARGIN);
    mpz_pow_ui(p, p, d);

    roots = (mpz_t *)malloc(d * sizeof(mpz_t));
    for (i = 0; i < d; i++)
        mpz_init(roots[i]);

    do {
        mpz_nextprime(p, p);
        nr = mpz_poly_roots_mpz(roots, f, p, rstate);
        /* If f has no root mod p and degree <= 3, it is irreducible,
           since a degree 2 polynomial can only factor into 1+1 or 2,
           and a degree 3 polynomial can only factor into 1+2 or 3. */
        if (nr == 0 && d <= 3) {
            ret = 1;
            goto clear_and_exit;
        }
    } while (nr == 0);

    g.coeff = (mpz_t **)malloc((dg + 2) * sizeof(mpz_t *));
    g.NumRows = g.NumCols = dg + 1;
    for (i = 0; i <= dg + 1; i++) {
        g.coeff[i] = (mpz_t *)malloc((dg + 2) * sizeof(mpz_t));
        for (j = 0; j <= dg + 1; j++) {
            mpz_init(g.coeff[i][j]);
        }
    }

    mpz_init(det);
    mpz_init_set_ui(a, 1);
    mpz_init_set_ui(b, 1);
    mpz_init_set(r, roots[0]);
    for (i = 0; i <= dg + 1; i++) {
        for (j = 0; j <= dg + 1; j++) {
            mpz_set_ui(g.coeff[i][j], 0);
        }
    }

    for (i = 1; i <= dg + 1; i++) {
        for (j = 1; j <= dg + 1; j++) {
            if (i == 1) {
                if (j == 1) {
                    mpz_set(g.coeff[j][i], p);
                } else {
                    mpz_neg(g.coeff[j][i], r);
                    mpz_mul(r, r, roots[0]);
                }
            } else
                mpz_set_ui(g.coeff[j][i], i == j);
        }
    }

    LLL(det, g, nullptr, a, b);

    for (j = 1; j <= dg + 1; j++) {
        /* the coefficients of vector j are in g.coeff[j][i], 1 <= i <= dg + 1
         */
        mpz_abs(a, g.coeff[j][1]);
        for (i = 2; i <= dg + 1; i++)
            if (mpz_cmpabs(g.coeff[j][i], a) > 0)
                mpz_abs(a, g.coeff[j][i]);
        /* now a = max (|g.coeff[j][i]|, 1 <= i <= dg+1) */
        if (j == 1 || mpz_cmpabs(a, b) < 0)
            mpz_set(b, a);
    }
    /* now b is the smallest infinity norm */
    if (mpz_sizeinbase(b, 2) < normf + MARGIN / 2)
        ret = 0;
    else
        ret = 1;

    for (i = 0; i <= dg + 1; i++) {
        for (j = 0; j <= dg + 1; j++)
            mpz_clear(g.coeff[i][j]);
        free(g.coeff[i]);
    }
    free(g.coeff);

    mpz_clear(det);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(r);
clear_and_exit:
    for (i = 0; i < d; i++)
        mpz_clear(roots[i]);
    free(roots);

    gmp_randclear(rstate);
    mpz_clear(p);

    return ret;
}

/* factoring polynomials */

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l)
{
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l)
{
    for (int i = 0; i < l->alloc; i++)
        mpz_poly_clear(l->factors[i]->f);
    free(l->factors);
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l)
{
    /* There's a design choice here. We may elect to free everything.
     * Instead, we'll simply mark everything as zero, but keep all
     * allocated space.
     */
    for (int i = 0; i < l->alloc; i++)
        l->factors[i]->f->deg = -1;
    l->size = 0;
}

static void mpz_poly_factor_list_prepare_write(mpz_poly_factor_list_ptr l,
                                               int index)
{
    if (index >= l->alloc) {
        l->alloc = index + 4 + l->alloc / 4;
        checked_realloc(l->factors, l->alloc);
        /* We need to set something. A zero polynomial has NULL storage
         * area, so that will do (realloc behaves as needed).  */
        for (int i = l->size; i < l->alloc; i++) {
            mpz_poly_init(l->factors[i]->f, -1);
            l->factors[i]->m = 0;
        }
    }
    if (l->size <= index)
        l->size = index + 1;
}

void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly_srcptr f,
                               int m)
{
    mpz_poly_factor_list_prepare_write(l, l->size);
    mpz_poly_set(l->factors[l->size - 1]->f, f);
    l->factors[l->size - 1]->m = m;
}

void mpz_poly_factor_list_fprintf(FILE * fp, mpz_poly_factor_list_srcptr l)
{
    for (int i = 0; i < l->size; i++) {
        char * res;
        int const rc = mpz_poly_asprintf(&res, l->factors[i]->f);
        ASSERT_ALWAYS(rc >= 0);
        if (i)
            fprintf(fp, "*");
        fprintf(fp, "(%s)^%d", res, l->factors[i]->m);
        free(res);
    }
    fprintf(fp, "\n");
}

void mpz_poly_factor_list_accumulate(mpz_poly_ptr f,
                                     mpz_poly_factor_list_srcptr l)
{
    cxx_mpz_poly T;
    mpz_poly_set_ui(f, 1);
    for (int i = 0; i < l->size; i++) {
        mpz_poly_pow_ui(T, l->factors[i]->f, l->factors[i]->m);
        mpz_poly_mul(f, f, T);
    }
}
/* Squarefree factorization */

/* This auxiliary function almost does the sqf. It fills
 * lf->factors[stride*i] (i from 1 to deg(f)) with the factors with
 * multiplicity i in f.  lf->factors[0] is filled with the product whose
 * multiplicity is a multiple of the field characteristic.  returns max
 * multiplicity stored (multiplied by the stride value).
 *
 * (stride has an importance for the recursive call, for p small. E.g. on
 * f=g^p, we get called with g=f^(1/p) and stride=p).
 *
 * We make no effort to check that lf is clean on input, which is to be
 * guaranteed by the caller (e.g. sufficiently many polynomials, all
 * equal to 1 -- or 0, which can easily be identified as something
 * unset).
 *
 * Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_sqf_inner(mpz_poly_factor_list_ptr lf,
                                     mpz_poly_srcptr f, int stride,
                                     mpz_srcptr p)
{
    int r = 0;

    mpz_poly g, mi, mi1;
    mpz_poly t0, t1, T, tmp;
    mpz_poly_init(g, f->deg);
    mpz_poly_init(mi, f->deg);
    mpz_poly_init(mi1, f->deg);
    mpz_poly_init(t0, f->deg);
    mpz_poly_init(t1, f->deg);
    mpz_poly_init(T, f->deg);
    mpz_poly_init(tmp, f->deg);

    mpz_poly_derivative(t0, f);
    mpz_poly_gcd_mpz(g, f, t0, p);
    mpz_poly_divexact(mi, f, g, p);
    /* mi is f/gcd(f,f') == all repeated prime factors of f whose
     * multiplicity isn't a multiple of the field characteristic.
     */

    mpz_poly_set_xi(T, 0);
    for (int i = 1; mi->deg > 0; i++) {
        /* note the weird argument ordering */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, mi, nullptr, i, p);
        mpz_poly_divexact(t1, f, t0, p);
        /* t1 = all polynomials in mi taken out from f with multiplicity i */
        mpz_poly_gcd_mpz(mi1, t1, mi, p);
        /* mi1 = almost like mi, but since factors with multiplicity i
         * are no longer in t1, there's absent from mi1 too. Whence
         * mi/mi1 is exactly the product of factors of multiplicity 1.
         */
        mpz_poly_factor_list_prepare_write(lf, i * stride);
        /* Use tmp so that we don't absurdly keep storage within
         * lf->factors */
        mpz_poly_divexact(tmp, mi, mi1, p);
        mpz_poly_set(lf->factors[i * stride]->f, tmp);
        /* multiplicity field still unused at this point */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, lf->factors[i * stride]->f, nullptr,
                                      i, p);
        mpz_poly_mul(T, T, t0);
        mpz_poly_mod_mpz(T, T, p, nullptr);
        mpz_poly_swap(mi, mi1);
        r = i * stride;
    }

    mpz_poly_factor_list_prepare_write(lf, 0);
    mpz_poly_divexact(lf->factors[0]->f, f, T, p);

    mpz_poly_clear(g);
    mpz_poly_clear(tmp);
    mpz_poly_clear(mi);
    mpz_poly_clear(mi1);
    mpz_poly_clear(t0);
    mpz_poly_clear(t1);
    mpz_poly_clear(T);
    return r;
}

/* Fills lf->factors[i] (i from 1 to deg(f)) with product of the factors
 * with multiplicity i in f, and to 1 if there are none. lf->factors[0]
 * is set to 1.  return the largest multiplicity stored. */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0,
                        mpz_srcptr p)
{
    /* factoring 0 doesn't make sense, really */
    ASSERT(f0->deg >= 0);

    ASSERT_ALWAYS(mpz_cmp_ui(p, 0) > 0);

    /* We'll call mpz_poly_factor_sqf_inner, possibly several times if
     * we are in small characteristic.
     */
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    mpz_poly_makemonic_mod_mpz(f, f0, p);
    ASSERT(mpz_cmp_ui(mpz_poly_lc(f), 1) == 0);

    int m = 0;
    // coverity[zero_return]
    int const pu = mpz_get_ui(p); // see below
    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    for (int stride = 1;; stride *= pu) {
        int const r = mpz_poly_factor_sqf_inner(lf, f, stride, p);
        if (r > m)
            m = r;
        if (lf->factors[0]->f->deg == 0) {
            // if p is LAAAARGE, then of course we'll never have a linear
            // polynomial out of sqf_inner, thus we'll break early here.
            break;
        }
        /* divide coefficients */
        for (int i = 0; i <= lf->factors[0]->f->deg; i++) {
            if (i % pu == 0) {
                mpz_set(f->_coeff[i / pu], lf->factors[0]->f->_coeff[i]);
            } else {
                ASSERT(mpz_cmp_ui(lf->factors[0]->f->_coeff[i], 0) == 0);
            }
        }
        f->deg = lf->factors[0]->f->deg / pu;
        mpz_poly_set_xi(lf->factors[0]->f, 0);
    }
    /* Now make sure that all factors in the factor list are non-zero */
    for (int i = 0; i < lf->size; i++) {
        if (lf->factors[i]->f->deg < 0) {
            mpz_poly_set_xi(lf->factors[i]->f, 0);
        }
    }
    mpz_poly_clear(f);
    return m;
}

/* This performs distinct degree factorization */
/* Input polynomial must be squarefree -- otherwise repeated factors
 * probably won't show up in the factor list, or maybe at the wrong place
 * as parasites. */
/* Returns max degree of factors found (i.e. if the largest factors we
 * have are two factors of degree 7, we return 7, not 14). */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_ddf_inner(mpz_poly_factor_list_ptr lf,
                                     mpz_poly_srcptr f0, mpz_srcptr p,
                                     int only_check_irreducible)
{
    mpz_poly g, gmx, x, tmp;
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    int i;

    /* factoring 0 doesn't make sense, really */
    ASSERT(f0->deg >= 0);

    mpz_poly_makemonic_mod_mpz(f, f0, p);
    ASSERT(mpz_cmp_ui(mpz_poly_lc(f), 1) == 0);

    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    mpz_poly_init(g, 2 * f->deg - 1);
    mpz_poly_init(gmx, 2 * f->deg - 1);
    mpz_poly_init(x, 1);
    mpz_poly_init(tmp, f->deg);

    mpz_poly_set_xi(x, 1);
    mpz_poly_set_xi(g, 1);

    for (i = 1; i <= f->deg; ++i) {
        if (2 * i > f->deg) {
            /* Then we know that the remaining f is irreducible.  */
            mpz_poly_factor_list_prepare_write(lf, f->deg);
            for (; i < f->deg; i++) {
                /* multiplicity field still unused at this point */
                mpz_poly_set_xi(lf->factors[i]->f, 0);
            }
            /* multiplicity field still unused at this point */
            mpz_poly_swap(lf->factors[f->deg]->f, f);
            break;
        }

        /* g <- g^p mod fp */
        mpz_poly_pow_mod_f_mod_mpz(g, g, f, p, p);

        /* subtract x */
        mpz_poly_sub(gmx, g, x);
        mpz_poly_mod_mpz(gmx, gmx, p, nullptr);

        /* lf[i] <- gcd (f, x^(p^i)-x) */
        mpz_poly_factor_list_prepare_write(lf, i);
        /* multiplicity field still unused at this point */

        /* see remark in _sqf regarding the relevance of tmp for storage */
        mpz_poly_gcd_mpz(tmp, f, gmx, p);
        mpz_poly_divexact(f, f, tmp, p);
        mpz_poly_set(lf->factors[i]->f, tmp);

        /* Note for a mere irreducibility test: the length of the loop in
         * the irreducible case would still be deg(f)/2, and the penalty
         * caused by storing factors can be neglected.
         */
        if (only_check_irreducible && lf->factors[i]->f->deg > 0)
            break;

        if (f->deg == 0)
            break;
    }

    mpz_poly_clear(g);
    mpz_poly_clear(x);
    mpz_poly_clear(gmx);
    mpz_poly_clear(tmp);
    mpz_poly_clear(f);

    return i;
}

int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0,
                        mpz_srcptr p)
{
    return mpz_poly_factor_ddf_inner(lf, f0, p, 0);
}

/* Note that this also works for non squarefree polynomials -- the factor
 * list returned by mpz_poly_factor_ddf will be rubbish, but the m ==
 * f->deg test will tell the truth. */
int mpz_poly_is_irreducible(mpz_poly_srcptr f, mpz_srcptr p)
{
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    int const m = mpz_poly_factor_ddf_inner(lf, f, p, 1);
    mpz_poly_factor_list_clear(lf);
    return m == f->deg;
}

/* Naive factorization of polynomials over GF(2). We assume we're always
 * dealing with polynomials of degree at most 10, so we're beter off
 * simply enumerating potential factors...
 */

/*
 * Add 1 to f. If the constant term is equal to 1, set this term to 0 and
 *  propagate the addition of 1 to the next coefficient, and so on.
 *
 * f: the polynomial on which the addition is computed, the modifications are
 *  made on f.
 */
static void mpz_poly_add_one_in_F2(mpz_poly_ptr f)
{
    ASSERT(f->deg >= 1);

    int i = 0;
    while (mpz_cmp_ui(f->_coeff[i], 1) == 0) {
        mpz_poly_setcoeff_si(f, i, 0);
        i++;
        if (i > f->deg) {
            break;
        }
    }
    mpz_poly_setcoeff_si(f, i, 1);
}

/*
 * Factorize naively a mpz_poly mod 2.
 *
 * Return the number of factors found.
 */
static int mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly_srcptr f,
                            mpz_srcptr p)
{
    ASSERT(mpz_cmp_ui(p, 2) == 0);

    // make a copy of f.
    mpz_poly fcopy;
    mpz_poly_init(fcopy, f->deg);
    mpz_poly_set(fcopy, f);

    // reduce all the coefficient mod p.
    mpz_t coeff;
    mpz_init(coeff);
    for (int i = 0; i <= f->deg; i++) {
        mpz_set(coeff, mpz_poly_coeff_const(f, i));
        mpz_mod(coeff, coeff, p);
        mpz_poly_setcoeff(fcopy, i, coeff);
    }

    // Purge list.
    mpz_poly_factor_list_flush(list);

    // If deg(f) in F2 is less than 1, we have the factor.
    if (fcopy->deg < 1) {
        mpz_clear(coeff);
        mpz_poly_clear(fcopy);
        ASSERT(list->size == 0);
        return list->size;
    } else if (fcopy->deg == 1) {
        mpz_clear(coeff);
        mpz_poly_factor_list_push(list, fcopy, 1);
        mpz_poly_clear(fcopy);
        ASSERT(list->size == 1);
        return list->size;
    }

    if (mpz_poly_is_irreducible(f, p)) {
        // If f is irreducible mod 2, fcopy is the factorisation of f mod 2.
        mpz_poly_factor_list_push(list, fcopy, 1);
    } else {
        // Create the first possible factor.
        mpz_poly tmp;
        mpz_poly_init(tmp, 1);
        mpz_poly_setcoeff_int64(tmp, 1, 1);

        // Enumerate all the possible factor of f mod 2.
        while (tmp->deg <= fcopy->deg) {
            // tmp is a possible factor.
            if (mpz_poly_is_irreducible(tmp, p)) {
                mpz_poly q;
                mpz_poly_init(q, 0);
                mpz_poly r;
                mpz_poly_init(r, 0);
                // Euclidean division of fcopy
                mpz_poly_div_qr_mod_mpz(q, r, fcopy, tmp, p);
                // Power of the possible factor.
                unsigned int m = 0;
                // While fcopy is divisible by tmp.
                while (r->deg == -1) {
                    // Increase the power of tmp.
                    m++;
                    mpz_poly_set(fcopy, q);
                    if (fcopy->deg == 0 || fcopy->deg == -1) {
                        // No other possible factor.
                        break;
                    }
                    mpz_poly_div_qr_mod_mpz(q, r, fcopy, tmp, p);
                }
                if (m != 0) {
                    // Push tmp^m as a factor of f mod 2.
                    mpz_poly_factor_list_push(list, tmp, m);
                }
                mpz_poly_clear(q);
                mpz_poly_clear(r);
            }
            // Go to the next possible polynomial in F2.
            mpz_poly_add_one_in_F2(tmp);
        }
        mpz_poly_clear(tmp);
    }

#ifndef NDBEBUG
    // Verify if the factorisation is good.
    mpz_poly_mod_mpz(fcopy, f, p, nullptr);

    cxx_mpz_poly fmul;
    mpz_poly_factor_list_accumulate(fmul, list);
    mpz_poly_mod_mpz(fmul, fmul, p, nullptr);

    ASSERT(mpz_poly_cmp(fcopy, fmul) == 0);

#endif // NDBEBUG

    mpz_clear(coeff);
    mpz_poly_clear(fcopy);

    return list->size;
}

/*
 * this tries to split between squares and non-squares -- it's the most
 * basic split of course, and the other splits merely randomize on top of
 * this. Note that this building block must be changed for characteristic
 * two
 *
 * returns non-zero if a non-trivial or split is obtained.
 *
 * k is the degree of factors we are looking for.
 *
 * We split in two parts:
 *
 *  - factors whose roots are squares in GF(p^k).
 *  - factors whose roots are non squares in GF(p^k).
 *
 * Note that we do not find the factor X this way; this is to be done by
 * the caller.
 *
 * Obviously, we include some shift, and hope that eventually there is a
 * shift that works. Large characteristic is generally happy with some
 * translation shift. Small characteristic may need more general shifts.
 *
 * Coefficients of f0 need not be reduced mod p.
 * Coefficients of g[0] and g[1] are reduced mod p.
 */
static void mpz_poly_factor_edf_pre(mpz_poly g[2], mpz_poly_srcptr f, int k,
                                    mpz_srcptr p)
{
    int nontrivial = 0;
    mpz_poly_set_xi(g[0], 0);
    mpz_poly_set_xi(g[1], 0);

    ASSERT_ALWAYS(mpz_cmp_ui(p, 0) > 0);

    ASSERT_ALWAYS(f->deg > k);

    mpz_poly xplusa;
    mpz_poly_init(xplusa, 1);

    mpz_t half_pk;
    mpz_init(half_pk);
    mpz_pow_ui(half_pk, p, k);
    mpz_fdiv_q_ui(half_pk, half_pk, 2); /* (p^k-1)/2 */

    for (unsigned long a = 0; !nontrivial; a++) {
        /* we want to iterate on monic polynomials of degree <= k-1. */
        /* In order to bear in mind what happens in large enough
         * characteristic, we'll name these polynomials xplusa, although
         * it does not have to be x+a (and it can't be restricted to only
         * x+a if the characteristic is small -- that does not give
         * enough legroom).
         */
        if (mpz_fits_ulong_p(p)) {
            // coverity[zero_return]
            unsigned long const pz = mpz_get_ui(p);
            if (a == 0) {
                /* special case, really */
                mpz_poly_set_xi(xplusa, 1);
            } else {
                /* write a in base p, and add 1 */
                int i = 0;
                for (unsigned long tmp = a; tmp; i++, tmp /= pz) {
                    mpz_poly_setcoeff_ui(xplusa, i, tmp % pz);
                }
                mpz_poly_setcoeff_ui(xplusa, i, 1);
            }
        } else {
            /* take the polynomial x+a */
            mpz_poly_set_xi(xplusa, 1);
            mpz_poly_setcoeff_ui(xplusa, 0, a);
        }

        mpz_poly_pow_mod_f_mod_mpz(g[0], xplusa, f, half_pk, p);

        mpz_poly_add_ui(g[1], g[0], 1); /* (x+a)^((p^k-1)/2) + 1 */
        mpz_poly_sub_ui(g[0], g[0], 1); /* (x+a)^((p^k-1)/2) - 1 */

        mpz_poly_mod_mpz(g[0], g[0], p, nullptr);
        mpz_poly_mod_mpz(g[1], g[1], p, nullptr);

        mpz_poly_gcd_mpz(g[0], g[0], f, p);
        mpz_poly_gcd_mpz(g[1], g[1], f, p);

        if (g[0]->deg + g[1]->deg < f->deg) {
            /* oh, we're lucky. x+a is a factor ! */
            int const s = g[0]->deg > g[1]->deg;
            /* multiply g[s] by (x+a) */
            mpz_poly_mul_mod_f_mod_mpz(g[s], g[s], xplusa, f, p, nullptr,
                                       nullptr);
        }
        ASSERT(g[0]->deg + g[1]->deg == f->deg);
        ASSERT(g[0]->deg % k == 0);
        ASSERT(g[1]->deg % k == 0);

        nontrivial += g[0]->deg != 0 && g[0]->deg != f->deg;
        nontrivial += g[1]->deg != 0 && g[1]->deg != f->deg;
    }
    mpz_clear(half_pk);
    mpz_poly_clear(xplusa);
}

/* This factors f, and for each factor q found, store q in lf.
 * Return the number of distinct factors found. */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_edf_inner(mpz_poly_factor_list_ptr lf,
                                     mpz_poly_srcptr f, int k, mpz_srcptr p,
                                     gmp_randstate_t rstate)
{
    if (f->deg == k) {
        mpz_poly_factor_list_push(lf, f, 1);
        return 1;
    }
    if (f->deg == 0) {
        return 0;
    }

    mpz_poly h[2];

    mpz_poly_init(h[0], f->deg);
    mpz_poly_init(h[1], f->deg);

    mpz_poly_factor_edf_pre(h, f, k, p);

    int n = 0;

    n += mpz_poly_factor_edf_inner(lf, h[0], k, p, rstate);
    mpz_poly_clear(h[0]);

    n += mpz_poly_factor_edf_inner(lf, h[1], k, p, rstate);
    mpz_poly_clear(h[1]);

    return n;
}

// returns f0->deg / d
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int k,
                        mpz_srcptr p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list_flush(lf);

    if (mpz_cmp_ui(p, 2) == 0) {
        /* we need some other code for edf. Currently we have very naive
         * code, but that's good enough for small degree. */
        return mpz_poly_factor2(lf, f, p);
    }

    int const v = mpz_poly_valuation(f);
    if (v) {
        /* Since our input is square-free, then we expect v==1.
         * Furthermore, k prescribes the extension field where the
         * expected roots live, thus for 0 it should really be 1. */
        ASSERT_ALWAYS(v == 1 && k == 1);
        mpz_poly_factor_list_prepare_write(lf, lf->size);
        mpz_poly_set_xi(lf->factors[lf->size - 1]->f, 1);

        mpz_poly f1;
        mpz_poly_init(f1, f->deg - 1);
        mpz_poly_div_xi(f1, f, v);
        int const n = 1 + mpz_poly_factor_edf_inner(lf, f1, k, p, rstate);
        mpz_poly_clear(f1);
        return n;
    }

    int const ret = mpz_poly_factor_edf_inner(lf, f, k, p, rstate);
    return ret;
}

typedef int (*sortfunc_t)(void const *, void const *);

static int mpz_poly_with_m_cmp(mpz_poly_with_m const * a,
                               mpz_poly_with_m const * b)
{
    int const r = mpz_poly_cmp((*a)->f, (*b)->f);
    if (r)
        return r;
    return ((*a)->m > (*b)->m) - ((*b)->m > (*a)->m);
}

/* putting it all together */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly_srcptr f, mpz_srcptr p,
                    gmp_randstate_t rstate)
{
    mpz_poly_factor_list sqfs, ddfs, edfs;
    mpz_poly_factor_list_init(sqfs);
    mpz_poly_factor_list_init(ddfs);
    mpz_poly_factor_list_init(edfs);

    mpz_poly_factor_list_flush(lf);

    int const maxmult = mpz_poly_factor_sqf(sqfs, f, p);
    for (int m = 1; m <= maxmult; m++) {
        ASSERT_ALWAYS(sqfs->factors[m]->f->deg >= 0);
        if (sqfs->factors[m]->f->deg == 0)
            continue;
        int const maxdeg = mpz_poly_factor_ddf(ddfs, sqfs->factors[m]->f, p);
        for (int k = 1; k <= maxdeg; k++) {
            ASSERT_ALWAYS(ddfs->factors[k]->f->deg >= 0);
            if (ddfs->factors[k]->f->deg == 0)
                continue;
            mpz_poly_factor_edf(edfs, ddfs->factors[k]->f, k, p, rstate);
            for (int j = 0; j < edfs->size; j++) {
                /* register this factor, with multiplicity m */
                /* cheat a bit... */
                mpz_poly_factor_list_prepare_write(lf, lf->size);
                mpz_poly_swap(lf->factors[lf->size - 1]->f,
                              edfs->factors[j]->f);
                mpz_poly_makemonic_mod_mpz(lf->factors[lf->size - 1]->f,
                                           lf->factors[lf->size - 1]->f, p);
                lf->factors[lf->size - 1]->m = m;
            }
        }
    }
    mpz_poly_factor_list_clear(edfs);
    mpz_poly_factor_list_clear(ddfs);
    mpz_poly_factor_list_clear(sqfs);

    /* sort factors by degree and lexicographically */
    qsort(lf->factors, lf->size, sizeof(mpz_poly_with_m),
          (sortfunc_t)&mpz_poly_with_m_cmp);

    return lf->size;
}

std::vector<std::pair<cxx_mpz_poly, int>>
mpz_poly_factor(mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    mpz_poly_factor(lf, f, p, rstate);
    std::vector<std::pair<cxx_mpz_poly, int>> res;
    res.reserve(lf->size);
    for (int i = 0; i < lf->size; i++) {
        res.emplace_back(lf->factors[i]->f, lf->factors[i]->m);
    }
    mpz_poly_factor_list_clear(lf);
    return res;
}

static void
mpz_poly_factor_list_set(mpz_poly_factor_list_ptr lf,
                         std::vector<std::pair<cxx_mpz_poly, int>> const & xlf)
{
    mpz_poly_factor_list_flush(lf);
    for (auto const & fm: xlf)
        mpz_poly_factor_list_push(lf, fm.first, fm.second);
}

int mpz_poly_number_of_real_roots(mpz_poly_srcptr f)
{
    /* This is coded in usp.c, with an interface pretty different from
     * what we have here.
     *
     * 0.0 in usp.c means: find a bound by yourself.
     */
    return mpz_poly_number_of_real_roots_extra(f, 0.0, nullptr);
}

/* return the product of P lf[a:b].
 * Store in dst[a:b] the coproducts, defined by coproduct[i] = P / lf[i].
 */

static cxx_mpz_poly
coproduct_tree(std::vector<cxx_mpz_poly> & dst,
               std::vector<std::pair<cxx_mpz_poly, int>> const & lf,
               mpz_srcptr modulus, mpz_srcptr invm, size_t a, size_t b,
               bool terminal = false)
{
    ASSERT_ALWAYS(b >= a);
    if (b - a == 0) {
        return 1;
    } else if (b - a == 1) {
        dst[a] = 1;
        return lf[a].first;
    }
    size_t const c = (a + b) / 2;
    cxx_mpz_poly Pleft = coproduct_tree(dst, lf, modulus, invm, a, c);
    cxx_mpz_poly Pright = coproduct_tree(dst, lf, modulus, invm, c, b);
    for (size_t i = a; i < c; i++) {
        mpz_poly_mul(dst[i], dst[i], Pright);
        if (modulus)
            mpz_poly_mod_mpz(dst[i], dst[i], modulus, invm);
    }
    for (size_t i = c; i < b; i++) {
        mpz_poly_mul(dst[i], dst[i], Pleft);
        if (modulus)
            mpz_poly_mod_mpz(dst[i], dst[i], modulus, invm);
    }
    if (terminal)
        return 1;
    cxx_mpz_poly P;
    mpz_poly_mul(P, Pleft, Pright);
    if (modulus)
        mpz_poly_mod_mpz(P, P, modulus, invm);
    return P;
}

static std::vector<cxx_mpz_poly>
coproduct_tree(std::vector<std::pair<cxx_mpz_poly, int>> const & lf,
               mpz_srcptr modulus = nullptr, mpz_srcptr invm = nullptr)
{
    std::vector<cxx_mpz_poly> dst(lf.size());
    coproduct_tree(dst, lf, modulus, invm, 0, lf.size(), true);
    return dst;
}

static cxx_mpz_poly prod(std::vector<std::pair<cxx_mpz_poly, int>> const & lf,
                         mpz_srcptr modulus, mpz_srcptr invm, size_t a,
                         size_t b)
{
    ASSERT_ALWAYS(b >= a);
    if (b - a == 0) {
        return 1;
    } else if (b - a == 1) {
        return lf[a].first;
    }
    size_t const c = (a + b) / 2;
    cxx_mpz_poly Pleft = prod(lf, modulus, invm, a, c);
    cxx_mpz_poly Pright = prod(lf, modulus, invm, c, b);
    cxx_mpz_poly P;
    mpz_poly_mul(P, Pleft, Pright);
    if (modulus)
        mpz_poly_mod_mpz(P, P, modulus, invm);
    return P;
}
cxx_mpz_poly prod(std::vector<std::pair<cxx_mpz_poly, int>> const & lf,
                  mpz_srcptr modulus, mpz_srcptr invm)
{
    return prod(lf, modulus, invm, 0, lf.size());
}

static int
mpz_poly_factor_list_lift(std::vector<std::pair<cxx_mpz_poly, int>> & lf,
        std::vector<cxx_mpz_poly> & cofactor_inverses,
        mpz_poly_srcptr f, mpz_srcptr ell, mpz_srcptr ell2)
{
    {
        /* sanity check */
        cxx_mpz ell_ell;
        mpz_mul(ell_ell, ell, ell);
        ASSERT_ALWAYS(mpz_cmp(ell_ell, ell2) >= 0);

        for (auto const & [ u, e ]: lf) {
            if (e != 1) {
                fprintf(stderr, "Ramified ell not supported\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    auto coprod = coproduct_tree(lf, ell);

    if (cofactor_inverses.empty()) {
        /* Initialize the auxiliary data with the inverse of (f/u) mod u,
         * for each factor */
        for (size_t i = 0; i < lf.size(); i++) {
            cxx_mpz_poly d, a, b;

            /* get v = product of other factors */
            cxx_mpz_poly const & v = coprod[i];
            cxx_mpz_poly const & u = lf[i].first;

            /* get a*u + b*v = 1 */
            mpz_poly_xgcd_mpz(d, u, v, a, b, ell);

            cofactor_inverses.push_back(b);
        }
    } else {
        /* cofactor_inverses contains inverses that are not modulo ell
         * but rather modulo sqrt(ell). It's easy to update.
         * If u0 = lf[i].first is a factor mod ell and v0 = coprod[i] is
         * (f/u0) mod ell, then we want to update b_ to the solution b0 of
         * 1 - b * v0 == 0 mod (u0, ell), which is
         *
         * b_ + ell*b__ = b_ + (1/v0) * (1 - b_ * v0) mod (u0, ell)
         *
         * so b__ = b_ * (1 - b_ * v0) mod (u0, ell).
         *
         * Note that b_ is already such that b_*v0 = 1 mod (u0,
         * sqrt(ell))
         */
        for (size_t i = 0; i < lf.size(); i++) {
            cxx_mpz_poly const & u = lf[i].first;
            cxx_mpz_poly const & v = coprod[i];
            cxx_mpz_poly & b = cofactor_inverses[i];
            cxx_mpz_poly t;
            mpz_poly_mul_mod_f_mod_mpz(t, b, v, u, ell, nullptr, nullptr);
            mpz_poly_sub_ui(t, t, 1);
            mpz_poly_mul_mod_f_mod_mpz(t, b, t, u, ell, nullptr, nullptr);
            mpz_poly_sub_mod_mpz(b, b, t, ell);
        }
    }

    /* compute f - product(everyone) mod ell^2 */
    cxx_mpz_poly f1 = prod(lf, ell2);
    mpz_poly_sub_mod_mpz(f1, f, f1, ell2);
    mpz_poly_divexact_mpz(f1, f1, ell);

    /* Lift all factors.
     *
     * Let u0 = lf[i].first be a factor mod ell, and b0 = cofactor_inverses[i]
     * the inverse of (f/u0) mod (u0, ell). We also have v0 = coprod[i] =
     * (f/u0).
     *
     * We need two steps of Newton lifting.
     *
     * First, we want to lift u0 to the solution u of f - (u * v) = 0 mod
     * ell^2, which is written as u = u0 + ell*u1, with
     *  u = u0 + 1/v0 * (f - u0 * v0)
     * u1 = 1/v0 * f1
     *    = b0 * f1 mod ell
     *
     * Furthermore we can reduce u1 mod u0 as well, since adjusting by
     * ell*k*u0 can be compensated by an adjustment by ell*k*v0 of the
     * other term.
     *
     * In a second step, once we have collected all new factors, we want
     * to update bi. Note that this requires not only next(ui) but also
     * next(vi), which we actually get from the next call to
     * coproduct_tree(). Therefore it's done in the next iteration, after
     * the call to that function.
     *
     */

    for (size_t i = 0; i < lf.size(); i++) {
        cxx_mpz_poly & u = lf[i].first;
        cxx_mpz_poly const & b = cofactor_inverses[i];

        cxx_mpz_poly t;
        mpz_poly_mul_mod_f_mod_mpz(t, b, f1, u, ell, nullptr, nullptr);
        mpz_poly_mul_mpz(t, t, ell);
        mpz_poly_add(u, u, t);
    }

    return 1;
}

std::vector<std::pair<cxx_mpz_poly, int>>
mpz_poly_factor_and_lift_padically(mpz_poly_srcptr f, mpz_srcptr ell, int prec,
                                   gmp_randstate_t rstate)
{
    // this is a false positive
    // coverity[exception_thrown]
    ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_lc(f), 1) == 0);

    auto xfac = mpz_poly_factor(f, ell, rstate);

    /* For each irreducible divisor u of f, we want to precompute the
     * inverse of (f/u) modulo u. It will be needed in the lifting
     * phase.
     */
    std::vector<cxx_mpz_poly> cofactor_inverses;

    std::vector<int> precs;
    for (int p = prec; p != 1; p -= p / 2)
        precs.push_back(p);

    int k = 1;
    cxx_mpz ell_k_minus_1 = 1; /* always ell^(k-1) */
    cxx_mpz ell_k = ell;       /* always ell^k */

    for (size_t i = precs.size(); i--;) {
        cxx_mpz ell_next;
        int const p = precs[i];
        if (p == 2 * k) {
            mpz_mul(ell_next, ell_k, ell_k);
            mpz_mul(ell_k_minus_1, ell_k, ell_k_minus_1);
        } else if (p == 2 * k - 1) {
            mpz_mul(ell_next, ell_k, ell_k_minus_1);
            mpz_mul(ell_k_minus_1, ell_k_minus_1, ell_k_minus_1);
        } else {
            ASSERT_ALWAYS(0);
        }
        mpz_poly_factor_list_lift(xfac, cofactor_inverses, f, ell_k, ell_next);
        mpz_swap(ell_next, ell_k);
        k = p;
    }

    return xfac;
}

int mpz_poly_factor_and_lift_padically(mpz_poly_factor_list_ptr fac,
                                       mpz_poly_srcptr f, mpz_srcptr ell,
                                       int prec, gmp_randstate_t rstate)
{
    mpz_poly_factor_list_set(
        fac, mpz_poly_factor_and_lift_padically(f, ell, prec, rstate));

    return 1;
}

std::string cxx_mpz_poly::print_poly(std::string const & var) const
{
    std::ostringstream os;
    if (x->deg < 0)
        os << "0";
    for (int i = 0; i <= x->deg; i++) {
        int const r = mpz_cmp_ui(x->_coeff[i], 0);
        if (r == 0)
            continue;
        if (r > 0 && os.str().size())
            os << "+";
        if (i == 0) {
            os << x->_coeff[i];
        } else {
            if (mpz_cmp_ui(x->_coeff[i], -1) == 0) {
                os << "-";
            } else if (mpz_cmp_ui(x->_coeff[i], 1) != 0) {
                os << x->_coeff[i] << "*";
            }
            os << var;
            if (i > 1)
                os << "^" << i;
        }
    }
    return os.str();
}

/* this proxy is used for pretty printing in gdb */
std::string cxx_mpz_poly::print_poly() const
{
    return print_poly("x");
}


/* functions for Joux--Lercier and Generalized Joux--Lercier */
/**
 * \brief set the (mpz_t) coefficients of f from ::counter
 * \param[out] f               polynomial
 * \param[out] max_abs_coeffs  largest absolute value of coefficients of f
 * \param[out] next_counter    the next counter to try after the present one,
 * larger than counter+1 sometimes
 * \param[in]  counter         counter storing the coefficients of f
 * \param[in]  bound           max absolute value of coefficients of f
 *
 * \return 1 if the poly is valid (content = 1, not a duplicate, and +/-1 is not
 * a root) (in fact the value returned is != 0) and in this case, max_abs_coeffs
 * is set to the largest absolute value of the coefficients of f
 * \return 0 if the poly corresponding to the counter is a duplicate, and in
 * that case, max_abs_coeffs is set to an error_code (poly duplicate, +/- 1
 * root, content != 1...)
 *
 * Given ::counter, compute the coefficients of ::f.
 * Assume that the leading coefficient of f is > 0, the next one is >= 0, and
 * the last one is not 0 1 <= f[deg] <= bound 0 <= f[deg-1] <= bound -bound <=
 * f[i] <= bound for 0 < i < deg-1 no other test is made: you need to check
 * after if the poly is square-free, irreducible, etc.
 */
int mpz_poly_setcoeffs_counter(mpz_poly_ptr f, int * max_abs_coeffs,
                               unsigned long * next_counter, int deg,
                               unsigned long counter, unsigned int bound)
{
    unsigned int i;
    unsigned long idx = counter;
    int j;
    int *fint, ok = 1;
    // unsigned long idx_next_poly_ok = idx;
    int error_code = 0;       // because I want to know what is the pb
#define POLY_EQUIV_INV_X -1   // the poly is the same as (+/-)f(1/x)
#define POLY_EQUIV_MINUS_X -2 // the poly is the same as (+/-)f(-x)
#define POLY_EQUIV_MINUS_COEFFS -3
#define POLY_ROOT_ONE -4
#define POLY_ROOT_MINUS_ONE -5
#define POLY_CONTENT -6

    fint = (int *)malloc((deg + 1) * sizeof(int));
    // by default the next counter is counter+1
    *next_counter = counter + 1;

    /* compute polynomial f and fint of index idx */
    /* assume that f was already initialized with mpz_poly_init(f, deg) */
    f->deg = deg;

    /* we take 1 <= f[deg] <= bound */
    fint[deg] = 1 + (idx % bound);
    idx = idx / bound;
    /* we take 0 <= f[deg-1] <= bound */
    fint[deg - 1] = idx % (bound + 1);
    idx = idx / (bound + 1);
    for (i = deg - 2; i > 0; i--) {
        /* we take -bound <= f[i] <= bound */
        fint[i] = (idx % (2 * bound + 1)) - bound;
        idx = idx / (2 * bound + 1);
    }

    /* we take -bound <= f[0] <= bound, f[0] <> 0,
       which makes 2*bound possible values */
    ASSERT(idx < 2 * bound);
    fint[0] = (idx < bound) ? idx - bound : idx - (bound - 1);

    /* since f and the reversed polynomial are equivalent, we can assume
       |f[deg]| < |f[0]| or (|f[deg]| = |f[0]| and |f[deg-1]| < |f[1]|) or ...
     */

    /* other point of view
    ok = 1;
    j = 0;
    while ((abs(fint[deg-j]) == abs(fint[j])) && (2*j < deg)){
      j++;
    }
    ok = ((j*2 >= deg) || (abs(fint[deg-j]) < abs(fint[j])));
    */
    ok = 1;
    /* out of the 85176 raw polynomials of degree 4 for bound=6,
       the following test discards 41106, i.e., about 48% */
    for (int i = 0; 2 * i < deg; i++) {
        if (abs(fint[deg - i]) != abs(fint[i])) {
            /* we want |f[deg-i]| < |f[i]| */
            ok = abs(fint[deg - i]) < abs(fint[i]);
            break;
        }
    }
    if (!ok) {
        error_code =
            POLY_EQUIV_INV_X; // the poly is an equivalent to (+/-)x^d*f(1/x)
    }
    if (abs(fint[deg]) >= abs(fint[0]) && (bound > 1)) {
        // next valid poly: the next one is an increment of the leading
        // coefficient. but either leading coeff  = |constant coeff|, and
        // incrementing by one will lead to leading coeff > |constant coeff|,
        // and this is not a valid poly, or we already are in the case: leading
        // coeff > |constant coeff|. In both cases, since the leading
        // coefficient is > 0, it means: set the leading coeff to 1 and
        // increment the next coefficient.
        *next_counter = counter - (counter % bound) + bound;
        // then, what if f[0] = +/- 1?
        if (abs(fint[0]) == 1) {
            // means that setting ld=1 might not be enough: what if f[d-1] >
            // |f[1]| now ? in that case, counter++ is not enough because f[d]=1
            // is the only possibility, so the next counter we are looking for
            // is with f[d-1]++
            if (abs(fint[deg - 1] + 1) > abs(fint[1])) {
                unsigned long mod_j, mod_k;
                // setting ld=1 then incrementing the second high deg coeff is
                // not enough
                idx = *next_counter;
                // set the second high deg coeff to 0
                *next_counter =
                    idx - (idx % (bound * (bound + 1))) + (idx % bound);
                // and increment the third high deg coeff
                *next_counter += bound * (bound + 1);
                j = 2;
                mod_j = bound * (bound + 1) * (2 * bound + 1);
                mod_k = mod_j * (2 * bound + 1);
                while ((deg - j > j) && (fint[j - 1] == 0)) {
                    // setting the deg-j+1 coeff to 0 and incrementing the next
                    // deg-j coeff migh not be enough
                    if (abs(fint[deg - j] + 1) > abs(fint[j])) {
                        if ((fint[deg - j] + 1) < fint[j]) { // negative values
                            // set fint[deg-j] to fint[j]
                            idx = *next_counter;
                            *next_counter = idx - (idx % mod_j) +
                                            (-abs(fint[j]) + bound) * mod_j;
                        } else { //(fint[deg-j]+1) is too large anyway: set it
                                 //to its minimal value which is -|fint[j]| and
                                 //do fint[d-j-1]++
                            // set fint[deg-j] to -|fint[j]| and
                            // fint[deg-j-1]++
                            idx = *next_counter;
                            *next_counter = idx - (idx % mod_j) +
                                            (-abs(fint[j]) + bound) * mod_j +
                                            mod_k;
                        }
                    }
                    j++;
                    mod_j = mod_k;
                    mod_k *= (2 * bound + 1);
                } // end while
            }
        }
    } // end of the computation of the next valid non-duplicate counter

    /* Since f(x) is equivalent to (-1)^deg*f(-x), if f[deg-1] = 0, then the
       largest i = deg-3, deg-5, ..., such that f[i] <> 0 should have f[i] > 0.
       Out of the 44070 remaining polynomials for degree 4 and bound = 6, this
       test discards 3276, i.e., about 7.4% */
    if (ok && fint[deg - 1] == 0) {
        for (int i = deg - 3; i >= 0; i -= 2) {
            if (fint[i]) {
                ok = fint[i] > 0;
                break;
            }
        }
    }
    if ((!ok) && (error_code == 0)) {
        error_code =
            POLY_EQUIV_MINUS_X; // the poly is an equivalent to (-1)^deg*f(-x)
    }

    /* If |f[i]| = |f[deg-i]| for all i, then [f[deg], f[deg-1], ..., f[1],
       f[0]] is equivalent to [s*f[0], s*t*f[1], s*t^2*f[2], ...], where s =
       sign(f[0]), and s*t^i is the sign of f[i] where i is the smallest odd
       index > 0 such that f[i] <> 0. Out of the 40794 remaining polynomials for
       degree 4 and bound = 6, this test discards 468, i.e., about 1.1% */
    if (ok && fint[deg] == abs(fint[0])) /* f[deg] > 0 */
    {
        int s = (fint[0] > 0) ? 1 : -1, t = 0, i;
        for (i = 1; i <= deg; i++) {
            if (2 * i < deg && abs(fint[deg - i]) != abs(fint[i]))
                break;
            if (t == 0 && fint[i] != 0 && (i & 1))
                t = (fint[i] > 0) ? s : -s;
        }
        if (2 * i >= deg) /* |f[i]| = |f[deg-i]| for all i */
        {
            /* if t=0, then all odd coefficients are zero, but since
               |f[deg-i]| = |f[i]| this can only occur for deg even, but then
               we can set t to 1, since t only affects the odd coefficients */
            if (t == 0)
                t = 1;
            for (i = 0; i <= deg; i++) {
                if (fint[deg - i] != s * fint[i]) {
                    ok = fint[deg - i] > s * fint[i];
                    break;
                }
                s = s * t;
            }
        }
    }
    if ((!ok) && (error_code == 0)) {
        error_code = POLY_EQUIV_MINUS_COEFFS;
    }

    /* Out of the 40326 remaining polynomials for degree 4 and bound = 6, this
       test discards 3480, i.e., about 8.6% */
    if (ok) {
        /* check if +-1 is a root of f (very common): compute the sum of the
         * coefficients, and the alterned sum: it should be != 0 */
        // evaluating at 1 means sum the coefficients
        int sum_even_coeffs = 0;
        int sum_odd_coeffs = 0;
        for (j = 0; j <= deg; j += 2) {
            sum_even_coeffs += fint[j];
        }
        for (j = 1; j <= deg; j += 2) {
            sum_odd_coeffs += fint[j];
        }
        ok = (sum_even_coeffs + sum_odd_coeffs) != 0;
        if (!ok) {
            error_code = POLY_ROOT_ONE;
        } else { // eval at -1: alterning sum of coefficients
            ok = (sum_even_coeffs - sum_odd_coeffs) != 0;
            if (!ok) {
                error_code = POLY_ROOT_MINUS_ONE;
            }
        }
    }

    /* set mpz_poly */
    f->deg = deg;
    for (j = 0; j <= deg; j++)
        mpz_set_si(f->_coeff[j], fint[j]);

    /* Out of the 36846 remaining polynomials for degree 4 and bound = 6, this
       test discards 1640, i.e., about 4.5% */
    if (ok) {
        /* content test */
        mpz_t content;
        mpz_init(content);
        mpz_poly_content(content, f);
        ok = mpz_cmp_ui(content, 1) == 0; /* duplicate with f/t */
        if (ok == 0) {
            error_code = POLY_CONTENT;
        }
        mpz_clear(content);
    }

    // compute the max absolute value of coefficients,
    // usefull for computing the lognorm after
    if (ok) {
        *max_abs_coeffs = fint[0]; // this is positive anyway
        for (j = 0; j <= deg; j++) {
            if ((fint[j] < 0) && (*max_abs_coeffs < -fint[j])) {
                *max_abs_coeffs = -fint[j];
            } else { // fint[i] >= 0
                *max_abs_coeffs = std::max(*max_abs_coeffs, fint[j]);
            }
        }
    } else {
        *max_abs_coeffs = error_code;
    }

    free(fint);
    return ok;

#undef POLY_EQUIV_INV_X
#undef POLY_EQUIV_MINUS_X
#undef POLY_EQUIV_MINUS_COEFFS
#undef POLY_ROOT_ONE
#undef POLY_ROOT_MINUS_ONE
#undef POLY_CONTENT
}

/**
 return the total number of polys f that satisfy:
 deg(f) = deg exactly (leading coefficient != 0)
 |f_i| <= bound: the coefficients are bounded by bound
 f_deg > 0: leading coeff strictly positive
 f_{deg-1} > 0
 f_0 != 0: constant coeff non-zero

This corresponds to the number of polynomials that can be enumerated with the
function mpz_poly_setcoeffs_counter
*/
unsigned long mpz_poly_cardinality(int deg, unsigned int bound)
{
    /* we take 1 <= f[d] <= bound */
    /* we take 0 <= f[d-1] <= bound */
    /* we take -bound <= f[i] <= bound for d-2 >= i > 0 */
    /* we take -bound <= f[0] < bound, f[0] <> 0,
       which makes 2*bound possible values */
    unsigned long number_polys = 0;
    int i;
    if (deg >= 3) {
        number_polys = bound * (bound + 1);
        for (i = deg - 2; i > 0; i--) {
            number_polys *= 2 * bound + 1;
        }
        number_polys *= 2 * bound;
    }
    return number_polys;
}

/**
 * return the counter in basis ::bound corresponding to f
 *
 */
unsigned long mpz_poly_getcounter(mpz_poly_ptr f, unsigned int bound)
{
    unsigned int counter;
    int i;
    // leading coeff: 1 <= f->_coeff[deg] <= bound
    counter = mpz_get_ui(f->_coeff[f->deg]);
    counter *= bound;
    // next leading coeff: 0 <= f->_coeff[deg-1] <= bound
    counter += mpz_get_ui(f->_coeff[f->deg - 1]);
    counter *= (bound + 1);
    // next coeffs: -bound <= f->_coeff[i] <= bound
    for (i = f->deg - 2; i > 0; i--) {
        counter += mpz_get_si(f->_coeff[i]) + bound;
        counter *= 2 * bound + 1;
    }
    if (mpz_sgn(f->_coeff[0]) < 0) {
        counter += mpz_get_si(f->_coeff[0]) + bound;
    } else { // f->_coeff[0] > 0
        counter += mpz_get_ui(f->_coeff[0]) + bound - 1;
    }
    return counter;
}

void mpz_poly_setcoeffs_counter_print_error_code(int error_code)
{
    switch (error_code) {
    case -1:
        printf("-1: f <-> +/- x^d f(1/x) \n");
        break;
    case -2:
        printf("-2: f <-> +/- f(-x) \n");
        break;
    case -3:
        printf("-3: f <-> another one with permuted/sign coeffs\n");
        break;
    case -4:
        printf("-4: f(1) = 0\n");
        break;
    case -5:
        printf("-5: f(-1) = 0\n");
        break;
    case -6:
        printf("-6: content(f) > 1\n");
        break;
    default:
        printf(" 0: undetermined error.\n");
    }
}

struct mpz_poly_parser_traits {
    std::string x;
    mpz_poly_parser_traits(std::string const & x)
        : x(x)
    {
    }
    struct unexpected_literal
        : public cado_expression_parser_details::parse_error {
        std::string msg;
        unexpected_literal(std::string const & v)
            : msg(std::string("unexpected literal " + v))
        {
        }
        char const * what() const noexcept override { return msg.c_str(); }
    };
    static constexpr int const accept_literals = 1;
    using type = cxx_mpz_poly;
    using number_type = cxx_mpz;
    static void add(cxx_mpz_poly & c, cxx_mpz_poly const & a,
                    cxx_mpz_poly const & b)
    {
        mpz_poly_add(c, a, b);
    }
    static void sub(cxx_mpz_poly & c, cxx_mpz_poly const & a,
                    cxx_mpz_poly const & b)
    {
        mpz_poly_sub(c, a, b);
    }
    static void neg(cxx_mpz_poly & c, cxx_mpz_poly const & a)
    {
        mpz_poly_neg(c, a);
    }
    static void mul(cxx_mpz_poly & c, cxx_mpz_poly const & a,
                    cxx_mpz_poly const & b)
    {
        mpz_poly_mul(c, a, b);
    }
    static void pow(cxx_mpz_poly & c, cxx_mpz_poly const & a,
                       unsigned long e)
    {
        mpz_poly_pow_ui(c, a, e);
    }
    static void swap(cxx_mpz_poly & a, cxx_mpz_poly & b)
    {
        mpz_poly_swap(a, b);
    }
    static void set(cxx_mpz_poly & a, cxx_mpz const & z)
    {
        mpz_poly_set_mpz(a, z);
    }
    void set_literal_power(cxx_mpz_poly & a, std::string const & v,
                           unsigned long e) const
    {
        if (v == x)
            mpz_poly_set_xi(a, e);
        else
            throw unexpected_literal(v);
    }
};

std::istream & operator>>(std::istream & in,
                          cado::named_proxy<cxx_mpz_poly &> const & F)
{
    std::string line;
    for (;; in.get()) {
        int const c = in.peek();
        if (in.eof() || !isspace(c))
            break;
    }
    if (!getline(in, line))
        return in;
    std::istringstream is(line);

    using poly_parser = cado_expression_parser<mpz_poly_parser_traits>;
    poly_parser P(F.x());
    P.tokenize(is);

    try {
        F.c = P.parse();
    } catch (cado_expression_parser_details::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

std::ostream &
operator<<(std::ostream & o,
           cado::named_proxy<cxx_mpz_poly const &> const & F)
{
    return o << F.c.print_poly(std::string(F.x()));
}

std::ostream & operator<<(std::ostream & os, mpz_poly_coeff_list const & P)
{
    if (P.P.degree() < 0)
        return os << "0";
    for (int i = 0; i <= P.P.degree(); i++) {
        if (i)
            os << P.sep;
        os << P.P->_coeff[i];
    }
    return os;
}

int mpz_poly_set_from_expression(mpz_poly_ptr f, char const * value)
{
    cxx_mpz_poly tmp;
    if (!(std::istringstream(value) >> tmp)) {
        return 0;
    }
    mpz_poly_set(f, tmp);
    return 1;
}

/* Pseudo-reduce a plain polynomial p modulo a non-monic polynomial F.
   The result is of type mpz_polymodF P, and satisfies:
   P->p = lc(F)^P->v * p mod F.

   */
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_poly_reducemodF(
    mpz_polymodF_ptr P, mpz_poly_srcptr p, mpz_poly_srcptr F) const
{
    int v = 0;

    ASSERT_ALWAYS(P->p != p);

    mpz_poly_set(P->p, p);
    P->v = 0;

    if (p->deg < F->deg)
        return;

    int const d = F->deg;

    while (P->p->deg >= d) {
        int const k = P->p->deg;
        int i;

        /* We compute F[d]*p - p[k]*F. In case F[d] divides p[k], we can simply
           compute p - p[k]/F[d]*F. However this will happen rarely with
           Kleinjung's polynomial selection, since lc(F) is large. */

        /* FIXME: in msieve, Jason Papadopoulos reduces by F[d]^d*F(x/F[d])
           instead of F(x). This might avoid one of the for-loops below. */

        // temporary hack: account for the possibility that we're indeed
        // using f_hat instead of f.
        if (mpz_cmp_ui(F->_coeff[d], 1) != 0) {
            v++; /* we consider p/F[d]^v */
#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
            for (i = 0; i < k; ++i)
                mpz_mul(P->p->_coeff[i], P->p->_coeff[i], F->_coeff[d]);
        }

#ifdef HAVE_OPENMP
#pragma omp                                                                    \
    parallel for if (!std::is_same<inf, mpz_poly_notparallel_info>::value)
#endif
        for (i = 0; i < d; ++i)
            mpz_submul(P->p->_coeff[k - d + i], P->p->_coeff[k], F->_coeff[i]);

        mpz_poly_cleandeg(P->p, k - 1);
    }

    P->v = v;
}

/* Set Q=P1*P2 (mod F). Warning: Q might equal P1 (or P2). */
void mpz_polymodF_mul(mpz_polymodF_ptr Q, mpz_polymodF_srcptr P1,
                      mpz_polymodF_srcptr P2, mpz_poly_srcptr F)
{
    mpz_poly_notparallel_info().mpz_polymodF_mul(Q, P1, P2, F);
}
template <typename inf>
void mpz_poly_parallel_interface<inf>::mpz_polymodF_mul(mpz_polymodF_ptr Q,
                                                        mpz_polymodF_srcptr P1,
                                                        mpz_polymodF_srcptr P2,
                                                        mpz_poly_srcptr F) const
{
    mpz_poly prd;
    int v;

    /* beware: if P1 and P2 are zero, P1->p->deg + P2->p->deg = -2 */
    mpz_poly_init(prd, (P1->p->deg == -1) ? -1 : P1->p->deg + P2->p->deg);

    ASSERT_ALWAYS(mpz_poly_normalized_p(P1->p));
    ASSERT_ALWAYS(mpz_poly_normalized_p(P2->p));

    mpz_poly_mul(prd, P1->p, P2->p);
    v = P1->v + P2->v;

    mpz_poly_reducemodF(Q, prd, F);
    Q->v += v;
    mpz_poly_clear(prd);
}

void mpz_polymodF_set_ui(mpz_polymodF_ptr P, unsigned long x)
{
    P->v = 0;
    mpz_poly_realloc(P->p, 1);
    mpz_set_ui(P->p->_coeff[0], x);
    mpz_poly_cleandeg(P->p, 0);
}

void mpz_polymodF_set(mpz_polymodF_ptr P, mpz_polymodF_srcptr Q)
{
    P->v = Q->v;
    mpz_poly_set(P->p, Q->p);
}

void mpz_polymodF_set_from_ab(mpz_polymodF_ptr P, mpz_srcptr a, mpz_srcptr b)
{
    P->v = 0;
    if (mpz_cmp_ui(b, 0) == 0) {
        mpz_poly_realloc(P->p, 1);
        mpz_set(P->p->_coeff[0], a);
        mpz_poly_cleandeg(P->p, 0);
    } else {
        mpz_poly_realloc(P->p, 2);
        mpz_set(P->p->_coeff[0], a);
        mpz_neg(P->p->_coeff[1], b);
        mpz_poly_cleandeg(P->p, 1);
    }
}

void mpz_polymodF_init(mpz_polymodF_ptr P, int x)
{
    P->v = 0;
    mpz_poly_init(P->p, x);
}

void mpz_polymodF_clear(mpz_polymodF_ptr P)
{
    mpz_poly_clear(P->p);
}

void mpz_polymodF_swap(mpz_polymodF_ptr P, mpz_polymodF_ptr Q)
{
    std::swap(P->v, Q->v);
    mpz_poly_swap(P->p, Q->p);
}

struct product_tree {
    size_t i0, i1;
    cxx_mpz_poly A;
    cxx_mpz_poly Ad;
    std::shared_ptr<product_tree> children[2];
};

static std::shared_ptr<product_tree>
polyfromroots_and_derivative(std::vector<cxx_mpz> const & points, size_t i0 = 0,
                             size_t i1 = SIZE_MAX)
{
    if (i1 >= points.size())
        i1 = points.size();
    ASSERT_ALWAYS(i0 <= i1);
    if (i1 == i0)
        return nullptr;
    auto ret = std::make_shared<product_tree>();
    ret->i0 = i0;
    ret->i1 = i1;
    if (i1 == i0 + 1) {
        mpz_poly_set_xi(ret->A, 1);
        mpz_poly_sub_mpz(ret->A, ret->A, points[i0]);
        mpz_poly_set_xi(ret->Ad, 0);
    } else {
        size_t j = (i0 + i1) / 2;
        ret->children[0] = polyfromroots_and_derivative(points, i0, j);
        ret->children[1] = polyfromroots_and_derivative(points, j, i1);

        mpz_poly_mul(ret->A, ret->children[0]->A, ret->children[1]->A);
        cxx_mpz_poly tmp;
        mpz_poly_mul(ret->Ad, ret->children[0]->A, ret->children[1]->Ad);
        mpz_poly_mul(tmp, ret->children[0]->Ad, ret->children[1]->A);
        mpz_poly_add(ret->Ad, ret->Ad, tmp);
    }
    return ret;
}

static void multievaluate_derivative(std::vector<cxx_mpz> & ret,
                                     std::shared_ptr<product_tree> const & T,
                                     cxx_mpz_poly const & P)
{
    cxx_mpz_poly tmp;
    mpz_poly_mod(tmp, P, T->A);
    if (T->i1 == T->i0 + 1) {
        ret.emplace_back(tmp.coeff(0));
    } else {
        multievaluate_derivative(ret, T->children[0], tmp);
        multievaluate_derivative(ret, T->children[1], tmp);
    }
}

static void interpolate_inner(cxx_mpz_poly & num, cxx_mpz & den,
                              std::shared_ptr<product_tree> const & T,
                              std::vector<cxx_mpz> const & multipliers,
                              std::vector<cxx_mpz> const & evaluations)
{
    if (T->i1 == T->i0 + 1) {
        num = T->Ad;
        mpz_poly_mul_mpz(num, num, evaluations[T->i0]);
        den = multipliers[T->i0];
    } else {
        cxx_mpz_poly cnum[2];
        cxx_mpz cden[2];
        interpolate_inner(cnum[0], cden[0], T->children[0], multipliers,
                          evaluations);
        interpolate_inner(cnum[1], cden[1], T->children[1], multipliers,
                          evaluations);
        mpz_poly_mul(cnum[0], cnum[0], T->children[1]->A);
        mpz_poly_mul(cnum[1], cnum[1], T->children[0]->A);
        mpz_poly_mul_mpz(cnum[0], cnum[0], cden[1]);
        mpz_poly_mul_mpz(cnum[1], cnum[1], cden[0]);
        mpz_poly_add(num, cnum[0], cnum[1]);
        mpz_mul(den, cden[0], cden[1]);
    }
}

int mpz_poly_interpolate(mpz_poly_ptr resultant,
                         std::vector<cxx_mpz> const & points,
                         std::vector<cxx_mpz> const & evaluations)
{
    auto T = polyfromroots_and_derivative(points);
    std::vector<cxx_mpz> multipliers;
    multievaluate_derivative(multipliers, T, T->Ad);

    cxx_mpz_poly num;
    cxx_mpz den;
    interpolate_inner(num, den, T, multipliers, evaluations);

    cxx_mpz tmp;
    mpz_poly_content(tmp, num);
    mpz_gcd(tmp, tmp, den);
    mpz_poly_divexact_mpz(num, num, tmp);
    mpz_divexact(den, den, tmp);

    mpz_poly_set(resultant, num);
    if (mpz_sgn(den) < 0) {
        mpz_poly_neg(resultant, resultant);
        mpz_abs(den, den);
    }
    return mpz_cmp_ui(den, 1) == 0;
}

template struct mpz_poly_parallel_interface<mpz_poly_notparallel_info>;
template struct mpz_poly_parallel_interface<mpz_poly_parallel_info>;

/* D <- discriminant (f+k*g), which has degree d */
void mpz_poly_discriminant_of_linear_combination(mpz_poly_ptr D,
                                                 mpz_poly_srcptr f0,
                                                 mpz_poly_srcptr g)
{
    uint32_t **M, pivot;

    cxx_mpz_poly f = f0;

    int const d = f->deg;

    ASSERT_ALWAYS(d <= 9);

    /* we first put in D[i] the value of disc(f + i*g) for 0 <= i <= d,
       thus if disc(f + k*g) = a[d]*k^d + ... + a[0], then
       D[0] = a[0]
       D[1] = a[0] + a[1] + ... + a[d]
       ...
       D[d] = a[0] + a[1]*d + ... + a[d]*d^d */

    mpz_poly_discriminant(mpz_poly_coeff(D, 0), f);
    for (int i = 1; i <= d; i++) {
        /* add g */
        mpz_poly_rotation_ui(f, f, g, 1, 0);
        mpz_poly_discriminant(mpz_poly_coeff(D, i), f);
    }
    mpz_poly_cleandeg(D, d);

    /* initialize matrix coefficients */
    M = (uint32_t **)malloc((d + 1) * sizeof(uint32_t *));
    for (int i = 0; i <= d; i++)
        M[i] = (uint32_t *)malloc((d + 1) * sizeof(uint32_t));

    /* Set M to a vandermonde matrix (M[i][j] = i**(j-1)) */
    for (int i = 0; i <= d; i++) {
        M[i][0] = 1;
        for (int j = 1; j <= d; j++)
            M[i][j] = i * M[i][j - 1];
    }

    /* current_D is M * target_D */

    for (int j = 0; j < d; j++) {
        /* invariant: D[i] = M[i][0] * a[0] + ... + M[i][d] * a[d]
           with M[i][k] = 0 for k < j and k < i */
        for (int i = j + 1; i <= d; i++) {
            /* eliminate M[i][j] */
            pivot = M[i][j] / M[j][j];
            mpz_submul_ui(mpz_poly_coeff(D, i), mpz_poly_coeff_const(D, j),
                          pivot);
            for (int k = j; k <= d; k++)
                M[i][k] -= pivot * M[j][k];
        }
    }

    /* now we have an upper triangular matrix */
    for (int j = d; j > 0; j--) {
        for (int k = j + 1; k <= d; k++)
            mpz_submul_ui(mpz_poly_coeff(D, j), mpz_poly_coeff_const(D, k),
                          M[j][k]);
        ASSERT_ALWAYS(mpz_divisible_ui_p(mpz_poly_coeff_const(D, j), M[j][j]));
        mpz_divexact_ui(mpz_poly_coeff(D, j), mpz_poly_coeff_const(D, j),
                        M[j][j]);
    }

    for (int i = 0; i <= d; i++)
        free(M[i]);
    free(M);
}
