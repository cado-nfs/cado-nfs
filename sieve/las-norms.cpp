#include "cado.h" // IWYU pragma: keep

#include <cinttypes>              // for PRId64, PRIu32
#include <climits>                // for UCHAR_MAX
#include <cstdlib>                // for free, malloc, abs
#include <cstring>                // for memset, size_t, NULL
#include <algorithm>              // for min, max
#include <cmath>                  // for fabs, log2, sqrt, pow, trunc, ceil
#include <cstdint>                // for int64_t, uint32_t
#include <cstdarg>             // IWYU pragma: keep
#include <iomanip>                // for operator<<, setprecision
#include <list>                   // for _List_const_iterator, list
#include <sstream>                // IWYU pragma: keep
#include <utility>                // for swap, pair
#include <gmp.h> // IWYU pragma: keep // for gmp_vfprintf, mpz_srcptr, ...
#include "cxx_mpz.hpp"
#include "las-norms.hpp"
#include "fb-types.h"             // for sublat_t
#include "las-config.h"           // for LOG_BUCKET_REGION, LOGNORM_GUARD_BITS
#include "las-siever-config.hpp"  // for siever_config::side_config, siever_...
#include "las-todo-entry.hpp"     // for las_todo_entry
#include "rho.h"        // dickman_rho_local
#include "verbose.h"    // verbose_output_print
#include "macros.h"

using namespace std;

/***********************************************************************/
/* {{{ some utility stuff */
/****************************************************************************
 * Tricky arithmetic functions to compute some values around ~log2.
 * Log2(n) is computed by the mantissa of n with IEEE 754 representation.
 ****************************************************************************/

/* lg2 is an approximation to log2, that is reasonably fast to compute.
 *
 * Let L be the function defined as follows. 
 * Let z = 2^a * (1+x) for a \in Z, and 0<=x<1. Let L(z) = a + x.
 *
 * We have 0 <= log2(z) - L(z) < 0.0861, for any real z > 0.
 *
 * lg2_raw(z) below returns (L(z) + 0x3FF) * 2^20
 *
 * lg2(z,a,s) below returns (lg2_raw(z) - a) * s.
 *
 * In particular, for a = 0x3FF00000 - G*2^20/s, we have:
 *
 *      lg2_raw(z,a,s/2^20) = L(z)*s+G
 *
 * Note that 0 <= log2(z)*s - L(z)*s < 0.0861*s
 * (IOW, the larger the scale, the larger the error we make here).
 */
static inline double lg2_raw (double i) {
    /* This is claimed to take 3 cycles only. I'm yet to be convinced...
    */
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    __asm__ __volatile__ (
            "psrlq $0x20,  %0    \n"
            "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even though
                                        it is not packed! */
            : "+&x" (i));            /* Really need + here! i is not modified in C code! */
    return i;
#else
    /* Same function, but in x86 gcc needs to convert the input i from a
     * xmm register to a classical register. No other way than use memory.
     * So this function needs at least 6 to 8 cycles more than the previous,
     * which uses ~3 cycles.

     * NOTE: tg declaration is mandatory: it's the place where gcc uses
     * memory to do the conversion. Without it, a warning appears but the
     * code is wrong!
     *
     * TODO: check that. I think we're playing dirty games here, which
     * will blow up sooner or later.
     */
    void *tg = &i;
    return (double)(*(uint64_t *)tg >> 0x20);
#endif
}

/* This is "morally" 2^x -- except that we need to take into account the
 * fact that lg2 is a bit fuzzy.
 */
static double lg2_reciprocal(double x) {
    int x0 = x;
    return exp2(x0) * (1 + x-x0);
}


static inline double lg2 (double i, double offset, double scale) {
    /* see comment above */
    return (lg2_raw(i) - offset) * scale;
}

/* n0 is expected to be a constant */
static inline void memset_with_writeahead(void *s, int c, size_t n, size_t n0)
{
#ifdef __GNUC__
    if (__builtin_constant_p(n) || !__builtin_constant_p(n0)) {
        memset(s, c, n);
        return;
    }
#endif
    if (n <= n0)
        memset(s, c, n0);
    else
        memset(s, c, n);
}


/*}}}*/
/* {{{ some utility stuff which I believe is now covered by double_poly, in
 * fact. Should check that.  */

/* This computes u = scale^deg(f) * f(x/scale)
 * not the same as double_poly_scale, deleted in commit
 * c480fe82174a9de96e1cd35b2317fdf0de3678ab
 *
 * Note that this is a local function only, let's keep it here...
 */
static inline void
double_poly_reverse_scale(double_poly_ptr u, double_poly_srcptr f, const double scale)
{
    double_poly_realloc(u, f->deg + 1);
    if (f->deg < 0) {
        double_poly_set_zero(u);
        return;
    }
    int d = f->deg;
    u->coeff[d] = f->coeff[d];
    for(double h = scale ; d-- ; h *= scale) u->coeff[d] = f->coeff[d] * h;
    double_poly_cleandeg(u, f->deg);
}
/* }}} */


/***********************************************************************/
/* {{{ max absolute norms for a homogeneous polynomial along a line.
 *
 * TODO: We need to do this over a circle as well, because that should be
 * the proper way to do the sieve_info_update_norm_data_Jmax function.
 * Presumably it would suffice to parameterize this by tan(t/2), so that
 * we are led to this function again.
 */
/* {{{ get_maxnorm_aux (for x in (0,s)) */
/* return max |g(x)| for x in (0, s) where s can be negative,
   and g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double get_maxnorm_aux (double_poly_srcptr poly, double s)
{
  double_poly derivative;
  const int d = poly->deg;

  ASSERT_ALWAYS(d >= 0);

  if (d == 0)
    return fabs (poly->coeff[0]);

  double *roots = (double*) malloc (poly->deg * sizeof (double));
  FATAL_ERROR_CHECK(roots == NULL, "malloc failed");

  /* Compute the derivative of polynomial */
  double_poly_init (derivative, d - 1);
  double_poly_derivative (derivative, poly);

  /* Look for extrema of the polynomial, i.e., for roots of the derivative */
  const unsigned int nr_roots = double_poly_compute_roots(roots, derivative, s);

  /* now abscissae of all extrema of poly are 0, roots[0], ...,
     roots[nr_roots-1], s */
  double gmax = fabs (poly->coeff[0]);
  for (unsigned int k = 0; k <= nr_roots; k++)
    {
      double x = (k < nr_roots) ? roots[k] : s;
      double va = fabs (double_poly_eval (poly, x));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  double_poly_clear(derivative);
  return gmax;
}
/* }}} */

/* {{{ get_maxnorm_aux_pm (for x in (-s,s)) */
/* Like get_maxnorm_aux(), but for interval [-s, s] */
static double
get_maxnorm_aux_pm (double_poly_srcptr poly, double s)
{
  double norm1 = get_maxnorm_aux(poly, s);
  double norm2 = get_maxnorm_aux(poly, -s);
  return (norm2 > norm1) ? norm2 : norm1;
}
/* }}} */

/* {{{ get_maxnorm_rectangular */
/* returns the maximal value of |F(x,y)| for -X <= x <= X, 0 <= y <= Y,
 * where F(x,y) is a homogeneous polynomial of degree d.
 * Let F(x,y) = f(x/y)*y^d, and F(x,y) = rev(F)(y,x).
 * Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
 * of the rectangle, i.e.:
 * (a) either on F(X,y) for -Y <= y <= Y (right-hand-side border, and the
 *     mirrored image of the left-hand-side border); We want the maximum of
 *     rev(F)(y,X) = rev(f)(y/X) * X^d for -Y <= j <= Y;
 *     this is rev(f)(t) * X^d for -Y/X <= t <= Y/X.
 * (b) either on F(x,Y) for -X <= x <= X (top border)
 *     = f(x/Y) * Y^d; this is f(t) * Y^d for -X/Y <= t <= X/J.
 * (d) or on F(x,0) for -X <= x <= X (lower border, on the abscissa), but this
 *     maximum is f[d]*X^d, and is attained in (a).
 */
double
get_maxnorm_rectangular (double_poly_srcptr src_poly, const double X,
			 const double Y)
{
  const unsigned int d = src_poly->deg;
  double norm, max_norm;

  /* Make copy of polynomial as we need to revert the coefficients */
  double_poly poly;
  double_poly_init (poly, d);
  double_poly_set (poly, src_poly);

  /* (b) determine the maximum of |f(x)| * Y^d for -X/Y <= x <= X/Y */
  max_norm = get_maxnorm_aux_pm (poly, X/Y) * pow(Y, (double)d);

  /* (a) determine the maximum of |g(y)| for -1 <= y <= 1, with g(y) = F(s,y) */
  double_poly_revert(poly, poly);
  norm = get_maxnorm_aux_pm (poly, Y/X) * pow(X, (double)d);
  if (norm > max_norm)
    max_norm = norm;

  double_poly_clear(poly);

  return max_norm;
}
/* }}} */
/* }}} */

/***********************************************************************/

/*
 * These are two initializations of the algebraic and rational norms.
 * These 2 initializations compute F(i, const j)=f(i) for each line.
 * f(i)=log2(abs(sum[k=0...d] Ak i^k)+1.)*scale+LOGNORM_GUARD_BITS = [LOGNORM_GUARD_BITS...254],
 * where Ak are the coefficients of the polynomial.
 *
 * The classical initialization is slow; the only optimization is done
 * by the fast computation of the log2. The associated error guarantees
 * the precision of the results +/- 1.
 *
 * The "smart" version is faster (more than 10* for the rational side,
 * almost 100* for the algebraic side).
 *
 * Both versions are implemented further down in this file.
 */

/* {{{ On the ::fill function in lognorm_base derived classes. */
/* This function is used to initialize lognorms
 * (=log2(F(i,j)*scale+LOGNORM_GUARD_BITS)
   for the bucket_region S[] number J.
   It's a wrapper; except for trivial degree = 0, it extracts the interesting
   parameters of the lognorm structure and calls the right function.

   - For degree 0, S[] is initialized by a memset: always exact.
   - A special ultra fast init function is used for degree = 1; it could be
     considered as exact (the maximal error is always -/+ 1 on S[]).
   - For smart = 0 and others degrees, the exact values F(i,j) are
     computed with a fast log2. The maximal error is always -/+ 1 on S[].
     It's obviously slow.

   - For smart != 0 and degree > 1, we first approximate the polynomial
     by piecewise linear functions, which are correct up to the
     prescribed multiplicative factor.
*/

/* }}} */

lognorm_base::lognorm_base(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J)
    : logI(logI), J(J)
    /*{{{*/
{
    int64_t H[4] = { Q.a0, Q.b0, Q.a1, Q.b1 };

    /* Update floating point version of polynomial. They will be used in
     * get_maxnorm_rectangular(). */

    mpz_poly_homography (fij, cpoly->pols[side], H);

    /* This is checked for in choose_sieve_area. Homographies with a
     * degree drop are always discarded */
    ASSERT_ALWAYS(fij->deg == cpoly->pols[side]->deg);

    if (Q.doing.side == side) {
        ASSERT_ALWAYS(mpz_poly_divisible_mpz(fij, Q.doing.p));
        mpz_poly_divexact_mpz(fij, fij, Q.doing.p);
    }
    double_poly_set_mpz_poly(fijd, fij);
    // Take sublat into account: multiply all coefs by m^deg.
    // We do it only for the floating point version, that is used to
    // compute a bound on the norms, and in the norm_init phase.
    if (Q.sublat.m > 0)
        double_poly_mul_double(fijd, fijd, pow(Q.sublat.m, fijd->deg));

    int I = 1 << logI;

    maxlog2 = log2(get_maxnorm_rectangular (fijd, (double)(I/2), (double)J));

    /* we know that |F(a,b)| < 2^(maxlog2) or |F(a,b)/q| < 2^(maxlog2)
       depending on the special-q side. */
    /* we used to increase artificially maxlog2, purportedly "to allow larger values of J". I don't see why it's useful. */
    // maxlog2 += 2.0;

    /* we want to map 0 <= x < maxlog2 to LOGNORM_GUARD_BITS <= y < UCHAR_MAX,
       thus y = LOGNORM_GUARD_BITS + x * (UCHAR_MAX-LOGNORM_GUARD_BITS)/maxlog2. */
    scale = (UCHAR_MAX - LOGNORM_GUARD_BITS) / maxlog2;

    /* We require that scale is of the form (int) * 0.025, so that only a small
       number of different factor base slicings can occur. */
    /* Note: if we replace 1/40 by a power-of-two fraction of unity, that
     * could enable some tricks with norms initialization */
    scale = (int)(scale * 40) * 0.025;

    verbose_output_start_batch();
    verbose_output_print (0, 2,
            "# Side %d: log2(maxnorm)=%1.2f scale=%1.2f, logbase=%1.6f",
            side, maxlog2, scale, exp2 (1. / scale));

    /* we want to select relations with a cofactor of less than r bits */
    double max_lambda = (maxlog2 - LOGNORM_GUARD_BITS / scale) / sc.sides[side].lpb;
    double lambda = sc.sides[side].lambda;
    double r = maxlog2 - LOGNORM_GUARD_BITS / scale;

    /* when lambda = 0 (automatic), we take mfb/lpb + 0.3, which is
       experimentally close to optimal in terms of seconds per relation
       (+ 0.2 might be even better on the rational side) */
    if (lambda == 0.0)
        lambda = 0.3 + (double) sc.sides[side].mfb /
            (double) sc.sides[side].lpb ;

    r = std::min(r, lambda * sc.sides[side].lpb);

    /* other option is to ditch this +0.3 and define r (and hence the
     * bound) from mfb directly. Maybe a constant offse would be better
     * than a multiplicative one...
     * r = std::min(r, lambda ? lambda * sc.sides[side].lpb : sc.sides[side].mfb);
     */

    bound = (unsigned char) (r * scale + LOGNORM_GUARD_BITS);

    verbose_output_print (0, 2, " bound=%u\n", bound);
    if (lambda > max_lambda)
        verbose_output_print (0, 2, "# Warning, lambda>%.1f on side %d does "
                "not make sense (capped to limit)\n", max_lambda, side);

    verbose_output_end_batch();
}/*}}}*/

void lognorm_base::norm(mpz_ptr x, int i, unsigned int j) const {
    mpz_poly_homogeneous_eval_siui(x, fij, i, j);
}
unsigned char lognorm_base::lognorm(int i, unsigned int j) const {
    cxx_mpz x;
    norm(x, i, j);
    return log2(mpz_get_d(x)) * scale + LOGNORM_GUARD_BITS;
}

    /* common definitions -- for the moment it's a macro, eventually I
     * expect it's gonna be something else, probably simply replicated
     * code, who knows... */
#define LOGNORM_FILL_COMMON_DEFS()				         \
    unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);\
    unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);\
    unsigned int regions_per_line = 1 << log_regions_per_line;           \
    unsigned int region_rank_in_line = N & (regions_per_line - 1);       \
    unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    \
    unsigned int j1 = j0 + (1 << log_lines_per_region);                  \
    int I = 1 << logI;						         \
    int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;           \
    int i1 = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));                   \
    do {} while (0)

#define LOGNORM_COMMON_HANDLE_ORIGIN() do {				\
    bool has_haxis = !j0;                                               \
    bool has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);   \
    bool has_origin = has_haxis && has_vaxis;                           \
    if (UNLIKELY(has_origin)) {						\
	/* compute only the norm for i = 1. Everybody else is 255. */	\
        memset(S, 255, i1-i0);						\
        if (has_origin) {						\
            double norm = (log2(fabs(fijd->coeff[fijd->deg]))) * scale;	\
            S[1 - i0] = LOGNORM_GUARD_BITS + (unsigned char) (norm);			\
        }								\
        /* And now make sure we start at the next line */		\
        S+=I;								\
	j0++;								\
    }									\
} while (0)

/***********************************************************************/

/* {{{ reference slow code for computing lognorms */
lognorm_reference::lognorm_reference(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J) : lognorm_base(sc, cpoly, side, Q, logI, J)/*{{{*/
{
    /* Knowing the norm on the rational side is bounded by 2^(2^k), compute
     * lognorms approximations for k bits of exponent + NORM_BITS-k bits
     * of mantissa.
     * Do the same for the algebraic side (with the corresponding bound for
     * the norms.
     */
    int k = (int) ceil (log2 (maxlog2));
    int K = 1 << k;
    ASSERT_ALWAYS(lognorm_reference::NORM_BITS >= k);
    int l = lognorm_reference::NORM_BITS - k;
    int L = 1 << l;

    /* extract k bits from the exponent, and l bits from the mantissa */
    double h = 1.0 / (double) L;
    double e,m;
    int i,j;
    for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
    {
        /* e = 2^i for 0 <= i < 2^k */
        for (m = 1.0, j = 0; j < L; j++, m += h)
        {
            /* m = 1 + j/2^l for 0 <= j < 2^l */
            double norm = m * e;
            /* Warning: since sdata->maxlog2 does not usually correspond to
               a power a two, and we consider full binades here, we have to
               take care that values > sdata->maxlog2 do not wrap around to 0 */
            norm = log2 (norm);
            if (norm >= maxlog2)
                lognorm_table[(i << l) + j] = 255;
            else
            {
                norm = norm * scale;
                lognorm_table[(i << l) + j] = LOGNORM_GUARD_BITS + (unsigned char) norm;
            }
        }
    }
}/*}}}*/
/* {{{ void lognorm_fill_rat_reference */
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 *
 * This is adapted from the reference version which was slaughtered in
 * commit 80df430a. It is slow, but at least readable.
 *
 * nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
static void lognorm_fill_rat_reference(
        unsigned char *S,
        uint32_t N,
        int logI,
        double scale,
        cxx_double_poly const & fijd,
        double maxlog2,
        const unsigned char * L) /* L = lognorm_table */
{
    LOGNORM_FILL_COMMON_DEFS();
    LOGNORM_COMMON_HANDLE_ORIGIN();

    double u0 = fijd->coeff[0];
    double u1 = fijd->coeff[1];

    int l = lognorm_reference::NORM_BITS - (int) ceil(log2(maxlog2));

    for(unsigned int j = j0 ; j < j1 ; j++) {
	double z = u0 * j + u1 * i0;
	for (int i = i0; i < i1; i++) {
            uint64_t y;
            /* clang doesn't seem to like this one */
#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && !defined(__clang__)
            /* This does: y = *(uint64_t*)&z */
            asm("":"=x"(y):"0"(z));
#else
            union { uint64_t y; double z; } foo;
            foo.z=z;
            y=foo.y;
#endif
            /* we first get rid of the sign bit, then unshift the
             * exponent.  */
            y = ((y<<1) - (UINT64_C(0x3FF)<<53)) >> (53-l);
            unsigned char n = L[y];
	    ASSERT(n > 0);
	    *S++ = n;
	    z += u1;
	}
    }
}
/* }}} */
/* {{{ void lognorm_fill_alg_reference */
/* Exact initialisation of F(i,j) with degre >= 2 (not mandatory). Slow.
   Internal function, only with simple types, for unit/integration testing. */
static void lognorm_fill_alg_reference (unsigned char *S, uint32_t N, int logI, double scale, cxx_double_poly const & fijd)
{
    LOGNORM_FILL_COMMON_DEFS();
    LOGNORM_COMMON_HANDLE_ORIGIN();

    double modscale = scale/0x100000;
    const double offset = 0x3FF00000 - LOGNORM_GUARD_BITS / modscale;

    cxx_double_poly u;
    for (unsigned int j = j0 ; j < j1 ; j++) {
        double_poly_reverse_scale(u, fijd, j);
        for(int i = i0; i < i0 + I; i++) {
            *S++ = lg2(fabs(double_poly_eval (u, i)), offset, modscale);
        }
    }
}
/* }}} */
void lognorm_reference::fill(unsigned char * S, int N) const/*{{{*/
{
    if (fijd->deg == 1)
        lognorm_fill_rat_reference(S, N, logI, scale, fijd, maxlog2, lognorm_table);
    else
        lognorm_fill_alg_reference(S, N, logI, scale, fijd);
}
/*}}}*/
/* }}} */

/***********************************************************************/

/* {{{ faster code */
lognorm_smart::lognorm_smart(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J) : lognorm_base(sc, cpoly, side, Q, logI, J)/*{{{*/
{
    /* See init_degree_one_norms_bucket_region_smart for the explanation of
     * this table */
    for (int i = 0; i < 257 ; i++)
        cexp2[i] = lg2_reciprocal((i-LOGNORM_GUARD_BITS)/scale);

    if (fijd->deg > 1) {
        piecewise_linear_approximator A(fijd, log(2)/scale);
        int I = 1 << logI;
        G = A.logapprox(-(I/2), I/2);
    }
}/*}}}*/

static inline double compute_y(double G, double offset, double modscale) {
    double res = lg2 ((G) + 1., offset, modscale);
    return res;
}

/* {{{ void lognorm_fill_rat_smart. */
/* Initialize lognorms of the bucket region S[] number N, for F(i,j) with
 * degree = 1.
 */
static void lognorm_fill_rat_smart_inner (unsigned char *S, int i0, int i1, unsigned int j0, unsigned int j1, double scale, cxx_double_poly const & fijd, const double *cexp2)
{
    double modscale = scale / 0x100000;
    double offset = 0x3FF00000 - LOGNORM_GUARD_BITS / modscale;
    /* The lg2 function wants  a positive value.
     * Here, in the classical rational initialization, the degree of the
     * used polynomial F(i,j) is hardcoded to 1. So, it's possible to have
     * F(i,j) = 0, even if i, j and the coefficients of F(i,j) are integers
     * and != 0.  So, I add 1.0 on all G values.  It's not useful to do a
     * fabs(G) here because the code uses always COMPUTE_Y(G) with G >= 0.
     */
// #define COMPUTE_Y(G) lg2 ((G) + 1., offset, modscale)
#define COMPUTE_Y(G) compute_y(G, offset, modscale)
    /* COMPUTE_Y(z) returns LOGNORM_GUARD_BITS + L(z) * scale. Recall that
     * we have 0 <= log2(z) - L(z) < 0.0861, for any real z > 0.
     */

    double u0 = fijd->coeff[0];
    double u1 = fijd->coeff[1];

    if (UNLIKELY(u1 == 0)) {
        /* constant approximating functions do occur, see bug 21684. This
         * is triggered in particular by quadratic polynomials, which we
         * do get with nonlinear polynomial selection.
         */
        for (unsigned int j = j0; j < j1; j++) {
            double g = fabs(u0 * j);
            uint8_t y = COMPUTE_Y(g);
            size_t di = i1 - i0;
            memset_with_writeahead(S, y, di, MEMSET_MIN);
            S += di;
        }
        return;
    }

    using std::signbit;
    double u1_abs = fabs(u1);
    double inv_u1_abs = 1/u1_abs;

    for (unsigned int j = j0; j < j1; j++) {
	int i = i0;
	double g = fabs(u0 * j + u1 * i0);
	uint8_t y;
	double root = -u0 * j / u1;

	bool root_ahead = false;

	/*
	 * Caveat: the sign of g is not significant when g is almost zero
	 * -- then we'd like to use the sign of u1 as indicating the sign
	 *  (after all, g+u1 is the value for the next i).
	 *
	 * To check for "g is almost zero", see bug #16388 and commit
	 * 5f682c6. We use this test, which is admittedly a bit surprising.
	 */
	if (LIKELY(g * (1ULL << 51) >= fabs(u0 * j))) {
	    y = COMPUTE_Y(g);
	    root_ahead = root >= i0;
	} else {
	    y = LOGNORM_GUARD_BITS;
	}

	if (root_ahead) {
            /* One of two possible cases here:
             * g > 0, real root ahead, so decreasing logs, and u1 < 0
             * g < 0, real root ahead, so decreasing logs, and u1 > 0
             */

	    /* Handle sequence of decreasing logs, compute stop points
	     * using incremental updates on the ix values */
            /* Each time we've rounded (down) the log value, since we
             * assume that logs are decreasing, we expect that the exact
             * spot where the log actually takes the very rounded value
             * we've computed is some steps away. Compute that, and
             * advance there.
             *
             * To do that, cexp2[y] incorporates a reciprocal to lg2abs.
             */
            for (; i < i1; y--) {
                double ix = root - cexp2[y] * inv_u1_abs;
                /* conversion rounding is towards zero, while we want it
                 * to be towards +infinity. */
                ix += (ix >= 0);
                if (UNLIKELY(ix >= (double) i1))
		    ix = i1;
		size_t di = (int) ix - i;	/* The cast matters ! */
		i = ix;
		if (!di)
                    break;
		memset_with_writeahead(S, y, di, MEMSET_MIN);
		S += di;
	    }
	    if (i >= i1)
		continue;

	    /* Now compute until y really changes signs. We don't expect
	     * that this will be needed more than very few times.  */
	    g = u0 * j + u1 * i;
            for( ; i < i1 && signbit(g) != signbit(u1) ; i++, g+=u1)
                *S++ = COMPUTE_Y(fabs(g));
	    if (i >= i1)
		continue;

	    /* change sign */
	    g = fabs(g);
	    y = COMPUTE_Y(g);
	}

        /* Invariant: g = fabs(u0 * j + u1 * i), and y = COMPUTE_Y(g) */

	/* One of two possible cases here:
	 * g < 0, real root behind, so increasing logs, and u1 < 0
	 * g > 0, real root behind, so increasing logs, and u1 > 0
         */

        /* If we are in the steep region where the log increases very
         * fast, we can't go for the general loop just now. We have to
         * advance a bit, at least until the wild jumps have calmed down
         * a bit.
         */
        for (; i < i1;) {
            *S++ = y;
            i++;
            g += u1_abs;
            uint8_t oy = y;
            y = COMPUTE_Y(g);
            if (y == oy)
                break;
        }
        if (i >= i1)
            continue;

        /* Now that logs are increasing, we want ix to be the position
         * where the logs reaches the value y+1, of course.
         */
	for (; i < i1; y++) {
	    double ix = root + cexp2[(unsigned int) y + 1] * inv_u1_abs;
            /* conversion rounding is towards zero, while we want it
             * to be towards +infinity. */
            ix += (ix >= 0);
            // ASSERT ((int) ix >= i); // see bug 21518
            if (LIKELY(ix >= (double) i)) {
                /* It is most likely that we'll only see this branch. Yet
                 * bug 21518 seems to trigger a nasty corner case */
                if (UNLIKELY(ix >= i1))
                    ix = i1;
                size_t di = (int) ix - i;	/* The cast matters ! */
                i = ix;
                memset_with_writeahead(S, y, di, MEMSET_MIN);
                S += di;
            }
	}
    }
#undef COMPUTE_Y
}

static void lognorm_fill_rat_smart (unsigned char *S, uint32_t N, int logI, double scale, cxx_double_poly const & fijd, const double *cexp2)
{
    LOGNORM_FILL_COMMON_DEFS();
    LOGNORM_COMMON_HANDLE_ORIGIN();
    lognorm_fill_rat_smart_inner(S, i0, i1, j0, j1, scale, fijd, cexp2);
}
/* }}} */

static void lognorm_fill_alg_smart (unsigned char *S, uint32_t N, int logI, double scale, cxx_double_poly const & fijd, piecewise_linear_function const & G, const double *cexp2) /* {{{ */
{
    LOGNORM_FILL_COMMON_DEFS();
    LOGNORM_COMMON_HANDLE_ORIGIN();
    for(unsigned int j = j0 ; j < j1 ; j++) {
        /* G approximates F(x,1).
         *
         * We have F(i,j) = j^deg(F) * F(i/j,1), so that given a set of linear
         * functions {g} which match F(x,1) on the intervals {[r0,r1[},
         * selecting those which intersect [-I/j, I/j[. We want to
         * evaluate the functions g(x/j)*j^d on the intervals [j*r0,j*r1[.
         *
         * Note that when g=u+vx,
         *      g(x/j)*j^d = u*j^d+v*j^(d-1)x
         *                 = u*j^(d-1) * j + v*j^(d-1) * x
         */
        auto it = G.endpoints.begin();
        double r0 = *it++;
        cxx_double_poly uv_pol(1);
        uv_pol->deg = 1;
        // int i = i0;
        double mj1 = j;
        ASSERT(fijd->deg > 1);
        for(int d = fijd->deg ; --d > 1 ; mj1 *= j);
        for(std::pair<double, double> uv : G.equations) {
            double r1 = *it++;
            int Gi0 = j*r0 + (j*r0 >= 0);
            int Gi1 = j*r1 + (j*r1 >= 0);
            Gi0 = std::max(Gi0, i0);
            Gi1 = std::min(Gi1, i1);
            // ASSERT(Gi0 >= i);
            if (Gi1 > Gi0) {
                // ASSERT(Gi0 == i);
                uv_pol->coeff[0] = uv.first * mj1;
                uv_pol->coeff[1] = uv.second * mj1;
                lognorm_fill_rat_smart_inner(S, Gi0, Gi1, j, j+1, scale, uv_pol, cexp2);
                S += Gi1 - Gi0;
                // i = Gi1;
            }
            r0 = r1;
        }
    }
}
/* }}} */

void lognorm_smart::fill(unsigned char * S, int N) const/*{{{*/
{
    if (fijd->deg > 1)
        lognorm_fill_alg_smart(S, N, logI, scale, fijd, G, cexp2);
    else
        lognorm_fill_rat_smart(S, N, logI, scale, fijd, cexp2);
}
/*}}}*/

/* }}} */


/* TODO: some of the sieve_range_adjust stuff can quite probably be
 * refactored on top of lognorm_base */

/* returns the maximal value of |F(X*cos(t),Y*sin(t))|,
   where F(x,y) is a homogeneous polynomial of degree d, 0 <= t <= pi.
   Let cos(t) = u/sqrt(1+u^2) and sin(t) = 1/sqrt(1+u^2) for -Inf < u < Inf,
   we have F^2(X*cos(t),Y*sin(t)) = F^2(X*u,Y)/(1+u^2)^d where d = deg(F),
   thus its derivative with respect to u is (up to a factor):
   X * f'(u*X/Y) * (1+u^2) - d * Y * u * f(u*X/Y).
*/
static double
get_maxnorm_circular (double_poly_srcptr src_poly, const double X,
		      const double Y)
{
  const unsigned int d = src_poly->deg;
  double_poly poly;
  double *roots, x, y, v, max_norm, t;
  unsigned int nr, i;

  double_poly_init (poly, d + 1);

  /* first compute X * f'(u*X/Y) * (1+u^2) */
  poly->coeff[0] = 0.0;
  poly->coeff[1] = 0.0;
  for (i = 1; i <= d; i++)
    {
      t = pow (X / Y, (double) i - 1.0);
      /* the following will set coefficients of degree 2, 3, ..., d+1 */
      poly->coeff[i + 1] = X * (double) i * src_poly->coeff[i] * t;
      /* the following will add to coefficients of degree 0, 1, ..., d-1 */
      poly->coeff[i - 1] += X * (double) i * src_poly->coeff[i] * t;
    }

  /* now subtract d * Y * u * f(u*X/Y) */
  for (i = 0; i <= d; i++)
    {
      t = src_poly->coeff[i] * pow (X / Y, (double) i);
      poly->coeff[i + 1] -= (double) d * Y * t;
    }
  double_poly_cleandeg(poly, d + 1);

  roots = (double*) malloc ((d + 2) * sizeof (double));
  nr = double_poly_compute_all_roots (roots, poly);

  /* evaluate at y=0 */
  max_norm = fabs (src_poly->coeff[d] * pow (X, (double) d));
  for (i = 0; i < nr; i++)
    {
      double u = roots[i];
      x = X * u / sqrt (1.0 + u * u);
      y = Y / sqrt (1.0 + u * u);
      v = double_poly_eval (src_poly, x / y) * pow (y, (double) d);
      v = fabs (v);
      if (v > max_norm)
	max_norm = v;
    }

  free (roots);
  double_poly_clear(poly);

  return max_norm;
}

/* {{{ various strategies to adjust the sieve area */

void sieve_range_adjust::prepare_fijd()/*{{{*/
{
    int64_t H[4] = { Q.a0, Q.b0, Q.a1, Q.b1 };
    /* We need to get the floating point polynomials. Yes, it will be
     * done several times in the computation, but that's a trivial
     * computation anyway.
     */
    int nsides = cpoly->nb_polys;
    fijd.assign(nsides, {});
    for (int side = 0; side < nsides; side++) {
        cxx_mpz_poly fz;
        mpz_poly_homography (fz, cpoly->pols[side], H);
        if (Q.doing.side == side) {
            ASSERT_ALWAYS(mpz_poly_divisible_mpz(fz, Q.doing.p));
            mpz_poly_divexact_mpz(fz, fz, Q.doing.p);
        }
        double_poly_set_mpz_poly(fijd[side], fz);
    }
}/*}}}*/

/* sieve_range_adjust::sieve_info_update_norm_data_Jmax {{{
 *
 * choose the largest possible J by simply bounding the Fij and Gij
 * polynomials
 * 
 * The image in the a,b-plane of the sieve region might be slanted at an
 * angle against the abscissa, thus even though it has the correct
 * skewness, it might include larger a,b-coordinates than a non-slanted
 * image would.
 * 
 * We compute the maximum norm that would occur if we had a perfect
 * lattice basis in the sense that its image forms a rectangle -A/2 <= a <
 * A/2, 0 <= b < B, with A/2/B = skew, and A*B = I*J*q (assuming J=I/2
 * here).  Thus we have B = A/2/skew, A*A/2/skew = I*I/2*q, A =
 * I*sqrt(q*skew).  The optimal maximum norm is then the maximum of
 * |F(a,b)| in this rectangle, divided by q on the special-q side.
 * 
 * Then we reduce J so that the maximum norm in the image of the actual
 * sieve region is no larger than this optimal maximum, times some
 * constant fudge factor.
 */
int
sieve_range_adjust::sieve_info_update_norm_data_Jmax (bool keep_logI)
{
  // The following parameter controls the scaling on the norm.
  // Relevant values are between 1.0 and 3.0. A higher value means we
  // select higher values of J, and therefore we find more relations, but
  // this increases the time per relation.
  // The value 2.0 seems to be a good compromise. Setting 1.5 reduces the
  // time per relation and number of relations by about 1.5% on a typical
  // RSA704 benchmark.
  const double fudge_factor = 2.0;
  if (!keep_logI)
      logI = (logA+1)/2;
  const double I = (double) (1 << logI);
  const double q = mpz_get_d(Q.doing.p);
  const double skew = cpoly->skew;
  const double A = (1 << logI)*sqrt(q*skew);
  const double B = (1 << (logA - logI))*sqrt(q/skew);
  double Jmax = (1 << (logA - logI));

  prepare_fijd();
  int nsides = cpoly->nb_polys;

  for (int side = 0; side < nsides; side++)
    {
      /* Compute the best possible maximum norm, i.e., assuming a nice
         circular sieve region in the a,b-plane */
      double_poly dpoly;
      double_poly_init (dpoly, cpoly->pols[side]->deg);
      double_poly_set_mpz_poly (dpoly, cpoly->pols[side]);
      if (Q.sublat.m > 0)
          double_poly_mul_double(dpoly, dpoly, pow(Q.sublat.m, cpoly->pols[side]->deg));
      double maxnorm = get_maxnorm_circular (dpoly, fudge_factor*A/2.,
              fudge_factor*B);
      double_poly_clear (dpoly);
      if (side == Q.doing.side)
        maxnorm /= q;

      /* in the (i,j)-plane, the sieving region is rectangular */
      double v = get_maxnorm_rectangular (fijd[side], I/2, Jmax);

      if (v > maxnorm)
      { /* use dichotomy to determine largest Jmax */
          double a, b, c;
          a = 0.0;
          b = Jmax;
          while (trunc (a) != trunc (b))
          {
              c = (a + b) * 0.5;
              v = get_maxnorm_rectangular (fijd[side], I/2, c);
              if (v < maxnorm)
                  a = c;
              else
                  b = c;
          }
          Jmax = trunc (a) + 1; /* +1 since we don't sieve for j = Jmax */
      }
    }

  J = (unsigned int) Jmax;

  return round_to_full_bucket_regions(__func__);
}//}}}

/* {{{ a few helpers (otherwise I'll manage to get things wrong
 * eventually)
 */
sieve_range_adjust::vec<double> operator*(sieve_range_adjust::vec<double> const& a, sieve_range_adjust::mat<int> const& m) {
    return sieve_range_adjust::vec<double>(a[0]*m(0,0)+a[1]*m(1,0),a[0]*m(0,1)+a[1]*m(1,1));
}

qlattice_basis operator*(sieve_range_adjust::mat<int> const& m, qlattice_basis const& Q) {
    qlattice_basis R = Q;       // copy all non-matrix fields.
    R.a0 = m(0,0) * Q.a0 + m(0,1) * Q.a1;
    R.a1 = m(1,0) * Q.a0 + m(1,1) * Q.a1;
    R.b0 = m(0,0) * Q.b0 + m(0,1) * Q.b1;
    R.b1 = m(1,0) * Q.b0 + m(1,1) * Q.b1;
    return R;
}
//}}}

/* estimate_yield_in_sieve_area {{{ 
 *
 * Approximate the integral of
 * dickman_rho(log2(F(x,y))/lpb)*dickman_rho(log2(G(x,y))/lpb), which
 * supposedly gives an idea of the expected yield over this sieve area
 *
 * The integral is taken over a range of size 2^A by integrating over
 * 2^(2N-1) points.
 *
 * The range is taken as one half of the 0-centered range whose width and
 * height are 2^(ceil(A/2)-squeeze)) times 2^(floor(A/2)+squeeze+1) --
 * this larger range has size 2^(A+1).
 * 
 * Integration points are chosen as centers of *rectangles* (not
 * squares) which are proportional to the sieve area.
 *
 * The shuffle[] argument can be used to specify an alternate basis
 */
double sieve_range_adjust::estimate_yield_in_sieve_area(mat<int> const& shuffle, int squeeze, int N)
{
    int nsides = cpoly->nb_polys;
    int nx = 1 << (N - squeeze);
    int ny = 1 << (N + squeeze);

    /* We'll integrate over [-X/2,X/2] times [0,Y/2], which has to have
     * size 2^A */
    double X = 1UL << ((logA-logA/2) - squeeze);
    double Y = 1UL << (logA/2     + squeeze + 1);
    /* In other words, X is I, and Jmax is Y/2. We can see that X*Y/2 =
     * 2^A as desired */

    /* Beware, we're really using (nx+1)*(ny/2+1) points, but weighted. */

    double weightsum = 0;

    double sum = 0;
    for(int i = -nx/2 ; i <= nx/2 ; i++) {
        double x = X/nx * i;
        /* We're doing half of the computation on the y axis, since
         * it's symmetric anyway */
        for(int j = 0 ; j <= ny/2 ; j++) {
            double y = Y/ny * j;
            vec<double> xys = vec<double>(x,y) * shuffle;

            double weight = 1;
            if (i == -nx/2 || i == nx/2) weight /= 2;
            if (j == 0 || j == ny/2) weight /= 2;
            verbose_output_print(0, 4, "# %d %d (%.2f) %.1f %.1f", i, j, weight, xys[0], xys[1]);

            double sprod = 0;
            double p0 = 1;
            for(int side = 0 ; side < nsides ; side++) {
                double z = double_poly_eval_homogeneous(fijd[side], xys[0], xys[1]);
                double a = log2(fabs(z));
                double d = dickman_rho_local(a/conf.sides[side].lpb, fabs(z));
                verbose_output_print(0, 4, " %d %e %e", side, z, d);
                if (side == 0)
                    p0 = d;
                else
                    sprod += p0 * d;
            }
            verbose_output_print(0, 4, " %e\n", sprod);

            weightsum += weight;
            sum += weight*sprod;
        }
    }
    sum /= weightsum;
    sum *= 1UL << logA;
    return sum;
}//}}}

int sieve_range_adjust::adjust_with_estimated_yield()/*{{{*/
{
    prepare_fijd(); // side-effect of the above

    /* List a few candidate matrices which can be used for distorting the
     * sieving range.
     *
     * Because we are including the "swapped" matrices here, we will not
     * investigate negative values of the squeeze parameter.
     *
     * It is important, though, that matrices come here in swapped pairs,
     * since we use that to avoid part of the computation (for squeeze==0,
     * the range is square when A is even, so swapping makes no sense).
     */
    mat<int> shuffle_matrices[] = {
        mat<int>( 1, 0, 0, 1 ),
        mat<int>( 0, 1, 1, 0 ),

        mat<int>( 1, 1, 1, 0 ),
        mat<int>( 1, 0, 1, 1 ),

        mat<int>( 1, 1, 0, -1 ),
        mat<int>( 0, -1, 1, 1 ),

        mat<int>( 1, -1, 0, 1 ),
        mat<int>( 0, 1, 1, -1 ),

        mat<int>( 1, -1, 1, 0 ),
        mat<int>( 1, 0, 1, -1 ),

        /* We're also adding matrices with twos, although we're not really
         * convinced it's worth it. It depends on the polynomial, anyway.
         *
         * On the hsnfs dlp-1024, only 15% of the special-q's among a test
         * of 1000 find a better estimated yield with the matrices below than
         * without.
         */
        mat<int>( 1, 2, 0, 1 ),
        mat<int>( 0, 1, 1, 2 ),

        mat<int>( 1, -2, -1, 1 ),
        mat<int>( -1, 1, 1, -2 ),

        mat<int>( 2, 1, 1, 1 ),
        mat<int>( 1, 1, 2, 1 ),

        mat<int>( -2, 1, 1, 0 ),
        mat<int>( 1, 0, -2, 1 ),

        mat<int>( 1, 0, 2, 1 ),
        mat<int>( 2, 1, 1, 0 ),

        mat<int>( 1, 2, 1, 1 ),
        mat<int>( 1, 1, 1, 2 ),

        mat<int>( 2, -1, 1, -1 ),
        mat<int>( 1, -1, 2, -1 ),

        mat<int>( -1, 2, 0, 1 ),
        mat<int>( 0, 1, -1, 2 ),
    };
    const int nmatrices = sizeof(shuffle_matrices)/sizeof(shuffle_matrices[0]);
#if 0
    /*
     * The following magma code can be used to generate the list of matrices
     * above. I'm only slightly editing the result so that the identity
     * matrix comes first.
MM:=[Matrix(2,2,[a,b,c,d]):a,b,c,d in [-2..2]];
MM:=[M:M in MM|IsUnit(Determinant(M))];
d:=DiagonalMatrix([1,-1]);
s:=Matrix(2,2,[0,1,1,0]);
npos:=func<v|#[a:a in Eltseq(v)|a gt 0]>;
bestrep:=func<s|x[i] where _,i is Max([npos(v):v in x]) where x is Setseq(s)>;
prepr:=func<m|[Eltseq(m),Eltseq(s*m)]>;
B:=[bestrep(a):a in {{a*b*c*x:a in {1,-1},b in {1,d},c in {1,s}}:x in MM}];
  &cat [prepr(x):x in B];
     */
#endif

    double best_sum = 0;
    int best_r = -1;
    int best_squeeze = -1;

    /* We integrate on 2^(2*N-1) points (well, morally 2^(2N), but we halve
     * that by homogeneity) */
    int N = 5;

    double reference = estimate_yield_in_sieve_area(shuffle_matrices[0], 0, N);
    for(int squeeze = ADJUST_STRATEGY2_MIN_SQUEEZE ; squeeze <= ADJUST_STRATEGY2_MAX_SQUEEZE ; squeeze++) {
        for(int r = 0 ; r < nmatrices ; r++) {
            if (squeeze == 0 && (r & 1)) continue;
            mat<int> const & Sr(shuffle_matrices[r]);
            double sum = estimate_yield_in_sieve_area(Sr, squeeze, N);
            if (sum > best_sum) {
                best_r = r;
                best_squeeze = squeeze;
                best_sum = sum;
            }
            verbose_output_print(0, 4, "# estimated yield for rectangle #%d,%d: %e\n", r, squeeze, sum);
        }
    }

    mat<int> const& shuffle (shuffle_matrices[best_r]);

    logI = ((logA-logA/2) - best_squeeze);
    Q = shuffle * Q;
    J = 1 << (logA/2    + best_squeeze);

    std::ostringstream adjust_message;
    adjust_message << "adjusting by ";
    shuffle.print_me(adjust_message);
    adjust_message << " (+" << std::fixed << std::setprecision(2) <<
            100.0*(best_sum/reference-1) << "%)";

    return round_to_full_bucket_regions(__func__, adjust_message.str());
}/*}}}*/

/* return 0 if we should discard that special-q because the rounded
 * region in the (a,b)-plane is flat to the point of having height 0.
 * For diagnostic, we set this->J to the unrounded value (rounding would
 * give 0) and then we return "false".
 *
 * The current check for discarding is whether we do fill one bucket
 * region or not. If we don't even achieve that, we should of course
 * discard.
 *
 * Now for efficiency reasons, the ``minimum reasonable'' number of
 * bucket regions should be more than that.
 */
int sieve_range_adjust::sieve_info_adjust_IJ()/*{{{*/
{
    using namespace std;
    /* compare skewed max-norms: let u0 = [a0, b0] and u1 = [a1, b1],
     * and u'0 = [a0/sqrt(s), sqrt(s)*b0],
     * u'1 = [a1/sqrt(s), sqrt(s)*b1] where s is the skewness.
     * Assume |u'0| <= |u'1|, so that |u'0|^2 <= |u'0|.|u'1|.

     * We know from Gauss reduction that u'0 and u'1 form an angle of at
     * least pi/3, thus since their determinant is q, we have
     * q = |u'0|*|u'1|*sin(angle(u'0,u'1))>=|u'0|*|u'1|*sqrt(3)/2

     * So that
     *  |u'0|^2 <= |u'0|.|u'1| <= 2*q/sqrt(3)

     * Define B := sqrt(q)*sqrt(2/sqrt(3)), then we have:
     *      |a0|/sqrt(s) <= B
     *  and |b0|*sqrt(s) <= B

     * If we take J <= I/2*min(sqrt(s)*B/|a1|, B/sqrt(s)/|b1|), then for
     * any j <= J we have
     *     |a1|/sqrt(s)*J <= I/2*B
     * and |b1|*sqrt(s)*J <= I/2*B,

     * thus for any i,j with |i|<=I/2 and 0<=j<J, we have:
     *     |a|/sqrt(s) = |a0*i+a1*j|/sqrt(s) <= sqrt(q)*I*sqrt(2/sqrt(3))
     * and |b|*sqrt(s) = |b0*i+b1*j|*sqrt(s) <= sqrt(q)*I*sqrt(2/sqrt(3)).
     */

    /*
     * So the strategy above is cado-nfs's legacy strategy for choosing
     * J, based on I.
     *
     * We now design another strategy. Let A be the target sieve area
     * size. We would like to reach the conditions:
     *
     *     for any i,j with |i|<=I/2 and 0<=j<J:
     *     |a|/sqrt(s) = |a0*i+a1*j|/sqrt(s) <= sqrt(q)*sqrt(A)*sqrt(2/sqrt(3))
     *     |b|*sqrt(s) = |b0*i+b1*j|*sqrt(s) <= sqrt(q)*sqrt(A)*sqrt(2/sqrt(3))
     *     I*J = A
     *
     * The geometrical reasoning given above would be the same. However,
     * to achieve I*J=A, we will have to bump I and J slightly.
     *
     * So we're going to start by the previous approach, and adjust it in
     * a second step.
     */
    const double skew = cpoly->skew;
    const double rt_skew = sqrt(skew);
    verbose_output_vfprint(0, 3, gmp_vfprintf,
            "# Called sieve_info_adjust_IJ((a0=%" PRId64 "; b0=%" PRId64
            "; a1=%" PRId64 "; b1=%" PRId64 "), p=%Zd, skew=%f)\n",
            Q.a0, Q.b0, Q.a1, Q.b1,
            (mpz_srcptr) Q.doing.p, skew);
    if (Q.skewed_norm0(skew) > Q.skewed_norm1(skew)) {
        /* exchange u0 and u1, thus I and J */
        swap(Q.a0, Q.a1);
        swap(Q.b0, Q.b1);
    }
    double maxab1 = MAX(abs(Q.a1) / rt_skew, abs(Q.b1) * rt_skew);
    double B = sqrt (2.0 * mpz_get_d(Q.doing.p) / sqrt (3.0));

    uint32_t I = 1UL << ((logA+1)/2);
    J = 1UL << ((logA-1)/2);
    J = (uint32_t) (B / maxab1 * (double) J);

    /* make sure J does not exceed I/2 */
    J = MIN((uint32_t) J, I >> 1);
    logI = (logA + 1) / 2;

    return round_to_full_bucket_regions(__func__);
}/*}}}*/

int sieve_range_adjust::round_to_full_bucket_regions(const char * origin, std::string const & message)/*{{{*/
{
    /* Compute number of i-lines per bucket region, must be integer */
    uint32_t i = 1U << MAX(LOG_BUCKET_REGION - logI, 0);

    // we no longer need to do that.
    // i *= nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */

    /* Bug 15617: if we round up, we are not true to our promises */
    uint32_t nJ = (J / i) * i; /* Round down to multiple of i */

    if (message.empty()) {
        verbose_output_print(0, 3, "# %s(): logI=%d J=%" PRIu32 "\n", origin, logI, nJ);
    } else {
        verbose_output_print(0, 3, "# %s(): logI=%d J=%" PRIu32 " [%s]\n", origin, logI, nJ, message.c_str());
    }
    /* XXX No rounding if we intend to abort */
    if (nJ > 0) J = nJ;
    return nJ > 0;
}/*}}}*/

/*}}}*/

 
uint32_t sieve_range_adjust::get_minimum_J()
{
    /* no longer dependent on the number of threads */
    return 1 << MAX(1, (LOG_BUCKET_REGION - logI));
}

void sieve_range_adjust::set_minimum_J_anyway()
{
    J = get_minimum_J();
}
