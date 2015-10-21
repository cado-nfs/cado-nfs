#include "cado.h"
#include <string.h>
#include <limits.h>
#include <cmath>               /* ceil signbit */
#include <pthread.h>

#ifdef HAVE_SSE41
#include <smmintrin.h>
#elif defined(HAVE_SSSE3)
#include <tmmintrin.h>
#elif defined(HAVE_SSE3)
#include <pmmintrin.h>
#elif defined(HAVE_SSE2)
#include <emmintrin.h>
#endif

#include "portability.h"

#include "las-config.h"
#include "las-types.h"
#include "las-debug.h"
#include "las-norms.h"
#include "utils.h"
#include "verbose.h"

using namespace std;

static long lg_page;

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET)

static inline uint64_t cputicks()
{
  uint64_t r;
  __asm__ __volatile__ ("rdtsc\n shlq $0x20, %%rdx\n orq %%rdx, %%rax\n" : "=a"(r) :: "%rdx", "cc");
  return r;
}

size_t min_stos = 0x8000; /* 32 KB = max size for Intel L0 cache */
size_t max_cache = 0x800000; /* 8 MB = average max size for biggest cache */

void tune_las_memset()
{
    /* larger than the biggest cache on all x86 processors */
    double xput[3];
    uint64_t ticks;
    size_t nmax = 0x8000000;
    void * S = malloc_aligned(nmax + 0x40, 0x1000);
    FATAL_ERROR_CHECK(!S, "malloc failed");

    /* First tuning: comparison between the speed of rep stosq and
     * movaps. All the memsets are here really small; the speed of the
     * memsets are completly different if the adress of the memset is L0
     * cache aligned or not.
     * So, to compute the real speed, we compute 0x40 memsets
     * for each test, with a L0 aligned cache adress + {0...0x3F} ; only
     * the best result is kept.
     */

    /* test for n, 1.25n, 1.5n, 1.75n, 2n, and then replace n by 2n and
     * iterate -- so in a sense, we're testing on a moderately geometric
     * progression of the size n, the ratio being roughly 2^(1/4)=1.189
     */
    for(size_t n = 0x80, base_n = 0x80, stops=0 ;
            n != nmax ;
            n = (++stops&3) ? (n + (base_n>>2)) : (base_n <<= 1))
    {
        /* measure bytes pre cpu tick for method write128 (aka movaps) */
        min_stos = SIZE_MAX;
        ticks = UINT64_MAX;
        ASSERT_ALWAYS(0x40 + n >= 0x20);
        for (unsigned int k = 0; k < 3; k++) {
            uint64_t c;
            c = -cputicks();
            for (size_t shift = 0x40; shift--; )
                las_memset (pointer_arith(S, shift), 0, n);
            c += cputicks();
            if (c < ticks) ticks = c;
        }
        xput[0] = (double) 0x40 * n / ticks;

        /* measure bytes pre cpu tick for method rep_stosq */
        min_stos = n;
        ticks = UINT64_MAX;
        for (unsigned int k = 0; k < 4; k++) {
            uint64_t c;
            c = -cputicks();
            for (size_t shift = 0x40; shift--; )
                las_memset (pointer_arith(S, shift), 0, n);
            c += cputicks();
            if (c < ticks) ticks = c;
        }
        xput[1] = (double) 0x40 * n / ticks;

        verbose_output_print(0, 4, "n=%zu(%zx), write128=%.2f, rep stosq=%.2f\n", n, n, xput[0], xput[1]);

        if (xput[1] > xput[0]) break;
    }
    verbose_output_print(0, 2, "# movaps / rep-stosq cutoff: %zu(0x%zx)\n", min_stos, min_stos);

    /* Second tuning: comparison between the speed of rep stosq and
     * movntps. All the memsets are here really big, because movntps is
     * always slower than rep stosq when the size of the memset is
     * smaller than the biggest cache.
     * So, the memsets are always L0 cache aligned here.
     */

    /* measure bytes pre cpu tick for method direct128 (aka movntps).
     * We're only taking one measure here, as we assume that it is
     * relatively independent of the size.
     */
    {
        size_t n = nmax;
        max_cache = nmax;       /* force direct128 */
        ticks = UINT64_MAX;
        for (unsigned int k = 0; k < 4; k++) {
            uint64_t c;
            c = -cputicks();
            las_memset(S, 0, n);
            c += cputicks();
            if (c < ticks) ticks = c;
        }
        xput[2] = (double) n / ticks;
    }

    max_cache = SIZE_MAX;       /* force rep-stosq */
    /* Now do the same size progression as before, but in reverse order.  */
    for(size_t n = nmax, base_n = nmax, stops=0 ;
            n >= 0x80 ;
            n = (++stops&3) ? (n - (base_n>>3)) : (base_n >>= 1))
    {
        /* measure bytes pre cpu tick for method rep_stosq */
        ticks = UINT64_MAX;
        for (unsigned int k = 0; k < 4; k++) {
            uint64_t c;
            c = -cputicks();
            las_memset (S, 0, n);
            c += cputicks();
            if (c < ticks) ticks = c;
        }
        xput[1] = (double) n / ticks;

        verbose_output_print(0, 4, "n=%zu(%zx), rep stosq=%.2f, direct128=%.2f\n", n, n, xput[1], xput[2]);

        if (xput[1] > xput[2])
            break;

        max_cache = n;
    }
    if (max_cache == SIZE_MAX) {
        verbose_output_print(0, 2, "# rep-stosq / movntps cutoff: never\n");
    } else {
        verbose_output_print(0, 2, "# rep-stosq / movntps cutoff: %zu(0x%zx)\n",
                max_cache, max_cache);
    }
    free_aligned(S);
}

#else /* defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET) */

void tune_las_memset()
{
    verbose_output_print(0, 2, "# las_memset: noasm build, special code disabled\n");
}

#endif /* defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && defined(LAS_MEMSET) */
/* End of X86-64 optimized memset pack */


/****************************************************************************
 * Tricky arithmetic functions to compute some values around ~log2.
 * Log2(n) is computed by the mantissa of n with IEEE 754 representation.
 ****************************************************************************/

/* Input: i, double. i > 0 needed! If not, the result has no sense at all.
   Output: o , trunc(o) == trunc(log2(i)) && o <= log2(i) < o + 0.0861.
   Careful: o ~= log2(i) iif add = 0x3FF00000 & scale = 1/0x100000.
   Add & scale are need to compute o'=f(log2(i)) where f is an affine function.
*/
static inline double lg2 (double i, double add, double scale) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psrlq $0x20,  %0    \n"
	    "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
	    : "+&x" (i));            /* Really need + here! i is not modified in C code! */
  return (i - add) * scale;
#else
  /* Same function, but in x86 gcc needs to transfer the input i from a
     xmm register to a classical register. No other way than use memory.
     So this function needs at least 6 to 8 cycles more than the previous,
     which uses ~3 cycles.
     NOTE: tg declaration is mandatory: it's the place where gcc use memory
     to do the transfert. Without it, a warning appears but the code is false!
  */
  void *tg = &i;
  return (((double)(*(uint64_t *)tg >> 0x20) - add) * scale);
#endif
}

/* Same than previous with a fabs on i before the computation of the log2.
 */
static inline double lg2abs (double i, double add, double scale) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psllq $0x01,  %0    \n"
            "psrlq $0x21,  %0    \n"
	    "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
	    : "+&x" (i));            /* Really need + here! i is not modified in C code! */
  return (i - add) * scale;
#else
  void *tg = &i;
  return (((double)((*(uint64_t *)tg << 1) >> 0x21) - add) * scale);
#endif
}

/* Same than previous with the result is duplicated 8 times
   in a "false" double in SSE version, i.e. in a xmm register, or in a 64 bits
   general register in non SSE version, and written at addr[decal].
   The SSE version has 4 interests: no memory use, no General Register use,
   no xmm -> GR conversion, no * 0x0101010101010101 (3 to 4 cycles) but only
   2 P-instructions (1 cycle).
*/
MAYBE_UNUSED static inline void w64lg2abs(double i, double add, double scale, uint8_t *addr, ssize_t decal) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psllq $0x01,  %0                 \n"
            "psrlq $0x21,  %0                 \n"
	    "cvtdq2pd      %0,              %0\n"
	    : "+&x" (i));
  i = (i - add) * scale;
  __asm__ __volatile__ (
	    "cvttpd2dq     %0,       %0       \n" /* 0000 0000 0000 000Y */
	    "punpcklbw     %0,       %0       \n" /* 0000 0000 0000 00YY */
	    "pshuflw    $0x00,       %0,    %0\n" /* 0000 0000 YYYY YYYY */
	    : "+&x" (i));
  *(double *)&addr[decal] = i;
#else
  void *tg = &i;
  *(uint64_t *)&addr[decal] = UINT64_C(0x0101010101010101) * (uint64_t) (((double)((*(uint64_t *)tg << 1) >> 0x21) - add) * scale);
#endif
}

#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
/* This function is for the SSE2 init algebraics.
   I prefer use X86 ASM directly and avoid intrinsics because the trick of
   cvtdq2pd (I insert 2 doubles values and do a cvtdq2pd on them in order to
   compute their log2).
   This function does a double abs of the __m128d i, and computes the double
   log2 of the result.
   Careful: if i == 0, the result is not predictible.
*/
MAYBE_UNUSED static inline void w128lg2abs(__m128d i, const __m128d add, const __m128d scale, uint8_t *addr, const ssize_t decal) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n" /* i = fast_abs(i) */
	   "shufps    $0xED,    %0,    %0\n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+&x"(i));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  __asm__ __volatile__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "punpcklbw %0,       %0       \n" /* 0000 0000 00XX 00YY */
	   "punpcklwd %0,       %0       \n" /* 0000 XXXX 0000 YYYY */
	   "pshufd    $0xA0,    %0,    %0\n" /* XXXX XXXX YYYY YYYY */
	   : "+&x"(i));
  *(__m128d *)&addr[decal] = i; /* addr and decal are 16 bytes aligned: MOVAPD */
}

/* This function is for the SSE2 init algebraics, but for the exact initialization.
   Same than previous but return a SSE2 register with only 16 lowest bits are computed
   as the 2 results (2 8-bits values).
   CAREFUL! This function returns the computed value in the data i.
   the i value is modified by this function !
*/
static inline __m128i _mm_lg2abs(__m128d *i, const __m128d add, const __m128d scale) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n"
	   "shufps    $0xED,    %0,    %0\n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+&x"(*i));
  *i = _mm_mul_pd(_mm_sub_pd(*i, add), scale);
  __asm__ __volatile__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "packuswb  %0,       %0       \n" /* 0000 0000 0000 00XY */
	   : "+&x"(*i));
  return *(__m128i *) i;
}
#endif  /* defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) */

static inline double compute_f (const unsigned int d, const double *u, const double h) {
  size_t k = (size_t) d;
  double f = u[k];
  switch (k) {
  default: do { f = f * h + u[--k]; } while (k > 9);
  case 9: f = f * h + u[8];
  case 8: f = f * h + u[7];
  case 7: f = f * h + u[6];
  case 6: f = f * h + u[5];
  case 5: f = f * h + u[4];
  case 4: f = f * h + u[3];
  case 3: f = f * h + u[2];
  case 2: f = f * h + u[1];
  case 1: f = f * h + u[0];
  case 0: break;
  }
  return f;
}

/* This function computes the roots of the polynome F(i,1) = ln2(f(i))
   and the roots of the first and second derivates of F.
   These roots are useful to begin the smart initialization of the normalization.
   cf las-config.h, SMART_NORM.

   This function is internal. Don't use it. Use the wrapper below.

   Inputs :
   F(i,1): degree and coeff(icients)
   max_abs_root = Absolute maximum value of a root (before and after, the root is useless).
   precision = minimal need precision for each root (usually 1/J = 2/I).
   Output:
   nroots = number of roots.
   roots = the roots of F, F', F", and pseudo root 0.0 added in F roots.

   NEED:
   1. roots size must >= 1 + degree  + max (0, degree - 1) + max (0, 2 * degree - 2)
   so degree * 4 + 1.
   2. p & Roots: the biggest polynome is F", degree = 2 * (degree - 1).
   So, p size is 2 * degree - 1, and Roots size is 2 * degree - 2.
*/

/* Only for qsort. MacOSX doesn't accept a function is declared into another. Odd... */
static int cmp_root (const void *a, const void *b) {
  return (((root_ptr) a)->value > ((root_ptr)b)->value) ? 1 : -1;
}

/* The norm initialisation code needs the page size; we set it here as this
   function is not thread-safe, to avoid race conditions */
static void set_lg_page()
{
  if (lg_page == 0)
    lg_page = pagesize();
}

void init_norms_roots_internal (unsigned int degree, double *coeff, double max_abs_root, double precision, unsigned int *nroots, root_ptr roots)
{

  const double_poly_t f = {{ (int) degree, coeff }};
  double_poly_t df, ddf, f_ddf, df_df, d2f;
  mpz_t           p[(degree << 1) + 1];
  root_struct Roots[degree << 1];
  unsigned int n, cumul_nroots;
  size_t k;

  set_lg_page();
  for (k = degree << 1, mpz_init (p[k]); k--; mpz_init (p[k]), root_struct_init (&(Roots[k])));
  for (k = degree + 1; k--; mpz_set_d (p[k], coeff[k]));

  /* Pseudo root 0.0 is inserted first as a root of F */
  roots[0].derivate = 0;
  roots[0].value = 0.;
  cumul_nroots = 1;

  if (degree) {
    /* The roots of F are inserted in roots */
    n = numberOfRealRoots (p, degree, max_abs_root, 0, Roots);
    for (k = n; k--; ) {
        roots[k + cumul_nroots].derivate = 0;
        roots[k + cumul_nroots].value = rootRefine (&(Roots[k]), p, degree, precision);
    }
    cumul_nroots += n;

    /* Computation of F' */
    double_poly_init (df, MAX(0,((int)degree - 1)));
    double_poly_derivative (df, f);

    /* The roots of F' are inserted in roots */
    for (k = df->deg + 1; k--; mpz_set_d (p[k], df->coeff[k]));
    n = numberOfRealRoots (p, df->deg, max_abs_root, 0, Roots);
    for (k = n; k--; ) {
        roots[k + cumul_nroots].derivate = 1;
        roots[k + cumul_nroots].value = rootRefine (&(Roots[k]), p, df->deg, precision);
    }
    cumul_nroots += n;

    /* Computation of F" */
    double_poly_init (df_df, df->deg + df->deg);
    double_poly_init (ddf, MAX(0,((int)df->deg - 1)));
    double_poly_init (f_ddf, f->deg + ddf->deg);
    double_poly_init (d2f, MAX(f_ddf->deg, df_df->deg));
    double_poly_product (df_df, df, df);
    double_poly_derivative (ddf, df);
    double_poly_product (f_ddf, f, ddf);
    double_poly_subtract (d2f, f_ddf, df_df);

    /* The roots of F" are inserted in roots */
    for (k = d2f->deg + 1; k--; mpz_set_d (p[k], d2f->coeff[k]));
    n = numberOfRealRoots (p, d2f->deg, max_abs_root, 0, Roots);
    for (k = n; k--; ) {
        roots[k + cumul_nroots].derivate = 2;
        roots[k + cumul_nroots].value = rootRefine (&(Roots[k]), p, d2f->deg, precision);
    }
    cumul_nroots += n;

    /* Clear double poly */
    double_poly_clear(df);
    double_poly_clear(ddf);
    double_poly_clear(f_ddf);
    double_poly_clear(df_df);
    double_poly_clear(d2f);

    /* roots must be sort */
    qsort (roots, cumul_nroots, sizeof(root_t), cmp_root);
  }
  *nroots = cumul_nroots;

  for (k = degree << 1, mpz_clear (p[k]); k--; mpz_clear (p[k]), root_struct_clear (&(Roots[k])));
}

/* A wrapper for the function above */
void init_norms_roots (sieve_info_ptr si, unsigned int side)
{
  init_norms_roots_internal (si->cpoly->pols[side]->deg, si->sides[side]->fijd,
			     (double) ((si->I + 16) >> 1), 1. / (double) ((si->I) >> 1),
			     &(si->sides[side]->nroots), si->sides[side]->roots);
}

/* Initialize lognorms of the bucket region S[] number J, for F(i,j) with
 * degree = 1.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
/* Internal function, only with simple types, for unit/integration testing */
void init_degree_one_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, double u0, double u1, double *cexp2) {

  /* The computation of the log2 by inttruncfastlog2 needs a value >= 1.
     Here, in the classical rational initialization, the degree of the used
     polynome F(i,j) is hardcoded to 1. So, it's possible to have F(i,j) = 0,
     even if i, j and the coefficients of F(i,j) are integers and != 0.
     So, I add 1.0 on all G values.
     It's not useful to do a fabs(G) here because the code uses always COMPUTE_Y(G)
     with G >= 0. */
#define COMPUTE_Y(G) lg2 ((G) + 1., add, scale)
  /* For internal debug: print all */
  // #define DEBUG_INIT_RAT 1
  double add, d0_init, g, rac, u0J, d0, d1, i;
  size_t ts;
  uint32_t endJ;
  int int_i;
  unsigned int inc;
  uint8_t oy, y;

  ASSERT_ALWAYS (u1 != 0.);

  const int Idiv2 = (int) I >> 1;
  const double Idiv2_double = (double) Idiv2, Idiv2_double_minus_one = Idiv2_double - 1., invu1 = 1./u1;

  endJ = LOG_BUCKET_REGION - cado_ctz(I);
  J <<= endJ;
  endJ = (1U << endJ) + J;

#ifdef DEBUG_INIT_RAT
  fprintf (stderr, "Begin: j=%u, u0=%.20e u1=%.20e, scale=%.20e, rac=%.20e\n", J, u0, u1, scale, u0 * J * (-invu1));
#endif

  scale *= 1./0x100000;
  add = 0x3FF00000 - GUARD / scale;
  u0J = u0 * J;
  d0_init = cexp2[((unsigned int)GUARD) - 1U];
  for (; J < endJ ; J++, u0J += u0) {
    int_i = -Idiv2;
    g = u0J + u1 * int_i;
    rac = u0J * (-invu1);
    d0 = d0_init;
    d1 = rac - d0 * rac;
    /* g sign is mandatory at the beginning of the line intialization.
       The sign of g is not significant if fabs(g) * 1ULL<<51 < fabs(u0J)
       In this case, the sign of g is the sign of u1, because g+u1
       is the next g value.
    */
    if (LIKELY(fabs(g) * (1ULL << 51) >= fabs(u0J)))
      if (signbit(g)) {
	g = -g;
	y = COMPUTE_Y(g);
	if (rac >= -Idiv2) goto cas3; else goto cas2;
      }
      else {
	y = COMPUTE_Y(g);
	if (rac >= -Idiv2) goto cas1; else goto cas4;
      }
    else {
      y = GUARD;
      if (signbit(u1)) goto cas2; else goto cas4;
    }
  cas1:
    /* In this case, we exit from the loop when ts == 0, at the exception
       of the first iteration. In this special case, old_i = -Idiv2 and
       int_i = trunc (i), where i=[inverse of the function g](trunc(y)) and
       y=g(old_i).
       So, it's possible if y is very near trunc(y), old_i == int_i, so ts == 0.
       We have to iterate at least one time to avoid this case => this is the
       use of inc here. */
    for (i = rac + cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A1.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A1 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN)) {
	if (!(ts + inc)) goto np1;
	memset(S, y, MEMSET_MIN);
      }
      else memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
      inc = 0;
    }
  np1:
    g = u0J + u1 * int_i;
    if (UNLIKELY(trunc(rac) >= Idiv2_double_minus_one)) {
      for ( ; int_i < Idiv2; int_i++) {
	y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g += u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g += u1) {
      y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = COMPUTE_Y(g);
  cas2:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i >= Idiv2) goto nextj;
      oy = y;
      g -= u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac - cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A4.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A4 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN))
	memset (S, y, MEMSET_MIN);
      else
	memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
    }

    /* Now, the same from cas1 but log2(-g): CAREFUL, not the same formula */
  cas3:
    for (i = rac - cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B1.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B1 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN)) {
	if (!(ts + inc)) goto np2;
	memset(S, y, MEMSET_MIN);
      }
      else memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
      inc = 0;
    }
  np2:
    g = -(u0J + u1 * int_i);
    if (UNLIKELY(trunc(rac) >= Idiv2_double_minus_one)) {
      for ( ; int_i < Idiv2; int_i++) {
	y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g -= u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g -= u1) {
      y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = COMPUTE_Y(g);
  cas4:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i >= Idiv2) goto nextj;
      oy = y;
      g += u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac + cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B4.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B4 : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN))
	memset(S, y, MEMSET_MIN);
      else
	memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
    }
  nextj:
    for (;0;); /* gcc needs something after a label */
  }
}

static inline void
poly_scale_double (double  *u, const double *t, unsigned int d, const double h)
{
  u[d] = t[d];
  for (double hpow = h; d--; hpow *= h) u[d] = t[d] * hpow;
}

#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
static inline void
poly_scale_m128d (__m128d  *u, const double *t, unsigned int d, const double h)
{
  u[d] = _mm_set1_pd (t[d]);
  for (double hpow = h; d--; hpow *= h) u[d] = _mm_set1_pd (t[d] * hpow);
}
#endif  /* defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) */

/* Exact initialisation of F(i,j) with degre >= 2 (not mandatory). Slow.
   Internal function, only with simple types, for unit/integration testing. */
void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd)
{
  unsigned char *beginS = S;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - cado_ctz (I);
  scale *= 1./0x100000;
  const double add = 0x3FF00000 - GUARD / scale;
  J <<= endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;

#if !defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) || !defined(HAVE_SSSE3) /* Light optimization, only log2 */

  for (; J < endJ; J++) {
    const unsigned char *endS = S + I;
    double f, h, u[d+1];
    poly_scale_double (u, fijd, d, (double) J);
    h = (double) (-(int32_t) (I >> 1));
    do {
      f = compute_f (d, u, h);
      h += 1.;
      *S++ = (unsigned char) lg2abs (f, add, scale);
    } while (S < endS);
  }

#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM && HAVE_SSSE3 : optimized part. Stupid but fast code. */

#define BEGIN ".balign 8\n 0:\n add $0x10,%[S]\n"
#define FU(A) "movapd %[u" #A "],%[f]\n"
#define U8 "mulpd %[h],%[f]\n addpd %[u8],%[f]\n"
#define U7 "mulpd %[h],%[f]\n addpd %[u7],%[f]\n"
#define U6 "mulpd %[h],%[f]\n addpd %[u6],%[f]\n"
#define U5 "mulpd %[h],%[f]\n addpd %[u5],%[f]\n"
#define U4 "mulpd %[h],%[f]\n addpd %[u4],%[f]\n"
#define U3 "mulpd %[h],%[f]\n addpd %[u3],%[f]\n"
#define U2 "mulpd %[h],%[f]\n addpd %[u2],%[f]\n"
#define U1 "mulpd %[h],%[f]\n addpd %[u1],%[f]\n"
#define U0 "mulpd %[h],%[f]\n addpd %[u0],%[f]\n"
#define LG2ABS							       \
  "psllq $1,%[f]\n"						       \
    "psrlq $1,%[f]\n"						       \
    "shufps $0xED,%[f],%[f]\n"					       \
    "cvtdq2pd %[f],%[f]\n"					       \
    "subpd %[_add],%[f]\n"					       \
    "mulpd %[_scale],%[f]\n"					       \
    "cvttpd2dq %[f],%[f]\n"					       \
    "packssdw %[f],%[f]\n"					       \
    "packuswb %[f],%[f]\n"					       \
    "palignr $2,%[cumul],%[f]\n"				       \
    "movapd %[f],%[cumul]\n"					       \
    "addpd %[_two],%[h]\n"
#define END					\
  "cmp %[endS],%[S]\n"				\
    "movapd %[f],-0x10(%[S])\n"			\
    "jl 0b\n"

  for (; J < endJ; J++) {
    const __m128d _two = _mm_set1_pd(2.), _scale = _mm_set1_pd (scale), _add = _mm_set1_pd (add);
    __m128d f, h, u[d+1];
    __m128i cumul;
    unsigned char *endS = S + I;
    poly_scale_m128d (u, fijd, d, (double) J);
    h = _mm_set_pd ((double) (1 - (int32_t) (I >> 1)), (double) (- (int32_t) (I >> 1)));

    /* These ASM & switch are really ugly. But it's the ONLY way to
       be sure all the line of S is set only with registers. */
    switch (d) {
    case 2:
      __asm__ __volatile__ ( BEGIN
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]));
      break;
    case 3:
      __asm__ __volatile__ ( BEGIN
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]));
      break;
    case 4:
      __asm__ __volatile__ ( BEGIN
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]));
      break;
    case 5:
      __asm__ __volatile__ ( BEGIN
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]));
      break;
    case 6:
      __asm__ __volatile__ ( BEGIN
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]));
      break;
    case 7:
      __asm__ __volatile__ ( BEGIN
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]));
      break;
    case 8:
      __asm__ __volatile__ ( BEGIN
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]), [u8]"x"(u[8]));
      break;
    case 9:
      __asm__ __volatile__ ( BEGIN
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_add]"x"(_add),
	[u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]), [u8]"x"(u[8]), [u9]"x"(u[9]));
      break;
    default:
      do {
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&f, _add, _scale), cumul, 2);
	*(__m128i *) S = cumul;
	S += 16;
      } while (S < endS);
    } /* End of the switch */
  } /* End of the line */
#undef BEGIN
#undef FU
#undef U0
#undef U1
#undef U2
#undef U3
#undef U4
#undef U5
#undef U6
#undef U7
#undef U8
#undef LG2ABS
#undef END
#endif /* End of HAVE_GCC_STYLE_AMD64_INLINE_ASM */
  /* Special ultra rare case. The correction of log2(F(0,0)) is false, because
     the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) beginS[I>>1] = GUARD;
}

/* Approximates [S[x],S[y][ by the segment (x,fx),(y,fy).
   (y - x) / (fy - fx) should be really larger than one for speed.
   CAREFUL : Need y > x. */
static inline void Fill_S (unsigned char *S, int x, double fx, int y, double fy) {
  int next_f = (int) fx, step_f = (int) fy - next_f;

  if (LIKELY (step_f)) {
    double m = (double) (y - x) / (fy - fx), next_m = (double) (x + 1);
    if (step_f < 0.) {
      m = fabs(m);
      step_f = -1;
      next_m += (fx - floor(fx)) * m;
    }
    else {
      step_f = 1;
      next_m += (ceil(fx) - fx) * m;
    }
    for (;;) {
      int next_x = (int) next_m;
      if (UNLIKELY (next_x >= y)) break;
      memset (S + x, next_f, (size_t) (next_x - x));
      x = next_x;
      next_f += step_f;
      next_m += m;
    }
  }
  memset (S + x, next_f, (size_t) (y - x));
}

/* This functions sets all the segments of the contiguous unset values of the line S
   with F(i, const j) by a polygonal approximation on these segments.
   Derecursivate optimal version.
*/
static inline void poly_approx_on_S (unsigned char *S, const unsigned int degree, const double *coeff, const double my_scale MAYBE_UNUSED, const double scale MAYBE_UNUSED, const double add MAYBE_UNUSED, const unsigned int nsg, const sg_t *sg) {
#define SIZE_STACK 256 /* In fact, on RSA704 benchmark, the max is... 10 */
  typedef struct x_fx_s {
    int x;
    double fx; } x_fx_t;
  size_t current_stack, max_stack = SIZE_STACK;
  x_fx_t begin, current, *stack = (x_fx_t *) malloc_check (sizeof(*stack) *  max_stack);

  /* At least 3 segments: -Idiv2, neighbourhood of 0., Idiv2 - 1 */
  ASSERT(nsg >= 3);
  for (size_t parse_sg = 1; parse_sg < nsg; ++parse_sg) {
    ASSERT (sg[parse_sg - 1].end < sg[parse_sg].begin);
    begin.x = sg[parse_sg - 1].end;
    begin.fx = sg[parse_sg - 1].f_end;
    current.x = sg[parse_sg].begin;
    current.fx = sg[parse_sg].f_begin;
    current_stack = 0;
    for (;;) {
      if (LIKELY((unsigned int) (current.x - begin.x) >= SMART_NORM_LENGTH)) {
	x_fx_t possible_new;
	possible_new.x = (begin.x + current.x) >> 1;
	possible_new.fx = compute_f (degree, coeff, (double) possible_new.x);
	possible_new.fx = lg2abs (possible_new.fx, add, scale);
	if (UNLIKELY(fabs (possible_new.fx + possible_new.fx - current.fx - begin.fx)) > (2 * SMART_NORM_DISTANCE)) {
	  stack[current_stack++] = current;
	  current = possible_new;
	  if (UNLIKELY(current_stack == max_stack)) {
	    max_stack += (max_stack >> 1);
	    if ((stack = (x_fx_t *) realloc (stack, sizeof(*stack) * max_stack)) == NULL) {
	      fprintf (stderr, "Error, realloc of %zu bytes failed\n", sizeof(*stack) * max_stack);
	      abort ();
	    }
	  }
	  continue;
	}
      }
      Fill_S (S, begin.x, begin.fx, current.x, current.fx);
      if (!current_stack) break;
      begin = current;
      current = stack[--current_stack];
    }
  }
  free (stack);
}

/* Smart initialization of the normalization. Only useful for degree >= 2.
   Cf las-config.h, SMART_INIT for the algorithm.
   No SSE version: it's completly unreadable, and the gain is not really interesting.
   Internal function, only with simple types, for unit/integration testing */
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double original_scale, unsigned int d, double *fijd, unsigned int nroots, root_ptr roots)
{
  ASSERT (d >= 2);
  /* F, F' and F" roots needs stability for their neighbourhood ?
     F roots, sure; F', sure not; F"... maybe. */
  const unsigned char stability_for_derivate[3] =
    { SMART_NORM_STABILITY, 0, SMART_NORM_STABILITY };
  const ssize_t Idiv2 = (ssize_t) (I >> 1);
  sg_t sg[d * 4 + 3]; /* For (F, F', F") roots, -Idiv2, Idiv2 - 1, 0 */
  size_t nsg;

  unsigned char *beginS;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - cado_ctz (I);
  double scale = original_scale * (1./0x100000);
  const double add = ((double) 0x3FF00000) - GUARD / scale;
  J <<= endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;
  S += Idiv2; /* CAREFUL! Here, *S is the middle of the first line we have to compute! */
  beginS = S;

  for (; J < endJ; S += I, J++) {
    double u[d+1], g, hl, hr;
    unsigned char fg, f1, f2;
    unsigned int cptf2id = 0, cpt;
    ssize_t ih;
    {
      ASSERT_ALWAYS(lg_page != 0);
      for (ssize_t k = I; (k -= lg_page) >= 0; __builtin_prefetch (S + k, 1));
    }
    poly_scale_double (u, fijd, d, (double) J);

    /* Insertion of point (-Idiv2, F(-Idiv2)) in sg[0] : an artificial one-point segment. */
    g = compute_f (d, u, (double) -Idiv2);
    g = lg2abs (g, add, scale);
    sg[0].begin   = sg[0].end = -(int)Idiv2;
    sg[0].f_begin = sg[0].f_end = g;
    nsg = 1;
    S[-Idiv2] = (uint8_t) g;

    for (size_t r = 0; r < nroots; r++) {
      hl = floor((double)J * roots[r].value);
      ih = (ssize_t) hl;
      /* Need stability for this root ? */
      if (UNLIKELY(!stability_for_derivate[roots[r].derivate])) {
	/* No. It's an one-point segment. */
	/* Is it not in the interesting region or in the last segment ? Yes -> next root */
	if (ih <= -Idiv2 || ih >= Idiv2 - 1 || ih <= sg[nsg - 1].end) continue;
	/* Ok, we insert this one-point segment */
	sg[nsg].begin = sg[nsg].end = ih;
	g = compute_f (d, u, hl);
	g = lg2abs (g, add, scale);
	sg[nsg].f_begin = sg[nsg].f_end = g;
	S[ih] = (unsigned char) g;
      }
      else {
	/* It's a real (non an one-point) segment. */
	/* Is the root really far after the interesting zone ? Yes -> end of roots */
	if (ih > Idiv2 - 1 + SMART_NORM_INFLUENCE) break;

	/* Is the root really far before the interesting zone or is this segment is
	   included (in all sides) in the previous segment ? Yes -> next root */
	if (ih < -Idiv2 - SMART_NORM_INFLUENCE ||
	    ih <= sg[nsg - 1].end - SMART_NORM_STABILITY) continue;

	/* OK, we have the right side to compute, and maybe the left side */
	hr = hl;

	/* Special case for left side. Is the root :
	   1. Near and before the interesting zone OR
	   2. in the previous segment AND [the beginning of the previous segment
   	   is far enough of the root OR this previous segment is the first, so
              its beginning is the beginning of the interesting zone] ? */
	if (UNLIKELY(ih <= -Idiv2 || (sg[nsg - 1].end >= ih &&
				      (nsg == 1 || sg[nsg - 1].begin + SMART_NORM_STABILITY <= ih)))) {
	  /* OK, this segment and the previous could be fusionned on this left side.
	     It's not 100% true in fact, but really very probable. */
	  sg[nsg].begin   = sg[nsg - 1].begin;
	  sg[nsg].f_begin = sg[nsg - 1].f_begin;
	  fg = (ih >= -Idiv2 && ih < Idiv2) ? S[ih] : 0;
	} else {

	  /* Here we compute the left side. -Idiv2 < ih <= Idiv2 - 1 + SMART_NORM_INFLUENCE */
	  /* First, we compute the root itself */
	  g = compute_f (d, u, hl);
	  g = lg2abs (g, add, scale);
	  fg = (unsigned char) g;
	  if (LIKELY(ih < Idiv2)) S[ih] = fg;    /* Right guard for the write */

	  /* The loop on the left side */
	  for (f1 = fg, f2 = 0, cpt = 0; ; ) {
	    --ih; hl -= 1.;
	    if (UNLIKELY(ih <= -Idiv2)) {        /* Left guard: in fact ==, not <= */
	      ASSERT(ih == -Idiv2);
	      sg[nsg].begin = sg[0].begin;
	      sg[nsg].f_begin = sg[0].f_begin;
	      break;
	    }
	    g = compute_f (d, u, hl);
	    g = lg2abs (g, add, scale);
	    f1 = (unsigned char) g;
	    if (LIKELY(ih < Idiv2)) S[ih] = f1;  /* Right guard for the write */
	    if (LIKELY(f1 == f2)) {
	      if (++cptf2id >= SMART_NORM_STABILITY) goto end_left;
	    } else {
	      cptf2id = 0;
	      f2 = f1;
	    }
	    if (UNLIKELY(++cpt >= SMART_NORM_INFLUENCE)) {
	    end_left:
	      sg[nsg].begin = ih;
	      sg[nsg].f_begin = g;
	      break;
	    }
	  }
	}
	/* If the left side is totally out of the interesting region ? -> next root */
	if (sg[nsg].begin > Idiv2 - 1) continue;

	/* Now, the right side. */
	ih = (ssize_t) hr; /* -Idiv2-SMART_NORM_INFLUENCE <= ih <= Idiv2-1+SMART_NORM_INFLUENCE */

	/* Immediate right guard for -Idiv2-SMART_NORM_INFLUENCE <= ih < Idiv2-1 */
	if (UNLIKELY(ih >= Idiv2 - 1)) {
	  ih = Idiv2 - 1;
	  hr = (double) ih;
	  g = compute_f (d, u, hr);
	  g = lg2abs (g, add, scale);
	} else {
	  /* The loop on the right side */
	  for (f1 = fg, f2 = 0, cpt = 0;;) {
	    ++ih; hr += 1.;
	    g = compute_f (d, u, hr);
	    g = lg2abs (g, add, scale);
	    f1 = (unsigned char) g;
	    if (LIKELY(ih > -Idiv2)) S[ih] = f1; /* Left guard for the write */
	    if (UNLIKELY(ih >= Idiv2 - 1)) {     /* Right guard: in fact ==, not >= */
	      ASSERT(ih == Idiv2 - 1);
	      break;
	    }
	    if (LIKELY (f1 == f2)) {
	      if (++cptf2id >= SMART_NORM_STABILITY) break;
	    } else {
	      cptf2id = 0;
	      f2 = f1;
	    }
	    if (++cpt >= SMART_NORM_INFLUENCE) break;
	  }
	  /* If the right side is totally out of the interesting region ? -> next root */
	  if (ih <= -Idiv2) continue;
	}
	sg[nsg].end = ih;
	sg[nsg].f_end = g;
	S[ih] = (unsigned char) g;
      }
      /* We have to do the possible fusions */
      while (nsg && sg[nsg - 1].end + 1 >= sg[nsg].begin) {
	if (sg[nsg - 1].begin > sg[nsg].begin) {
	  sg[nsg - 1].begin = sg[nsg].begin; sg[nsg - 1].f_begin = sg[nsg].f_begin; }
	if (sg[nsg - 1].end   < sg[nsg].end  ) {
	  sg[nsg - 1].end   = sg[nsg].end;   sg[nsg - 1].f_end   = sg[nsg].f_end;   }
	--nsg;
      }
      ++nsg;
    } /* End of the roots */

    /* The last sg end must be Idiv2 - 1 */
    if (LIKELY(sg[nsg - 1].end != Idiv2 - 1)) {
      /* We have to add a segment, except if the last ends at Idiv2 - 2 */
      g = compute_f (d, u, (double) (Idiv2 - 1));
      g = lg2abs (g, add, scale);
      S[Idiv2 - 1] = (unsigned char) g;
      if (LIKELY(sg[nsg - 1].end != Idiv2 - 2)) {
	sg[nsg].begin = Idiv2 - 1;
	sg[nsg].f_begin = g;
	++nsg;
      }
      sg[nsg - 1].end = Idiv2 - 1;
      sg[nsg - 1].f_end = g;
    }

#if 0
    for (size_t k = 0; k < nroots; k++)
      fprintf (stderr, "Racine %zu (derivate=%u): %e - %u\n", k, roots[k].derivate, roots[k].value * J, J);
    for (size_t k = 0; k < nsg; k++)
      fprintf (stderr, "Segments %zu: (%d, %e) - (%d, %e)\n", k, sg[k].begin, sg[k].f_begin, sg[k].end, sg[k].f_end);
#endif
    /* Here, we have to set all the unset S values by polygonal approximation on all
       sg[x].end...sg[x+1].begin */
    poly_approx_on_S (S, d, u, original_scale, scale, add, nsg, sg);
  } /* End of the line J; */

  /* Special ultra rare case. The correction of log2(F(0,0)) could be false,
     because the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) *beginS = GUARD;
} /* True end of the function (finally!) */

/* This function is used to initialize lognorms (=log2(F(i,j)*scale+GUARD)
   for the bucket_region S[] number J.
   It's a wrapper; except for trivial degree = 0, it extracts the interesting
   parameters of the complex structure si and calls the right function.

   - For degree 0, S[] is initialized by a memset: always exact.
   - A special ultra fast init function is used for degree = 1; it's could be
     considered as exact (the maximal error is always -/+ 1 on S[]).
   - For smart = 0 and others degrees, the exact values F(i,j) are computed with
     a fast log2. The maximal error is always -/+ 1 on S[]. It's obvious slow.

   - For smart != 0 and others degrees, cf las-config.h for the smart init algo.
   It's ~10 times faster than the exact init.

   This smart algo needs :
   -> The roots of d^2(F(i,1)/d(i)^2 must be in si->sides[side]->roots;
   -> si->sides[side]->nroots is the number of roots + 1;
   -> si->roots[si->sides[side]->nroots - 1] must be = 0.0 : it's a "pseudo" root
      in order to correct the neigborhood of F(0, const j).
   Of course, if si->sides[side]->nroots = 0, no correction is done.
*/
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info_ptr si, unsigned int side, unsigned int smart)
{
  unsigned int degree = si->sides[side]->fij->deg;
  switch (degree) {
  case 0 :
    memset (S, (int) (log2(1.+fabs(si->sides[side]->fijd[0])) * si->sides[side]->scale) + GUARD, 1U << LOG_BUCKET_REGION);
    break;
  case 1 :
    init_degree_one_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, si->sides[side]->fijd[0], si->sides[side]->fijd[1], si->sides[side]->cexp2);
    break;
  default:
    if (smart)
      init_smart_degree_X_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, degree, si->sides[side]->fijd, si->sides[side]->nroots, si->sides[side]->roots);
    else
      init_exact_degree_X_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, degree, si->sides[side]->fijd);
    break;
  }
}

/* return max |g(x)| for x in (0, s) where s can be negative,
   and g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double_poly_srcptr poly, double s)
{
  double_poly_t derivative;
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

/* Like get_maxnorm_aux(), but for interval [-s, s] */
static double
get_maxnorm_aux_pm (double_poly_srcptr poly, double s)
{
  double norm1 = get_maxnorm_aux(poly, s);
  double norm2 = get_maxnorm_aux(poly, -s);
  return (norm2 > norm1) ? norm2 : norm1;
}

/* returns the maximal norm of |F(x,y)| for -X <= x <= X, 0 <= y <= Y,
   where F(x,y) is a homogeneous polynomial of degree d.
   Let F(x,y) = f(x/y)*y^d, and F(x,y) = rev(F)(y,x).
   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(X,y) for -Y <= y <= Y (right-hand-side border, and the
       mirrored image of the left-hand-side border); We want the maximum of
       rev(F)(y,X) = rev(f)(y/X) * X^d for -Y <= j <= Y;
       this is rev(f)(t) * X^d for -Y/X <= t <= Y/X.
   (b) either on F(x,Y) for -X <= x <= X (top border)
       = f(x/Y) * Y^d; this is f(t) * Y^d for -X/Y <= t <= X/J.
   (d) or on F(x,0) for -X <= x <= X (lower border, on the abscissa), but this
       maximum is f[d]*X^d, and is attained in (a).
*/
double
get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y)
{
  const unsigned int d = src_poly->deg;
  double norm, max_norm;

  /* Make copy of polynomial as we need to revert the coefficients */
  double_poly_t poly;
  double_poly_init (poly, d);
  double_poly_set (poly, src_poly);

  /* (b) determine the maximum of |f(x)| * Y^d for -X/Y <= x <= X/Y */
  max_norm = get_maxnorm_aux_pm (poly, X/Y) * pow(Y, (double)d);

  /* (a) determine the maximum of |g(y)| for -1 <= y <= 1, with g(y) = F(s,y) */
  double_poly_revert(poly);
  norm = get_maxnorm_aux_pm (poly, Y/X) * pow(X, (double)d);
  if (norm > max_norm)
    max_norm = norm;

  double_poly_clear(poly);

  return max_norm;
}

/*
   Allocates:
     si->sides[side]->fij
     si->sides[side]->fijd
*/

void sieve_info_init_norm_data(sieve_info_ptr si)
{
  for (size_t side = si->cpoly->nb_polys; side--; ) {
    unsigned int d = si->cpoly->pols[side]->deg;
    sieve_side_info_ptr s = si->sides[side];
    mpz_poly_init (s->fij, d);
    s->fijd = (double *) malloc_aligned((d + 1) * sizeof(*(s->fijd)), 16);
    FATAL_ERROR_CHECK(s->fijd == NULL, "malloc failed");
    s->roots = (root_ptr) malloc_aligned((d * 4 + 1) * sizeof(*(s->roots)), 16);
    FATAL_ERROR_CHECK(s->roots == NULL, "malloc failed");
    s->nroots = 0;
  }
}

void sieve_info_clear_norm_data(sieve_info_ptr si)
{
  for(size_t side = si->cpoly->nb_polys ; side-- ; ) {
    sieve_side_info_ptr s = si->sides[side];
    mpz_poly_clear (s->fij);
    free(s->fijd);
    free(s->roots);
    s->nroots = 0;
  }
}

/* return largest possible J by simply bounding the Fij and Gij polynomials

  The image in the a,b-plane of the sieve region might be slanted at an angle
  against the abscissa, thus even though it has the correct skewness, it
  might include larger a,b-coordinates than a non-slanted image would.

  We compute the maximum norm that would occur if we had a perfect lattice
  basis in the sense that its image forms a rectangle A/2 <= a < A/2,
  0 <= b < B, with A/2/B = skew, and A*B = I*J*q (assuming J=I/2 here).
  Thus we have B = A/2/skew, A*A/2/skew = I*I/2*q, A = I*sqrt(q*skew).
  The optimal maximum norm is then the maximum of |F(a,b)| in this rectangle,
  divided by q on the special-q side.

  Then we reduce J so that the maximum norm in the image of the actual sieve
  region is no larger than this optimal maximum, times some constant fudge
  factor.
*/
static unsigned int
sieve_info_update_norm_data_Jmax (sieve_info_ptr si)
{
  // The following parameter controls the scaling on the norm.
  // Relevant values are between 1.0 and 3.0. A higher value means we
  // select higher values of J, and therefore we find more relations, but
  // this increases the time per relations.
  // The value 2.0 seems to be a good compromise. Setting 1.5 reduces the
  // time per relation and number of relations by about 1.5% on a typical
  // RSA704 benchmark.
  const double fudge_factor = 2.0;
  const double I = (double) (si->I);
  const double q = mpz_get_d(si->doing->p);
  const double skew = si->cpoly->skew;
  const double A = I*sqrt(q*skew);
  const double B = A/2./skew;
  double Jmax = I/2.;

  for (int side = 0; side < si->cpoly->nb_polys; side++)
    {
      sieve_side_info_ptr s = si->sides[side];
      mpz_poly_ptr ps = si->cpoly->pols[side];

      /* Compute the best possible maximum norm, i.e., assuming a nice
         rectangular sieve region in the a,b-plane */
      double_poly_t dpoly;
      double_poly_init (dpoly, ps->deg);
      double_poly_set_mpz_poly (dpoly, ps);
      double maxnorm = get_maxnorm_alg (dpoly, fudge_factor*A/2.,
              fudge_factor*B);
      double_poly_clear (dpoly);
      if (side == si->doing->side)
        maxnorm /= q;

      double_poly_t F;
      F->deg = ps->deg;
      F->coeff = s->fijd;

      double v = get_maxnorm_alg (F, I/2, Jmax);

      if (v > maxnorm)
        { /* use dichotomy to determine largest Jmax */
          double a, b, c;
          a = 0.0;
          b = Jmax;
          while (trunc (a) != trunc (b))
            {
              c = (a + b) * 0.5;
              v = get_maxnorm_alg (F, I/2, c);
              if (v < maxnorm)
                a = c;
              else
                b = c;
            }
          Jmax = trunc (a) + 1; /* +1 since we don't sieve for j = Jmax */
        }
    }

  return (unsigned int) Jmax;
}


/* return 0 if we should discard that special-q because the rounded
 * region in the (a,b)-plane is flat to the point of having height 0.
 * For diagnostic, we set si->J to the unrounded value (rounding would
 * give 0) and then we return "false".
 *
 * The current check for discarding is whether we do fill one bucket or
 * not. If we don't even achieve that, we should of course discard.
 *
 * Now for efficiency reasons, the ``minimum reasonable'' number of
 * buckets should be more than that.
 */
int sieve_info_adjust_IJ(sieve_info_ptr si, int nb_threads)/*{{{*/
{
    /* compare skewed max-norms: let u = [a0, b0] and v = [a1, b1],
       and u' = [a0, s*b0], v' = [a1, s*b1] where s is the skewness.
       Assume |u'| <= |v'|.
       We know from Gauss reduction that u' and v' form an angle of at
       least pi/3, thus since their determinant is q*s, we have
       |u'|^2 <= |u'| * |v'| <= q*s/sin(pi/3) = 2*q*s/sqrt(3).
       Define B := sqrt(2/sqrt(3)*q/s), then |a0| <= s*B and |b0| <= B.

       If we take J <= I/2*min(s*B/|a1|, B/|b1|), then for any j <= J:
       |a1|*J <= I/2*s*B and |b1|*J <= I/2*B, thus
       |a| = |a0*i+a1*j| <= s*B*I and |b| <= |b0*i+b1*j| <= B*I.
    */
    const double skewness = si->cpoly->skew;
    const int verbose = 0;
    if (verbose) {
        printf("# Called sieve_info_adjust_IJ((a0=%" PRId64 "; b0=%" PRId64
               "; a1=%" PRId64 "; b1=%" PRId64 "), p=%lu, skew=%f, nb_threads=%d)\n",
               si->qbasis.a0, si->qbasis.b0, si->qbasis.a1, si->qbasis.b1,
               mpz_get_ui(si->doing->p), skewness, nb_threads);
    }
    double maxab1, maxab0;
    maxab1 = si->qbasis.b1 * skewness;
    maxab1 = maxab1 * maxab1 + si->qbasis.a1 * si->qbasis.a1;
    maxab0 = si->qbasis.b0 * skewness;
    maxab0 = maxab0 * maxab0 + si->qbasis.a0 * si->qbasis.a0;
    if (maxab0 > maxab1) { /* exchange u and v, thus I and J */
        int64_t oa[2] = { si->qbasis.a0, si->qbasis.a1 };
        int64_t ob[2] = { si->qbasis.b0, si->qbasis.b1 };
        si->qbasis.a0 = oa[1]; si->qbasis.a1 = oa[0];
        si->qbasis.b0 = ob[1]; si->qbasis.b1 = ob[0];
    }
    maxab1 = MAX(labs(si->qbasis.a1), labs(si->qbasis.b1) * skewness);
    /* make sure J does not exceed I/2 */
    /* FIXME: We should not have to compute this B a second time. It
     * appears in sieve_info_init_norm_data already */
    double B = sqrt (2.0 * mpz_get_d(si->doing->p) / (skewness * sqrt (3.0)));
    if (maxab1 >= B * skewness)
        si->J = (uint32_t) (B * skewness / maxab1 * (double) (si->I >> 1));
    else
        si->J = si->I >> 1;

    /* Make sure the bucket region size divides the sieve region size,
       partly covered bucket regions may lead to problems when
       reconstructing p from half-empty buckets. */
    /* Compute number of i-lines per bucket region, must be integer */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->conf->logI);
    uint32_t i = 1U << (LOG_BUCKET_REGION - si->conf->logI);
    i *= nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */

    /* Bug 15617: if we round up, we are not true to our promises */
    uint32_t nJ = (si->J / i) * i; /* Round down to multiple of i */

    if (verbose) printf("# %s(): Final J=%" PRIu32 "\n", __func__, nJ);
    /* XXX No rounding if we intend to abort */
    if (nJ > 0) si->J = nJ;
    return nJ > 0;
}/*}}}*/


/* this function initializes the scaling factors and report bounds on the
   rational and algebraic sides */
/*
   Updates:
     si->sides[side]->logmax
     si->sides[side]->scale
     si->sides[side]->cexp2[]
     si->sides[side]->bound

*/

void
sieve_info_update_norm_data (sieve_info_ptr si, int nb_threads)
{
  int64_t H[4] = { si->qbasis.a0, si->qbasis.b0, si->qbasis.a1, si->qbasis.b1 };

  double step, begin;
  double r, maxlog2;
  double_poly_t poly;

  /* Update floating point version of both polynomials. They will be used in
   * get_maxnorm_alg(). */
  for (int side = 0; side < si->cpoly->nb_polys; side++) {
      sieve_side_info_ptr s = si->sides[side];
      mpz_poly_ptr ps = si->cpoly->pols[side];
      mpz_poly_homography (s->fij, ps, H);
      /* On the special-q side, divide all the coefficients of the
         transformed polynomial by q */
      if (si->conf->side == side) {
          for (int i = 0; i <= ps->deg; i++) {
              ASSERT_ALWAYS(mpz_divisible_p(s->fij->coeff[i], si->doing->p));
              mpz_divexact(s->fij->coeff[i], s->fij->coeff[i], si->doing->p);
          }
      }
      for (int k = 0; k <= ps->deg; k++)
          s->fijd[k] = mpz_get_d (s->fij->coeff[k]);
  }


  for (int side = 0; side < si->cpoly->nb_polys; ++side) {
    sieve_side_info_ptr sideptr = si->sides[side];

    /* Compute the roots of the polynomial F(i,1) and the roots of its
       inflexion points d^2(F(i,1))/d(i)^2.  These roots and 0.0 are need to
       be corrected the norms initialization on their neighboroods in
       init_smart_degree_X_norms_bucket_region_internal.  */

#ifdef SMART_NORM
    if (si->cpoly->pols[side]->deg >= 2)
      init_norms_roots (si, side);
#endif

    /* Compute the maximum norm of the polynomial over the sieve region.
       The polynomial coefficient in fijd are already divided by q
       on the special-q side. */
    poly->deg = si->cpoly->pols[side]->deg;
    poly->coeff = si->sides[side]->fijd;
    sideptr->logmax = log2(get_maxnorm_alg (poly, (double)si->I/2, (double)si->I/2));

    /* we know that |F(a,b)| < 2^(logmax) or |F(a,b)/q| < 2^(logmax)
       depending on sqside. 0 */
    /* we increase artificially 'logmax', to allow larger values of J */
    sideptr->logmax += 2.0;
    maxlog2 = sideptr->logmax;

    /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
       thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2.
       We require that scale is of the form (int) * 0.025, so that only a small
       number of different factor base slicings can occur. */
    sideptr->scale = (int)(((double) UCHAR_MAX - GUARD) / maxlog2 * 40.)*0.025;
    verbose_output_print (0, 1,
        "# Side %d: log2(maxnorm)=%1.2f scale=%1.2f, logbase=%1.6f",
        side, maxlog2, sideptr->scale, exp2 (1. / sideptr->scale));
    step = 1. / sideptr->scale;
    begin = -step * GUARD;
    for (unsigned int inc = 0; inc < 257; begin += step)
      sideptr->cexp2[inc++] = exp2(begin);
    /* we want to select relations with a cofactor of less than r bits */
    {
      double max_lambda = (maxlog2 - GUARD / sideptr->scale) /
        si->conf->sides[side]->lpb;
      double lambda = si->conf->sides[side]->lambda;
      if (lambda == 0) {
        r = MIN(si->conf->sides[side]->mfb, maxlog2 - GUARD / sideptr->scale);
      } else {
        r = MIN(lambda * (double) si->conf->sides[side]->lpb,
            maxlog2 - GUARD / sideptr->scale);
      }
      sideptr->bound = (unsigned char) (r * sideptr->scale + GUARD);
      verbose_output_print (0, 1, " bound=%u\n", sideptr->bound);
      if (lambda > max_lambda)
        verbose_output_print (0, 1, "# Warning, lambda>%.1f on side %d does "
            "not make sense (capped to limit)\n", max_lambda, side);
    }
  }

  /* improve bound on J if possible */
  unsigned int Jmax;
  Jmax = sieve_info_update_norm_data_Jmax (si);
  if (Jmax > si->J)
    {
      /* see sieve_info_adjust_IJ */
      ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->conf->logI);
      uint32_t i = 1U << (LOG_BUCKET_REGION - si->conf->logI);
      i *= nb_threads;
      si->J = (Jmax / i) * i; /* cannot be zero since the previous value
                                 of si->J was already a multiple of i,
                                 and this new value is larger */
    }
}
