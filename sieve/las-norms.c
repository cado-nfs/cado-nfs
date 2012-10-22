#define __STDC_LIMIT_MACROS
#include <string.h>
#include <limits.h>
#include <math.h>               /* ceil */
#ifdef SSE_NORM_INIT
#include <emmintrin.h>
#endif

#include "las-config.h"
#include "las-debug.h"
#include "las-norms.h"
#include "utils.h"
#include "mpz_poly.h"

/************************** sieve info stuff *********************************/

/* initialize array C[0..255]: C[i] is zero whenever the log-norm i
   is considered as a potential report, and 127 otherwise.
   The large prime bound is L = 2^l.
*/
static void
sieve_info_init_lognorm (unsigned char *C, unsigned char threshold,
                         unsigned long B MAYBE_UNUSED,
                         unsigned long l MAYBE_UNUSED,
                         double scale MAYBE_UNUSED)
{
  unsigned long k;

  for (k = 0; k < 256; k++)
    C[k] = (k <= threshold) ? 0 : 127;

#ifdef COFACTOR_TRICK
  {
    unsigned char k0, k1;
    double lost = 6.0; /* maximal number of bits lost due to prime powers */

    /* for L < R < B^2, a cofactor R cannot be L-smooth, since then it should
       have at least two prime factors in [B,L], and then should be >= B^2.
       Note:
       - the lognorms S[] are GUARD larger than the real values, thus the
         interval [L,R^2] corresponds to [k0+GUARD,k1+GUARD], but to take into
         account the roundoff errors we consider [k0+2*GUARD,k1]
       - 'lost' takes into account the error because we do not sieve
         prime powers; it applies only to the lower bound, and has to be scaled
       - the additional guard value takes into account additional errors
    */
    k0 = (unsigned char) (((double) l + lost) * scale) + 2 * GUARD + 2;
    k1 = (unsigned char) (2.0 * log2 ((double) B) * scale) - 2;
    for (k = k0; k <= k1 && k < 256; k++)
      C[k] = 127;
  }
#endif
}

/* {{{ initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
 * lognorms approximations for k bits of exponent + NORM_BITS-k bits
 * of mantissa.
 * Do the same for the algebraic side (with the corresponding bound for
 * the norms.
 */
void
init_norms (sieve_info_ptr si)
{
  for(int s = 0 ; s < 2 ; s++) {
      sieve_side_info_ptr sdata = si->sides[s];

      unsigned char *S = sdata->lognorm_table;

      int k = (int) ceil (log2 (sdata->logmax));
      int K = 1 << k;
      ASSERT_ALWAYS(NORM_BITS >= k);
      int l = NORM_BITS - k;
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
              /* Warning: since sdata->logmax does not usually correspond to
                 a power a two, and we consider full binades here, we have to
                 take care that values > sdata->logmax do not wrap around to 0 */
              norm = log2 (norm);
              if (norm >= sdata->logmax)
                  S[(i << l) + j] = 255;
              else
              {
                  norm = norm * sdata->scale;
                  S[(i << l) + j] = GUARD + (unsigned char) norm;
              }
          }
      }
  }
}

/* }}} */

/* {{{ initialize norms for bucket regions */
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void init_rat_norms_bucket_region(unsigned char *S,
                                 /* no condition array here ! */
                                 const int N,
                                 sieve_info_ptr si)
{
    sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];

    uint64_t mask = (1 << NORM_BITS) - 1;
    union { double z; uint64_t x; } zx[1];
    int l = NORM_BITS - (int) ceil(log2(rat->logmax));
    /* The degree d is always 1 */

    const unsigned char * L = rat->lognorm_table;

    int32_t I = si->I;
    int halfI = I / 2;
    int logI = si->logI;

    int i MAYBE_UNUSED;

    unsigned int j0 = N << (LOG_BUCKET_REGION - logI);
    unsigned int j1 = j0 + (1 << (LOG_BUCKET_REGION - logI));

    /* parity means:
     * 1 : (i,j) actually represents (real_i,real_j) = (2i+1,2j)
     * 2 : (i,j) actually represents (real_i,real_j) = (2i,2j+1)
     * 3 : (i,j) actually represents (real_i,real_j) = (2i+1,2j+1)
     */
    /* G_q(ri,rj) = g1 * (a0*ri+a1*rj) + g0 * (b0*ri+b1*rj)
       = (g1*a0+g0*b0) * ri + (g1*a1+g0*b1) * rj
       = gi * i + gj * rj
       = 2 * gi * i + 2 * gj * j + some offset gc;
     */

    double u[2] = {
        si->sides[RATIONAL_SIDE]->fijd[0],  // gj
        si->sides[RATIONAL_SIDE]->fijd[1],  // gi
    };

    /* Note that if ratq holds, then fijd has q divided out already */

    /* bucket_region is a multiple of I. Let's find the starting j
     * corresponding to N and the last j.
     */
#ifdef UGLY_HACK
    memset(S, 255, 1 << LOG_BUCKET_REGION);
#else
    unsigned int j = j0;
    /* if j = 0, it will be the first value */
    if (j == 0) {
	// compute only the norm for i = 1. Everybody else is 255.
        memset(S, 255, I);
	double norm = (log2(fabs(u[1]))) * rat->scale;
	S[halfI + 1] = GUARD + (unsigned char) (norm);
        S+=I;
	++j;
    }
    for( ; j < j1 ; j++) {
	zx->z = u[0] * (double) j - u[1] * (double) halfI;
	__asm__("### Begin rational norm loop\n");
/* LAZY_NORMS does not seem to be so interesting on the rational side,
 * since wrong values generate a lot of additional work thereafter */
#if 0   
        uint64_t y;
        unsigned char n, oldn = 0;
        const int normstride=128; // must be a power of 2 dividing I.
        double gii = gi * normstride;
        for (i = 0; i < (int) si->I; i+=normstride) {
            y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
            n = rat->S[y & mask];
            ASSERT (n > 0);
            if (i > 0 && oldn != n) {
                // the lognorm has changed: recompute more precisely the
                // previous stride.
                zx->z -= gii;
                S -= normstride;
                for (int ii = 0; ii < normstride; ++ii) {
                    y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
                    n = rat->S[y & mask];
                    ASSERT (n > 0);
                    zx->z += gi;
                    *S++ = n;
                }
                y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
                n = rat->S[y & mask];
                ASSERT (n > 0);
            }
            for (int ii = 0; ii < normstride; ++ii)
                *S++ = n; 
            zx->z += gii;
            oldn = n;
        }
#else
#ifndef SSE_NORM_INIT
	uint64_t y;
	unsigned char n;
	for (i = 0; i < (int) I; i++) {
	    /* the double precision number 1.0 has high bit 0 (sign),
	       then 11-bit biased exponent 1023, and 52-bit mantissa 0 */
	    /* the magic constant here is simply 1023*2^52, where
	       1023 is the exponent bias in binary64 */
	    y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
	    n = L[y & mask];
	    ASSERT(n > 0);
	    *S++ = n;
	    zx->z += u[1];
	}
#else
	union {
	    __v2df dble;
	    __v2di intg;
	    struct {
		uint64_t y0;
		uint64_t y1;
	    } intpair;
	} y_vec;
        if ((j & 1) == 1) {
            __v2df gi_vec = { 2*u[1], 2*u[1] };
            __v2di mask_vec = { mask, mask };
            __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
                (uint64_t) 0x3FF0000000000000 };
           
            // in spite of the appearance, only the low word gives the shift
            // count. The high word is ignored.
            __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
            __v2df z_vec = { zx->z, zx->z+u[1] };
             
            for (i = 0; i < halfI; ++i) {
                y_vec.dble = z_vec;
                y_vec.intg -= cst_vec;
                y_vec.intg = _mm_srl_epi64(y_vec.intg, shift_value);
                y_vec.intg &= mask_vec;
                //            *S++ = L[((uint32_t *)(&y_vec.intg))[0]];
                //             *S++ = L[((uint32_t *)(&y_vec.intg))[2]];
                *S++ = L[y_vec.intpair.y0];
                *S++ = L[y_vec.intpair.y1];
                z_vec += gi_vec;
            }
        } else { // skip when i and j are both even
            __v2df gi_vec = { 4*u[1], 4*u[1] };
            __v2di mask_vec = { mask, mask };
            __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
                (uint64_t) 0x3FF0000000000000 };

            // in spite of the appearance, only the low word gives the shift
            // count. The high word is ignored.
            __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
            __v2df z_vec = { zx->z + u[1], zx->z+3*u[1] };

            for (i = 0; i < halfI; i+=2) {
                y_vec.dble = z_vec;
                y_vec.intg -= cst_vec;
                y_vec.intg = _mm_srl_epi64(y_vec.intg, shift_value);
                y_vec.intg &= mask_vec;
                //            *S++ = L[((uint32_t *)(&y_vec.intg))[0]];
                //             *S++ = L[((uint32_t *)(&y_vec.intg))[2]];
                *S++ = 255;
                *S++ = L[y_vec.intpair.y0];
                *S++ = 255;
                *S++ = L[y_vec.intpair.y1];
                z_vec += gi_vec;
            }
        }
#endif
#endif
	__asm__("### End rational norm loop\n");
    }
#endif
}


/* {{{ some utility stuff */
#ifdef SSE_NORM_INIT
static inline void init_fpoly_v2df(__v2df * F, const double *f, const int deg)
{
    int i;
    for (i = 0; i <= deg; ++i) {
	__v2df tmp = { f[i], f[i] };
	F[i] = tmp;
    }
}

static inline __v2df fpoly_eval_v2df(const __v2df * f, const __v2df x, int d)
{
    __v2df r;
    r = f[d--];
    for (; d >= 0; --d)
	r = r * x + f[d];
    return r;
}
#endif

static inline void fpoly_scale(double * u, const double * t, int d, double h)
{
    u[d] = t[d];        /* This one never changes. */
    double hpow = h;
    u[d-1] = t[d-1] * hpow;
    for (int k = 2; k <= d; k++) {
        hpow *= h;
        u[d - k] = t[d - k] * hpow;
    }
}
/* }}} */

/* Initialize lognorms on the algebraic side for the bucket
 * number N.
 * Only the survivors of the other sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 *
 * ``the other side'' is checked via xB[*xS] > 127.
 */
int
init_alg_norms_bucket_region(unsigned char *S,
			     const unsigned char *xS MAYBE_UNUSED, const int N,
			     sieve_info_ptr si)
{
    sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];

    uint64_t mask = (1 << NORM_BITS) - 1;
    union { double z; uint64_t x; } zx[1];
    int l = NORM_BITS - (int) ceil(log2(alg->logmax));
    int d = si->cpoly->alg->degree;

    const unsigned char * L = alg->lognorm_table;

#ifndef LAZY_NORMS
    /* Used for checking the bound on the other side */
    const unsigned char * xB = si->sides[RATIONAL_SIDE]->Bound;
    const unsigned char * xS0 MAYBE_UNUSED = xS;
#endif

    int32_t I = si->I;
    double halfI = I / 2;
    int logI = si->logI;
    unsigned int j0 = N << (LOG_BUCKET_REGION - logI);
    unsigned int j1 = j0 + (1 << (LOG_BUCKET_REGION - logI));

    int report = 0;

    double u[d+1];

    /* bucket_region is a multiple of I. Let's find the starting j
     * corresponding to N and the last j.
     */
    for (unsigned int j = j0 ; j < j1 ; j++) {
#ifdef LAZY_NORMS
        unsigned int jmask = VERT_NORM_STRIDE - 1;
        if (j > j0 && ((j & jmask) != 0)) {
            // copy norms from the row below.
            unsigned char *old_S;
            old_S = S - si->I;
            for (unsigned int ii = 0; ii < si->I; ++ii) {
                *S++ = *old_S++;
            }
            continue;
        }
#endif
        /* scale by rj^(d-k) the coefficients of fij */
        fpoly_scale(u, alg->fijd, d, j);

#ifdef LAZY_NORMS
        const int normstride=NORM_STRIDE;
        const double fpnormstride=(double)NORM_STRIDE;
        for (int i = -halfI; i < halfI; i += fpnormstride) {
            unsigned char n;
            zx->z = fpoly_eval (u, d, i);
            /* 4607182418800017408 = 1023*2^52 */
            uint64_t y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
            report++;
            n = L[y & mask];
            ASSERT (n > 0);
            for (int ii = 0; ii < normstride; ++ii)
                *S++ = n;
        }
#else  /* full norm computation */


#ifdef SSE_NORM_INIT
	__v2df u_vec[d+1];
	init_fpoly_v2df(u_vec, u, d);
	double previous_ri = 0.0;
	unsigned char * previous_S= NULL;
	int cpt = 0;
#endif
	/* now compute norms */
	for (double ri = -halfI ; ri < halfI; ri += 1.0) {
#ifdef TRACE_K
	    if (trace_on_spot_Nx(N, xS - xS0)) {
		fprintf(stderr, "# init_alg_norms_bucket_region: "
			"xS[%zd] = %u, -log(prob) = %u\n",
			xS - xS0,
			(unsigned int) *xS,
			(unsigned int) xB[*xS]);
	    }
#endif
	    if (xB[*xS] < 127) {
//  The SSE version seems to be slower on the algebraic side...
//  Let's forget it for the moment.
//  TODO: try with single precision.
//#ifndef SSE_NORM_INIT
#if 1
		zx->z = fpoly_eval(u, d, ri);
		/* the constant is 1023*2^52 */
		uint64_t y = (zx->x - UINT64_C(0x3FF0000000000000)) >> (52 - l);
		report++;
		unsigned char n = L[y & mask];
		ASSERT(n > 0);
		*S = n;
#else
		__v2di mask_vec = { mask, mask };
		__v2di cst_vec = {
                    UINT64_C(0x3FF0000000000000),
		    UINT64_C(0x3FF0000000000000),
		};
		__v2di shift_value =
		    { (uint64_t) (52 - l), (uint64_t) (52 - l) };
		report++;
		if (cpt == 0) {
		    previous_ri = ri;
		    previous_S = S;
		    cpt++;
		} else {
		    cpt--;
		    __v2df i_vec = { previous_ri, ri };
		    union {
			__v2df dble;
			__v2di intg;
			struct {
			    uint64_t y0;
			    uint64_t y1;
			} intpair;
		    } fi_vec;
		    fi_vec.dble = fpoly_eval_v2df(u_vec, i_vec, d);
		    fi_vec.intg -= cst_vec;
		    fi_vec.intg = _mm_srl_epi64(fi_vec.intg, shift_value);
		    fi_vec.intg &= mask_vec;
		    *previous_S = L[fi_vec.intpair.y0];
		    *S = L[fi_vec.intpair.y1];
		}
#endif
	    } else {
		*S = UCHAR_MAX;
	    }
	    S++;
	    xS++;
	}
#ifdef SSE_NORM_INIT
	if (cpt == 1) {		// odd number of computations
	    zx->z = fpoly_eval(u, d, previous_ri);
	    uint64_t y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
	    *previous_S = L[y & mask];
	}
#endif
#endif  /* LAZY_NORMS */
    }
    return report;
}

/* }}} */

/* return max |g(x)| for 0 <= x <= s,
   where g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double *g, const unsigned int d, double s)
{
  unsigned int k, l, sign_change, new_sign_change;
  double **dg;    /* derivatives of g */
  double a, va, b, vb;
  double *roots, gmax;

  dg = (double**) malloc (d * sizeof (double*));
  FATAL_ERROR_CHECK(dg == NULL, "malloc failed");
  dg[0] = g;
  for (k = 1; k < d; k++) { /* dg[k] is the k-th derivative, thus has
                             degree d-k, i.e., d-k+1 coefficients */
    dg[k] = (double*) malloc ((d - k + 1) * sizeof (double));
    FATAL_ERROR_CHECK(dg[k] == NULL, "malloc failed");
  }
  roots = (double*) malloc (d * sizeof (double));
  FATAL_ERROR_CHECK(roots == NULL, "malloc failed");
  for (k = 1; k < d; k++)
    for (l = 0; l <= d - k; l++)
      dg[k][l] = (l + 1) * dg[k - 1][l + 1];
  /* now dg[d-1][0]+x*dg[d-1][1] is the (d-1)-th derivative: it can have at
     most one sign change, iff dg[d-1][0] and dg[d-1][0]+dg[d-1][1] have
     different signs */
  if (dg[d-1][0] * (dg[d-1][0] + dg[d-1][1]) < 0)
    {
      sign_change = 1;
      roots[0] = - dg[d-1][0] / dg[d-1][1]; /* root of (d-1)-th derivative */
    }
  else
    sign_change = 0;
  roots[sign_change] = s; /* end of interval */
  for (k = d - 1; k-- > 1;)
    {
      /* invariant: sign_change is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_change-1], and roots[sign_change] = s. */
      a = 0.0;
      va = dg[k][0]; /* value of dg[k] at x=0 */
      new_sign_change = 0;
      for (l = 0; l <= sign_change; l++)
        {
          b = roots[l]; /* root of dg[k+1], or end of interval */
          vb = fpoly_eval (dg[k], d - k, b);
          if (va * vb < 0) /* root in interval */
            roots[new_sign_change++] = fpoly_dichotomy (dg[k], d - k,
                                                        a, b, va, 20);
          a = b;
          va = vb;
        }
      roots[new_sign_change] = s; /* end of interval */
      sign_change = new_sign_change;
    }
  /* now all extrema of g are 0, roots[0], ..., roots[sign_change] = s */
  gmax = fabs (g[0]);
  for (k = 0; k <= sign_change; k++)
    {
      va = fabs (fpoly_eval (g, d, roots[k]));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  for (k = 1; k < d; k++)
    free (dg[k]);
  free (dg);
  return gmax;
}

/* returns the maximal value of log2 |F(a,b)/q| for
   a = a0 * i + a1 * j, b = b0 * i + b1 * j and q >= q0,
   -I/2 <= i <= I/2, 0 <= j <= I/2*min(s*B/|a1|,B/|b1|)
   where B >= sqrt(2*q/s/sqrt(3)) for all special-q in the current range
   (s is the skewness, and B = si->B, see lattice.tex).

   Since |a0| <= s*B and |b0| <= B, then
   |a0 * i + a1 * j| <= s*B*I and |b0 * i + b1 * j| <= B*I,
   thus it suffices to compute M = max |F(x,y)| in the rectangle
   -s <= x <= s, 0 <= y <= 1, and to multiply M by (B*I)^deg(F).

   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(s,y) for 0 <= y <= 1
   (b) either on F(x,1) for -s <= x <= s
   (c) either on F(-s,y) for 0 <= y <= 1
   (d) or on F(x,0) for -s <= x <= s, but this maximum is f[d]*s^d,
       and is attained in (a) or (c).
*/
static double
get_maxnorm (cado_poly cpoly, sieve_info_ptr si, double q0d)
{
  unsigned int d = cpoly->alg->degree, k;
  double *fd; /* double-precision coefficients of f */
  double norm, max_norm, pows, tmp;

  fd = (double*) malloc ((d + 1) * sizeof (double));
  FATAL_ERROR_CHECK(fd == NULL, "malloc failed");
  for (k = 0; k <= d; k++)
    fd[k] = mpz_get_d (cpoly->alg->f[k]);

  /* (b1) determine the maximum of |f(x)| for 0 <= x <= s */
  max_norm = get_maxnorm_aux (fd, d, cpoly->skew);

  /* (b2) determine the maximum of |f(-x)| for 0 <= x <= s */
  norm = get_maxnorm_aux (fd, d, -cpoly->skew);
  if (norm > max_norm)
    max_norm = norm;

  for (pows = 1.0, k = 0; k <= d; k++)
    {
      fd[k] *= pows;
      pows *= cpoly->skew;
    }
  /* swap coefficients; if d is odd, we need to go up to k = floor(d/2) */
  for (k = 0; k <= d / 2; k++)
    {
      tmp = fd[k];
      fd[k] = fd[d - k];
      fd[d - k] = tmp;
    }

  /* (a) determine the maximum of |g(y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, 1.0);
  if (norm > max_norm)
    max_norm = norm;

  /* (c) determine the maximum of |g(-y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, -1.0);
  if (norm > max_norm)
    max_norm = norm;

  free (fd);

  /* multiply by (B*I)^d and divide by q0 if sieving on alg side */
  tmp = max_norm * pow (si->B * (double) si->I, (double) d);
  if (!si->ratq)
      tmp /= q0d;
  return log2(tmp);
}

void sieve_info_init_norm_data(sieve_info_ptr si, mpz_srcptr q0)
{
    for (int side = 0; side < 2; side++) {
        int d = si->cpoly->pols[side]->degree;
        si->sides[side]->fij = (mpz_t *) malloc((d + 1) * sizeof(mpz_t));
        FATAL_ERROR_CHECK(si->sides[side]->fij == NULL, "malloc failed");
        si->sides[side]->fijd = (double *) malloc((d + 1) * sizeof(double));
        FATAL_ERROR_CHECK(si->sides[side]->fijd == NULL, "malloc failed");
        for (int k = 0; k <= d; k++)
            mpz_init(si->sides[side]->fij[k]);
    }



    double q0d = mpz_get_d(q0);

  double r, scale;
  unsigned char alg_bound, rat_bound;
  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];

  /* initialize bounds for the norm computation, see lattice.tex */
  si->B = sqrt (2.0 * q0d / (si->cpoly->skew * sqrt (3.0)));
  alg->logmax = get_maxnorm (si->cpoly, si, q0d); /* log2(max norm) */

  /* We want some margin, (see below), so that we can set 255 to discard
   * non-survivors.*/
  r = si->cpoly->alg->lambda * (double) si->cpoly->alg->lpb;
  scale = alg->logmax + si->cpoly->alg->lambda * (double) si->cpoly->alg->lpb;

  fprintf (si->output, "# Alg. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           alg->logmax, exp2 (scale / LOG_MAX));
  // second guard, due to the 255 trick!
  alg->scale = scale = (LOG_MAX - GUARD) / scale;
  alg_bound = (unsigned char) (r *  scale) + GUARD;
  fprintf (si->output, " bound=%u\n", alg_bound);
  sieve_info_init_lognorm (alg->Bound, alg_bound, si->cpoly->alg->lim, si->cpoly->alg->lpb,
                           scale);

  /* similar bound on the rational size: |a| <= s*I*B and |b| <= I*B */
  double x = fabs (mpz_get_d (si->cpoly->rat->f[1])) * si->cpoly->skew
        + fabs (mpz_get_d (si->cpoly->rat->f[0]));
  x *= si->B * (double) si->I;
  if (si->ratq)
      x /= q0d;
  rat->logmax = log2 (x);

  scale = rat->logmax + si->cpoly->rat->lambda * (double) si->cpoly->rat->lpb;
  fprintf (si->output, "# Rat. side: log2(maxnorm)=%1.2f logbase=%1.6f",
          rat->logmax, exp2 (scale / LOG_MAX ));
  /* we subtract again GUARD to avoid that non-reports overlap the report
     region due to roundoff errors */
  rat->scale = scale = (LOG_MAX - GUARD) / scale;
  rat_bound = (unsigned char) (r *  scale) + GUARD;
  fprintf (si->output, " bound=%u\n", rat_bound);
  sieve_info_init_lognorm (rat->Bound, rat_bound, si->cpoly->rat->lim, si->cpoly->rat->lpb,
                           rat->scale);
}


void sieve_info_clear_norm_data(sieve_info_ptr si)
{
    for (int side = 0; side < 2; side++) {
        sieve_side_info_ptr s = si->sides[side];
        cado_poly_side_ptr ps = si->cpoly->pols[side];
        for (int k = 0; k <= ps->degree; k++)
            mpz_clear(s->fij[k]);
        free(s->fij);
        free(s->fijd);
    }
}

void sieve_info_update_norm_data(sieve_info_ptr si)
{
    int32_t H[4] = { si->a0, si->b0, si->a1, si->b1 };
    /* Update floating point version of algebraic poly (do both, while
     * we're at it...) */
    for (int side = 0; side < 2; side++) {
        sieve_side_info_ptr s = si->sides[side];
        cado_poly_side_ptr ps = si->cpoly->pols[side];
        mp_poly_homography(s->fij, ps->f, ps->degree, H);
        double invq = 1.0;
        if (si->ratq == (side == RATIONAL_SIDE))
            invq /= mpz_get_d(si->q);
        for (int k = 0; k <= ps->degree; k++)
            s->fijd[k] = mpz_get_d(s->fij[k]) * invq;
    }
}
