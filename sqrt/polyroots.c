/*
 * This file has been extracted from Jason Papadopoulos's msieve-1.42,
 * for inclusion in cado-nfs. It has been very marginally edited (notably
 * making it C99 by using uint32_t and complex type, thereby saving some
 * work).
 *
 * The original version carries the following note:
 *
 */

/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 4/3/09
--------------------------------------------------------------------*/
#include "cado.h"
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdlib.h>

#define MAX_ROOTFINDER_DEGREE   10

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define xxxMAIN

#ifdef  MAIN
#include <stdlib.h>
#include <stdio.h>
#endif

#include "polyroots.h"

/* This rootfinder uses a simplified version of the all-complex
   Jenkins-Traub algorithm.
   
   [ET for cado-nfs : the note below has been dropped because we dropped
   higher precision stuff --]
   
   then switches to Newton's method
   in extended precision to polish the roots that are found. In
   practice we cannot use only Newton's method, because this
   behaves essentially randomly unless given a very good initial
   root approximation. Jenkins-Traub is much more complex, and 
   converges just about unconditionally.

   The J-T code was kindly converted by Brian Gladman from the
   awful Fortran code of Algorithm 419 from comm. acm, vol. 15, 
   no. 0. I've removed a lot of the gold-plating of the full 
   implementation that is unlikely to be necessary for NFS 
   polynomials, and converted the paired-array-of-doubles code
   to much simpler array-of-complex code */

/* Jenkins-Traub rootfinder */

static const double are = DBL_EPSILON;
static const double mre = 2.0 * M_SQRT2 * DBL_EPSILON;

typedef struct {
	double complex poly[MAX_ROOTFINDER_DEGREE + 1]; 
	double complex poly_aux[MAX_ROOTFINDER_DEGREE + 1];
	double complex hpoly[MAX_ROOTFINDER_DEGREE + 1]; 
	double complex hpoly_aux[MAX_ROOTFINDER_DEGREE + 1];

	double complex angle, angle_inc;
	uint32_t hpoly_root_found;
	uint32_t degree;
} jt_t;

/*-----------------------------------------------------------------------*/
static double cauchy_bound(uint32_t n, const double complex p[]) {

	/* computes a lower bound on the moduli of the zeros of a
	   polynomial p(x). norms(x) is a polynomial whose i_th coefficient
	   is the modulus of the i_th coeffcient of p(x) but
	   whose constant term is negative. The lower bound is the 
	   (unique) positive root x of norms(x) */

	double x, xmax, f, dx, df;
	double norms[MAX_ROOTFINDER_DEGREE + 1];
	uint32_t i;

	for (i = 0; i < n; i++)
		norms[i] = cabs(p[i]);
	norms[i] = -cabs(p[i]);

	/* compute upper estimate of bound: assume all the
	   middle terms of norms(x) are zero */

	xmax = exp((log(-norms[n]) - log(norms[0])) / (double) n);

	/* if ignoring the nonlinear terms of norms(x) produces
	   a smaller root, use that instead */

	if (norms[n - 1] != 0.0) {
		x = -norms[n] / norms[n - 1];
		xmax = MIN(x, xmax);
	}

	/* chop the interval (0, x) until until x is about
	   to make norms(x) change sign */

	do {
		x = xmax;
		xmax = 0.1 * x;

		f = norms[0];
		for (i = 1; i <= n; i++)
			f = f * xmax + norms[i];
	} while (f > 0.0);

	/* do newton iteration until x converges to two decimal places */

	dx = x;
	while (fabs(dx / x) > 0.005) {
		df = 0;
		f = norms[0];
		for (i = 1; i <= n; i++) {
			df = df * x + f;
			f = f * x + norms[i];
		}
		dx = f / df;
		x -= dx;
	}

	return x;
}

/*-----------------------------------------------------------------------*/
static double complex poly_val(uint32_t n, double complex s, 
			const double complex p[], double complex q[]) {

	/* evaluates a polynomial p at s by the horner 
	   recurrence, placing the partial sums in q and 
	   returning the computed value */

	double complex pv;
	uint32_t i;

	pv = q[0] = p[0];
	for (i = 1; i <= n; i++)
		pv = q[i] = p[i] + pv * s;

	return pv;
}

/*-----------------------------------------------------------------------*/
static void next_hpoly(double complex correction, jt_t * w) {

	/* calculates the next shifted h polynomial */

	uint32_t i;
	double complex *poly_aux = w->poly_aux;
	double complex *hpoly_aux = w->hpoly_aux;
	double complex *hpoly = w->hpoly;

	if (w->hpoly_root_found == 0) {
		hpoly[0] = poly_aux[0];
		for (i = 1; i < w->degree; i++) {
			hpoly[i] = poly_aux[i] + hpoly_aux[i - 1] * correction;
		}
	}
	else {
		/* we are essentially at a root of h(x); remove 
		   it by deflating the polynomial. Calling code 
		   always expects h(x) to have the same degree, 
		   so the high-order coefficient becomes zero */

		hpoly[0] = 0;
		for (i = 1; i < w->degree; i++)
			hpoly[i] = hpoly_aux[i - 1];
	}
}

/*-----------------------------------------------------------------------*/
static double complex next_correction(double complex pval, double complex curr_root,
				jt_t * w) {

	/* computes -pval / hpoly(curr_root)
	   sets flag to true if hval is essentially zero. */

	double complex *hpoly = w->hpoly;
	double complex *hpoly_aux = w->hpoly_aux;
	double complex hval = poly_val(w->degree - 1, curr_root, hpoly, hpoly_aux);

	if (cabs(hval) <= 10.0 * DBL_EPSILON * cabs(hpoly[w->degree - 1])) {
		w->hpoly_root_found = 1;
		return 0;
	}
	else {
		w->hpoly_root_found = 0;
		return -pval/hval;
	}
}

/*-----------------------------------------------------------------------*/
#define STAGE3_ITER 10

static uint32_t stage3(double complex *root, jt_t *w) {

	/* carries out the third stage iteration,
	   returns 1 if iteration converges */

	double mp, ms, tp;
	uint32_t i, j;
	double complex pval;
	double complex correction;
	double complex curr_root = *root;
	double complex *poly = w->poly;
	double complex *poly_aux = w->poly_aux;

	for (i = 0; i < STAGE3_ITER; i++) {

		/* evaluate poly at current root value */

		pval = poly_val(w->degree, curr_root, poly, poly_aux);

		/* calculate bound on the error in evaluating the polynomial 
		   by the horner recurrence */

		mp = cabs(pval);
		ms = cabs(curr_root);
		tp = cabs(poly_aux[0]) * mre / (DBL_EPSILON + mre);
		for (j = 0; j <= w->degree; j++)
			tp = tp * ms + cabs(poly_aux[j]);
		tp = tp * (DBL_EPSILON + mre) - mp * mre;

		if (mp <= 20.0 * tp) {
			/* polynomial value is smaller in value 
			   than a bound on the error in evaluating p, 
			   terminate the iteration */
			*root = curr_root;
			return 1;
		}

		/* calculate next h polynomial */

		correction = next_correction(pval, curr_root, w);
		next_hpoly(correction, w);

		/* use the next h polynomial to calculate the next
		   root estimate, using the current root estimate */

		correction = next_correction(pval, curr_root, w);
		curr_root = curr_root+correction;
	}

	return 0;
}

/*-----------------------------------------------------------------------*/
static uint32_t stage2(uint32_t stage2_iter, double complex *root, jt_t *w) {

	uint32_t i;
	double complex curr_root;
	double complex correction;
	double complex pval;

	/* calculate first correction */

	curr_root = *root;
	pval = poly_val(w->degree, curr_root, w->poly, w->poly_aux);
	correction = next_correction(pval, curr_root, w);

	for (i = 0; i < stage2_iter; i++) {

		/* compute next h polynomial and new correction;
		   note that the fixed-shift iteration never changes
		   the value of curr_root, only the h polynomial */

		next_hpoly(correction, w);
		correction = next_correction(pval, curr_root, w);

		if (w->hpoly_root_found == 1)
			break;
	}

	/* attempt stage 3 with the final h polynomial and
	   final correction */

	*root = curr_root+correction;
	return stage3(root, w);
}

/*-----------------------------------------------------------------------*/
#define STAGE1_ITER 5

static void stage1(uint32_t n, double complex p[], double complex h[]) {

	uint32_t i, j;

	/* the initial h polynomial is a scaled version of the
	   derivative of the input polynomial p(x) */

	for (i = 0; i < n; i++)
		h[i] = p[i] * (double) (n - i) / n;

	/* compute a series of no-shift h polynomials */

	for (i = 0; i < STAGE1_ITER; i++) {

		/* if the constant term is essentially zero, 
		   shift the h coefficients */

		if (cabs(h[n-1]) <= 10.0 * DBL_EPSILON * cabs(p[n-1])) {

			for (j = n - 1; j; j--)
				h[j] = h[j-1];
			h[j] = 0;
		}
		else {

			double complex tmp = -p[n]/h[n-1];
			for (j = n - 1; j; j--)
				h[j] = p[j] + h[j-1] * tmp;
			h[j] = p[0];
		}
	}
}

/*-----------------------------------------------------------------------*/
static int find_one_root(double complex *root, jt_t *w) {

	uint32_t i, j, k;
	double bound;
	double complex hpoly_start[MAX_ROOTFINDER_DEGREE + 1];

	/* find linear roots immediately */

	if (w->degree <= 1) {
		*root = -w->poly[1]/w->poly[0];
		return 1;
	}

	/* calculate a lower bound on the modulus of the zeros */

	bound = cauchy_bound(w->degree, w->poly);

	/* stage 1 sets up the initial h polynomial only */

	stage1(w->degree, w->poly, hpoly_start);

	/* try the fixed-shift sequence twice */

	for (i = 0; i < 2; i++) {

		/* inner loop to select a shift */

		for (j = 1; j < 10; j++) {

			/* start point is chosen with modulus 'bound'
			   and a pseudo-random angle. In practice
			   we don't want to repeat previous work,
			   so the starting angle is rotated a fixed 
			   amount (94 degrees) from the previous 
			   start point */

			w->angle = w->angle*w->angle_inc;
			*root = bound*w->angle;

			/* do the second stage, with a varying
			   number of iterations.
			   
			   Note that every starting point uses the same
			   h polynomial. This is a change from all other
			   cpoly() versions, including the original 1972 
			   fortran, which uses a global h array that is
			   not reinitialized when a new start point is
			   chosen (I think erroneously) */

			for (k = 0; k < w->degree; k++)
				w->hpoly[k] = hpoly_start[k];

			if (stage2(10 * j, root, w) == 1)
				return 1;
		}
	}

	return 0;
}

/*-----------------------------------------------------------------------*/
static uint32_t jenkins_traub(double complex poly[], 
			uint32_t degree, double complex roots[]) {

	/* main Jenkins-Traub driver; returns number 
	   of roots found */

	uint32_t i; 
	uint32_t roots_found;
	jt_t w;

	/* remove any zeros at the origin */

	for (i = degree, roots_found = 0; i; i--, roots_found++) {
		if (poly[i] != 0.0)
			break;

		roots[roots_found] = 0;
	}
	w.degree = i;

	/* initialize */

	for (i = 0; i <= w.degree; i++)
		w.poly[i] = poly[i];

	w.angle = M_SQRT1_2-I*M_SQRT1_2;
	w.angle_inc = cexp(I * 94 * M_PI / 180);

	/* loop to find roots */

	for (; roots_found < degree; roots_found++) {

		if (find_one_root(roots + roots_found, &w) == 0)
			break;

		/* deflate the polynomial */

		(w.degree)--;
		for (i = 0; i <= w.degree; i++)
			w.poly[i] = w.poly_aux[i];
	}

	return roots_found;
}
/*-----------------------------------------------------------------------*/

#define NEWTON_ITER 10

static uint32_t polish_root(long double complex *poly, uint32_t degree,
			long double complex x, long double complex *root,
			double eps)
{
	uint32_t i = 0;

	for (i = 0; i < NEWTON_ITER; i++) {
		uint32_t j = degree;
		long double complex f = poly[j];
		long double complex df = 0;
		long double complex dx;
		long double abs_x, abs_dx;

		for (j--; (int32_t)j >= 0; j--) {
			df = x * df + f;
			f = x * f + poly[j];
		}
		dx = f / df;
		x = x - dx;

		abs_x = cabsl(x);
		abs_dx = cabsl(dx);

		if (abs_dx <= eps * abs_x) {
                    *root = x;
                    return 0;
                }
	}

	*root = x;
	return 1;
}
/*------------------------------------------------------------------*/
uint32_t poly_roots_longdouble(double *poly, uint32_t degree, long double complex *eroots) {

	uint32_t i;
	double complex rev_dccoeffs[MAX_ROOTFINDER_DEGREE + 1];
	long double complex ldccoeffs[MAX_ROOTFINDER_DEGREE + 1];

	if (degree == 1) {
		eroots[0] = -poly[0]/poly[1];
		return 0;
	}

	for (i = 0; i <= degree; i++) {
		rev_dccoeffs[degree - i] = poly[i];
		ldccoeffs[i] = poly[i];
	}

	/* find the roots to a relative error close to the
	   double-precision limit */

        double complex * droots = malloc(degree * sizeof(double complex));

	if (jenkins_traub(rev_dccoeffs, degree, droots) != degree)
		return 1;

	/* polish each root */

        int rc = 0;

	for (i = 0; i < degree; i++) {
		if (polish_root(ldccoeffs, degree,
                            droots[i], eroots + i, 1e-30) != 0)
                    rc++;
        }
        free(droots);

        /* change roots with very small imaginary part to
           be explicitly real roots */
	for (i = 0; i < degree; i++) {
		if (fabsl(cimagl(eroots[i])) <= 1e-30l*fabsl(creall(eroots[i])))
			eroots[i] = creall(eroots[i]);
	}
	return rc;
}

uint32_t poly_roots_double(double *poly, uint32_t degree, double complex *roots) {
        long double complex * eroots = malloc(degree * sizeof(long double complex));
        uint32_t res = poly_roots_longdouble(poly, degree, eroots);
        for (uint32_t i = 0; i < degree; i++) {
            roots[i]=eroots[i];
        }
        free(eroots);
        return res;
}



#ifdef  MAIN
int main(int argc, char * argv[])
{
    int degree = argc - 2;
    printf("degree %d\n", degree);
    double coeffs[MAX_ROOTFINDER_DEGREE];
    for(int i = 0 ; i <= degree ; i++) {
        coeffs[i] = atof(argv[1+i]);
    }
    double complex roots[MAX_ROOTFINDER_DEGREE];
    poly_roots_double(coeffs, degree, roots);
    for(int i = 0 ; i < degree ; i++) {
        printf("%f+i*%f\n", creal(roots[i]), cimag(roots[i]));
    }
}
#endif
