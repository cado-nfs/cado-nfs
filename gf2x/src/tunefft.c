/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

/* Program to tune the FFT multiplication over GF(2). */

/* This program is called by "make tune-fft". You may also want to call
 * it independently.
 *
 * You are advised to run "make tune-lowlevel" and "make tune-toom"
 * beforehand, or otherwise the tuning results that you obtain with this
 * program are likely to have limited significance.
 *
 * _IF_ you want to run this program independently of the "make
 * tune-fft" scenario, the instructions are as follows.

   1) beforehand, tune the low-level and mid-sized multiplication routines
      with ``make tune-lowlevel'' and ``make tune-toom''

   2) build the tunefft binary using ``make tunefft''
   3) run tunefft, giving as argument the maximum word size, for example
      ./tunefft 8000000 will tune multiplication of polynomials up to
      degree 512000000 on a 64-bit machine.

      The results are printed on stdout by default, but is possible and
      perhaps preferrable to save them to another file using the -o
      option.

      Giving a multiplicative step value with the -s argument makes the
      tuning quicker, but less accurate.

      The output file tunefft.res can then be used as input for the
      modify-thresholds program, like in the following:

       cat ../gf2x-thresholds.h > tuned_thresholds.h
       ./modify-thresholds -o tuned_thresholds.h < tunefft.res
       mkdir -p ../already_tuned/tuned/
       mv tuned_thresholds.h ../already_tuned/tuned/gf2x-thresholds.h
       rm -f ../gf2x-thresholds.h
       ln -sf already_tuned/tuned/gf2x-thresholds.h ../
*/

#define _DEFAULT_SOURCE /* _BSD_SOURCE is deprecated */
#define _POSIX_C_SOURCE 200112L /* solaris needs >= 199506L for ctime_r */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>		/* for LONG_MAX */
#include <float.h>		/* for DBL_MAX */
#include <time.h>
#include <sys/utsname.h>	/* for uname */
#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
#include "fft/gf2x-fft-impl-utils.h"
#include "timing.h"
#include "tuning-common.h"

/* This version of tunefft uses the midpoint of each stair */

/* Must be at least >= GF2X_TERNARY_FFT_MINIMUM_SIZE (28), but it saves time to set larger */
#define GF2X_MUL_FFT_BEGIN_TUNE 1000

#define STEPMAX 50

/* Return the largest m >= n such that a product of m*m words with an FFT of
   length K leads to the same size of pointwise products as with n words.
   More precisely if f(n) = K/3 * ceil(6*n*GF2X_WORDSIZE/K^2) then m is the
   largest integer such that
   ceil(f(n)/GF2X_WORDSIZE) = ceil(f(m)/GF2X_WORDSIZE).
   In fact we take for m the largest integer such that g(n) = g(m) with
   g(n) = ceil(6*n*GF2X_WORDSIZE/K^2).
*/

static inline void wantzero(int x)
{
    if (x != 0) {
        abort();
    }
}

long end_of_stair(long n, long K)
{
    long g = 6 * n * GF2X_WORDSIZE;

    g = 1 + (g - 1) / K;	/* ceil(g/K) */
    g = 1 + (g - 1) / K;	/* ceil(g/K^2) */
    g = g * K * K;
    return g / (6 * GF2X_WORDSIZE);
}

long next_step(long n, long K)
{
    // Return the minimum so steps are not too large
    return MIN(end_of_stair(n, K), MAX(n+1, n*mulstep));

}

void usage(int rc)
{
    FILE * f = rc ? stderr : stdout;

    fprintf(f, "Usage: tunefft [options] [min_word_size] max_word_size\n");
    fprintf(f, "Allowed options:\n");
    fprintf(f, "\t-s <multiplicative step>\n");
    fprintf(f, "\t-o <output file> (default is stdout or file desc. 3)\n");
    exit(rc);
}

int main(int argc, char *argv[])
{
    long minn, maxn, mid, n, n2, ns, i;
    long besti;			/* 0 for TC, 1, 2, ... for FFT(K0*3^(bestK-1)) */
    long bestK;
    long K, K0 = 3;		/* try K0, 3*K0, 9*K0 */
    double T[4];		/* T[0] is for TC, T[1] for K0, T[2] for 3*K0, T[3] for 9*K0 */
    double t1[4], t2[4];
    unsigned long *a, *b, *c, *t, *u, *v;
    int nsz = 0;
    int tc_takes_too_long = 0;
    const char * reference = "TC";

    maxn = 1000000;		// default
    minn = GF2X_MUL_FFT_BEGIN_TUNE / 2 + 1;

    char * progname = argc ? argv[0] : "";

    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        int r;

        if (strcmp(argv[0], "--help") == 0) {
            usage(0);
        }

        if (strcmp(argv[0], "--no-toom") == 0) {
            tc_takes_too_long = 1;
            reference = "F1(K0)";
            continue;
        }

        r = handle_tuning_mulstep(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;
        r = handle_tuning_outfile(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;

        if (strcmp(argv[0], "-k0") == 0) {
            argc--,argv++;
            if (! argc) usage(1);
            K0 = atoi(argv[0]);
            continue;
        }

        if (nsz == 0) {
            maxn = atoi(argv[0]);
            nsz++;
            continue;
        }
        if (nsz == 1) {
            minn = maxn;
            maxn = atoi(argv[0]);
            nsz++;
            continue;
        }
        usage(1);
    }

    if (nsz == 0)
        usage(1);

    set_tuning_output();

    {
	char date[40];
	time_t t;
	size_t u;
	struct utsname buf;
	time(&t);
	ctime_r(&t, date);
	u = strlen(date);
	for (; u && isspace(date[u - 1]); date[--u] = '\0');
	uname(&buf);

        /* strip the dirname */
        char * ptr = strrchr(progname, '/');
        if (ptr) {
            ptr++;
        } else {
            ptr = progname;
        }

        fprintf(rp, "info-fft \"%s -s %.2f %ld run on %s on %s ; based on %s\"\n",
                ptr,mulstep,maxn,buf.nodename,date,GF2X_TOOM_TUNING_INFO);
    }

    printf("Tuning FFT multiplication to wordsize %ld\n\n", maxn);

    a = (unsigned long *) malloc(maxn * sizeof(unsigned long));
    b = (unsigned long *) malloc(maxn * sizeof(unsigned long));
    c = (unsigned long *) malloc(2 * maxn * sizeof(unsigned long));
    u = (unsigned long *) malloc(2 * maxn * sizeof(unsigned long));
    v = (unsigned long *) malloc(2 * maxn * sizeof(unsigned long));
    t = (unsigned long *) malloc(gf2x_toomspace(maxn) * sizeof(unsigned long));

    random_wordstring(a, maxn);
    random_wordstring(b, maxn);

/* Skip n if (2*n < GF2X_MUL_FFT_BEGIN_TUNE) as this is too small for the FFT */


    for (n = minn; n <= maxn;) {
	n2 = next_step(n, 3 * K0);	// End of interval
	if (n2 > maxn)		// Only go as far
	    n2 = maxn;		// as maxn.
	mid = (n + n2) / 2;	// Mid-point
	printf("%ld..%ld ", n, n2);
	fflush(stdout);

        if (tc_takes_too_long) {
            T[0] = DBL_MAX;
        } else {
            TIME(T[0], gf2x_mul_toom(u, a, b, mid, t));	// Time Toom-Cook
            printf("TC:%1.1e ", T[0]);
        }
	fflush(stdout);
	besti = 0;
	bestK = 1;
        K = K0;
        i = 1;
ugly_label:
	for ( ; i <= 3; i++, K *= 3) {
	    TIME(t1[i], wantzero(gf2x_mul_fft(c, a, mid, b, mid, K) == 0));
            if (tc_takes_too_long) {
                memcpy(u, c, 2 * maxn * sizeof(unsigned long));
            }
            check(a, mid, b, mid, reference, u, "F1", c);
            if (K >= GF2X_WORDSIZE) {
                TIME(t2[i], wantzero(gf2x_mul_fft(v, a, mid, b, mid, -K) == 0));
                check(a, mid, b, mid, "F1", c, "F2", v);
            } else {
                t2[i] = DBL_MAX;
            }
	    if (t1[i] < t2[i]) {
		T[i] = t1[i];
		printf("F1(%ld):%1.1e ", K, T[i]);
	    } else {
		T[i] = t2[i];
		printf("F2(%ld):%1.1e ", K, T[i]);
	    }
	    fflush(stdout);
	    if (T[i] < T[besti]) {
		besti = i;
		bestK = (t2[i] > t1[i]) ? K : -K;	/* -K for FFT2(|K|) */
	    }
	}

	if (T[3] < T[1] && T[3] < T[2]) {
            if (besti) {
                if (besti == 1)
                    abort();
                besti--;
            }
	    K0 *= 3;
            /* K just stays as it was */
            i = 3;
            T[1] = T[2];
            T[2] = T[3];
            goto ugly_label;
            /* Notice that we can't loop forever here. If we have T[3] <
             * T[2], this will ensure T[2] < T[1] at the next turn,
             * thereby forcing the other case not to happen */
        } else if (T[1] < T[2] && T[1] < T[3] && K0 > 3) {
	    K0 /= 3;
        }

        /* OK, this stair is done */

	if (bestK == 1)
	    printf("TC");
	else {
	    if (bestK > 0)
		printf("F1(%ld)", bestK);
	    else
		printf("F2(%ld)", -bestK);
	}
	printf("\n");
	fflush(stdout);

        if (T[0] >= 4 * T[besti] && !tc_takes_too_long) {
            printf("TC is taking too long, disabling for next sizes\n");
            tc_takes_too_long = 1;
            reference = "F1(K0)";
        }

	/* go to next size */
	ns = n;
	n = next_step(n, 3 * K0);	/* middle value of K */
	if (n > n2)
	    n = n2;		/* end of last stair if K0 increased */
	n++;
	if (n < mid) {		/* redo the last stair if K0 decreased */
	    n = ns;
        } else {
            fprintf(rp, "fft %ld %ld\n",
                    ns == minn ? 1 : ns, ns == minn ? 1 : bestK);
        }
    }

    free(a);
    free(b);
    free(c);
    free(t);
    free(u);
    free(v);

    return 0;
}
