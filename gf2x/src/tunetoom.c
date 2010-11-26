/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
/* Program to tune Toom-Cook multiplication over GF(2). */

/* How to use this program:

   1) beforehand, tune the low-level multiplication routines
      with ``make tune-lowlevel''
   2) build the tunetoom binary using ``make tunetoom''
   3) run tunetoom, giving as argument the maximum word size, for example
      ./tunetoom 2000 will tune multiplication of polynomials up to
      degree 128000 on a 64-bit machine. For higher degrees the FFT
      is probably faster - see tunefft.

      tunetoom works in two phases - first the "balanced" routines with
      equal-sized inputs are tuned, then the "unbalanced" routines where
      one input is about twice as large as the other are tuned.

      The results are printed on stdout by default, but is possible and
      perhaps preferrable to save them to another file using the -o
      option.

      Steps are equal-sized by 1, but with the -s option, it is possible
      to force multiplicative steps, e.g. -s 1.01 for 1% steps.

      Summary:
       ./tunetoom -s 1.05 2048 -o tunetoom.res
        -> does 5% steps up to size 2048, prepare output in tunetoom.res
       ./tunetoom 2048 -o tunetoom.res
        -> does 1 by 1 steps up to size 2048, prepare output in tunetoom.res

      The output file tunetoom.res can then be used as input for the
      update-thresholds program, like in the following:

       cat ../gf2x-thresholds.h > tuned_thresholds.h
       ./update-thresholds -o tuned_thresholds.h < tunetoom.res
       mkdir -p ../already_tuned/tuned/
       mv tuned_thresholds.h ../already_tuned/tuned/gf2x-thresholds.h
       rm -f ../gf2x-thresholds.h
       ln -sf already_tuned/tuned/gf2x-thresholds.h ../

   4) compile and run tunefft to tune FFT multiplication
      (see instructions in tunefft.c).
*/

#define _BSD_SOURCE
#define _POSIX_C_SOURCE 200112L /* solaris needs >= 199506L for ctime_r */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>		/* for LONG_MAX */
#include <string.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <sys/utsname.h>        /* for uname */
#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-small.h"
#include "timing.h"
#include "tuning-common.h"

#define MINI_GF2X_MUL_TOOM_THRESHOLD 	17
#define MINI_GF2X_MUL_TOOMW_THRESHOLD	8
#define MINI_GF2X_MUL_TOOM4_THRESHOLD	30
#define MINI_GF2X_MUL_TOOMU_THRESHOLD	33

#define BESTMIN (GF2X_MUL_KARA_THRESHOLD-1)
#define BESTMINU (GF2X_MUL_TOOMU_THRESHOLD-1)

const char * gf2x_toom_select_string[] = {
    [GF2X_SELECT_KARA] = "TC2",
    [GF2X_SELECT_TC3]  = "TC3",
    [GF2X_SELECT_TC3W] = "TC3W",
    [GF2X_SELECT_TC4]  = "TC4",
};

const char * gf2x_utoom_select_string[] = {
    [GF2X_SELECT_UNB_DFLT] = "default",
    [GF2X_SELECT_UNB_TC3U]  = "TC3U",
};

FILE *rp;

void tunetoom(long tablesz)
{
    long high, n;
    int k;
    double T3[1], TK[1], TW[1], T4[1];
    double mint;
    unsigned long *a, *b, *c, *d, *t;

    high = tablesz;
    if (high < BESTMIN)
	high = BESTMIN;

    if (high > GF2X_TOOM_TUNING_LIMIT) {
	fprintf(stderr,
	    "Increase constant GF2X_TOOM_TUNING_LIMIT in thresholds.h to %ld\n",
	     high);
	exit(1);
    }

    a = (unsigned long *) malloc(high * sizeof(unsigned long));
    b = (unsigned long *) malloc(high * sizeof(unsigned long));
    c = (unsigned long *) malloc(2 * high * sizeof(unsigned long));
    d = (unsigned long *) malloc(2 * high * sizeof(unsigned long));
    t = (unsigned long *) malloc(gf2x_toomspace(high) * sizeof(unsigned long));

    for (n = BESTMIN + 1; n <= high; ) {
	srandom(1);
	TK[0] = T3[0] = TW[0] = T4[0] = 0.0;
	printf("%ld ", n);
	fflush(stdout);
	random_wordstring(a, n);
	random_wordstring(b, n);
	if (n >= GF2X_MUL_KARA_THRESHOLD)
	    TIME(TK[0], gf2x_mul_kara(c, a, b, n, t));
	if (n >= MINI_GF2X_MUL_TOOM_THRESHOLD) {
	    TIME(T3[0], gf2x_mul_tc3(d, a, b, n, t));
	    check(a, n, b, n, "Kara", c, "TC3", d);
	}
	if (n >= MINI_GF2X_MUL_TOOMW_THRESHOLD) {
	    TIME(TW[0], gf2x_mul_tc3w(d, a, b, n, t));
	    check(a, n, b, n, "Kara", c, "TC3W", d);
	}
	if (n >= MINI_GF2X_MUL_TOOM4_THRESHOLD) {
	    TIME(T4[0], gf2x_mul_tc4(d, a, b, n, t));
	    check(a, n, b, n, "Kara", c, "TC4", d);
	}
	printf("TC2:%1.2e TC3:%1.2e TC3W:%1.2e TC4:%1.2e ",
	       TK[0], T3[0], TW[0], T4[0]);
	mint = TK[0];
	k = GF2X_SELECT_KARA;
	if ((T3[0] < mint) && (n >= MINI_GF2X_MUL_TOOM_THRESHOLD)) {
	    mint = T3[0];
	    k = GF2X_SELECT_TC3;
	}
	if ((TW[0] < mint) && (n >= MINI_GF2X_MUL_TOOMW_THRESHOLD)) {
	    mint = TW[0];
	    k = GF2X_SELECT_TC3W;
	}
	if ((T4[0] < mint) && (n >= MINI_GF2X_MUL_TOOM4_THRESHOLD)) {
	    mint = T4[0];
	    k = GF2X_SELECT_TC4;
	}
	printf("best:%1.2e %s\n", mint, gf2x_toom_select_string[k]);
        fprintf(rp, "toom %ld %d\n", n, k);
	fflush(stdout);
        long nn = MAX(n * mulstep, n + 1);
        for( ; n < nn && n <= high ; n++) {
            best_tab[n - 1] = k;
        }
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(t);

    return;
}

/* Forms c := a*b where b has size sb, a has size sa = (sb+1)/2 (words),
   representing polynomials over GF(2).  c needs space for sa+sb words.
   Needs space 2*sa + gf2x_toomspace(sa) words in stk[0] ...

   The code is essentially the same as in HalfGCD.c   */

static void gf2x_mul21(unsigned long *c, const unsigned long *b, long sb,
		  const unsigned long *a, unsigned long *stk)
{
    long i, j;
    long sa = (sb + 1) / 2;
    long sc = sa + sb;
    unsigned long *v;
    v = stk;
    stk += 2 * sa;
    for (i = 0; i < sc; i++)
	c[i] = 0;
    do {
	if (sa == 0)
	    break;

	if (sa == 1) {
	    c[sb] ^= gf2x_addmul_1_n(c, c, b, sb, a[0]);
	    break;
	}

	for (i = 0; i + sa <= sb; i += sa) {
	    gf2x_mul_toom(v, a, b + i, sa, stk);	// Generic Toom-Cook mult.
	    for (j = 0; j < 2 * sa; j++)
		c[i + j] ^= v[j];
	}

	{
	    const unsigned long *t;
	    t = a;
	    a = b + i;
	    b = t;
	}
	{
	    long t;
	    t = sa;
	    sa = sb - i;
	    sb = t;
	}
	c = c + i;
    }
    while (1);
}

void checku(const unsigned long *a, const unsigned long *b, long n)
{
    long i;
    for (i = 0; i < n; i++) {
	if (a[i] == b[i]) continue;
        fprintf(stderr,
                "Error detected: mul_toom3u and "
                "default give different results\n");
        printf("index %ld\n", i);
        exit(1);
    }
}

void tuneutoom(long tabsz)
{
    long high;
    int k;
    double T3[1], TK[1];
    double mint;
    unsigned long *a, *b, *c, *d, *t;

    high = tabsz;
    if (high < BESTMINU)
	high = BESTMINU;

    if (high > GF2X_TOOM_TUNING_LIMIT) {
	fprintf(stderr,
                "Increase constant GF2X_TOOM_TUNING_LIMIT in thresholds.c to %ld\n",
	     high);
	exit(1);
    }

    long sa = high;
    long sb = (sa + 1) / 2;

    long sp1 = gf2x_toomuspace(sa);	// space for mul_toom3u
    long sp2 = gf2x_toomspace(sb) + 2 * sb;	// space for mul21
    long sp = (sp1 > sp2) ? sp1 : sp2;

    a = (unsigned long *) malloc(sa * sizeof(unsigned long));
    b = (unsigned long *) malloc(sb * sizeof(unsigned long));
    c = (unsigned long *) malloc(3 * sb * sizeof(unsigned long));
    d = (unsigned long *) malloc(3 * sb * sizeof(unsigned long));
    t = (unsigned long *) malloc(sp * sizeof(unsigned long));


    for (sa = BESTMINU + 1; sa <= high; ) {
	sb = (sa + 1) / 2;
	random_wordstring(a, sa);
	random_wordstring(b, sb);
	TK[0] = T3[0] = 0.0;
	printf("%ld ", sa);
	fflush(stdout);
	TIME(TK[0], gf2x_mul21(c, a, sa, b, t));
        if (0) {
	// if (sa >= MINI_GF2X_MUL_TOOMU_THRESHOLD) {
	    TIME(T3[0], gf2x_mul_tc3u(d, a, sa, b, t));
	    checku(c, d, sa + sb);
	}
	printf("default:%1.2e TC3U:%1.2e ", TK[0], T3[0]);
	mint = TK[0];
	k = GF2X_SELECT_UNB_DFLT;
	if ((T3[0] < mint) && (sa >= MINI_GF2X_MUL_TOOMU_THRESHOLD)) {
	    mint = T3[0];
	    k = GF2X_SELECT_UNB_TC3U;
	}
	printf("best:%1.2e %s\n", mint, gf2x_utoom_select_string[k]);
	fflush(stdout);
        fprintf(rp, "utoom %ld %d\n", sa, k);
        long nn = MAX(sa * mulstep, sa + 1);
        for( ; sa < nn && sa <= high ; sa++) {
            best_utab[sa - 1] = k;
        }
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(t);

    return;
}

void usage(int rc)
{
    FILE * f = rc ? stderr : stdout;

    fprintf(f, "Usage: tunetoom [options] table-size1 [table-size2]\n");
    fprintf(f, " where %d <= table-size1 <= %d\n",
            BESTMIN, GF2X_TOOM_TUNING_LIMIT);
    fprintf(f, " and  %d <= table-size2 <= %d\n",
            BESTMINU, GF2X_TOOM_TUNING_LIMIT);
    fprintf(f, "Allowed options:\n");
    fprintf(f, "\t-s <multiplicative step>\n");
    fprintf(f, "\t-o <output file> (default is stdout or file desc. 3)\n");
    exit(rc);
}

int main(int argc, char *argv[])
{
    long tabsz1 = 0, tabsz2 = 0;

    int nsz = 0;

    char * progname = argc ? argv[0] : "";

    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        int r;
        if (strcmp(argv[0], "--help") == 0) {
            usage(0);
        }
        r = handle_tuning_mulstep(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;
        r = handle_tuning_outfile(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;

        if (nsz == 0) {
            tabsz1 = tabsz2 = atoi(argv[0]);
            nsz++;
            continue;
        }
        if (nsz == 1) {
            tabsz2 = atoi(argv[0]);
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
        for(;u && isspace(date[u-1]);date[--u]='\0');
        uname(&buf);

        /* strip the dirname */
        char * ptr = strrchr(progname, '/');
        if (ptr) {
            ptr++;
        } else {
            ptr = progname;
        }

        fprintf(rp, "info-toom \"%s -s %.2f %ld %ld run on %s on %s\"\n",
                ptr,mulstep,tabsz1,tabsz2,buf.nodename,date);
    }

    tunetoom(tabsz1);		// Tune balanced routines
    tuneutoom(tabsz2);		// Tune unbalanced routines

    fflush(stdout);

    return 0;
}

/* vim: set sw=4 sta et: */
