/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2014, 2015
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
#define _ISOC99_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <math.h>		// for random

#include "gf2x.h"
#include "gf2x/gf2x-config.h"
#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-thresholds.h"

#include "test-tools.h"
#include "check-and-bench-common.h"

void usage_and_die(char **argv)
{
    fprintf(stderr, "usage: %s nmin nmax [-s step|-f factor] [-o order]\n", argv[0]);
    fprintf(stderr, "  where nmin and nmax are numbers of 64-words\n");
    exit(1);
}

/* we are talking numbers of words here */
void bench(int nwords, int dim, int quick, long supplied)
{
    size_t n32bitwords = nwords;
#if GF2X_WORDSIZE == 64
    n32bitwords *= 2;
#endif

    int d = quick ? 1 : dim;

    unsigned long ** f = (unsigned long**) malloc(d*d*sizeof(unsigned long*));
    unsigned long ** g = (unsigned long**) malloc(d*d*sizeof(unsigned long*));
    unsigned long ** h = (unsigned long**) malloc(d*d*sizeof(unsigned long*));

    for(int i = 0 ; i < d ; i++)
        for(int j = 0 ; j < d ; j++) {
            f[i*d+j] = (unsigned long *) malloc(nwords * sizeof(unsigned long));
            g[i*d+j] = (unsigned long *) malloc(nwords * sizeof(unsigned long));
            h[i*d+j] = (unsigned long *) malloc(2 * nwords * sizeof(unsigned long));
        }

    uint32_t start=0x32104567UL;
    uint32_t ratio=0x76540123UL;

    uint32_t v = start;
    for(int i = 0 ; i < d ; i++)
        for(int j = 0 ; j < d ; j++) {
            v = fill(f[i*d+j], n32bitwords, v, ratio);
        }
    for(int i = 0 ; i < d ; i++)
        for(int j = 0 ; j < d ; j++) {
            v = fill(g[i*d+j], n32bitwords, v, ratio);
        }


    int rep;

    time_total = 0;
    time_dft = 0;
    time_ift = 0;
    time_conv = 0;

    set_extra_arg_from_n32bitwords(n32bitwords, supplied);

    long order = 0;
    for(rep = 0 ; time_total < 0.5 ; rep++) {
        order = ENGINE_mul(h, f, nwords, g, nwords, d);
    }

    double mmul_ft = ((double) dim/d) * ((double) dim/d);
    double mmul_cv = ((double) dim/d) * ((double) dim/d) * ((double) dim/d);
    printf("%d %f %f %f %f %ld\n",
            nwords,
            quick ? NAN : ((double) time_total / rep),
            mmul_ft * time_dft / 2.0 / rep,
            mmul_ft * time_ift / rep,
            mmul_cv * time_conv/ rep,
            order
          );
    fflush(stdout);

    for(int i = 0 ; i < d ; i++)
        for(int j = 0 ; j < d ; j++) {
            free(f[i*d+j]);
            free(g[i*d+j]);
            free(h[i*d+j]);
        }
    free(h);
    free(g);
    free(f);
}

/* usage:
 *   ./bench nmin nmax [step]
 */

int main(int argc, char **argv)
{
    int nmin = 0, nmax = 0;
    int step = 1;
    double factor = 1.0;
    int additive = 1;
    int i;
    long supplied = 0;
    int d = 1;
    int quick = 0;

    for(i = 1 ; i < argc ; i++) {
        if (strcmp(argv[i], "-h") == 0) usage_and_die(argv);
        if (strcmp(argv[i], "-") == 0) usage_and_die(argv);
        if (strcmp(argv[i], "--help") == 0) usage_and_die(argv);

        if (i < argc - 1 && strcmp(argv[i], "-f") == 0) {
	    additive = 0;
	    factor = atof(argv[++i]);
            continue;
        }
        if (i < argc - 1 && strcmp(argv[i], "-s") == 0) {
	    additive = 1;
	    step = atoi(argv[++i]);
            continue;
        }
        if (i < argc - 1 && strcmp(argv[i], "-d") == 0) {
	    d = atoi(argv[++i]);
            continue;
        }
        if (strcmp(argv[i], "-q") == 0) {
            quick=1;
            continue;
        }
        if (i < argc - 1 && strcmp(argv[i], "-o") == 0) {
            supplied = atol(argv[++i]);
            continue;
        }
        if (!nmin) {
            nmin = atoi(argv[i]);
            if (!nmin) usage_and_die(argv);
            continue;
        }
        if (!nmax) {
            nmax = atoi(argv[i]);
            if (!nmax) usage_and_die(argv);
            continue;
        }
        usage_and_die(argv);
    }


    if (nmin < 33) {
	fprintf(stderr,
		"Warning: nmin adjusted to the minimal value allowed: 33.\n");
	nmin = 33;
    }
    if (nmin > nmax) {
	fprintf(stderr, "Error: nmin should be less or equal to nmax\n");
	usage_and_die(argv);
    }

    for(int n = nmin ; n < nmax ; ) {
        bench(n, d, quick, supplied);
        if (additive)
            n += step;
        else {
            int new_n = (int) (((double) n) * factor);
            if (new_n <= n)
                n++;
            else
                n = new_n;
        }
    }

    return 0;
}
