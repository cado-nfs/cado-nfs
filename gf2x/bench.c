/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015, 2016
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

#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <stdlib.h>
#include "gf2x.h"

/* Usage: ./bench <limit size> */

/* return runtime in seconds */
double runtime(void)
{
  return (double) clock () / CLOCKS_PER_SEC;
}

int main(int argc, char * argv[])
{
    unsigned long i;
    unsigned long n0 = 1;
    unsigned long n = 200;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    if (argc > 2) {
        n0 = n;
        n = atoi(argv[2]);
    }

    setbuf(stdout,NULL);

    unsigned long * a, * b, * c;

    a = malloc(n * sizeof(unsigned long));
    b = malloc(n * sizeof(unsigned long));
    c = malloc(2 * n * sizeof(unsigned long));

    for(i = 0 ; i < n ; i++) {
        a[i] = rand();
        b[i] = rand();
    }

    unsigned long jmax = 100000;

    for(i = n0 ; i < n ; i<<=1) {
        /* warm up */
        gf2x_mul(c,a,i,b,i);
        gf2x_mul(c,a,i,b,i);
        gf2x_mul(c,a,i,b,i);
        double t = runtime();
        unsigned long j;
        double d;
#if 0
        for(j = 0 ; j < jmax ; j++) {
            gf2x_mul(c,a,i,b,i);
        }
        d = runtime() - t;

        if (d > 1 && jmax > 8) {
            jmax -= jmax / 3;
        }
#else
        for(j = 0 ; (d=runtime()-t) < 1 && (j < 2000 || d < 0.2) ; j++) {
            gf2x_mul(c,a,i,b,i);
        }
#endif
        printf("%u %e\t%u %.2f\n", i, d/j,j,d);

    }

    free(a);
    free(b);
    free(c);
}

