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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"

#include "test-tools.h"


int main(int argc, char * argv[])
{
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: ./check <size1> [<size2>]\n");
        fprintf(stderr, "Output of this program must be"
                " checked against precomputed tables\n");
        exit(1);
    }

    int n1 = atoi(argv[1]);
    int n2 = n1;
    if (argc >= 3) {
        n2 = atoi(argv[2]);
    }
    int n3 = n1 + n2;

    /* n1 and n2 must be understood as a number of 32-bit words */

    int N1 = n1;
    int N2 = n2;

#if GF2X_WORDSIZE == 64
    N1 = N1 / 2 + (N1 & 1);
    N2 = N2 / 2 + (N2 & 1);
#endif

    int N3 = N1 + N2;


    unsigned long * a = malloc(N1 * sizeof(unsigned long));
    unsigned long * b = malloc(N2 * sizeof(unsigned long));
    unsigned long * c = malloc(N3 * sizeof(unsigned long));

    uint32_t start=0x32104567UL;
    uint32_t ratio=0x76540123UL;

    start |= 1UL;

    uint32_t v;
    uint32_t check0, check1;

    v = start;
    v = fill(a, n1, v, ratio);
    v = fill(b, n2, v, ratio);
    gf2x_mul(c,a,N1,b,N2);
    check0 = gf2x_crc32(c, n3, 0);

    v = start;
    v = fill(c, n1, v, ratio);
    gf2x_mul(c,c,N1,b,N2);
    check1 = gf2x_crc32(c, n3, 0);

    if (check0 != check1) {
        printf("aliasing test failed\n");
        exit(1);
    }

    v = start;
    v = fill(a, n1, v, ratio);
    v = fill(c, n2, v, ratio);
    gf2x_mul(c,a,N1,c,N2);
    check1 = gf2x_crc32(c, n3, 0);

    if (check0 != check1) {
        printf("aliasing test failed\n");
        exit(1);
    }

    printf("%d %d %08" PRIx32 "\n", n1, n2, check0);

    free(a);
    free(b);
    free(c);

    return 0;
}

