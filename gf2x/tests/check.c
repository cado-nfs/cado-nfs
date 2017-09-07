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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdint.h>
#include <inttypes.h>

#include "gf2x/gf2x-config.h"
#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-thresholds.h"

#include "test-tools.h"
#include "check-and-bench-common.h"

int print_crc = 0;

/* n32bitwords is a number of 32-bit words ; k is a matrix dimension */
int test(int n32bitwords, int d)
{
    int nwords = n32bitwords;
#if GF2X_WORDSIZE == 64
    nwords = nwords / 2 + (nwords & 1);
#endif

    /* now nwords is a number of unsigned longs */

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

    if (!print_crc) {
        printf("w := %u;\n", GF2X_WORDSIZE);
#ifdef  ENGINE_CANTOR
        printf("b := %u;\n", CANTOR_BASE_FIELD_SIZE);
#endif
        printf("f := [\n");
        for(int i = 0 ; i < d ; i++)
            for(int j = 0 ; j < d ; j++) {
                if (i || j) printf(", ");
                printf("[");
                for (int k = 0; k < nwords - 1; ++k)
                    printf("%lu, ", f[i*d+j][k]);
                printf("%lu\n]\n", f[i*d+j][nwords - 1]);
            }
        printf("];\n");
        printf("g := [\n");
        for(int i = 0 ; i < d ; i++)
            for(int j = 0 ; j < d ; j++) {
                if (i || j) printf(", ");
                printf("[");
                for (int k = 0; k < nwords - 1; ++k)
                    printf("%lu, ", g[i*d+j][k]);
                printf("%lu\n]\n", g[i*d+j][nwords - 1]);
            }
        printf("];\n");
    }

    //for (int i = 0; i < 10; ++i) 
    ENGINE_mul(h,f,nwords,g,nwords,d);

    if (!print_crc) {
        printf("fg := [\n");
        for(int i = 0 ; i < d ; i++)
            for(int j = 0 ; j < d ; j++) {
                if (i || j) printf(", ");
                printf("[");
                for (int k = 0; k < 2*nwords - 1; ++k)
                    printf("%lu, ", h[i*d+j][k]);
                printf("%lu\n]\n", h[i*d+j][2*nwords - 1]);
            }
        printf("];\n");
    } else {
        uint32_t check0 = 0;
        for(int i = 0 ; i < d ; i++)
            for(int j = 0 ; j < d ; j++)
                check0 = crc32(h[(d-1-i)*d+(d-1-j)], 2*n32bitwords, check0);
        printf("%d %d %08" PRIx32 "\n", n32bitwords, d, check0);
    }

    for(int i = 0 ; i < d ; i++)
        for(int j = 0 ; j < d ; j++) {
            free(f[i*d+j]);
            free(g[i*d+j]);
            free(h[i*d+j]);
        }
    free(h);
    free(g);
    free(f);
    return 0;
}
void usage_and_die(char **argv)
{
    fprintf(stderr, "usage: %s --n <n32bitwords> --k <k>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv)
{
    int n32bitwords = -1;
    int k = 1;

    for(int i = 1 ; i < argc ; i++) {
        if (strcmp(argv[i], "-h") == 0) usage_and_die(argv);
        if (strcmp(argv[i], "-") == 0) usage_and_die(argv);
        if (strcmp(argv[i], "--help") == 0) usage_and_die(argv);
        if (strcmp(argv[i], "--crc") == 0) {
            print_crc=1;
            continue;
        }

        if (i < argc - 1 && strcmp(argv[i], "--n") == 0) {
            n32bitwords = atoi(argv[++i]);
            continue;
        }
        if (i < argc - 1 && strcmp(argv[i], "--k") == 0) {
            k = atoi(argv[++i]);
            continue;
        }
        usage_and_die(argv);
    }

    if (n32bitwords < 0)
        usage_and_die(argv);

    set_extra_arg_from_n32bitwords(n32bitwords, 0);

    test(n32bitwords, k);

    return 0;
}
