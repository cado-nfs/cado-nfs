/* This file is part of the gf2x library.

   Copyright 2016, 2023
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
#include "gf2x.h"
#include "gf2x-thresholds.h"

#include SOURCEFILE
#define PAD(X,Y) X ## Y
#define FUNCTION_NAME(X) PAD(gf2x_mul, X)


unsigned long bbb[TESTING_SIZE];
unsigned long aaa[TESTING_SIZE];
unsigned long ccc[2*TESTING_SIZE];
unsigned long rrr[2*TESTING_SIZE];

int main()
{
    int errors=0;
    for(int i = 0 ; i < TESTING_SIZE * GF2X_WORDSIZE ; i++) {
        memset(aaa, 0, TESTING_SIZE*sizeof(unsigned long));
        aaa[i / GF2X_WORDSIZE] = 1UL << (i % GF2X_WORDSIZE);
        for(int j = 0 ; j < TESTING_SIZE * GF2X_WORDSIZE ; j++) {
            memset(bbb, 0, TESTING_SIZE*sizeof(unsigned long));
            bbb[j / GF2X_WORDSIZE] = 1UL << (j % GF2X_WORDSIZE);
            memset(ccc, 0, 2*TESTING_SIZE*sizeof(unsigned long));
            memset(rrr, 0, 2*TESTING_SIZE*sizeof(unsigned long));
#if TESTING_SIZE > 1
            FUNCTION_NAME(TESTING_SIZE)(ccc, aaa, bbb);
#else
            FUNCTION_NAME(TESTING_SIZE)(ccc, aaa[0], bbb[0]);
#endif
            gf2x_mul(rrr, aaa, TESTING_SIZE, bbb, TESTING_SIZE);
            if (memcmp(ccc, rrr, 2*TESTING_SIZE*sizeof(unsigned long)) != 0) {
                fprintf(stderr, "Error t^%d * t^%d\n", i, j);
                for(int i = 0 ; i < 2*TESTING_SIZE ; i++) {
                    if (rrr[i] != ccc[i]) {
                    fprintf(stderr, GF2X_WORDSIZE == 64 ?
                            "[%02d] %016lx %016lx %016lx\n"
                            :
                            "[%02d] %08lx %08lx %08lx\n"
                            ,
                            i,
                            rrr[i],
                            ccc[i],
                            rrr[i] ^ ccc[i]);
                    }
                }
                errors++;
            }
        }
    }

    if (errors) {
        fprintf(stderr, "Found %d errors\n", errors);
        exit(EXIT_FAILURE);
    }
    return 0;
}


