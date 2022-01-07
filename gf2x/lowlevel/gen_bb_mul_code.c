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

/* Generates low-level multiplication routines over GF(2)[x]. */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char *argv[])
{
    unsigned long i;
    unsigned long w, K, CHOP, REM, fn;
    unsigned long MASK, mask2;

    if (argc != 3) {
	fprintf(stderr, "Usage: %s <wordsize> k\n", argv[0]);
	exit(1);
    }

    w = atoi(argv[1]);
    K = atoi(argv[2]);
    if (w <= 0 || (w & (w-1)) || w > 64 || K <= 0 || K >= w) {
        fprintf(stderr, "Usage: %s <wordsize> k\n", argv[0]);
        exit(1);
    }
    MASK = 1;
    if ((MASK << K) < MASK) {
        fprintf(stderr, "Usage: %s <wordsize> k\n", argv[0]);
        exit(1);
    }

    printf(
"/* This file is part of the gf2x library.\n"
"\n"
"   Copyright 2007, 2008, 2009, 2010, 2011, 2012, 2013\n"
"   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann\n"
"\n"
"   This program is free software; you can redistribute it and/or modify it\n"
"   under the terms of either:\n"
"    - If the archive contains a file named toom-gpl.c (not a trivial\n"
"    placeholder), the GNU General Public License as published by the Free\n"
"    Software Foundation; either version 3 of the License, or (at your\n"
"    option) any later version.\n"
"    - If the archive contains a file named toom-gpl.c which is a trivial\n"
"    placeholder, the GNU Lesser General Public License as published by\n"
"    the Free Software Foundation; either version 2.1 of the License, or\n"
"    (at your option) any later version.\n"
"\n"
"   This program is distributed in the hope that it will be useful, but WITHOUT\n"
"   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or\n"
"   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.\n"
"\n"
"   You should have received a copy of the GNU General Public License as\n"
"   well as the GNU Lesser General Public License along with this program;\n"
"   see the files COPYING and COPYING.LIB.  If not, write to the Free\n"
"   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA\n"
"   02110-1301, USA.\n"
"*/\n"
"\n");
    printf("#ifndef GF2X_MUL1_H_\n");
    printf("#define GF2X_MUL1_H_\n");
    printf("\n");
    printf("/* This file was generated automatically with\n");
    printf("   %s %lu %lu. Don't edit it! */\n",
	   argv[0], w, K);
    printf("\n");
    printf(
        "#include \"gf2x.h\"\n"
"/* All gf2x source files for lowlevel functions must include gf2x-small.h\n"
" * This is mandatory for the tuning mechanism. */\n"
        "#include \"gf2x/gf2x-small.h\"\n"
        "\n");

    for (fn = 0; fn < 3; fn++) {
	if (fn == 0) {
	    printf("GF2X_STORAGE_CLASS_mul1 void\n");
	    printf("gf2x_mul1 (unsigned long *c, unsigned long a, unsigned long b)\n");
	} else if (fn == 1) {
		// continue;	/* mul_1_n is no longer used */
	    printf("GF2X_STORAGE_CLASS_mul_1_n unsigned long\n");
	    printf
		("gf2x_mul_1_n (unsigned long *cp, const unsigned long *bp, long sb, unsigned long a)\n");
	} else if (fn == 2) {
	    printf("GF2X_STORAGE_CLASS_addmul_1_n unsigned long\n");
	    printf
		("gf2x_addmul_1_n (unsigned long *dp, const unsigned long *cp, const unsigned long* bp, long sb,\n");
	    printf("        unsigned long a)\n");
	}
	printf("{\n");
	if (fn > 0) {		/* Mul1/AddMul1 */
	    printf("   long i;\n");
	    printf("   unsigned long carry = 0, b;\n");
	}
	printf("   unsigned long hi, lo, tmp, A[%lu];\n", 1UL << K);
	printf("\n");
	printf("   A[0] = 0;\t\t");
	/* do not truncate: fix non-considered bit products afterwards */
	printf("A[1] = a;\n");
	for (i = 2u; i < (1u << K); i++) {
	    if (i % 2 == 0)
		printf("   A[%lu] = A[%lu] << 1;\t", i, i / 2);
	    else
		printf("A[%lu] = A[%lu] ^ a;\n", i, i - 1);
	}
	printf("\n");

	if (fn > 0) {
	    printf("   for (i = 0; i < sb; i++)\n");
	    printf("     {\n");
	    printf("       b = bp[i];\n");
	}
	REM = w - 2 * K;
	MASK = (1UL << K) - 1UL;
	for (CHOP = 0; CHOP + 2 * K < w; CHOP += 2 * K);
	CHOP = w - CHOP;
	/* now 1 <= CHOP <= 2 * K, with w - CHOP multiple of 2K */
	if (CHOP <= K)		/* one slice is enough for the upper part */
	    printf("   lo = A[b >> %lu];\n", w - CHOP);
	else			/* two slices for the upper part */
	    printf("   lo = (A[b >> %lu] << %lu) ^ A[(b >> %lu) & %lu];\n",
		   w - (CHOP - K), K,
		   w - CHOP, MASK);
	printf("   hi = lo >> %lu;\n", w - 2 * K);
	for (i = w - CHOP - 2 * K; i >= 2 * K; i -= 2 * K) {
	    printf
		("   lo = (lo << %lu) ^ (A[(b >> %lu) & %lu] << %lu) ^ A[(b >> %lu) & %lu];\n",
		 2 * K, i + K, MASK, K, i, MASK);
	    printf("   hi = (hi << %lu) | (lo >> %lu);\n", 2 * K, REM);
	}
	/* special case for i=0 since a shift of 0 is undefined */
	printf
	    ("   lo = (lo << %lu) ^ (A[(b >> %lu) & %lu] << %lu) ^ A[b & %lu];\n",
	     2 * K, K, MASK, K, MASK);

	/* now perform the repair step for 2*K */
	MASK = 0;
	for (i = 0; i < w; i += 2 * K)
	    MASK |= 1UL << i;
	MASK = ~MASK;
	mask2 = MASK;
	printf("\n");
	for (i = 1; i < 2 * K; i++) {
	    /* bit w-i from a was not considered in blocks of
	       K bits from b for index j >= i */
	    /* Gaudry's no-branching optimization */
	    printf("   tmp = -((a >> %lu) & 1);\n", w - i);
            if (w == 32) {
                mask2 = (unsigned long) (uint32_t) mask2;
            }
	    printf("   tmp &= ((b & 0x%lx) >> %lu);\n", mask2, i);
	    printf("   hi = hi ^ tmp;\n");
	    mask2 = (mask2 << 1) & MASK;
	}
	printf("\n");
	if (fn == 0) {
	    printf("   c[0] = lo;\n");
	    printf("   c[1] = hi;\n");
	} else if (fn == 1) {
	    printf("   cp[i] = carry ^ lo;\n");
	    printf("   carry = hi;\n");
	} else if (fn == 2) {
	    printf("   dp[i] = cp[i] ^ (carry ^ lo);\n");
	    printf("   carry = hi;\n");
	}
	if (fn > 0) {
	    printf("    }\n");
	    printf("   return carry;\n");
	}
	printf("}\n\n");
    }

    printf("#endif\t/* GF2X_MUL1_H_ */\n");

    return 0;
}
