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

#include "gf2x/gf2x-config.h"
#include "gf2x/gf2x-thresholds.h"

#include "test-tools.h"

/* Fill the ulong* array a by n1 32-bit values in geometric progression.
 * Returns a value which is appropriate for chaining.
 */
uint32_t fill(unsigned long * a, int n1, uint32_t startvalue, uint32_t ratio)
{
    int i;
    uint32_t v = startvalue;
#if GF2X_WORDSIZE == 32
    for(i = 0 ; i < n1 ; i++) { a[i] = v; v *= ratio; }
#elif GF2X_WORDSIZE == 64
    int N1 = n1 / 2 + (n1 & 1);
    unsigned long w = 1;
    for(i = 0 ; i < N1 ; i++) {
        w = (unsigned long) v * (unsigned long) ratio;
        a[i] = w << 32 | (unsigned long) v;
        v = w * ratio;
    }
    if (n1 & 1) { v = w; a[N1-1] &= (1UL << 32) - 1; }
#else
#error "Config problem"
#endif
    return v;
}

#define RED(l,h) do {							\
        /* Compute crc mod x^32 + x^7 + x^6 + x^2 + 1 */		\
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;				\
        h  = h >> 30 ^ h >> 26 ^ h >> 25;				\
        /* h is at most 7 bits now. */					\
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;				\
        h = 0;								\
} while (0)

uint32_t crc32(unsigned long * c, int n3, uint32_t v)
{
    int i;
#if GF2X_WORDSIZE == 32
    for(i = 0 ; i < n3 ; i++) {
        uint32_t l,h;
        l = c[n3-1-i];
        h = v;
        RED(l,h);
        v = l;
    }
#elif GF2X_WORDSIZE == 64
    int N3 = n3 / 2 + (n3 & 1);
    i = 0;
    for(int j = 0 ; j < N3 ; j++) {
        uint32_t l, h;
        unsigned long cj = c[N3-1-j];

        if (2*(N3-1-j) + 1 < n3) {
            h = v;
            l = cj >> 32;
            RED(l,h);
            i++;
            v = l;
        }
        if (i == n3) break;

        if (2*(N3-1-j) < n3) {
            h = v;
            l = cj;
            RED(l,h);
            i++;
            v = l;
        }
        if (i == n3) break;
    }
#endif
    return v;
}

