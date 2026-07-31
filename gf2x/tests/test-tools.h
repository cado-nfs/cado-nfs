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

#ifndef TEST_TOOLS_H_
#define TEST_TOOLS_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* fill the area a[0..X[ with n1 32-bit integers, all equal to
 * startvalue*ratio^i (for to i<n1). This represent n1 unsigned longs on
 * 32-bits, and twice less on 64 bits). We assume little-endian.
 */
uint32_t fill(unsigned long * a, int n1, uint32_t startvalue, uint32_t ratio);

/* return the crc32 mod x^32 + x^7 + x^6 + x^2 + 1 of the memory area */
uint32_t gf2x_crc32(unsigned long * c, int n3, uint32_t v);

#ifdef __cplusplus
}
#endif

#endif	/* TEST_TOOLS_H_ */
