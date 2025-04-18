/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2015
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

#ifndef GF2X_MUL1_H_
#define GF2X_MUL1_H_

/* This file was generated automatically with
   ./gen_bb_mul_code 32 2. Don't edit it! */

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

GF2X_STORAGE_CLASS_mul1 void
gf2x_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
   unsigned long hi, lo, tmp, A[4];

   A[0] = 0;		A[1] = a;
   A[2] = A[1] << 1;	A[3] = A[2] ^ a;

   lo = (A[b >> 30] << 2) ^ A[(b >> 28) & 3];
   hi = lo >> 28;
   lo = (lo << 4) ^ (A[(b >> 26) & 3] << 2) ^ A[(b >> 24) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 22) & 3] << 2) ^ A[(b >> 20) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 18) & 3] << 2) ^ A[(b >> 16) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 14) & 3] << 2) ^ A[(b >> 12) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 10) & 3] << 2) ^ A[(b >> 8) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 6) & 3] << 2) ^ A[(b >> 4) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 2) & 3] << 2) ^ A[b & 3];

   tmp = -((a >> 31) & 1);
   tmp &= ((b & 0xeeeeeeee) >> 1);
   hi = hi ^ tmp;
   tmp = -((a >> 30) & 1);
   tmp &= ((b & 0xcccccccc) >> 2);
   hi = hi ^ tmp;
   tmp = -((a >> 29) & 1);
   tmp &= ((b & 0x88888888) >> 3);
   hi = hi ^ tmp;

   c[0] = lo;
   c[1] = hi;
}

GF2X_STORAGE_CLASS_mul_1_n unsigned long
gf2x_mul_1_n (unsigned long *cp, const unsigned long *bp, long sb, unsigned long a)
{
   long i;
   unsigned long carry = 0, b;
   unsigned long hi, lo, tmp, A[4];

   A[0] = 0;		A[1] = a;
   A[2] = A[1] << 1;	A[3] = A[2] ^ a;

   for (i = 0; i < sb; i++)
     {
       b = bp[i];
   lo = (A[b >> 30] << 2) ^ A[(b >> 28) & 3];
   hi = lo >> 28;
   lo = (lo << 4) ^ (A[(b >> 26) & 3] << 2) ^ A[(b >> 24) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 22) & 3] << 2) ^ A[(b >> 20) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 18) & 3] << 2) ^ A[(b >> 16) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 14) & 3] << 2) ^ A[(b >> 12) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 10) & 3] << 2) ^ A[(b >> 8) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 6) & 3] << 2) ^ A[(b >> 4) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 2) & 3] << 2) ^ A[b & 3];

   tmp = -((a >> 31) & 1);
   tmp &= ((b & 0xeeeeeeee) >> 1);
   hi = hi ^ tmp;
   tmp = -((a >> 30) & 1);
   tmp &= ((b & 0xcccccccc) >> 2);
   hi = hi ^ tmp;
   tmp = -((a >> 29) & 1);
   tmp &= ((b & 0x88888888) >> 3);
   hi = hi ^ tmp;

   cp[i] = carry ^ lo;
   carry = hi;
    }
   return carry;
}

GF2X_STORAGE_CLASS_addmul_1_n unsigned long
gf2x_addmul_1_n (unsigned long *dp, const unsigned long *cp, const unsigned long* bp, long sb,
        unsigned long a)
{
   long i;
   unsigned long carry = 0, b;
   unsigned long hi, lo, tmp, A[4];

   A[0] = 0;		A[1] = a;
   A[2] = A[1] << 1;	A[3] = A[2] ^ a;

   for (i = 0; i < sb; i++)
     {
       b = bp[i];
   lo = (A[b >> 30] << 2) ^ A[(b >> 28) & 3];
   hi = lo >> 28;
   lo = (lo << 4) ^ (A[(b >> 26) & 3] << 2) ^ A[(b >> 24) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 22) & 3] << 2) ^ A[(b >> 20) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 18) & 3] << 2) ^ A[(b >> 16) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 14) & 3] << 2) ^ A[(b >> 12) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 10) & 3] << 2) ^ A[(b >> 8) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 6) & 3] << 2) ^ A[(b >> 4) & 3];
   hi = (hi << 4) | (lo >> 28);
   lo = (lo << 4) ^ (A[(b >> 2) & 3] << 2) ^ A[b & 3];

   tmp = -((a >> 31) & 1);
   tmp &= ((b & 0xeeeeeeee) >> 1);
   hi = hi ^ tmp;
   tmp = -((a >> 30) & 1);
   tmp &= ((b & 0xcccccccc) >> 2);
   hi = hi ^ tmp;
   tmp = -((a >> 29) & 1);
   tmp &= ((b & 0x88888888) >> 3);
   hi = hi ^ tmp;

   dp[i] = cp[i] ^ (carry ^ lo);
   carry = hi;
    }
   return carry;
}

#endif	/* GF2X_MUL1_H_ */
