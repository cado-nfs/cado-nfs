/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015, 2019
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

#ifndef THRESHOLDS_H_
#define THRESHOLDS_H_

/* This file contains sensible defaults to start with, which are updated
 * by the tuning program */

/* If you read "placeholder" here, it means that the tuning program has
 * not been run (or has not completed) */
#define	GF2X_TOOM_TUNING_INFO		"placeholder"
#define GF2X_FFT_TUNING_INFO		"placeholder"

#define GF2X_WORDSIZE           32

#define GF2X_STORAGE_CLASS_mul1        static
#define GF2X_STORAGE_CLASS_mul_1_n     static
#define GF2X_STORAGE_CLASS_addmul_1_n  static
#define GF2X_STORAGE_CLASS_mul2        static
#define GF2X_STORAGE_CLASS_mul3        static
#define GF2X_STORAGE_CLASS_mul4        static
#define GF2X_STORAGE_CLASS_mul5        static
#define GF2X_STORAGE_CLASS_mul6        static
#define GF2X_STORAGE_CLASS_mul7        static
#define GF2X_STORAGE_CLASS_mul8        static
#define GF2X_STORAGE_CLASS_mul9        static

/* First size for which KARA is used. Essentially hard-coded, since the
 * sizes up to 9 words are already karatsuba, unrolled. The unrolled
 * routines handle their own temporary storage on the stack.
 */
#define GF2X_MUL_KARA_THRESHOLD	10

/* First size for which TC3W is used. It is assumed that TC3W is used
 * before TC3 kicks in. */
/* must be >= 8 */
#define GF2X_MUL_TOOMW_THRESHOLD		18

/* First size for which TC3 is used. This threshold is informative, it is
 * not used in the code */
/* must be >= 17 */
#define GF2X_MUL_TOOM_THRESHOLD		21

/* First size for which TC4 is used. This threshold is informative, it is
 * not used in the code */
/* must be >= 30 */
#define GF2X_MUL_TOOM4_THRESHOLD		235

/* Size above which TC4 is the only TC variant used. */
#define GF2X_MUL_TOOM4_ALWAYS_THRESHOLD		1600

/* First size for which TC3U is used */
/* must be >= 33 */
#define GF2X_MUL_TOOMU_THRESHOLD		49

/* Size above which TC3U is the only TCU variant used. */
#define GF2X_MUL_TOOMU_ALWAYS_THRESHOLD		2049

/* The default values here are appropriate for the tuning program.
 * Appropriate values are substituted later on. Note that the tuning
 * table is always created with this size, to that re-tuning is possible
 * at any time. */
#define	GF2X_TOOM_TUNING_LIMIT	2048
#define GF2X_BEST_TOOM_TABLE {}

#define	GF2X_BEST_UTOOM_TABLE {}

/* {n, K} means use FFT(|K|) up from n words, */
/* where |K|<3 stands for Toom-Cook 3, K < 0 means use FFT2 */

/* these thresholds are just plain crap, but probably better this than
 * nothing (which would amount to disabling fft altogether) */
#define GF2X_MUL_FFT_TABLE		{	\
	{ 1, 1 }, { 3516, -81 }, { 3587, 1 }, { 3733, -81 }, 	\
	{ 3808, 1 }, { 4614, 243 }, { 4921, 1 }, { 5020, -81 }, }

#endif	/* THRESHOLDS_H_ */

/* vim: set ft=c: */
