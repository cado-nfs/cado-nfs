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

#ifndef THRESHOLDS_H_
#define THRESHOLDS_H_

/* This file contains sensible defaults to start with, which are updated
 * by the tuning program */

/* If you read "placeholder" here, it means that the tuning program has
 * not been run (or has not completed) */
#define GF2X_TOOM_TUNING_INFO          "lt-tunetoom -s 1.05 2048 2048 run on cassoulet on Fri Nov 15 14:37:02 2013"
#define GF2X_FFT_TUNING_INFO		"lt-tunefft -s 1.02 8000000 run on cassoulet on Fri Nov 15 15:13:05 2013 ; based on lt-tunetoom -s 1.05 2048 2048 run on cassoulet on Fri Nov 15 14:37:02 2013"

#define GF2X_WORDSIZE           64

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
#define GF2X_MUL_FFT_TABLE		{	\
	{ 1, 1 }, { 9650, 243 }, { 9844, 1 }, { 16163, -243 }, 	\
	{ 17154, 1 }, { 17498, -243 }, { 17848, 1 }, { 18570, -243 }, 	\
	{ 18942, 1 }, { 19321, -243 }, { 19708, 1 }, { 21336, 729 }, 	\
	{ 21763, -729 }, { 22199, -243 }, { 31095, 729 }, { 31717, -729 }, 	\
	{ 33216, -243 }, { 35251, 729 }, { 36677, -243 }, { 38160, -729 }, 	\
	{ 38924, 729 }, { 40135, -729 }, { 58127, -243 }, { 59290, -729 }, 	\
	{ 60476, 2187 }, { 61686, -729 }, { 86387, 729 }, { 88115, -729 }, 	\
	{ 97289, -243 }, { 99235, -729 }, { 118600, -2187 }, { 124557, -729 }, 	\
	{ 134826, 2187 }, { 137523, -2187 }, { 149468, -729 }, { 158619, 2187 }, 	\
	{ 161792, -2187 }, { 224202, -729 }, { 233261, -2187 }, { 323847, 6561 }, 	\
	{ 336303, -2187 }, { 423493, -6561 }, { 448404, -2187 }, { 548049, 6561 }, 	\
	{ 560505, -6561 }, { 672606, -2187 }, { 742615, 6561 }, { 784706, -6561 }, 	\
	{ 896807, 19683 }, { 1008908, -6561 }, { 1793614, -19683 }, { 2017816, -6561 }, 	\
	{ 2690421, 19683 }, { 3026723, -19683 }, { 4035631, -6561 }, { 4708236, -19683 }, }

#endif	/* THRESHOLDS_H_ */

/* vim: set ft=c: */
