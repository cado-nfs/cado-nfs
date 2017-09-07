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
/* An implementation of Cantor's algorithm for multiplication of
   Copyright 2007 Pierrick Gaudry.
   Copyright 2008,2009,2010,2012 Emmanuel Thomé.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#ifndef GF2X_FFT_CANTOR_H_
#define GF2X_FFT_CANTOR_H_

#include <stdint.h>
#include <stdlib.h>

#include "gf2x.h"
#include "gf2x/gf2x-config-export.h"
#include "gf2x/gf2x-thresholds.h"       /* GF2X_WORDSIZE */

#ifdef __cplusplus
extern "C" {
#endif

// we must implement some sort of fallback on the fake_fft thing when the
// polynomials are too small, or something like this. The gf2x-fft
// interface requires that any size can be used.

struct gf2x_cantor_fft_info_struct {
    unsigned int k;
    size_t n;
};
typedef struct gf2x_cantor_fft_info_struct gf2x_cantor_fft_info_t[1];
typedef struct gf2x_cantor_fft_info_struct * gf2x_cantor_fft_info_ptr;
typedef const struct gf2x_cantor_fft_info_struct * gf2x_cantor_fft_info_srcptr;

#if !defined(GF2X_WORDSIZE) || (GF2X_WORDSIZE != 32 && GF2X_WORDSIZE != 64)
#error  "define GF2X_WORDSIZE to either 32 or 64"
#endif

#if !defined(CANTOR_BASE_FIELD_SIZE) || (CANTOR_BASE_FIELD_SIZE != 64 && CANTOR_BASE_FIELD_SIZE != 128)
#error  "Define CANTOR_BASE_FIELD_SIZE to 64 or 128"
#endif

#if CANTOR_BASE_FIELD_SIZE == 128
#if GF2X_WORDSIZE == 64
#include "mpfq/x86_64/mpfq_2_128.h"
#else
#include "mpfq/i386/mpfq_2_128.h"
#endif
typedef mpfq_2_128_elt gf2x_cantor_fft_base_field_elt;
#else
#if GF2X_WORDSIZE == 64
#include "mpfq/x86_64/mpfq_2_64.h"
#else
#include "mpfq/i386/mpfq_2_64.h"
#endif
typedef mpfq_2_64_elt gf2x_cantor_fft_base_field_elt;
#endif
typedef gf2x_cantor_fft_base_field_elt * gf2x_cantor_fft_ptr;
typedef const gf2x_cantor_fft_base_field_elt * gf2x_cantor_fft_srcptr;

/* sets up the appropriate data for multiplying polynomials having
 * numbers of coefficients bounded respectively by nF and nG. This does
 * not allocate a transform in itself: only the descriptor structure that
 * will be useful for allocating transforms with _alloc.
 */
extern void gf2x_cantor_fft_init(gf2x_cantor_fft_info_t p, size_t nF, size_t nG, ...);
/* paired with previous function. Zeroes out the initialization that was
 * done. */
extern void gf2x_cantor_fft_clear(gf2x_cantor_fft_info_t p GF2X_MAYBE_UNUSED);
/* allocates space for n transforms, according to the transform info at p */
extern gf2x_cantor_fft_ptr gf2x_cantor_fft_alloc(const gf2x_cantor_fft_info_t p, size_t n);
/* frees the n transforms allocated at x */
extern void gf2x_cantor_fft_free(
        const gf2x_cantor_fft_info_t p GF2X_MAYBE_UNUSED,
        gf2x_cantor_fft_ptr x,
        size_t n GF2X_MAYBE_UNUSED);
/* accesses the k-th transform stored at location x. */
extern gf2x_cantor_fft_ptr gf2x_cantor_fft_get(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, size_t k);
/* accesses the k-th transform stored at location x. */
extern gf2x_cantor_fft_srcptr gf2x_cantor_fft_get_const(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_srcptr x, size_t k);
/* initialize to zero the n transforms allocated at x */
extern void gf2x_cantor_fft_zero(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, size_t n);
/* performs a direct transform of F, having at most nF coefficients. */
extern void gf2x_cantor_fft_dft(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr x, const unsigned long * F, size_t nF);
/* composes the transforms of two polynomials. The result is by
 * definition the transform of their product. */
extern void gf2x_cantor_fft_compose(const gf2x_cantor_fft_info_t p,
		gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2);
/* compose the n transforms at x1[0..[ and x2[0..[, accumulates the sum
 * of the result to y[0..[ */
extern void gf2x_cantor_fft_addcompose_n(const gf2x_cantor_fft_info_t p,
		gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr * x1, gf2x_cantor_fft_srcptr * x2, size_t);
/* compose the n transforms at x1 and x2, accumulates the result to
 * y[0..[ */
extern void gf2x_cantor_fft_addcompose(const gf2x_cantor_fft_info_t p,
		gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2);
/* adds two transforms (the direct transform is linear). */
extern void gf2x_cantor_fft_add(const gf2x_cantor_fft_info_t p,
		gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x1, gf2x_cantor_fft_srcptr x2);
/* copy */
extern void gf2x_cantor_fft_cpy(const gf2x_cantor_fft_info_t p, gf2x_cantor_fft_ptr y, gf2x_cantor_fft_srcptr x);
/* performs an inverse transform. nH is the number of coefficients. */
extern void gf2x_cantor_fft_ift(const gf2x_cantor_fft_info_t p,
		unsigned long * H, size_t Hl, gf2x_cantor_fft_ptr h);
extern size_t gf2x_cantor_fft_size(gf2x_cantor_fft_info_srcptr p);
extern void gf2x_cantor_fft_init_similar(gf2x_cantor_fft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_cantor_fft_info_srcptr other);
extern int gf2x_cantor_fft_compatible(gf2x_cantor_fft_info_srcptr o1, gf2x_cantor_fft_info_srcptr o2);
static inline int gf2x_cantor_fft_recoverorder(gf2x_cantor_fft_info_srcptr o) { return o->k; }

#ifdef __cplusplus
}
#endif

#endif	/* GF2X_FFT_CANTOR_H_ */
