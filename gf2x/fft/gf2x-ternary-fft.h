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


#ifndef GF2X_TERNARY_FFT_H_
#define GF2X_TERNARY_FFT_H_

#ifdef __cplusplus
extern "C" {
#endif

struct gf2x_ternary_fft_info_s {
    size_t bits_a;  // number of bits of operand1
    size_t bits_b;  // number of bits of operand2
    size_t K;       // 0 indicates fallback.
    size_t M;
    unsigned long * tmp;
    size_t * perm;
    int split;  // boolean
};

typedef struct gf2x_ternary_fft_info_s gf2x_ternary_fft_info_t[1];
typedef struct gf2x_ternary_fft_info_s * gf2x_ternary_fft_info_ptr;
typedef const struct gf2x_ternary_fft_info_s * gf2x_ternary_fft_info_srcptr;

typedef unsigned long gf2x_ternary_fft_t;
typedef gf2x_ternary_fft_t * gf2x_ternary_fft_ptr;
typedef const gf2x_ternary_fft_t * gf2x_ternary_fft_srcptr;

extern size_t GF2X_EXPORTED gf2x_ternary_fft_size(gf2x_ternary_fft_info_srcptr o);
extern gf2x_ternary_fft_ptr GF2X_EXPORTED gf2x_ternary_fft_get(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr ptr, size_t k);
extern gf2x_ternary_fft_srcptr GF2X_EXPORTED gf2x_ternary_fft_get_const(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_srcptr ptr, size_t k);
extern void GF2X_EXPORTED gf2x_ternary_fft_zero(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr ptr, size_t n);
extern void GF2X_EXPORTED gf2x_ternary_fft_cpy(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr y, gf2x_ternary_fft_srcptr x);
extern gf2x_ternary_fft_ptr GF2X_EXPORTED gf2x_ternary_fft_alloc(gf2x_ternary_fft_info_srcptr o, size_t n);
extern void GF2X_EXPORTED gf2x_ternary_fft_free(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr ptr, size_t n);
extern void GF2X_EXPORTED gf2x_ternary_fft_dft(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr tr, const unsigned long * a, size_t bits_a);
extern void GF2X_EXPORTED gf2x_ternary_fft_compose(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr tc, gf2x_ternary_fft_srcptr ta, gf2x_ternary_fft_srcptr tb);
extern void GF2X_EXPORTED gf2x_ternary_fft_addcompose_n(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr tc, gf2x_ternary_fft_srcptr * ta, gf2x_ternary_fft_srcptr * tb, size_t n);
extern void GF2X_EXPORTED gf2x_ternary_fft_addcompose(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr tc, gf2x_ternary_fft_srcptr ta, gf2x_ternary_fft_srcptr tb);
extern void GF2X_EXPORTED gf2x_ternary_fft_add(gf2x_ternary_fft_info_srcptr o, gf2x_ternary_fft_ptr tc, gf2x_ternary_fft_srcptr ta, gf2x_ternary_fft_srcptr tb);
extern void GF2X_EXPORTED gf2x_ternary_fft_ift(gf2x_ternary_fft_info_srcptr o, unsigned long * c, size_t bits_c, gf2x_ternary_fft_ptr tr);
extern void GF2X_EXPORTED gf2x_ternary_fft_init(gf2x_ternary_fft_info_ptr o, size_t bits_a, size_t bits_b, ...);
extern void GF2X_EXPORTED gf2x_ternary_fft_init_similar(gf2x_ternary_fft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_ternary_fft_info_srcptr other);
extern int GF2X_EXPORTED gf2x_ternary_fft_compatible(gf2x_ternary_fft_info_srcptr o1, gf2x_ternary_fft_info_srcptr o2);
extern void GF2X_EXPORTED gf2x_ternary_fft_clear(gf2x_ternary_fft_info_ptr o);


static inline int gf2x_ternary_fft_recoverorder(gf2x_ternary_fft_info_srcptr o) { return ((o)->K * ((o)->split ? -1 : 1)); }

#ifdef __cplusplus
}
#endif

#endif	/* GF2X_TERNARY_FFT_H_ */
