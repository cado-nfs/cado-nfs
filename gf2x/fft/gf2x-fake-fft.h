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

#ifndef GF2X_FAKE_FFT_FFT_H_
#define GF2X_FAKE_FFT_FFT_H_

#include <stdlib.h>
#include <string.h>

#include "gf2x/gf2x-config-export.h"
#include "gf2x/gf2x-impl-export.h"
#include "gf2x/gf2x-thresholds.h"

/* This file is a placeholder for the typical requirements of an FFT
 * interface. Of course, there is nothing interesting being done here.
 * It's just an E-X-A-M-P-L-E. See also the .c file.
 *
 *
 * A unified adapter for the CXX code is in fft-adapter.hpp
 */

/* TODO
 *
 * There is some space wasted here. We're only copying the data when
 * doing the transforms, so we might as well keep the data as is. This
 * implies that it would perhaps make more sense not to have a separate
 * structure for this case -- although of course this could make the
 * calling interface more complicated.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

struct gf2x_fake_fft_info_struct {
    size_t n1,n2,n3;
    size_t acc;
    size_t size;
};
typedef struct gf2x_fake_fft_info_struct gf2x_fake_fft_info_t[1];
typedef struct gf2x_fake_fft_info_struct * gf2x_fake_fft_info_ptr;
typedef const struct gf2x_fake_fft_info_struct * gf2x_fake_fft_info_srcptr;

typedef unsigned long gf2x_fake_fft_t;
typedef gf2x_fake_fft_t * gf2x_fake_fft_ptr;
typedef const gf2x_fake_fft_t * gf2x_fake_fft_srcptr;

extern void GF2X_EXPORTED gf2x_fake_fft_init(gf2x_fake_fft_info_ptr p, size_t nF, size_t nG, ...); 
static inline void gf2x_fake_fft_clear(gf2x_fake_fft_info_ptr p GF2X_MAYBE_UNUSED) {}
static inline gf2x_fake_fft_ptr gf2x_fake_fft_alloc(gf2x_fake_fft_info_srcptr p, size_t n) {
    return (gf2x_fake_fft_ptr) malloc(n * p->size * sizeof(unsigned long));
}
static inline void gf2x_fake_fft_free(gf2x_fake_fft_info_srcptr p GF2X_MAYBE_UNUSED,
        gf2x_fake_fft_ptr x,
        size_t n GF2X_MAYBE_UNUSED)
{
    free(x);
}
static inline gf2x_fake_fft_ptr gf2x_fake_fft_get(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr x, size_t k) {
    return x + (k * p->size);
}
static inline gf2x_fake_fft_srcptr gf2x_fake_fft_get_const(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_srcptr x, size_t k) {
    return x + (k * p->size);
}
static inline void gf2x_fake_fft_zero(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr x, size_t n)
{
    memset(x, 0, n * p->size * sizeof(unsigned long));
}
extern void GF2X_EXPORTED gf2x_fake_fft_dft(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr x, const unsigned long * F, size_t nF);
extern void GF2X_EXPORTED gf2x_fake_fft_compose(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr y, gf2x_fake_fft_srcptr x1, gf2x_fake_fft_srcptr x2);
extern void GF2X_EXPORTED gf2x_fake_fft_addcompose_n(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr y, gf2x_fake_fft_srcptr *x1, gf2x_fake_fft_srcptr *x2, size_t n);
extern void GF2X_EXPORTED gf2x_fake_fft_addcompose(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr y, gf2x_fake_fft_srcptr x1, gf2x_fake_fft_srcptr x2);
extern void GF2X_EXPORTED gf2x_fake_fft_add(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr y, gf2x_fake_fft_srcptr x1, gf2x_fake_fft_srcptr x2);
extern void GF2X_EXPORTED gf2x_fake_fft_cpy(gf2x_fake_fft_info_srcptr p, gf2x_fake_fft_ptr y, gf2x_fake_fft_srcptr x);
extern void GF2X_EXPORTED gf2x_fake_fft_ift(gf2x_fake_fft_info_srcptr p, unsigned long * H, size_t Hl, gf2x_fake_fft_srcptr h);

extern size_t GF2X_EXPORTED gf2x_fake_fft_size(gf2x_fake_fft_info_srcptr p);

extern void GF2X_EXPORTED gf2x_fake_fft_init_similar(gf2x_fake_fft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_fake_fft_info_srcptr other);
extern int GF2X_EXPORTED gf2x_fake_fft_compatible(gf2x_fake_fft_info_srcptr o1, gf2x_fake_fft_info_srcptr o2);

static inline int gf2x_fake_fft_recoverorder(gf2x_fake_fft_info_srcptr o GF2X_MAYBE_UNUSED) { return 0; }
#ifdef __cplusplus
}
#endif

#endif	/* GF2X_FAKE_FFT_FFT_H_ */
/* vim: set sw=4 sta et: */
