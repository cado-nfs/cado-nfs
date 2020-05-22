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

#ifndef CHECK_AND_BENCH_COMMON_H_
#define CHECK_AND_BENCH_COMMON_H_

#include "gf2x/gf2x-config-export.h"
#include "gf2x/gf2x-impl-export.h"

extern int time_total;
extern int time_dft;
extern int time_conv;
extern int time_ift;

#if defined(ENGINE_CANTOR)
#include "gf2x-cantor-fft.h"

#define ENGINE_EXTRA_ARG_DEFAULT  0     /* does not make sense here (yet) */

#define ENGINE_info_t gf2x_cantor_fft_info_t
#define ENGINE_t gf2x_cantor_fft_elt
#define ENGINE_get_alloc_sizes gf2x_cantor_fft_info_get_alloc_sizes
#define ENGINE_srcptr gf2x_cantor_fft_srcptr
#define ENGINE_ptr gf2x_cantor_fft_ptr
#define ENGINE_init gf2x_cantor_fft_info_init
#define ENGINE_adjust gf2x_cantor_fft_info_adjust
#define ENGINE_alloc gf2x_cantor_fft_alloc
#define ENGINE_dft gf2x_cantor_fft_dft
#define ENGINE_zero gf2x_cantor_fft_zero
#define ENGINE_compose gf2x_cantor_fft_compose
#define ENGINE_addcompose gf2x_cantor_fft_addcompose
#define ENGINE_addcompose_n gf2x_cantor_fft_addcompose_n
#define ENGINE_free gf2x_cantor_fft_free
#define ENGINE_ift gf2x_cantor_fft_ift
#define ENGINE_get gf2x_cantor_fft_get
#define ENGINE_get_const gf2x_cantor_fft_get_const
#define ENGINE_clear gf2x_cantor_fft_info_clear
#define ENGINE_order gf2x_cantor_fft_info_order 
#define ENGINE_mul      test_gf2x_cantor_fft_mul
#elif defined(ENGINE_TERNARY)
#include "gf2x-ternary-fft.h"
#ifdef MULCOUNT
#error "MULCOUNT incompatible with gf2x_ternary_fft"
#endif

#define ENGINE_EXTRA_ARG_DEFAULT  0

#define ENGINE_info_t gf2x_ternary_fft_info_t
#define ENGINE_t gf2x_ternary_fft_elt
#define ENGINE_get_alloc_sizes gf2x_ternary_fft_info_get_alloc_sizes
#define ENGINE_srcptr gf2x_ternary_fft_srcptr
#define ENGINE_ptr gf2x_ternary_fft_ptr
#define ENGINE_init gf2x_ternary_fft_info_init
#define ENGINE_adjust gf2x_ternary_fft_info_adjust
#define ENGINE_alloc gf2x_ternary_fft_alloc
#define ENGINE_dft gf2x_ternary_fft_dft
#define ENGINE_zero gf2x_ternary_fft_zero
#define ENGINE_compose gf2x_ternary_fft_compose
#define ENGINE_addcompose gf2x_ternary_fft_addcompose
#define ENGINE_addcompose_n gf2x_ternary_fft_addcompose_n
#define ENGINE_free gf2x_ternary_fft_free
#define ENGINE_ift gf2x_ternary_fft_ift
#define ENGINE_get gf2x_ternary_fft_get
#define ENGINE_get_const gf2x_ternary_fft_get_const
#define ENGINE_clear gf2x_ternary_fft_info_clear
#define ENGINE_order gf2x_ternary_fft_info_order
#define ENGINE_mul      test_gf2x_ternary_fft_mul
#else
#error "Please define either ENGINE_CANTOR or ENGINE_TERNARY"
#endif

#ifdef __cplusplus
extern "C" {
#endif

long ENGINE_mul(unsigned long ** H, unsigned long ** F, size_t Fl, unsigned long ** G, size_t Gl, int n);

void set_extra_arg_from_n32bitwords(size_t n32 GF2X_MAYBE_UNUSED, long supplied);

#ifdef __cplusplus
}
#endif

#endif	/* CHECK_AND_BENCH_COMMON_H_ */
