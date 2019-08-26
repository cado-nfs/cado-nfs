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

#ifndef GF2X_FFT_H_
#define GF2X_FFT_H_

/* 64 may also be used */
#ifndef CANTOR_BASE_FIELD_SIZE
#define CANTOR_BASE_FIELD_SIZE 128
#endif

#include "gf2x.h"
#include "gf2x-fake-fft.h"
#include "gf2x-cantor-fft.h"
#include "gf2x-ternary-fft.h"

#ifdef __cplusplus
#include "gf2x-fft-adapter.hpp"
DEFINE_FFT_ADAPTER(gf2x_fake_fft)
DEFINE_FFT_ADAPTER(gf2x_cantor_fft)
DEFINE_FFT_ADAPTER(gf2x_ternary_fft)
#endif

#endif	/* GF2X_FFT_H_ */
