/* 
 * Copyright (C) 2009, 2011 William Hart
 * 
 * This file is part of FLINT.
 * 
 * FLINT is free software: you can redistribute it and/or modify it under the 
 * terms of the GNU Lesser General Public License (LGPL) as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>. */

#include "gmp.h"
#include "flint.h"
#include "fft.h"

/*
 * Set ``s = i1 + i2``, ``t = z1^i*(i1 - i2)`` modulo
 * ``B^limbs + 1`` where ``z1 = exp(Pi*I/n)`` corresponds to
 * multiplication by `2^w`. Requires `0 \leq i < n` where `nw =`
 * ``limbs*FLINT_BITS``.
 * 
 */
void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
		   mp_limb_t * i2, mp_size_t i, mp_size_t limbs,
		   flint_bitcnt_t w)
{
    mp_size_t y;
    flint_bitcnt_t b1;

    b1 = i * w;
    y = b1 / FLINT_BITS;
    b1 = b1 % FLINT_BITS;

    butterfly_lshB(s, t, i1, i2, limbs, 0, y);
    mpn_mul_2expmod_2expp1(t, t, limbs, b1);
}

/*
 * The radix 2 DIF FFT works as follows:
 * 
 * Input: ``[i0, i1, ..., i(m-1)]``, for `m = 2n` a power of `2`.
 * 
 * Output: ``[r0, r1, ..., r(m-1)]``\\ `` = FFT[i0, i1, ..., i(m-1)]``.
 * 
 * Algorithm:
 * 
 * `\bullet` Recursively compute ``[r0, r2, r4, ...., r(m-2)]``\
 *      ``= FFT[i0+i(m/2), i1+i(m/2+1), ..., i(m/2-1)+i(m-1)]``
 * 
 * `\bullet` Let ``[t0, t1, ..., t(m/2-1)]``\
 *      ``= [i0-i(m/2), i1-i(m/2+1), ..., i(m/2-1)-i(m-1)]``
 * 
 * `\bullet` Let ``[u0, u1, ..., u(m/2-1)]``\
 *      ``= [z1^0*t0, z1^1*t1, ..., z1^(m/2-1)*t(m/2-1)]``
 *      where ``z1 = exp(2*Pi*I/m)`` corresponds to multiplication
 *      by `2^w`.
 * 
 * `\bullet` Recursively compute ``[r1, r3, ..., r(m-1)]``\
 *      ``= FFT[u0, u1, ..., u(m/2-1)]``
 * 
 * The parameters are as follows:
 * 
 * `\bullet` ``2*n`` is the length of the input and output
 *      arrays
 * 
 * `\bullet` `w` is such that `2^w` is an `2n`-th root of unity
 *      in the ring `\mathbf{Z}/p\mathbf{Z}` that we are working in,
 *      i.e. `p = 2^{wn} + 1` (here `n` is divisible by
 *      ``GMP_LIMB_BITS``)
 * 
 * `\bullet` ``ii`` is the array of inputs (each input is an
 *      array of limbs of length ``wn/GMP_LIMB_BITS + 1`` (the
 *      extra limbs being a "carry limb"). Outputs are written
 *      in-place.
 * 
 * We require `nw` to be at least 64 and the two temporary space pointers to
 * point to blocks of size ``n*w + FLINT_BITS`` bits.
 * 
 */
void fft_radix2(mp_limb_t ** ii,
		mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1,
		mp_limb_t ** t2)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / GMP_LIMB_BITS;

    if (n == 1) {
	fft_butterfly(*t1, *t2, ii[0], ii[1], 0, limbs, w);

	SWAP_PTRS(ii[0], *t1);
	SWAP_PTRS(ii[1], *t2);

	return;
    }

    for (i = 0; i < n; i++) {
	fft_butterfly(*t1, *t2, ii[i], ii[n + i], i, limbs, w);

	SWAP_PTRS(ii[i], *t1);
	SWAP_PTRS(ii[n + i], *t2);
    }

    fft_radix2(ii, n / 2, 2 * w, t1, t2);
    fft_radix2(ii + n, n / 2, 2 * w, t1, t2);
}
