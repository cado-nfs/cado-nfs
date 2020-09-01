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
 * Let `w = 2k + 1`, `i = 2j + 1`. Set ``s = i1 + i2``, 
 * ``t = z1^i*(i1 - i2)`` modulo ``B^limbs + 1`` where 
 * ``z1^2 = exp(Pi*I/n)`` corresponds to multiplication by `2^w`. Requires 
 * `0 \leq i < 2n` where `nw =` ``limbs*FLINT_BITS``.
 * 
 * Here ``z1`` corresponds to multiplication by `2^k` then multiplication
 * by\\ ``(2^(3nw/4) - 2^(nw/4))``. We see ``z1^i`` corresponds to
 * multiplication by ``(2^(3nw/4) - 2^(nw/4))*2^(j+ik)``.
 * 
 * We first multiply by ``2^(j + ik + wn/4)`` then multiply by an
 * additional ``2^(nw/2)`` and subtract.
 * 
 */
void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t,
			 mp_limb_t * i1, mp_limb_t * i2, mp_size_t i,
			 mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)
{
    flint_bitcnt_t wn = limbs * FLINT_BITS;
    mp_limb_t cy = 0;
    mp_size_t j = i / 2, k = w / 2;
    mp_size_t y;
    flint_bitcnt_t b1;
    int negate = 0;

    b1 = j + wn / 4 + i * k;
    if (b1 >= wn) {
	negate = 1;
	b1 -= wn;
    }
    y = b1 / FLINT_BITS;
    b1 = b1 % FLINT_BITS;

    /* sumdiff and multiply by 2^{j + wn/4 + i*k} */
    butterfly_lshB(s, t, i1, i2, limbs, 0, y);
    mpn_mul_2expmod_2expp1(t, t, limbs, b1);

    /* multiply by 2^{wn/2} */
    y = limbs / 2;

    flint_mpn_copyi(temp + y, t, limbs - y);
    temp[limbs] = 0;
    if (y)
	cy = mpn_neg_n(temp, t + limbs - y, y);
    mpn_addmod_2expp1_1(temp + y, limbs - y, -t[limbs]);
    mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy);

    /* shift by an additional half limb (rare) */
    if (limbs & 1)
	mpn_mul_2expmod_2expp1(temp, temp, limbs, FLINT_BITS / 2);

    /* subtract */
    if (negate)
	mpn_sub_n(t, t, temp, limbs + 1);
    else
	mpn_sub_n(t, temp, t, limbs + 1);
}

/*
 * As per ``fft_truncate`` except that the transform is twice the usual 
 * length, i.e. length `4n` rather than `2n`. This is achieved by making use 
 * of twiddles by powers of a square root of 2, not powers of 2 in the first 
 * layer of the transform.  
 * 
 * We require `nw` to be at least 64 and the three temporary space pointers 
 * to point to blocks of size ``n*w + FLINT_BITS`` bits.
 * 
 */
void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
			mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
			mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / GMP_LIMB_BITS;

    if ((w & 1) == 0) {
	fft_truncate(ii, 2 * n, w / 2, t1, t2, trunc);
	return;
    }

    for (i = 0; i < trunc - 2 * n; i++) {
	fft_butterfly(*t1, *t2, ii[i], ii[2 * n + i], i / 2, limbs, w);

	SWAP_PTRS(ii[i], *t1);
	SWAP_PTRS(ii[i + 2 * n], *t2);

	i++;

	fft_butterfly_sqrt2(*t1, *t2, ii[i], ii[2 * n + i], i, limbs, w,
			    *temp);

	SWAP_PTRS(ii[i], *t1);
	SWAP_PTRS(ii[2 * n + i], *t2);
    }

    for (i = trunc - 2 * n; i < 2 * n; i++) {
	fft_adjust(ii[i + 2 * n], ii[i], i / 2, limbs, w);

	i++;

	fft_adjust_sqrt2(ii[i + 2 * n], ii[i], i, limbs, w, *temp);
    }

    fft_radix2(ii, n, w, t1, t2);
    fft_truncate1(ii + 2 * n, n, w, t1, t2, trunc - 2 * n);
}
