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
 * As for ``fft_radix2`` except that only the first ``trunc``
 * coefficients of the output are computed. The transform still needs all
 * `2n` input coefficients to be specified.
 * 
 */
void fft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
		   mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / FLINT_BITS;

    if (trunc == 2 * n)
	fft_radix2(ii, n, w, t1, t2);
    else if (trunc <= n) {
	for (i = 0; i < n; i++)
	    mpn_add_n(ii[i], ii[i], ii[i + n], limbs + 1);

	fft_truncate1(ii, n / 2, 2 * w, t1, t2, trunc);
    } else {
	for (i = 0; i < n; i++) {
	    fft_butterfly(*t1, *t2, ii[i], ii[n + i], i, limbs, w);

	    SWAP_PTRS(ii[i], *t1);
	    SWAP_PTRS(ii[n + i], *t2);
	}

	fft_radix2(ii, n / 2, 2 * w, t1, t2);
	fft_truncate1(ii + n, n / 2, 2 * w, t1, t2, trunc - n);
    }
}

/*
 * As for ``fft_radix2`` except that only the first ``trunc``
 * coefficients of the output are computed and the input is regarded as
 * having (implied) zero coefficients from coefficient ``trunc`` onwards.
 * The coefficients must exist as the algorithm needs to use this extra
 * space, but their value is irrelevant. The value of ``trunc`` must be
 * divisible by 2.
 * 
 */
void fft_truncate(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
		  mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / FLINT_BITS;

    if (trunc == 2 * n)
	fft_radix2(ii, n, w, t1, t2);
    else if (trunc <= n)
	fft_truncate(ii, n / 2, 2 * w, t1, t2, trunc);
    else {
	for (i = 0; i < trunc - n; i++) {
	    fft_butterfly(*t1, *t2, ii[i], ii[n + i], i, limbs, w);

	    SWAP_PTRS(ii[i], *t1);
	    SWAP_PTRS(ii[n + i], *t2);
	}

	for (; i < n; i++)
	    fft_adjust(ii[i + n], ii[i], i, limbs, w);

	fft_radix2(ii, n / 2, 2 * w, t1, t2);
	fft_truncate1(ii + n, n / 2, 2 * w, t1, t2, trunc - n);
    }
}
