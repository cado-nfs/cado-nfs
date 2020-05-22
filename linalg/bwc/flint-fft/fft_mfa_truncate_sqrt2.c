#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
/* flint uses unprotected openmp pragmas every so often */
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
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
#include "ulong_extras.h"
#include "fft.h"

/*
 * Set ``u = 2^b1*(s + t)``, ``v = 2^b2*(s - t)`` modulo 
 * ``B^limbs + 1``. This is used to compute 
 * ``u = 2^(ws*tw1)*(s + t)``,\\ ``v = 2^(w+ws*tw2)*(s - t)`` in the 
 * matrix fourier algorithm, i.e. effectively computing an ordinary butterfly 
 * with additional twiddles by ``z1^rc`` for row `r` and column `c` of the 
 * matrix of coefficients. Aliasing is not allowed.
 * 
 */
void fft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
			   mp_limb_t * s, mp_limb_t * t, mp_size_t limbs,
			   flint_bitcnt_t b1, flint_bitcnt_t b2)
{
    mp_limb_t nw = limbs * FLINT_BITS;
    mp_size_t x, y;
    int negate1 = 0;
    int negate2 = 0;

    if (b1 >= nw) {
	negate2 = 1;
	b1 -= nw;
    }
    x = b1 / FLINT_BITS;
    b1 = b1 % FLINT_BITS;

    if (b2 >= nw) {
	negate1 = 1;
	b2 -= nw;
    }
    y = b2 / FLINT_BITS;
    b2 = b2 % FLINT_BITS;

    butterfly_lshB(u, v, s, t, limbs, x, y);
    mpn_mul_2expmod_2expp1(u, u, limbs, b1);
    if (negate2)
	mpn_neg_n(u, u, limbs + 1);
    mpn_mul_2expmod_2expp1(v, v, limbs, b2);
    if (negate1)
	mpn_neg_n(v, v, limbs + 1);
}

/*
 * As for ``fft_radix2`` except that the coefficients are spaced by 
 * ``is`` in the array ``ii`` and an additional twist by ``z^c*i``
 * is applied to each coefficient where `i` starts at `r` and increases by
 * ``rs`` as one moves from one coefficient to the next. Here ``z`` 
 * corresponds to multiplication by ``2^ws``. 
 * 
 */
void fft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
			mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1,
			mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			mp_size_t c, mp_size_t rs)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / FLINT_BITS;

    if (n == 1) {
	mp_size_t tw1 = r * c;
	mp_size_t tw2 = tw1 + rs * c;
	fft_butterfly_twiddle(*t1, *t2, ii[0], ii[is], limbs, tw1 * ws,
			      tw2 * ws);

	SWAP_PTRS(ii[0], *t1);
	SWAP_PTRS(ii[is], *t2);

	return;
    }

    for (i = 0; i < n; i++) {
	fft_butterfly(*t1, *t2, ii[i * is], ii[(n + i) * is], i, limbs, w);

	SWAP_PTRS(ii[i * is], *t1);
	SWAP_PTRS(ii[(n + i) * is], *t2);
    }

    fft_radix2_twiddle(ii, is, n / 2, 2 * w, t1, t2, ws, r, c, 2 * rs);
    fft_radix2_twiddle(ii + n * is, is, n / 2, 2 * w, t1, t2, ws, r + rs, c,
		       2 * rs);
}

/*
 * As per ``fft_radix2_twiddle`` except that the transform is truncated 
 * as per\\ ``fft_truncate1``.
 * 
 */
void fft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
			   mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1,
			   mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			   mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w * n) / FLINT_BITS;

    if (trunc == 2 * n)
	fft_radix2_twiddle(ii, is, n, w, t1, t2, ws, r, c, rs);
    else if (trunc <= n) {
	for (i = 0; i < n; i++)
	    mpn_add_n(ii[i * is], ii[i * is], ii[(i + n) * is], limbs + 1);

	fft_truncate1_twiddle(ii, is, n / 2, 2 * w, t1, t2, ws, r, c, 2 * rs,
			      trunc);
    } else {
	for (i = 0; i < n; i++) {
	    fft_butterfly(*t1, *t2, ii[i * is], ii[(n + i) * is], i, limbs,
			  w);

	    SWAP_PTRS(ii[i * is], *t1);
	    SWAP_PTRS(ii[(n + i) * is], *t2);
	}

	fft_radix2_twiddle(ii, is, n / 2, 2 * w, t1, t2, ws, r, c, 2 * rs);
	fft_truncate1_twiddle(ii + n * is, is, n / 2, 2 * w,
			      t1, t2, ws, r + rs, c, 2 * rs, trunc - n);
    }
}

/*
 * This is as per the ``fft_truncate_sqrt2`` function except that the 
 * matrix fourier algorithm is used for the left and right FFTs. The total 
 * transform length is `4n` where ``n = 2^depth`` so that the left and
 * right transforms are both length `2n`. We require ``trunc > 2*n`` and
 * that ``trunc`` is divisible by ``2*n1`` (explained below).
 * 
 * The matrix fourier algorithm, which is applied to each transform of length
 * `2n`, works as follows. We set ``n1`` to a power of 2 about the square
 * root of `n`. The data is then thought of as a set of ``n2`` rows each
 * with ``n1`` columns (so that ``n1*n2 = 2n``). 
 * 
 * The length `2n` transform is then computed using a whole pile of short 
 * transforms. These comprise ``n1`` column transforms of length ``n2``
 * followed by some twiddles by roots of unity (namely ``z^rc`` where `r` 
 * is the row and `c` the column within the data) followed by ``n2``
 * row transforms of length ``n1``. Along the way the data needs to be
 * rearranged due to the fact that the short transforms output the data in 
 * binary reversed order compared with what is needed.
 * 
 * The matrix fourier algorithm provides better cache locality by decomposing
 * the long length `2n` transforms into many transforms of about the square 
 * root of the original length. 
 * 
 * For better cache locality the sqrt2 layer of the full length `4n` 
 * transform is folded in with the column FFTs performed as part of the first 
 * matrix fourier algorithm on the left half of the data.
 * 
 * The second half of the data requires a truncated version of the matrix
 * fourier algorithm. This is achieved by truncating to an exact multiple of 
 * the row length so that the row transforms are full length. Moreover, the 
 * column transforms will then be truncated transforms and their truncated 
 * length needs to be a multiple of 2. This explains the condition on 
 * ``trunc`` given above. 
 * 
 * To improve performance, the extra twiddles by roots of unity are combined
 * with the butterflies performed at the last layer of the column transforms.
 * 
 * We require `nw` to be at least 64 and the three temporary space pointers 
 * to point to blocks of size ``n*w + FLINT_BITS`` bits.
 * 
 */
void fft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
			    flint_bitcnt_t w, mp_limb_t ** t1,
			    mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1,
			    mp_size_t trunc)
{
    mp_size_t i, j, s;
    mp_size_t n2 = (2 * n) / n1;
    mp_size_t trunc2 = (trunc - 2 * n) / n1;
    mp_size_t limbs = (n * w) / FLINT_BITS;
    flint_bitcnt_t depth = 0;
    flint_bitcnt_t depth2 = 0;

    while ((UWORD(1) << depth) < n2)
	depth++;
    while ((UWORD(1) << depth2) < n1)
	depth2++;

    /* first half matrix fourier FFT : n2 rows, n1 cols */

    /* FFTs on columns */
    for (i = 0; i < n1; i++) {
	/* relevant part of first layer of full sqrt2 FFT */
	if (w & 1) {
	    for (j = i; j < trunc - 2 * n; j += n1) {
		if (j & 1)
		    fft_butterfly_sqrt2(*t1, *t2, ii[j], ii[2 * n + j], j,
					limbs, w, *temp);
		else
		    fft_butterfly(*t1, *t2, ii[j], ii[2 * n + j], j / 2,
				  limbs, w);

		SWAP_PTRS(ii[j], *t1);
		SWAP_PTRS(ii[2 * n + j], *t2);
	    }

	    for (; j < 2 * n; j += n1) {
		if (i & 1)
		    fft_adjust_sqrt2(ii[j + 2 * n], ii[j], j, limbs, w,
				     *temp);
		else
		    fft_adjust(ii[j + 2 * n], ii[j], j / 2, limbs, w);
	    }
	} else {
	    for (j = i; j < trunc - 2 * n; j += n1) {
		fft_butterfly(*t1, *t2, ii[j], ii[2 * n + j], j, limbs,
			      w / 2);

		SWAP_PTRS(ii[j], *t1);
		SWAP_PTRS(ii[2 * n + j], *t2);
	    }

	    for (; j < 2 * n; j += n1)
		fft_adjust(ii[j + 2 * n], ii[j], j, limbs, w / 2);
	}

	/* 
	 * FFT of length n2 on column i, applying z^{r*i} for rows going up
	 * in steps of 1 starting at row 0, where z => w bits */

	fft_radix2_twiddle(ii + i, n1, n2 / 2, w * n1, t1, t2, w, 0, i, 1);
	for (j = 0; j < n2; j++) {
	    mp_size_t s = n_revbin(j, depth);
	    if (j < s)
		SWAP_PTRS(ii[i + j * n1], ii[i + s * n1]);
	}
    }

    /* FFTs on rows */
    for (i = 0; i < n2; i++) {
	fft_radix2(ii + i * n1, n1 / 2, w * n2, t1, t2);
	for (j = 0; j < n1; j++) {
	    mp_size_t t = n_revbin(j, depth2);
	    if (j < t)
		SWAP_PTRS(ii[i * n1 + j], ii[i * n1 + t]);
	}
    }

    /* second half matrix fourier FFT : n2 rows, n1 cols */
    ii += 2 * n;

    /* FFTs on columns */
    for (i = 0; i < n1; i++) {
	/* 
	 * FFT of length n2 on column i, applying z^{r*i} for rows going up
	 * in steps of 1 starting at row 0, where z => w bits */

	fft_truncate1_twiddle(ii + i, n1, n2 / 2, w * n1, t1, t2, w, 0, i, 1,
			      trunc2);
	for (j = 0; j < n2; j++) {
	    mp_size_t s = n_revbin(j, depth);
	    if (j < s)
		SWAP_PTRS(ii[i + j * n1], ii[i + s * n1]);
	}
    }

    /* FFTs on relevant rows */
    for (s = 0; s < trunc2; s++) {
	i = n_revbin(s, depth);
	fft_radix2(ii + i * n1, n1 / 2, w * n2, t1, t2);

	for (j = 0; j < n1; j++) {
	    mp_size_t t = n_revbin(j, depth2);
	    if (j < t)
		SWAP_PTRS(ii[i * n1 + j], ii[i * n1 + t]);
	}
    }
}

/*
 * Just the outer layers of ``fft_mfa_truncate_sqrt2``.
 * 
 */
void fft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
				  flint_bitcnt_t w, mp_limb_t ** t1,
				  mp_limb_t ** t2, mp_limb_t ** temp,
				  mp_size_t n1, mp_size_t trunc)
{
    mp_size_t i, j;
    mp_size_t n2 = (2 * n) / n1;
    mp_size_t trunc2 = (trunc - 2 * n) / n1;
    mp_size_t limbs = (n * w) / FLINT_BITS;
    flint_bitcnt_t depth = 0;
    flint_bitcnt_t depth2 = 0;
    int k = 0;

    while ((UWORD(1) << depth) < n2)
	depth++;
    while ((UWORD(1) << depth2) < n1)
	depth2++;

    /* first half matrix fourier FFT : n2 rows, n1 cols */

    /* FFTs on columns */

#pragma omp parallel for private(i, j, k)
    for (i = 0; i < n1; i++) {
#if defined(_OPENMP)
	k = omp_get_thread_num();
#endif
	/* relevant part of first layer of full sqrt2 FFT */
	if (w & 1) {
	    for (j = i; j < trunc - 2 * n; j += n1) {
		if (j & 1)
		    fft_butterfly_sqrt2(t1[k], t2[k], ii[j], ii[2 * n + j], j,
					limbs, w, temp[k]);
		else
		    fft_butterfly(t1[k], t2[k], ii[j], ii[2 * n + j], j / 2,
				  limbs, w);

		SWAP_PTRS(ii[j], t1[k]);
		SWAP_PTRS(ii[2 * n + j], t2[k]);
	    }

	    for (; j < 2 * n; j += n1) {
		if (i & 1)
		    fft_adjust_sqrt2(ii[j + 2 * n], ii[j], j, limbs, w,
				     temp[k]);
		else
		    fft_adjust(ii[j + 2 * n], ii[j], j / 2, limbs, w);
	    }
	} else {
	    for (j = i; j < trunc - 2 * n; j += n1) {
		fft_butterfly(t1[k], t2[k], ii[j], ii[2 * n + j], j, limbs,
			      w / 2);

		SWAP_PTRS(ii[j], t1[k]);
		SWAP_PTRS(ii[2 * n + j], t2[k]);
	    }

	    for (; j < 2 * n; j += n1)
		fft_adjust(ii[j + 2 * n], ii[j], j, limbs, w / 2);
	}

	/* 
	 * FFT of length n2 on column i, applying z^{r*i} for rows going up
	 * in steps of 1 starting at row 0, where z => w bits */

	fft_radix2_twiddle(ii + i, n1, n2 / 2, w * n1, t1 + k, t2 + k, w, 0,
			   i, 1);
	for (j = 0; j < n2; j++) {
	    mp_size_t s = n_revbin(j, depth);
	    if (j < s)
		SWAP_PTRS(ii[i + j * n1], ii[i + s * n1]);
	}
    }

    /* second half matrix fourier FFT : n2 rows, n1 cols */
    ii += 2 * n;

    /* FFTs on columns */

#pragma omp parallel for private(i, j, k)
    for (i = 0; i < n1; i++) {
#if defined(_OPENMP)
	k = omp_get_thread_num();
#endif
	/* 
	 * FFT of length n2 on column i, applying z^{r*i} for rows going up
	 * in steps of 1 starting at row 0, where z => w bits */

	fft_truncate1_twiddle(ii + i, n1, n2 / 2, w * n1, t1 + k, t2 + k, w,
			      0, i, 1, trunc2);
	for (j = 0; j < n2; j++) {
	    mp_size_t s = n_revbin(j, depth);
	    if (j < s)
		SWAP_PTRS(ii[i + j * n1], ii[i + s * n1]);
	}
    }
}
