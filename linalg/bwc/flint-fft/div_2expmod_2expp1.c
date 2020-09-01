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
#include "longlong.h"

/* WARNING: relies on GCC's handling of >> as arithmetic shift right */

/*
 * Given ``i1`` a signed integer of ``limbs + 1`` limbs in twos
 * complement format reduced modulo ``B^limbs + 1`` up to some
 * overflow, compute ``t = i1/2^d`` modulo `p`. The result will not
 * necessarily be fully reduced. The number of bits ``d`` must be
 * nonnegative and less than ``FLINT_BITS``. Aliasing is permitted.
 * 
 * 
 * 
 */
void mpn_div_2expmod_2expp1(mp_limb_t * t, mp_limb_t * i1, mp_size_t limbs,
			    flint_bitcnt_t d)
{
    mp_limb_t lo;
    mp_limb_t *ptr;
    mp_limb_signed_t hi;

    if (d == 0) {
	if (t != i1)
	    flint_mpn_copyi(t, i1, limbs + 1);
    } else {
	hi = i1[limbs];
	lo = mpn_rshift(t, i1, limbs + 1, d);
	t[limbs] = (hi >> d);
	ptr = t + limbs - 1;
	sub_ddmmss(ptr[1], ptr[0], ptr[1], ptr[0], (mp_limb_t) 0, lo);
    }
}
