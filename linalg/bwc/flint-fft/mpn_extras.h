#ifndef MPN_EXTRAS_H_
#define MPN_EXTRAS_H_

#include "gmp.h"
#include "flint.h"

#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp,
				     int c, mp_bitcnt_t b, mp_ptr tp);

#endif				/* MPN_EXTRAS_H_ */
