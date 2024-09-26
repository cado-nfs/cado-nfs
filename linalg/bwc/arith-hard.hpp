#ifndef ARITH_HARD_HPP_
#define ARITH_HARD_HPP_

/* This provides a hard-coded interface to the underlying arithmetic.
 * This include file depends on compilation flags!
 */

#include "arith-modp.hpp"
#include "arith-mod2.hpp"

#if defined(ARITH_MOD2)
typedef arith_mod2::gf2<ARITH_SIMD_GROUPSIZE> arith_hard;
#elif defined(ARITH_MODP)
typedef arith_modp::gfp<ARITH_PRIME_WIDTH> arith_hard;
#else
#error "change hard-coding scheme"
#endif

#endif	/* ARITH_HARD_HPP_ */
