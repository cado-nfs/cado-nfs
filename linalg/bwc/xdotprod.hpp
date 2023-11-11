#ifndef XDOTPROD_HPP_
#define XDOTPROD_HPP_

#include <stdint.h>              // for uint32_t
#include "matmul_top.hpp"

/* This interface is relevant to both krylov and mksol, since it's used
 * both for checking and computing A files (obviously only krylov is
 * concerned by the latter aspect). */

void x_dotprod(arith_generic::elt * dst, uint32_t * xv, unsigned int m, unsigned int nx, mmt_vec const & v, int sign);

#endif	/* XDOTPROD_HPP_ */
