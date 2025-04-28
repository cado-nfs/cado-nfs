#ifndef CADO_XDOTPROD_HPP
#define CADO_XDOTPROD_HPP

#include <cstdint>              // for uint32_t

#include <vector>

#include "matmul_top.hpp"

/* This interface is relevant to both krylov and mksol, since it's used
 * both for checking and computing A files (obviously only krylov is
 * concerned by the latter aspect). */

void x_dotprod(arith_generic::elt * dst, std::vector<uint32_t> const & xv, unsigned int j0, unsigned int j1, unsigned int nx, mmt_vec const & v, int sign);

#endif	/* XDOTPROD_HPP_ */
