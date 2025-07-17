#ifndef CADO_LINGEN_QCODE_PRIME_HPP
#define CADO_LINGEN_QCODE_PRIME_HPP

#include <cstddef>

#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"

template<bool is_binary> struct bmstatus;

extern matpoly<false> bw_lingen_basecase(bmstatus<false> & bm, matpoly<false> & E);
extern void test_basecase(typename matpoly<false>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);

#endif	/* LINGEN_QCODE_PRIME_HPP_ */
