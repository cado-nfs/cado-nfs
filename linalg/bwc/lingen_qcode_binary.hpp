#ifndef CADO_LINGEN_QCODE_BINARY_HPP
#define CADO_LINGEN_QCODE_BINARY_HPP

#include <cstddef>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"

template<bool is_binary> struct bmstatus;

extern matpoly<true> bw_lingen_basecase(bmstatus<true> & bm, matpoly<true> & E);

extern void test_basecase(matpoly<true>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);

extern void test_basecase_bblas(matpoly<true>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);

#endif	/* LINGEN_QCODE_BINARY_HPP_ */
