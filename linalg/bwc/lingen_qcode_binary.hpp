#ifndef CADO_LINGEN_QCODE_BINARY_HPP
#define CADO_LINGEN_QCODE_BINARY_HPP

#include <cstddef>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"
struct bmstatus;

extern matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E);

extern void test_basecase(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);
extern void test_basecase_bblas(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);


#endif	/* LINGEN_QCODE_BINARY_HPP_ */
