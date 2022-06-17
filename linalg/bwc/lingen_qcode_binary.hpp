#ifndef LINGEN_QCODE_BINARY_HPP_
#define LINGEN_QCODE_BINARY_HPP_

#include <cstddef>
#include <gmp.h>
#include "lingen_matpoly_select.hpp"
struct bmstatus;

extern matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E);

extern void test_basecase(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);
extern void test_basecase_bblas(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);


#endif	/* LINGEN_QCODE_BINARY_HPP_ */
