#ifndef CADO_LINGEN_QCODE_PRIME_HPP
#define CADO_LINGEN_QCODE_PRIME_HPP

#include <cstddef>      // size_t

#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"
struct bmstatus;

extern matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E);
// extern int bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly const & E, std::vector<unsigned int> & delta);
extern void test_basecase(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate);


#endif	/* LINGEN_QCODE_PRIME_HPP_ */
