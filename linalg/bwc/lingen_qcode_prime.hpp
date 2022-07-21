#ifndef LINGEN_QCODE_PRIME_HPP_
#define LINGEN_QCODE_PRIME_HPP_

#include <cstddef>      // size_t
#include <gmp.h>        // gmp_randstate_t
#include "lingen_matpoly_select.hpp"
struct bmstatus;

extern matpoly bw_lingen_basecase(bmstatus & bm, matpoly & E);
// extern int bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly const & E, std::vector<unsigned int> & delta);
extern void test_basecase(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);


#endif	/* LINGEN_QCODE_PRIME_HPP_ */
