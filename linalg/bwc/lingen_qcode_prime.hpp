#ifndef LINGEN_QCODE_PRIME_HPP_
#define LINGEN_QCODE_PRIME_HPP_

#include <cstddef>
#include <vector>
#include <gmp.h>
#include "lingen_bmstatus.hpp"
#include "lingen_matpoly.hpp"

extern int
bw_lingen_basecase(bmstatus & bm, matpoly & pi, matpoly & E);
// extern int bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly const & E, std::vector<unsigned int> & delta);
extern void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);


#endif	/* LINGEN_QCODE_PRIME_HPP_ */
