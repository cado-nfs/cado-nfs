#ifndef LINGEN_QCODE_HPP_
#define LINGEN_QCODE_HPP_

#include "bw-common.h"
#include "lingen.hpp"
#include "lingen_matpoly.hpp"
#include <vector>
#include <tuple>

extern unsigned int expected_pi_length(bw_dimensions & d, unsigned int len = 0);
extern unsigned int expected_pi_length(bw_dimensions & d, std::vector<unsigned int> const & delta, unsigned int len);
extern unsigned int expected_pi_length_lowerbound(bw_dimensions & d, unsigned int len);
extern std::tuple<unsigned int, unsigned int> get_minmax_delta_on_solutions(bmstatus & bm, std::vector<unsigned int> const & delta);
extern unsigned int get_max_delta_on_solutions(bmstatus & bm, std::vector<unsigned int> const & delta);

extern int
bw_lingen_basecase(bmstatus & bm, matpoly & pi, matpoly & E, std::vector<unsigned int> & delta);
// extern int bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly const & E, std::vector<unsigned int> & delta);
extern void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);


#endif	/* LINGEN_QCODE_HPP_ */
