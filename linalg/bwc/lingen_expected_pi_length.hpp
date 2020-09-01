#ifndef LINGEN_EXPECTED_PI_LENGTH_HPP_
#define LINGEN_EXPECTED_PI_LENGTH_HPP_

#include <tuple>
#include <vector>
struct bw_dimensions;

extern std::tuple<unsigned int, unsigned int> get_minmax_delta(std::vector<unsigned int> const & delta);
// extern unsigned int get_min_delta(std::vector<unsigned int> const & delta);
// extern unsigned int get_max_delta(std::vector<unsigned int> const & delta);

extern unsigned int expected_pi_length(bw_dimensions & d, unsigned int len = 0);
extern unsigned int expected_pi_length(bw_dimensions & d, std::vector<unsigned int> const & delta, unsigned int len);
extern unsigned int expected_pi_length_lowerbound(bw_dimensions & d, unsigned int len);



#endif	/* LINGEN_EXPECTED_PI_LENGTH_HPP_ */
