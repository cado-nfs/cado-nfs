#ifndef CADO_LINGEN_EXPECTED_PI_LENGTH_HPP
#define CADO_LINGEN_EXPECTED_PI_LENGTH_HPP

#include <tuple>
#include <vector>

template<bool is_binary>
struct bw_dimensions;

extern std::tuple<unsigned int, unsigned int> get_minmax_delta(std::vector<unsigned int> const & delta);
// extern unsigned int get_min_delta(std::vector<unsigned int> const & delta);
// extern unsigned int get_max_delta(std::vector<unsigned int> const & delta);

template<bool is_binary>
extern unsigned int expected_pi_length(bw_dimensions<is_binary> & d, unsigned int len = 0);
template<bool is_binary>
extern unsigned int expected_pi_length(bw_dimensions<is_binary> & d, std::vector<unsigned int> const & delta, unsigned int len);
template<bool is_binary>
extern unsigned int expected_pi_length_lowerbound(bw_dimensions<is_binary> & d, unsigned int len);



#endif	/* LINGEN_EXPECTED_PI_LENGTH_HPP_ */
