#ifndef CADO_XVECTORS_HPP
#define CADO_XVECTORS_HPP

#include <cstdint>
#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "parallelizing_info.hpp"

std::vector<uint32_t> setup_x_random(unsigned int m, unsigned int nx,
                                           unsigned int nr,
                                           parallelizing_info_ptr pi,
                                           cxx_gmp_randstate & rstate,
                                           std::vector<unsigned int> const & forced = {});
void save_x(std::vector<uint32_t> const & xs, unsigned int m, unsigned int nx,
            parallelizing_info_ptr pi);

std::vector<uint32_t> load_x(unsigned int m, unsigned int & nx,
                                   parallelizing_info_ptr pi);
std::vector<uint32_t> set_x_fake(unsigned int m, unsigned int & nx,
                                       parallelizing_info_ptr pi);

#endif /* XVECTORS_HPP_ */
