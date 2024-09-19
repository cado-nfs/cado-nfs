#ifndef XVECTORS_HPP_
#define XVECTORS_HPP_

#include <gmp.h>
#include <cstdint>
#include <memory>
#include "parallelizing_info.hpp"

void setup_x_random(uint32_t * xs,
        unsigned int m, unsigned int nx, unsigned int nr,
        parallelizing_info_ptr pi, gmp_randstate_t rstate);
void save_x(const uint32_t * xs, unsigned int m, unsigned int nx,
        parallelizing_info_ptr pi);

std::unique_ptr<uint32_t[]>
load_x(unsigned int m, unsigned int & nx,
        parallelizing_info_ptr pi);
std::unique_ptr<uint32_t[]>
set_x_fake(unsigned int m, unsigned int & nx,
        parallelizing_info_ptr pi);

#endif	/* XVECTORS_HPP_ */
