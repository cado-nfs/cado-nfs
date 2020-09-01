#ifndef LINGEN_AVERAGE_MATSIZE_HPP_
#define LINGEN_AVERAGE_MATSIZE_HPP_

#ifndef SELECT_MPFQ_LAYER_u64k1
#include "mpfq_layer.h"
#else
#include "mpfq_fake.hpp"
#endif

double average_matsize(abdst_field ab, unsigned int m, unsigned int n, int ascii);

#endif	/* LINGEN_AVERAGE_MATSIZE_HPP_ */
