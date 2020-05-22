#ifndef PLINGEN_HPP_
#define PLINGEN_HPP_

/* This contains some definitions for lingen mod p.
 *
 * This code is to be understood as scotch tape around an old and quirky
 * implementation.
 */

#include <cstdlib>
#include <vector>

#ifndef SELECT_MPFQ_LAYER_u64k1
#include "mpfq_layer.h"
#else
#include "mpfq_fake.hpp"
#endif

#include "select_mpi.h"

/* TODO: Rename ! */
struct bw_dimensions {
    unsigned int m,n,nrhs;
    abfield ab;
};


#endif	/* PLINGEN_HPP_ */
