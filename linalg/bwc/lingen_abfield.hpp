#ifndef LINGEN_ABFIELD_HPP_
#define LINGEN_ABFIELD_HPP_

// This header file defines the abfield typedef, based on which
// SELECT_MPFQ_LAYER macro is defined. Note that in the precise case of
// u64k1

#ifndef SELECT_MPFQ_LAYER_u64k1
#include "mpfq_layer.h" // IWYU pragma: export
#else
#include "mpfq_fake.hpp" // IWYU pragma: export
#endif

#endif	/* LINGEN_ABFIELD_HPP_ */
