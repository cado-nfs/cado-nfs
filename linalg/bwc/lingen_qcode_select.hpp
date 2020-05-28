#ifndef LINGEN_QCODE_SELECT_HPP_
#define LINGEN_QCODE_SELECT_HPP_

#ifdef SELECT_MPFQ_LAYER_u64k1
#include "lingen_qcode_binary.hpp" // IWYU pragma: export
#else
#include "lingen_qcode_prime.hpp" // IWYU pragma: export
#endif

#endif	/* LINGEN_QCODE_SELECT_HPP_ */
