#ifndef CADO_ARITH_CONCRETE_BASE_HPP
#define CADO_ARITH_CONCRETE_BASE_HPP

#define ALIGNMENT_ON_ALL_BWC_VECTORS    64

/* See also MINIMUM_ITEMS_IN_BWC_CHUNKS in balancing.hpp */

/* see also
   tests/linalg/bwc/bwc-ptrace.sh
   tests/linalg/bwc/convert_magma.pl
   */

/* all concrete instantiations derive from this */

struct arith_concrete_base {
    struct elt {};
};

#endif	/* ARITH_CONCRETE_BASE_HPP_ */
