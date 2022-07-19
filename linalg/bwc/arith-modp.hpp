#ifndef ARITH_MODP_HPP_
#define ARITH_MODP_HPP_

/* arith-modp-main.hpp contains the base implementation. Faster code goes
 * in the specializations for critical functions.
 */
#include "arith-modp-main.hpp"

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
/* Now some specializations */

#include "arith-modp-specialization-p1.hpp"
#include "arith-modp-specialization-p2.hpp"
#include "arith-modp-specialization-p3.hpp"
#include "arith-modp-specialization-p4.hpp"
#include "arith-modp-specialization-p5.hpp"
#include "arith-modp-specialization-p6.hpp"
#include "arith-modp-specialization-p7.hpp"
#include "arith-modp-specialization-p8.hpp"

/* further specialization only seem to bring very marginal
 * improvements. This should probably go away. */

// this one is disabled.
// #include "arith-modp-carry-save.hpp"

#endif

namespace arith_modp {
/* expose only what we have in our public interface */
using details::gfp;
// using details::fast_type;
}

#endif
