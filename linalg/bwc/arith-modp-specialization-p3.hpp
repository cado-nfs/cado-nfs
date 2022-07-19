#ifndef ARITH_MODP_SPECIALIZATION_P3_HPP_
#define ARITH_MODP_SPECIALIZATION_P3_HPP_

#include "arith-modp-main.hpp"
#include "arith-modp-specializations.hpp"

namespace arith_modp {
namespace details {

/*  code for gfp<3, 1> */
#define ADDSUBMUL_CODE3(op, opc)                                               \
    FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)                                    \
    INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)                               \
    FINISH(op, opc, "%[z1]", "%[z2]", "%[z3]", r9, r10)                        \
      : "=&a"(foo)                                                             \
      , [z0] "+rm"(dst[0])                                                     \
      , [z1] "+rm"(dst[1])                                                     \
      , [z2] "+rm"(dst[2])                                                     \
      , [z3] "+rm"(dst[3])                                                     \
      : [s0] "0"(src[0])                                                       \
      , [s1] "rm"(src[1])                                                      \
      , [s2] "rm"(src[2])                                                      \
      , [mult] "rm"(x)                                                         \
      : "r8"                                                                   \
      , "r9"                                                                   \
      , "r10"                                                                  \
      , "rdx"

#define ADDSUB_CODE3(op, opc)                                                  \
    "" #op "q %q[s0], %q[d0]\n"                                                \
    "" #opc "q %q[s1], %q[d1]\n"                                               \
    "" #opc "q %q[s2], %q[d2]\n"                                               \
    "" #opc "q $0x0, %q[d3]\n"                                                 \
      : [d0] "+rm"(dst[0])                                                     \
      , [d1] "+rm"(dst[1])                                                     \
      , [d2] "+rm"(dst[2])                                                     \
      , [d3] "+rm"(dst[3])                                                     \
      : [s0] "r"(src[0])                                                       \
      , [s1] "r"(src[1])                                                       \
      , [s2] "r"(src[2])

EXPOSE_SPECIALIZATION(3);

}
}

#endif /* ARITH_MODP_SPECIALIZATION_P3_HPP_ */
