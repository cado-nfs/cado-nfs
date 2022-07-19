#ifndef ARITH_MODP_SPECIALIZATION_P7_HPP_
#define ARITH_MODP_SPECIALIZATION_P7_HPP_

#include "arith-modp-main.hpp"
#include "arith-modp-specializations.hpp"

namespace arith_modp {
namespace details {

/*  code for gfp<7, 1> */
#define ADDSUBMUL_CODE7(op, opc)                                               \
    FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)                                  \
    INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)                      \
    INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)                     \
    INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)                    \
    INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)                     \
    INNER_MUL(op, "0x20(%[z])", "0x30(%[s])", r8, r9, r10)                     \
    FINISH(op, opc, "0x28(%[z])", "0x30(%[z])", "0x38(%[z])", r9, r10)         \
      :                                                                        \
      : [z] "D"(&dst[0])                                                       \
        , [s] "S"(&src[0])                                                     \
        , [mult] "rm"(x)                                                       \
      : "r8"                                                                   \
      , "r9"                                                                   \
      , "r10"                                                                  \
      , "r11"                                                                  \
      , "rax"                                                                  \
      , "rdx"                                                                  \
      , "memory"

#define ADDSUB_CODE7(op, opc)                                                  \
    "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q %q[s6], 0x30(%[z])\n"				\
        "" #opc "q $0x0, 0x38(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5]),                                  \
                    [s6]"r"(src[6])                                   \
                : "memory"

EXPOSE_SPECIALIZATION(7);

}
}

#endif /* ARITH_MODP_SPECIALIZATION_P7_HPP_ */
