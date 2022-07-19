#ifndef ARITH_MODP_SPECIALIZATION_P4_HPP_
#define ARITH_MODP_SPECIALIZATION_P4_HPP_

#include "arith-modp-main.hpp"
#include "arith-modp-specializations.hpp"

namespace arith_modp {
namespace details {

/*  code for gfp<4, 1> */
/*
#define ADDSUBMUL_CODE4(op, opc)					\
            FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)			\
            INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)		\
            INNER_MUL(op, "%[z1]", "%[s3]", r9, r10, r11)		\
            FINISH(op, opc, "%[z2]", "%[z3]", "%[z4]", r10, r11)	\
            : "=&a"(foo),                                           \
                [z0]"+rm"(dst[0]),				\
                [z1]"+rm"(dst[1]),				\
                [z2]"+rm"(dst[2]),				\
                [z3]"+rm"(dst[3]),				\
                [z4]"+rm"(dst[4])					\
            :							\
                [s0]"0"(src[0]),					\
                [s1]"rm"(src[1]),					\
                [s2]"rm"(src[2]),					\
                [s3]"rm"(src[3]),					\
                [mult]"rm"(x)					\
            : "r8", "r9", "r10", "r11", "rdx"

#define xADDSUB_CODE4(op, opc)   \
"" #op  "q %q[s0], (%[z])\n"					\
"" #opc "q %q[s1], 0x8(%[z])\n"					\
"" #opc "q %q[s2], 0x10(%[z])\n"				\
"" #opc "q %q[s3], 0x18(%[z])\n"				\
"" #opc "q $0x0, 0x20(%[z])\n"					\
    :							\
    :							\
        [z]"r"(&dst[0]),				        \
        [s0]"r"(src[0]),					\
        [s1]"r"(src[1]),					\
        [s2]"r"(src[2]),					\
        [s3]"r"(src[3])                                   \
    : "memory"


            */
#define ADDSUB_CODE4(op, opc)                                                  \
    "" #op "q %q[s0], %q[d0]\n"                                                \
    "" #opc "q %q[s1], %q[d1]\n"                                               \
    "" #opc "q %q[s2], %q[d2]\n"                                               \
    "" #opc "q %q[s3], %q[d3]\n"                                               \
    "" #opc "q $0x0, %q[d4]\n"                                                 \
      : [d0] "+rm"(dst[0])                                                     \
      , [d1] "+rm"(dst[1])                                                     \
      , [d2] "+rm"(dst[2])                                                     \
      , [d3] "+rm"(dst[3])                                                     \
      , [d4] "+rm"(dst[4])                                                     \
      : [s0] "r"(src[0])                                                       \
      , [s1] "r"(src[1])                                                       \
      , [s2] "r"(src[2])                                                       \
      , [s3] "r"(src[3])

#define ADDSUBMUL_CODE4(op, opc)                                               \
    FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)                                  \
    INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)                      \
    INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)                     \
    FINISH(op, opc, "0x10(%[z])", "0x18(%[z])", "0x20(%[z])", r10, r11)        \
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

EXPOSE_SPECIALIZATION(4);

}
}

#endif /* ARITH_MODP_SPECIALIZATION_P4_HPP_ */
