#ifndef LAS_REDUCE_PLATTICE_REFERENCE2_ASM_HPP_
#define LAS_REDUCE_PLATTICE_REFERENCE2_ASM_HPP_

#include "las-plattice.hpp"
#include "reduce-plattice/plattice-proxy.hpp"

/* XXX This code is known to be buggy. It does not pass
 * test_reduce_plattice.
 * (It's not used beyond testing. And even then, this code should simply
 * go away)
 */
int
reference2_asm (plattice_proxy *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
{
    const int32_t hI = (int32_t) I;
    int32_t i0 = - (int32_t) p, i1 = (int32_t) r, j0, j1;

    /* Mac OS X 10.8 embarks a version of llvm which crashes on the code
     * below (could be that the constraints are exerting too much of the
     * compiler's behaviour).
     *
     * See tracker #16540
     */

#define RPA(LABEL)      \
    "addl %2, %0\n"						\
    "leal (%0,%2,4), %%edx\n"					\
    "addl %3, %1\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %0, %%eax\n"						\
    "jng " LABEL "\n"						\
    "addl %2, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovngl %%eax, %0\n"						\
    "cmovngl %%edx, %1\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"
#define RPB(LABEL)      \
    "addl %0, %2\n"						\
    "leal (%2,%0,4), %%edx\n"					\
    "addl %1, %3\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %2, %%eax\n"						\
    "jns " LABEL "\n"						\
    "addl %0, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovnsl %%eax, %2\n"						\
    "cmovnsl %%edx, %3\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"
#define RPC     \
    "cltd\n"							\
    "idivl %2\n"							\
    "imull %3, %%eax\n"						\
    "movl %%edx, %0\n"						\
    "subl %%eax, %1\n"
#define RPD "cltd\n"    \
    "idivl %0\n"							\
    "imull %1, %%eax\n"						\
    "movl %%edx, %2\n"						\
    "subl %%eax, %3\n"

    int32_t mhI;
    __asm__ __volatile__ (
            "xorl %1, %1\n"
            "cmpl %2, %5\n"
            "movl $0x1, %3\n"
            "jg 9f\n"
            "movl %5, %%eax\n"
            "negl %%eax\n"
            "movl %%eax, %4\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jl 0f\n"
            ".balign 8\n"
            "1:\n"
            "cmpl %0, %4\n"
            "jl 9f\n"
            RPB("3f")
            "2:\n"
            "cmpl %2, %5\n"
            "jg 9f\n"
            RPA("4f")
            "jmp 1b\n"
            ".balign 8\n"
            "3:\n"
            RPD
            "jmp 2b\n"
            ".balign 8\n"
            "4:\n"
            RPC
            "jmp 1b\n"
            ".balign 8\n"
            "0:\n"
            "movl %2, %%eax\n"
            "cltd\n"
            "idivl %0\n"
            "imull %1, %%eax\n"
            "subl %%eax, %3\n"
            "cmpl $0x19999999, %%edx\n"
            "movl %%edx, %2\n"
            "jle 2b\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "imull %3, %%eax\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jge 1b\n"
            "jmp 0b\n"
            "9:\n"
            : "+&r"(i0), "=&r"(j0), "+&r"(i1), "=&r"(j1),
        "=&rm"(mhI) : "rm"(hI) : "%rax", "%rdx", "cc");

#undef RPA
#undef RPB
#undef RPC
#undef RPD
    int64_t k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = -(int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}

#endif	/* LAS_REDUCE_PLATTICE_REFERENCE2_ASM_HPP_ */
