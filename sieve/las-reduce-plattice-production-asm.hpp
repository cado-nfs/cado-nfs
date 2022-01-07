#ifndef LAS_REDUCE_PLATTICE_PRODUCTION_ASM_HPP_
#define LAS_REDUCE_PLATTICE_PRODUCTION_ASM_HPP_

/* This is example code that is meant to go __into__ a class which has
 * several data members and function members defined. For maintenance,
 * all these pieces of code are put in separate files, but these must be
 * #included from within the class itself.
 */
#ifndef LAS_PLATTICE_H
#error "This prologue is only here to please clangd !!!"
#include <cstdint>
#include "macros.h"
struct mock_plattice_info {
        uint32_t mi0;  /* This encodes i0 = -mi0 */
        uint32_t j0;
        uint32_t i1;
        uint32_t j1;
        void reduce_with_vertical_vector(uint32_t I);
        bool needs_special_treatment(uint32_t) const;
#endif

void reduce_plattice_asm(uint32_t I) {
    if (needs_special_treatment(I)) {
        reduce_with_vertical_vector(I);
        return;
    }

    /* This is the main reduce_plattice loop */

    /* This version of the code tries to maintain the vectors
     *      u0 = (-i0,-j0)  and   u1 = (i1,j1)
     * inside 64-bit registers as much as possible.
     *
     * We know that -i0 == mi0 is always nonnegative, and -j0 is
     * definitely negative, so it is stored in two complement form.
     *
     * Occasionally, for control, we resort to a temporary register
     * to save either mi0 or i1. But most of the work looks at the
     * lower part of the registers that hold u0 and u1.
     *
     * We use a register marker to indicate the branch from which
     * we're exiting.
     *
     * We explored the possibility of collapsing the two loop
     * branches into one. It doesn't seem to bring any benefit.
     *
     * Another variant that was tried is, instead of comparing the
     * quotient to a bound to decide whether a division is relevant, and
     * otherwise do a simple subtractive algorithm, prefer a scheme where
     * the subtraction is done right away (and canceled by a cmov if
     * actually unneeded), and a division follows if we can see that this
     * wasn't enough. In fact, this doesn't bring a measurable benefit.
     *
     */

    int branch;
    uint64_t u0 = (((uint64_t)-j0)<<32) | (uint64_t) mi0;
    uint64_t u1 = (((uint64_t) j1)<<32) | (uint64_t) i1;

    // see 6.47.2.8 x86 Operand Modifiers regarding %k[] and such

#define SUBTRACTIVE_BLOCK_BRANCH_1_NO_HEAD			\
    "subq %[u1], %[u0]\n"

#define SUBTRACTIVE_BLOCK_BRANCH_1_NO_TAIL			\
    "movq %[u0], %%rdx\n"					\
    SUBTRACTIVE_BLOCK_BRANCH_1_NO_HEAD                          \
    "cmpl %k[u1], %%edx\n"

#define SUBTRACTIVE_BLOCK_BRANCH_1				\
    SUBTRACTIVE_BLOCK_BRANCH_1_NO_TAIL			        \
    "cmovb %%rdx, %[u0]\n"

#define SUBTRACTIVE_BLOCK_BRANCH_2_NO_HEAD			\
    "subq %[u0], %[u1]\n"

#define SUBTRACTIVE_BLOCK_BRANCH_2_NO_TAIL			\
    "movq %[u1], %%rdx\n"					\
    SUBTRACTIVE_BLOCK_BRANCH_2_NO_HEAD                          \
    "cmpl %k[u0], %%edx\n"

#define SUBTRACTIVE_BLOCK_BRANCH_2				\
    SUBTRACTIVE_BLOCK_BRANCH_2_NO_TAIL			        \
    "cmovb %%rdx, %[u1]\n"

#define xxxPREEMPTIVELY_SUBTRACT

    asm volatile(
            ".p2align 4,0x90\n"
            "# first branch\n"
            "0:\n"
            "xorl %[branch], %[branch]\n"
            /* if I > i1, then exit */
            "cmpl %k[u1], %[I]\n"
            "ja 9f\n"       /* exit */

#ifndef PREEMPTIVELY_SUBTRACT
            "leal (%k[u1], %k[u1], 2), %%edx\n"
            "shl $0x1, %%edx\n"
            "cmpl %k[u0], %%edx\n"
            "jbe 2f\n"      /* needs division */

            /* The optimal number of subtractive blocks may be a
             * matter of data-dependent and
             * microarchitecture-dependent optimization. For
             * reference, the distribution of quotients is 45% / 17%
             * / 9% / 6% / less.
             */
            SUBTRACTIVE_BLOCK_BRANCH_1_NO_HEAD
#else
            SUBTRACTIVE_BLOCK_BRANCH_1
#endif
            SUBTRACTIVE_BLOCK_BRANCH_1
            SUBTRACTIVE_BLOCK_BRANCH_1
            SUBTRACTIVE_BLOCK_BRANCH_1
            SUBTRACTIVE_BLOCK_BRANCH_1
#ifndef PREEMPTIVELY_SUBTRACT
            SUBTRACTIVE_BLOCK_BRANCH_1
#else
            SUBTRACTIVE_BLOCK_BRANCH_1_NO_TAIL
            // "cmovb %%rdx, %[u0]\n"

            // if b, then we're doing something trivial, meaning that
            // reduction is over at this point. Let's move on.
            // Otherwise (ae) we need to do a division now.
            "jae 2f\n"
            "movq %%rdx, %[u0]\n"
#endif

            "# second branch\n"
            "1:\n"
            "movl %[I], %[branch]\n"
            "cmpl %k[u0], %[I]\n"
            "ja 9f\n"       /* exit */

#ifndef PREEMPTIVELY_SUBTRACT
            "leal (%k[u0], %k[u0], 2), %%edx\n"
            "shl $0x1, %%edx\n"
            "cmpl %k[u1], %%edx\n"
            "jbe 3f\n"      /* needs division */

            SUBTRACTIVE_BLOCK_BRANCH_2_NO_HEAD
#else
            SUBTRACTIVE_BLOCK_BRANCH_2
#endif
            SUBTRACTIVE_BLOCK_BRANCH_2
            SUBTRACTIVE_BLOCK_BRANCH_2
            SUBTRACTIVE_BLOCK_BRANCH_2
            SUBTRACTIVE_BLOCK_BRANCH_2
#ifndef PREEMPTIVELY_SUBTRACT
            SUBTRACTIVE_BLOCK_BRANCH_2
#else
            SUBTRACTIVE_BLOCK_BRANCH_2_NO_TAIL
            // "cmovb %%rdx, %[u1]\n"
            "jae 3f\n"
            "movq %%rdx, %[u1]\n"
#endif
            "jmp 0b\n"

            "3:\n"
            "xorl %%edx, %%edx\n"
            "movl %k[u1], %%eax\n"
            "divl %k[u0]\n"
            /* edx is the remainder, eax is the quotient */
            /* we don't want u0 to change, but we do want to adjust
             * j1, and that's actually slightly tricky to do, since
             * we insist on doing this with an 32-bit multiply. */

            /* In these branches, we play shift / shift-with-feed
             * tricks to manipulate the high parts of the registers.
             * Note that the low part of u1 (= i1) is useless at this
             * point, since it will be replaced by the remainder edx
             */

            /* prepare old j1. We simply put it in the lower part of
             * u1 temporarily */
            "shr $0x20, %[u1]\n"

            "shl $0x20, %%rdx\n"    /* this will be i1. Move to high part */

            /* This copies the high part of u0 (= -j0) into
             * k[branch], which we can overwrite temporarily */
            "shld $0x20, %[u0], %q[branch]\n"

            /* store in k[u1] the result of j1 - (q * (-j0)) */
            "imull %k[branch], %%eax\n"
            "subl %%eax, %k[u1]\n" /* this is new j1 */

            /* we're done with j1, put it back in place, and feed in
             * with i1 which is still in the high part of rdx */
            "shld $0x20, %%rdx, %[u1]\n"   /* final u1 */
            "jmp 0b\n"

            "2:\n"
            "xorl %%edx, %%edx\n"
            "movl %k[u0], %%eax\n"
            "divl %k[u1]\n"
            /* edx is the remainder, eax is the quotient */
            "shr $0x20, %[u0]\n"
            "shl $0x20, %%rdx\n"    /* this will be mi0 */
            "shld $0x20, %[u1], %q[branch]\n"
            "imull %k[branch], %%eax\n"
            "subl %%eax, %k[u0]\n"
            "shld $0x20, %%rdx, %[u0]\n"
            "jmp 1b\n"

            "9:\n"
            :   [branch]"=&r"(branch),
        [u0]"+r"(u0),
        [u1]"+r"(u1)
            :   [I]"r"(I)
                 : "eax", "edx");
    mi0 = u0;
    i1 = u1;
    j0 = -(int32_t)(u0>>32);
    j1 = (int32_t)(u1>>32);
    if (branch) {
        if (mi0 == 0) {
            mi0 = i1;
            i1 = j0 ; j0 = j1 ; j1 = i1;
            i1 = 0;
            reduce_with_vertical_vector(I);
            return;
        }
        ASSERT(mi0 + i1 >= I);
        int a = (mi0 + i1 - I) / mi0;
        i1 -= a * mi0;
        j1 += a * j0;
    } else {
        if (i1 == 0) {
            // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
            j0 = j1 - j0;
            reduce_with_vertical_vector(I);
            return;
        }
        ASSERT(mi0 + i1 >= I);
        int a = (mi0 + i1 - I) / i1;
        mi0 -= a * i1;
        j0  += a * j1;
    }
}

#ifndef LAS_PLATTICE_H
#error "This epilogue is only here to please clangd !!!"
};
#endif

#endif	/* LAS_REDUCE_PLATTICE_PRODUCTION_ASM_HPP_ */
