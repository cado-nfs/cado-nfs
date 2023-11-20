#ifndef ARITH_MODP_SPECIALIZATIONS_HPP_
#define ARITH_MODP_SPECIALIZATIONS_HPP_

/* This header only defines a few macros that are used in various
 * arith-modp-specialization-pXXX.hpp headers
 */

/*  macros for assembly for further specializations */
#define FEED_IN_WITH_S0_IN_RAX(in1, r0, r1)                                    \
    /* status: s0 in rax */                                                    \
    "mulq   %[mult]\n"          /* rdx:rax = s0 * v */                         \
    "movq   %%rax, %%" #r0 "\n" /* lo contrib to d1 */                         \
    "movq   " in1 ", %%rax\n"   /* load s1          */                         \
    "movq   %%rdx, %%" #r1 "\n" /* hi contrib to d1 */
#define FEED_IN(in0, in1, r0, r1)                                              \
    "movq   " in0 ", %%rax\n" FEED_IN_WITH_S0_IN_RAX(in1, r0, r1)
#define INNER_MUL(op, out, in, r0, r1, r2)                                     \
    /* status: r1:r0 to be added to d_{i+1}:d_i, rax = s_{i+1} */              \
    "xorq   %%" #r2 ", %%" #r2 "\n"                                            \
    "mulq   %[mult]\n"              /* rdx:rax = s_{i+1} * v */                \
    "" #op "q %%" #r0 ", " out "\n" /* store d_i             */                \
    "adcq   %%rax, %%" #r1 "\n"     /* lo contrib to d_{i+1} */                \
    "adcq   %%rdx, %%" #r2 "\n"     /* hi contrib to d_{i+2} */                \
    "movq   " in ", %%rax\n"        /* load s_{i+2}          */
#define FINISH(op, opc, out0, out1, out2, r0, r1)                              \
    /* r1:r0 to be added to d_{i+1}:d_i ; rax = s_{i+2} */                     \
    "mulq   %[mult]\n"                                                         \
    "" #op "q   %%" #r0 ", " out0 "\n"                                         \
    "adcq   %%rax, %%" #r1 "\n"                                                \
    "adcq   $0x0, %%rdx\n"                                                     \
    "" #op "q   %%" #r1 ", " out1 "\n"                                         \
    "" #opc "q   %%rdx, " out2 "\n"

/*  this macro actually exposes the specialization in itself */
#define EXPOSE_SPECIALIZATION(n)                                               \
    template<>                                                                 \
    struct gfp<n> : public gfp_base<n, gfp<n>>                                 \
    {                                                                          \
        typedef gfp_base<n, gfp<n>> super;                                     \
        using typename super::elt;                                             \
        using typename super::elt_ur_for_add;                                  \
        template<typename... Args>                                             \
        gfp(Args&&... args)                                                    \
          : super(std::forward<Args>(args)...)                                 \
        {}                                                                     \
        using super::add_ur;                                                   \
        using super::sub_ur;                                                   \
        /* We need this so that the templates at the base level                \
         * participate in the resolution here.                                 \
         */                                                                    \
        using super::add;                                                      \
        using super::sub;                                                      \
        static inline void add(elt_ur_for_add& dst, elt const& src)            \
        {                                                                      \
            asm("# gfp<" #n ", 1>::add\n" ADDSUB_CODE##n(add, adc));           \
        }                                                                      \
        static inline void sub(elt_ur_for_add& dst, elt const& src)            \
        {                                                                      \
            asm("# gfp<" #n ", 1>::sub\n" ADDSUB_CODE##n(sub, sbb));           \
        }                                                                      \
        inline void add(elt_ur_for_add& dst,                            \
                               elt const& a,                                   \
                               elt const& b) const                                   \
        {                                                                      \
            set(dst, a);                                                 \
            add(dst, b);                                                       \
        }                                                                      \
        inline void sub(elt_ur_for_add& dst,                            \
                               elt const& a,                                   \
                               elt const& b) const                                   \
        {                                                                      \
            set(dst, a);                                                 \
            sub(dst, b);                                                       \
        }                                                                      \
        inline void add(elt_ur_for_add& dst, elt_ur_for_add const& src) const \
        {                                                                      \
            add_ur(dst, src);                                                  \
        }                                                                      \
        inline void sub(elt_ur_for_add& dst, elt_ur_for_add const& src) const \
        {                                                                      \
            sub_ur(dst, src);                                                  \
        }                                                                      \
        static inline void addmul_ui(elt_ur_for_add& dst,                      \
                                     elt const& src,                           \
                                     mp_limb_t x)                              \
        {                                                                      \
            mp_limb_t foo MAYBE_UNUSED;                                        \
            asm("# gfp<" #n ", 1>::addmul_ui\n" ADDSUBMUL_CODE##n(add, adc));  \
        }                                                                      \
        static inline void submul_ui(elt_ur_for_add& dst,                      \
                                     elt const& src,                           \
                                     mp_limb_t x)                              \
        {                                                                      \
            mp_limb_t foo MAYBE_UNUSED;                                        \
            asm("# gfp<" #n ", 1>::submul_ui\n" ADDSUBMUL_CODE##n(sub, sbb));  \
        }                                                                      \
    }

#endif /* ARITH_MODP_SPECIALIZATIONS_HPP_ */
