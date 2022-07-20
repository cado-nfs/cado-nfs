#ifndef ARITH_MODP_SPECIALIZATION_P2_HPP_
#define ARITH_MODP_SPECIALIZATION_P2_HPP_

#include "arith-modp-main.hpp"

namespace arith_modp {
namespace details {

/*  gfp<2,1> */
template<>
struct gfp<2> : public gfp_base<2, gfp<2>>
{
    typedef gfp_base<2, gfp<2>> super;
    using typename super::elt;
    using typename super::elt_ur_for_add;
    /* We need this so that the templates at the base level
     * participate in the resolution here.
     */
    using super::add;
    using super::sub;
    template<typename... Args>
    gfp(Args&&... args)
      : super(std::forward<Args>(args)...)
    {}
    static inline void add(elt_ur_for_add& dst, elt const& src)
    {
        asm("# gfp<2, 1>::add\n"
            "addq %q[s0], %q[d0]\n"
            "adcq %q[s1], %q[d1]\n"
            "adcq $0x0, %q[d2]\n"
            : [d0] "+rm"(dst[0]), [d1] "+rm"(dst[1]), [d2] "+rm"(dst[2])
            : [s0] "r"(src[0]), [s1] "r"(src[1]));
    }

    static inline void add(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        asm("# gfp<2, 1>::add\n"
            "addq %q[s0], %q[d0]\n"
            "adcq %q[s1], %q[d1]\n"
            "adcq %q[s2], %q[d2]\n"
            : [d0] "+rm"(dst[0]), [d1] "+rm"(dst[1]), [d2] "+rm"(dst[2])
            : [s0] "r"(src[0]), [s1] "r"(src[1]), [s2] "r"(src[2]));
    }

    static inline void sub(elt_ur_for_add& dst, elt const& src)
    {
        asm("# gfp<2, 1>::sub\n"
            "subq %q[s0], %q[d0]\n"
            "sbbq %q[s1], %q[d1]\n"
            "sbbq $0x0, %q[d2]\n"
            : [d0] "+rm"(dst[0]), [d1] "+rm"(dst[1]), [d2] "+rm"(dst[2])
            : [s0] "r"(src[0]), [s1] "r"(src[1]));
    }
    static inline void sub(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        asm("# gfp<2, 1>::sub\n"
            "subq %q[s0], %q[d0]\n"
            "sbbq %q[s1], %q[d1]\n"
            "sbbq %q[s2], %q[d2]\n"
            : [d0] "+rm"(dst[0]), [d1] "+rm"(dst[1]), [d2] "+rm"(dst[2])
            : [s0] "r"(src[0]), [s1] "r"(src[1]), [s2] "r"(src[2]));
    }
    inline void add(elt_ur_for_add& dst, elt const& a, elt const& b) const
    {
        set(dst, a);
        add(dst, b);
    }
    inline void sub(elt_ur_for_add& dst, elt const& a, elt const& b) const
    {
        set(dst, a);
        sub(dst, b);
    }

    static inline void addmul_ui(elt_ur_for_add& dst,
                                 elt const& src,
                                 mp_limb_t x)
    {
        mp_limb_t foo, bar;
        asm("# gfp<2, 1>::addmul_ui\n"
            "mulq   %[mult]\n"
            "addq   %%rax, %[z0]\n"
            "adcq   $0, %%rdx\n"
            "movq   %%rdx, %%rcx\n"
            "movq   %[s1], %%rax\n"
            "mulq   %[mult]\n"
            "addq   %%rcx, %%rax\n"
            "adcq   $0, %%rdx\n"
            "addq   %%rax, %[z1]\n"
            "adcq   $0, %%rdx\n"
            "addq   %%rdx, %[z2]\n"
            : "=&a"(foo),
              "=&d"(bar),
              [z0] "+rm"(dst[0]),
              [z1] "+rm"(dst[1]),
              [z2] "+rm"(dst[2])
            : [s0] "0"(src[0]), [s1] "rm"(src[1]), [mult] "rm"(x)
            : "rcx");
    }

    static inline void submul_ui(elt_ur_for_add& dst,
                                 elt const& src,
                                 mp_limb_t x)
    {
        mp_limb_t foo, bar;
        asm("# gfp<2, 1>::submul_ui\n"
            "mulq   %[mult]\n"
            "subq   %%rax, %[z0]\n"
            "adcq   $0, %%rdx\n"
            "movq   %%rdx, %%rcx\n"
            "movq   %[s1], %%rax\n"
            "mulq   %[mult]\n"
            "addq   %%rcx, %%rax\n"
            "adcq   $0, %%rdx\n"
            "subq   %%rax, %[z1]\n"
            "adcq   $0, %%rdx\n"
            "subq   %%rdx, %[z2]\n"
            : "=&a"(foo),
              "=&d"(bar),
              [z0] "+rm"(dst[0]),
              [z1] "+rm"(dst[1]),
              [z2] "+rm"(dst[2])
            : [s0] "0"(src[0]), [s1] "rm"(src[1]), [mult] "rm"(x)
            : "rcx");
    }
};

}
}

#endif /* ARITH_MODP_SPECIALIZATION_P2_HPP_ */
