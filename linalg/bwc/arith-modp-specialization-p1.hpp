#ifndef ARITH_MODP_SPECIALIZATION_P1_HPP_
#define ARITH_MODP_SPECIALIZATION_P1_HPP_

#include "arith-modp-main.hpp"

namespace arith_modp {
namespace details {

/* gfp<1,1> */
template<>
struct gfp<1> : public gfp_base<1, gfp<1>>
{
    typedef gfp_base<1, gfp<1>> super;
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

    static inline void add_ur(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        asm("# gfp<1, 1>::add\n"
            "addq %q2, %q0\n"
            "adcq %q3, %q1\n"
            : "+r"(dst[0]), "+r"(dst[1])
            : "rm"(src[0]), "rm"(src[1]));
    }

    static inline void add(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        add_ur(dst, src);
    }

    static inline void add(elt_ur_for_add& dst, elt const& src)
    {
        asm("# gfp<1, 1>::add\n"
            "addq %q2, %q0\n"
            "adcq $0x0, %q1\n"
            : "+r"(dst[0]), "+r"(dst[1])
            : "rm"(src[0]));
    }
    static inline void add(elt_ur_for_add& dst, elt const& a, elt const& b)
    {
        asm("# gfp<1, 1>::add\n"
            "movq %q2, %q0\n"
            "xorq %q1, %q1\n"
            "addq %q3, %q0\n"
            "adcq $0x0, %q1\n"
            : "=&r"(dst[0]), "=&r"(dst[1])
            : "rm"(a[0]), "rm"(b[0]));
    }

    static inline void sub(elt_ur_for_add& dst, elt const& src)
    {
        asm("# gfp<1, 1>::sub\n"
            "subq %q2, %q0\n"
            "sbbq $0x0, %q1\n"
            : "+r"(dst[0]), "+r"(dst[1])
            : "rm"(src[0]));
    }

    static inline void sub_ur(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        asm("# gfp<1, 1>::sub\n"
            "subq %q2, %q0\n"
            "sbbq %q3, %q1\n"
            : "+r"(dst[0]), "+r"(dst[1])
            : "rm"(src[0]), "rm"(src[1]));
    }
    static inline void sub(elt_ur_for_add& dst, elt_ur_for_add const& src)
    {
        sub_ur(dst, src);
    }

    static inline void sub(elt_ur_for_add& dst, elt const& a, elt const& b)
    {
        asm("# gfp<1, 1>::sub\n"
            "movq %q2, %q0\n"
            "xorq %q1, %q1\n"
            "subq %q3, %q0\n"
            "sbbq $0x0, %q1\n"
            : "=&r"(dst[0]), "=&r"(dst[1])
            : "rm"(a[0]), "rm"(b[0]));
    }

    static inline void addmul_ui(elt_ur_for_add& dst,
                                 elt const& src,
                                 mp_limb_t x)
    {
        mp_limb_t foo, bar;
        asm("# gfp<1, 1>::addmul_ui\n"
            "mulq   %[mult]\n"
            "addq   %%rax, %[z0]\n"
            "adcq   $0, %%rdx\n"
            "addq   %%rdx, %[z1]\n"
            : "=a"(foo), "=&d"(bar), [z0] "+rm"(dst[0]), [z1] "+rm"(dst[1])
            : "0"(src[0]), [mult] "r1m"(x));
    }
    static inline void submul_ui(elt_ur_for_add& dst,
                                 elt const& src,
                                 mp_limb_t x)
    {
        mp_limb_t foo, bar;
        asm("# gfp<1, 1>::submul_ui\n"
            "mulq   %[mult]\n"
            "subq   %%rax, %[z0]\n"
            "adcq   $0, %%rdx\n"
            "subq   %%rdx, %[z1]\n"
            : "=a"(foo), "=&d"(bar), [z0] "+rm"(dst[0]), [z1] "+rm"(dst[1])
            : "0"(src[0]), [mult] "r1m"(x));
    }
};

}
}

#endif /* ARITH_MODP_SPECIALIZATION_P1_HPP_ */
