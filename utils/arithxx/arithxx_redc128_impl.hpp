#ifndef CADO_UTILS_ARITHXX_REDC128_IMPL_HPP
#define CADO_UTILS_ARITHXX_REDC128_IMPL_HPP

#include <cstdint>

#include <vector>

#include "arithxx_redc128.hpp"

#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY                                                            \
    "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall "               \
    "abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif


template<typename layer>
void
arithxx_details::redc128<layer>::mul_ul(Residue & r, Residue const & a, uint64_t const b) const
{
    auto const & me = downcast();

    me.assertValid(a);

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    uint64_t u = me.invm;    // NOLINT(misc-const-correctness)
    uint64_t a0 = a.r[0]; // NOLINT(misc-const-correctness)
    __asm__ __VOLATILE(
            /* Product of low words */
            "mulq %[b]\n\t"         /* rdx:rax = a0*b <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t" /* t0:rax = a0*b <= (2^64-1)^2 */
            /* Compute u such that a0*b + u*m == 0 (mod 2^64), compute u*m, add
               to t0:rax */
            "imulq %[u], %%rax\n\t"
            "movq %%rax, %[u]\n\t" /* u <= 2^64-1 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"       /* rdx:rax = u*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"       /* if low word != 0, carry to high word */
            "movq %[u], %%rax\n\t" /* rax = u, independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t" /* t1:t0 = (a0*b + u*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"  /* rdx:rax = u*m1 */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0*b+u*m)/2^64 <= 2^126 - 2^62
            */
            /* Free slot in pipe 2 here */
            ABORT_IF_CY

            /* Product of low and high word */
            "mulq %[b]\n\t" /* rdx:rax = a1*b <= (2^63-2^32-1)*(2^64-1) */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
                                       <= ((2^126-1)*(2^64-1) + (2^64-1)*(2^126-1))
                                       / 2^64 < 2^127 - 2^63 - 1, thus no carry */
            ABORT_IF_CY
            /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
               <= ((m-1)*(2^64-1) + (2^64-1)*m) / 2^64
               = 2*m*(1-1/2^64) - 1*(1-1/2^64). May need to subtract m. */
            "movq %[t0], %%rax\n\t" /* See if result > m */
            "movq %[t1], %%rdx\n\t"
            "subq %[m0], %[t0]\n\t" /* Try subtracting m, see if there's a carry
            */
            "sbbq %[m1], %[t1]\n\t"
            "cmovc %%rax, %[t0]\n\t" /* Carry -> restore old result */
            "cmovc %%rdx, %[t1]\n\t"
            :
            [u] "+&r"(u), [t0] "=&r"(r.r[0]), [t1] "=&r"(r.r[1]), [a0] "+&a"(a0)
            : [a1] "g"(a.r[1]), [b] "rm"(b), [m0] "rm"(me.m[0]), [m1] "rm"(me.m[1])
               : "%rdx", "cc");
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */
    uint64_t pl, ph, t[2];

    /* m < 1/4 W^2,  a < m, b < W */

    /* Product of b and low word */
    u64arith_mul_1_1_2(&(t[0]), &(t[1]), a.r[0],
            b); /* t1:t0 = a0*b < W^2 */

    /* One REDC step */
    redc1(t, t); /* t1:t0 = (a0*b + k*m) / W < m + W < 1/4 W^2 + W */

    /* Product of b and high word  */
    u64arith_mul_1_1_2(&pl, &ph, a.r[1], b);    /* ph:pl < 1/4 W^2 */
    u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 1/2 W^2 + W */

    /* Result may be larger than m, but is < 2*m */
    u64arith_sub_2_2_ge(&(t[0]), &(t[1]), me.m[0], me.m[1]);

    r.r[0] = t[0];
    r.r[1] = t[1];
#endif
    me.assertValid(r);
}
/* This computes the plain representatives (not Montgomery) of c/a[i],
 * where the a[i] are small.
 *
 * It's really pretty much the same code as
 * arithxx_details::api<layer>::batchinv, just with a few redc1 calls
 * inserted.
 *
 * The level of general usefulness of this interface is not entirely
 * clear to me.
 */
template<typename layer>
auto
arithxx_details::redc128<layer>::batchinv_redc(std::vector<uint64_t> const & a, Integer const & c) const
-> std::vector<Integer>
{
    auto const & me = downcast();

    if (a.empty())
        return {};

    std::vector<Integer> r;
    r.reserve(a.size());

    /* Note that the code in modredc64.cpp is somewhat different, and
     * computes r[0] by a multiplication by the representative of one --
     * which is equivalent to a reduction. By not doing it, we 
     * implicitly assume that a[0] is reduced. Which makes sense if our
     * base assumption is that our Modulus is larger than 64 bits.
     */
    Residue R(me);
    R.r = a[0];
    r.emplace_back(R.r);

    /* beta' = 2^64, beta = 2^128 */
    for (size_t i = 1 ; i < a.size() ; i++) {
        mul_ul(R, R, a[i]);
        r.emplace_back(R.r);
        /* r[i] = beta'^{-i} \prod_{0 <= j <= i} a[j]
         *      = beta' * \prod_{0 <= j <= i} (a[j] / beta')
         *      = beta / beta' * \prod_{0 <= j <= i} (a[j] / beta')
         */

        /* Notice that unlike the modredc64 case, the data that we have
         * here, which is still beta' * \prod_{0 <= j <= i} (a[j] /
         * beta'), is no longer a representative of the product but a
         * representative of 1/beta' * \prod(a[j]/beta').
         * 
         * This is reflected later on.
         */
    }

    /* Computes R = beta^2/r[n-1] */
    if (!me.inv(R, R))
        return {};
    /* R = beta * beta' * prod (beta' / a[j])
     *   = beta^2 / beta' * prod (beta' / a[j])
     *   = beta * [beta / beta' * prod (beta' / a[j])]
     *   = beta^2 * beta'^{n-1} \prod_{0 <= j < n} a[j]^{-1} */

    /* if c==1, then Residue(me, 1) is the representation of 2^-128
     * modulo n, so multiplying by it is equivalent to calling redc1
     * twice.
     *
     * c is beta * c/beta (the montgomery representative of c/beta)
     */
    me.mul(R, R, c);

    /* R = beta beta'^{n-1} c \prod_{0 <= j < n} a[j]^{-1} */
    /* R = beta * [c / beta * beta / beta' * prod (beta' / a[j])]]
     */

    /* A last redc1 nicely compensates things. */
    me.redc1(R, R);
    /* R = beta beta'^{n-2} c \prod_{0 <= j < n} a[j]^{-1}
     *   = beta * [ c/beta'^2 * prod (beta' / a[j])]
     */

    for (size_t i = a.size() - 1; i > 0; i--) {
        /* Invariant: R = beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1} */

        me.mul(r[i], R.r, r[i - 1]);
        /* r[i] := R * r[i-1] / beta
                = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) *
           (1/beta'^{i-1} \prod_{0 <= j <= i-1} a[j]) / beta = c a[i]^{-1} */
        /* Note that again, we compensate here the extra beta/beta'
         * factor that is in r[i-1]: we have
         *r[i-1]= beta * [ 1 / beta' * \prod_{0 <= j <= i-1} (a[j] / beta')]
         *    R = beta * [ c/beta'^2 * prod (beta' / a[j])]
         * r[i] = beta * [ c/beta'^2 * 1/beta' * beta' / a[i]]
         *      = c / a[i]
         */

        /* And this is here to adjust for the loop invariant.
         */
        mul_ul(R, R, a[i]);
        /* R := R * a[i] / beta'
             = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) * a[i] / beta'
             = beta beta'^{i-2} c \prod_{0 <= j < i} a[j]^{-1},
           thus satisfying the invariant for i := i - 1 */
    }
    /* Here have R = beta beta'^{-1} c / a[0]. We need to convert the factor
       beta to a factor of beta', so that the beta' cancel. */
    me.redc1(R, R); /* R := beta * beta'^{-1} c / a[0] / beta',
                    with beta = beta'^2, this is c / a[0] */
    r[0] = R.r;
    /* Note that at this point, the r[i]'s are plain representatives of
     * the result, not Montgomery representatives. It's slightly
     * annoying, isn't it?
     *
     * We could possibly "fix" this by moving the redc1 call to before
     * the inversion, thereby creating some extra multiplicative offset.
     */

    return r;
}

#endif	/* UTILS_ARITHXX_REDC128_IMPL_HPP_ */
