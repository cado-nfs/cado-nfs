#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdlib> // for abort

#include <array>
#include <vector>

#include "macros.h"
#include "modredc126.hpp"
#include "u64arith.h"

#include "arithxx_common.hpp"

/* Only the .cpp source files that emit the non-inline symbols will
 * include this impl header file. So even though it does not look like
 * we're using it, in fact we are!  */
#include "arithxx_api_impl.hpp"      // IWYU pragma: keep
#include "arithxx_api128_impl.hpp"  // IWYU pragma: keep

// scan-headers: stop here

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
std::vector<arithxx_modredc126::Integer> arithxx_modredc126::Modulus::batchinv_redc(std::vector<uint64_t> const & a, Integer c) const
{
    auto const & me = *this;

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
    Residue R = Residue(me, Integer(a[0]));
    r.emplace_back(R.r);

    /* beta' = 2^64, beta = 2^128 */
    for (size_t i = 1 ; i < a.size() ; i++) {
        me.mul_ul(R, R, a[i]);
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
        me.mul_ul(R, R, a[i]);
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

arithxx_modredc126::Modulus::batch_Q_to_Fp_context::batch_Q_to_Fp_context(
        Integer const & num,
        Integer const & den)
    : remainder(num % den)
    , quotient((num - remainder).divexact(den))
    , D(den)
{
    ASSERT_ALWAYS(quotient.size_in_words() <= 1);
}


std::vector<uint64_t>
arithxx_modredc126::Modulus::batch_Q_to_Fp_context::operator()(
        std::vector<uint64_t> const & p, int k) const
{
    auto r = D.batchinv_redc(p, Integer(D.m) - remainder);
    if (r.empty())
        return {};

    std::vector<uint64_t> ri(p.size());

    /* It's a bit weird. den==D.m takes two words, but we're apparently
     * happy with one-word arithmetic.
     */
    for (size_t i = 0; i < p.size(); i++)
        ri[i] = u64arith_post_process_inverse(r[i][0], p[i],
                remainder[0], -D.invm, quotient[0], k);

    /* We used to have a "neg" flag as well. */
    return ri;
}

/* For each 0 <= i < n, compute r[i] = num/(den*2^k) mod p[i].
   den must be odd. If k > 0, then all p[i] must be odd.
   The memory pointed to be r and p must be non-overlapping.
   Returns 1 if successful. If any modular inverse does not exist,
   returns 0 and the contents of r are undefined. */
std::vector<uint64_t> arithxx_modredc126::Modulus::batch_Q_to_Fp(Integer const & num,
                            Integer const & den, int const k,
                            std::vector<uint64_t> const & p)
{
    return batch_Q_to_Fp_context(num, den)(p, k);
}

template struct arithxx_details::api<arithxx_modredc126>;
template struct arithxx_details::api128<arithxx_modredc126>;
