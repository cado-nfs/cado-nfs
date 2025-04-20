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

template<>
void
arithxx_details::api<arithxx_modredc126>::gcd(Integer & r, const Residue & A) const
{
    auto const & me = downcast();
    int sh;

    /* Since we do REDC arithmetic, we must have m odd */
    ASSERT_EXPENSIVE(me.m[0] % 2 != 0);

    if (me.is0(A)) {
        r = me.getmod();
        return;
    }

    std::array<uint64_t, 2> a = A.r, b = me.m;

    while (a[1] != 0 || a[0] != 0) {
        /* Make a odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros[(unsigned char)a[0]];
            u64arith_shrd(&(a[0]), a[1], a[0], sh);
            *(int64_t *)&(a[1]) >>= sh;
        } while (sh == 8);
#else
        if (a[0] == 0) /* ctzl does not like zero input */
        {
            a[0] = a[1];
            a[1] = ((int64_t)a[1] < 0L) ? (uint64_t)(-1L) : 0;
        }
        sh = u64arith_ctz(a[0]);
        u64arith_shrd(a.data(), a[1], a[0], sh);
        *(int64_t *)&(a[1]) >>= sh;
#endif

        /* Try to make the low two bits of b[0] zero */
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);
        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2(b.data(), b.data() + 1, a[0], a[1]);
        else
            u64arith_sub_2_2(b.data(), b.data() + 1, a[0], a[1]);

        if (b[0] == 0 && b[1] == 0) {
            if ((int64_t)a[1] < 0) {
                a[1] = -a[1];
                if (a[0] != 0)
                    a[1]--;
                a[0] = -a[0];
            }
            r = Integer(a);
            return;
        }

        /* Make b odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros[(unsigned char)b[0]];
            u64arith_shrd(&(b[0]), b[1], b[0], sh);
            *(int64_t *)&(b[1]) >>= sh;
        } while (sh == 8);
#else
        if (b[0] == 0) /* ctzl does not like zero input */
        {
            b[0] = b[1];
            b[1] = ((int64_t)b[1] < 0) ? (uint64_t)(-1) : 0;
        }
        sh = u64arith_ctz(b[0]);
        u64arith_shrd(b.data(), b[1], b[0], sh);
        *(int64_t *)&(b[1]) >>= sh;
#endif
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);

        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2(a.data(), a.data() + 1, b[0], b[1]);
        else
            u64arith_sub_2_2(a.data(), a.data() + 1, b[0], b[1]);
    }

    if ((int64_t)b[1] < 0) {
        b[1] = -b[1];
        if (b[0] != 0)
            b[1]--;
        b[0] = -b[0];
    }
    r = Integer(b);
}

template<>
int arithxx_details::api<arithxx_modredc126>::jacobi(Residue const & a_par) const
{
    auto const & me = downcast();
    Integer s;
    int t = 1;

    Integer a = me.get(a_par);
    Integer m = me.getmod();

    while (a != 0) {
        while ((a & 1) == 0) { /* TODO speedup */
            a >>= 1;
            if ((m & 7) == 3 || (m & 7) == 5)
                t = -t;
        }
        s = a; /* swap a and m */
        a = m;
        m = s;
        if ((a & 3) == 3 && (m & 3) == 3)
            t = -t;

        /* m is odd here */
        if (a >= m)
            a %= m;
    }
    if (m != 1)
        t = 0;

    return t;
}

template<>
bool arithxx_details::api<arithxx_modredc126>::inv(Residue & r, Residue const & A) const
{
    auto const & me = downcast();
    Integer a, u, v;
    int t;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue tmp(*this);

    set(tmp, A);
#endif

    me.assertValid(A);
    ASSERT_EXPENSIVE(m[0] % 2 != 0);

    if (me.is0(A))
        return false;

    Integer b = me.getmod();

    /* Let A = x*2^{2w}, so we want the Montgomery representation of 1/x,
       which is 2^{2w}/x. We start by getting a = x */

    /* We simply set a = x/2^{2w} and t=0. The result before correction
       will be 2^(2w+t)/x so we have to divide by t, which may be >64,
       so we may have to do one or more full and a variable width REDC. */
    /* TODO: If b[1] > 1, we could skip one of the two REDC */
    {
        Residue x(me);
        me.set(x, A);
        me.redc1(x, x);
        a = me.get(x);
    }
    /* Now a = x/2^w */
    t = -64;

    u = 1;
    v = 0; /* 0 is a valid pointer */

    /* make a odd */
    auto lsh = int(a.ctz());
    t += lsh;
    a >>= lsh;

    // Here a and b are odd, and a < b
    do {
        /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
        ASSERT_EXPENSIVE(a < b);
        ASSERT_EXPENSIVE((a & 1) == 1);
        ASSERT_EXPENSIVE((b & 1) == 1);
        ASSERT_EXPENSIVE((u & 1) == 1);
        ASSERT_EXPENSIVE((v & 1) == 0);

        do {
            b -= a;
            v += u;
            ASSERT_EXPENSIVE((b & 1) == 0);

            lsh = int(b.ctz());
            t += lsh;
            b >>= lsh;
            u <<= lsh;
        } while (a < b); /* ~50% branch taken :( */

        /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
        ASSERT_EXPENSIVE((a & 1) == 1);
        ASSERT_EXPENSIVE((b & 1) == 1);
        ASSERT_EXPENSIVE((u & 1) == 0);
        ASSERT_EXPENSIVE((v & 1) == 1);

        if (a == b)
            break;
        ASSERT_EXPENSIVE(a > b);

        /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
        do {
            a -= b;
            u += v;

            ASSERT_EXPENSIVE((a & 1) == 0);
            lsh = int(a.ctz());
            a >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (b < a); /* about 50% branch taken :( */
        /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
    } while (a != b);

    if (a != 1) /* Non-trivial GCD */
        return false;

    ASSERT_ALWAYS(t >= 0);

    /* Here, the inverse of a is u/2^t mod m. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t
       with impunity. */
    for (; t >= 64; t -= 64)
        me.redc1(u, u);

    if (t > 0) {
        uint64_t s[5], k;
        k = ((u.getWord(0) * me.invm) &
             ((UINT64_C(1) << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2(&s[0], &s[1], k, me.m[0]);
        /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
        u64arith_add_2_2(&s[0], &s[1], u.getWord(0), u.getWord(1));
        /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
        /* s[0] == 0 (mod 2^t) */
        ASSERT_EXPENSIVE((s[0] & ((1UL << t) - 1)) == 0);
        s[2] = 0;
        u64arith_mul_1_1_2(&(s[3]), &(s[4]), k, me.m[1]);
        u64arith_add_2_2(&(s[1]), &(s[2]), s[3], s[4]);

        /* Now shift s[2]:s[1]:s[0] right by t */
        u64arith_shrd(&(s[0]), s[1], s[0], t);
        u64arith_shrd(&(s[1]), s[2], s[1], t);

        u = Integer(std::array<uint64_t, 2> {s[0], s[1]});
        // t = 0;
    }

    r.r = u.get();
#ifdef WANT_ASSERT_EXPENSIVE
    mul(tmp, tmp, r);
    ASSERT_EXPENSIVE(is1(tmp));
#endif

    return true;
}

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
