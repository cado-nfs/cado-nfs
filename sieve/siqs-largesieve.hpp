#ifndef CADO_SIQS_LARGESIEVE_HPP
#define CADO_SIQS_LARGESIEVE_HPP

#include <cstdint>

#include <algorithm>
#include <limits>
#include <tuple>
#include <vector>

#include "bucket.hpp"
#include "las-arith.hpp"
#include "las-qlattice.hpp"
#include "macros.h"

/* Large sieve for SIQS.
 * For prime (or prime power) p, such that p > I, we want to compute all (i, j)
 * pairs such that f(rj+q*i) is divisible by p with i in [-I/2, I/2[ and j in
 * [0,2^#factors(sq)[.
 * The implementatation follows Kleinjung's paper "Quadratic Sieving" with
 * slight deviations.
 *
 * Assumes gcd(p, q) == 1 [ this case could be handle (see Remark 2.3 of
 * Kleinjung's paper) but in this case the test Q.is_coprime_to(e.p) should be
 * remove for fill_in_buckets_* functions ]
 *
 * Notations:
 *  q=prod(qk, 1 <= k <= n) the special-q
 *  R = {Rk} such that Rk root modulo qk, Rk = 0 mod q/qk and 0 <= Rk < q/2
 *  rj = add((-1)^(k-th bit of g)*Rk for 1 <= k <= n) where g is the gray code
 *  corresponding to j and 0 <= j < 2^n
 *  rp a fixed root mod p
 *
 *  (i, j) hit <=> rj+q*i = rp modulo p, with i in [-I/2, I/2[, j in [0,2^n[
 *             <=> i = (rp - rj)/q modulo p, with same condition
 *
 * Remarks:
 *  - we use [-I/2, I/2[ as sieving interval instead of
 *    [-rj/q - I/2, -rj/q + I/2[, as it simplifies the computation and does not
 *    change much because -rj/q is in ]-n/2, n/2].
 *  - it also allows to use Remark 2.7 more efficiently as all computation can
 *    be performed exactly modulo p and no rational nor double computation are
 *    needed.
 *
 * Precomputation:
 *  rho = { rhok := Rk/q }                          [ crt_data_modp attribute ]
 *
 * Let 2 <= n' <= n be an integer, n1 = floor(n/2), n2 = ceil(n/2).
 * Define T1, T2 as the following **sorted** arrays:
 *  T1 = {(sum_{k=1}^{n}{dk*Rk/q} mod p), l(d1,...,dn))
 *                  | dk = 1 if k > n1, dk in {-1,+1} otherwise }
 *  T2 = {(sum_{k=1}^{n}{dk*Rk/q} mod p), l(d1,...,dn))
 *                  | dk = 1 if k <= n1 or k > n', dk in {-1,+1} otherwise }
 *  with l is a value computed from {dk} that will be use to retrieve the
 *  corresponding j at the end (more details in section "j from label"
 *
 * Ti are of size 2^ni and computed in time O(2^ni) using the shift operator
 * defined by Kleinjung.                            [ shift and set_Ti methods ]
 *
 * Then, for each root rp and each possible value of {dk}_{k>n'} in
 * {-1,+1}^(n-n'), compute T'2 from T2
 *  T'2 = {(I/2 + rp/q + sum_{k=1}^{n}{dk*Rk/q} mod p), l(d1,...,dn))
 *                  | dk = 1 if k <= n1, with prescribed value if k > n' and
 *                      dk in {+/-1} otherwise }
 * Again, it is computed using the shift operator.  [ prepare_for_root method ]
 *
 * Finally, iterator over (t1, l1), (t2, l2) in T1 x T's such that t1+t2 is in
 * [0, I[ or [p, p+I[.
 * Note:
 *  t1+t2 = I/2 + rp/q - rj/q with j index of the gray code {(1+dk)/2}
 *  so i = t1+t2-I/2 is in [-I/2, I/2[ if t1+t2 (which is in [0, 2p[) is [0, I[
 *  or in [p, p+I[.
 * The fact that T1 and T'2 are sorted allows this operation to be done very
 * efficiently.
 *
 * # j from label
 *
 * Given a hit (t1, l1), (t2, l2) we want to be able to compute j from l1 and
 * l2. We could store {(1+dk)/2} in the labels, which would give us g = l2 || l1
 * the gray code corresponding to j. Then the k-th bit of j can be computed as
 * XOR(s-th bit of g for s >= k).
 * Instead we store in the label the partial XOR so that we are able to
 * reconstruct j directly as l1 XOR l2.
 */
class siqs_largesieve
{
  public:
    using T_elt = std::pair<uint32_t, uint32_t>;

  protected:
    uint32_t pp;
    uint32_t invq; /* 1/q modulo pp */
    /* crt_data_modp := [ Rk/q mod p for Rk in Q.crt_data_modq ] */
    std::vector<uint32_t> crt_data_modp;
    std::vector<uint32_t> roots; /* r/q modulo pp, for r roots modulo pp */
    size_t n; /* number of factors in the special-q */
    size_t n1, n2;
    std::vector<T_elt> T1, T2;
    slice_offset_t hint;

    /* Given a sorted vector (according to its first component) T, build a new
     * sorted vector (according to its first component) Tt defined as
     *      Tt = sorted(((t+r) % pp, l xor mask) for (t, l) in T)
     *
     *  Assumes 0 <= r < pp
     */
    void shift(
            std::vector<T_elt> & Tt,
            std::vector<T_elt> const & T,
            uint32_t const r,
            uint32_t const mask) const
    {
        Tt.reserve(T.size()); /* Tt will have the same size as T */

        auto over = [add = r - this->pp, mask](T_elt const & e) -> T_elt {
            return { e.first + add, e.second xor mask};
        };
        auto under = [add = r, mask](T_elt const & e) -> T_elt {
            return { e.first + add, e.second xor mask};
        };

        /* despite the name, it performs a binary search */
        auto zp = std::ranges::lower_bound(T, pp-r, {}, &T_elt::first);
        /* zp and the element after zp will need a reduction mod pp when r will
         * be added
         */

        /* Apply 'over' to [zp, end[ and store the result to the start of Tt */
        std::ranges::transform(zp, T.end(), std::back_inserter(Tt), over);
        /* Apply 'under' to [begin, zp[ and store the result in Tt */
        std::ranges::transform(T.begin(), zp, std::back_inserter(Tt), under);
    }

    void set_Ti(std::vector<T_elt> & T, size_t start, size_t end) const
    {
        T = {{0U, 0U}};
        std::vector<T_elt> Tp, Tm;
        uint32_t mask = ((uint32_t) 1U << (start+1U)) - 1U;

        auto comp = [](T_elt const & u, T_elt const & v) -> bool {
            return u.first < v.first;
        };

        for (size_t k = start; k < end; ++k, mask=(mask << 1U) + 1U) {
            size_t new_size = 2U*T.size();
            shift(Tp, T, crt_data_modp[k], mask);
            shift(Tm, T, pp-crt_data_modp[k], 0U);

            T.clear();
            T.reserve(new_size);
            std::ranges::merge(Tp, Tm, std::back_inserter(T), comp);
            Tp.clear();
            Tm.clear();
        }
    }

  public:
    template<class FB_ENTRY_TYPE>
    siqs_largesieve(
            siqs_special_q_data const & Q,
            FB_ENTRY_TYPE const & e,
            size_t n,
            slice_offset_t const hint)
        : pp(e.get_q())
        , n(n)
        , n1((n+1U)/2U)
        , n2(n-n1)
        , hint(hint)
    {
        ASSERT_ALWAYS(2U <= n && n <= Q.nfactors());
        e.compute_crt_data_modp(invq, crt_data_modp, Q, false);
        /* copy roots */
        if constexpr (std::is_same_v<FB_ENTRY_TYPE, fb_entry_general>)
        {
            for(unsigned char i=0U; i < e.nr_roots; ++i) {
                roots.push_back(e.roots[i].r);
            }
        } else {
            for (auto const & r: e.roots) {
                roots.push_back(r);
            }
        }
        /* transform them: r/q mod pp */
        for (auto & r: roots) {
            if (UNLIKELY(!(pp & 1U))) { /* p power of 2 */
                r = (r * invq) & (pp - 1U);
            } else {
                r = mulmodredc_u32<true>(r, invq, pp, e.invq);
            }
        }

        set_Ti(T1, 0U, n1);
        set_Ti(T2, n1, n1+n2);
    }

    void prepare_for_root(
            std::vector<T_elt> & T2s,
            uint32_t const root,
            int const logI,
            uint32_t const j_high) const
    {
        T2s.clear();
        uint32_t mask = ((uint32_t) 1U << n);
        bool is_nth_bit_set = j_high & mask;
        mask = mask - 1U; /* mask is 0...01..1 with n 1's */
        ASSERT_EXPENSIVE(!(j_high & mask));

        uint32_t s = (root + (1U << (logI-1))) % pp;
        uint32_t g = siqs_special_q_data::gray_code_from_j(j_high >> n);
        for (size_t i = n; i < crt_data_modp.size(); ++i, g>>=1) {
            if (g & 1U) {
                s = (s + crt_data_modp[i]) % pp;
            } else {
                s = (s + (pp - crt_data_modp[i])) % pp;
            }
        }
        ASSERT_EXPENSIVE(g == 0);
        /* We do not need to put j_high in the mask because it would force us to
         * remove it later later to compute the bucket index. But we still need
         * to xor the mask if the nth bit of j_high is 1 in order to compute the
         * correct lowest bits of j later.
         */
        shift(T2s, T2, s, (is_nth_bit_set ? mask : 0U));
    }

    uint32_t get_pp() const {
        return pp;
    }

    template <typename BA_t>
    friend void fill_in_buckets_siqs_compute_hits(
            siqs_largesieve & ple,
            BA_t & BA,
            int const logI,
            int const logB,
            uint32_t const j_high,
            slice_index_t const slice_index,
            std::vector<T_elt> & T2scratch,
            where_am_I & w);
};

#endif /* CADO_SIQS_LARGESIEVE_HPP */
