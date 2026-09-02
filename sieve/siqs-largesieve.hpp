#ifndef CADO_SIQS_LARGESIEVE_HPP
#define CADO_SIQS_LARGESIEVE_HPP

#include <cstdint>
#include <cstddef>

#include <array>
#include <concepts>
#include <span>
#include <utility>

#include "bucket.hpp"
#include "fb.hpp"
#include "fb-types.hpp"
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
 * Let 2 <= n' <= n be an integer, n1 = floor(n'/2), n2 = ceil(n'/2).
 * Define T1, T2 as the following **sorted** arrays:
 *  T1 = {(sum_{k=1}^{n1}{dk*Rk/q} mod p), l(d1,...,dn1)) | dk in {-1,+1}}
 *  T2 = {(sum_{k=n1+1}^{n'}{dk*Rk/q} mod p), l(dn1+1,...,dn')) | dk in {-1,+1}}
 *  with l is a value computed from {dk} that will be use to retrieve the
 *  corresponding j at the end (more details in section "j from label"
 *
 * Ti are of size 2^ni and computed in time O(2^ni) using the shift operator
 * defined by Kleinjung.                            [ shift and set_Ti methods ]
 *
 * Then, for a given root rp and given values {dk}_{n'<k<=n} in {-1,+1}^(n-n'),
 * compute T'2 from T2 as
 *  T'2 = {(s+t2, label from l2 and {dk}) | (t2, l2) in T2}
 *  with s = I/2 + rp/q + sum_{n'+1}^{n}{dk*Rk/q} mod p
 * Again, it is computed using the shift operator.  [ prepare_for_root method ]
 *
 * Finally, iterator over (t1, l1), (t2, l2) in T1 x T'2 such that t1+t2 is in
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
    struct T_elt
    {
        uint32_t v; /* value */
        uint32_t l; /* label */
    };

  protected:
    uint32_t pp;
    std::span<std::byte> mem;
    bool mem_owner;
    /* crt_data_modp := [ Rk/q mod p for Rk in Q.crt_data_modq ] */
    std::span<uint32_t> crt_data_modp;
    std::array<uint32_t, 2> roots; /* r/q modulo pp, for r roots modulo pp */
    unsigned int nroots;
    size_t n; /* number of non-fixed values of dk (should be less or equal to
               * crt_data_modp.size(); corresponds to n' in the comments)
               */
    std::span<T_elt> T1, T2;
    slice_offset_t hint;

    /* Given a sorted vector (as pointer + size) whose length is a power of 2
     * and a value, return a pointer to the first element not less than the
     * value, or pointer+size if no such value exist.
     *
     * Assumptions:
     *  - size is a power of 2
     *  - ptr[i].first <= ptr[j].first if i <= j
     *
     * Rationale for not using std::ranges::lower_bound.
     * The hope is that the compiler will make the inner loop branchless and
     * thus faster than the std::ranges::lower_bound which seems to implement a
     * classic binary search.
     */
    static T_elt const * binary_search(
            T_elt const * ptr,
            size_t size,
            uint32_t const value)
    {
        ASSERT_EXPENSIVE(std::has_single_bit(size));
        for (size >>= 1 ; size != 0; size >>= 1)
        {
            if (ptr[size].v < value) {
                ptr += size;
            }
        }
        return ptr + (ptr->v < value);
    }

    /* Given a sorted vector (according to its first component) T, build a new
     * sorted vector (according to its first component) Tt defined as
     *      Tt = sorted(((t+r) % pp, l xor mask) for (t, l) in T)
     *
     * Assumptions:
     *  - 0 <= r < pp
     *  - pminusr = pp - r
     *  - in and out point to an allocated memory zone of size at least 'size'
     *  - in and out do not overlap
     *
     * Only used with sizeof...(Mask) <= 1. Note that it may not be optimal if
     * sizeof...(Mask) > 1 as (... xor mask) may be compute multiple times.
     */
    template<typename ...Mask>
        requires (std::same_as<Mask, uint32_t> && ...)
    static void shift(
            T_elt * out,
            T_elt const * const input,
            std::size_t size,
            uint32_t const r,
            uint32_t const pminusr,
            Mask const ...mask)
    {
        const auto * const middle = binary_search(input, size, pminusr);

        T_elt const * const end = input + size;
        for (T_elt const * ptr = middle; ptr != end; ++ptr) {
            /* here v is >= p-r, so v+r >= p and we need to compute
             * v+r-p = v-(p-r)
             */
            *out++ = { ptr->v - pminusr, (ptr->l xor ... xor mask) };
        }
        for (T_elt const * ptr = input; ptr != middle; ++ptr) {
            /* here v is < p-r, so v+r < p and no reduction mod p is needed */
            *out++ = { ptr->v + r, (ptr->l xor ... xor mask) };
        }
    }

    /* Merge the two sorted arrays p1 and p2 into a sorted array out (sorted
     * according to its first component)
     *
     * Assumptions:
     *  - p1 and p2 point to an allocated memory zone of size at least 'size'
     *  - out point to an allocated memory zone of size at least 2*'size'
     */
    static void merge(
            T_elt * out,
            T_elt const * p1,
            T_elt const * p2,
            std::size_t size)
    {
        T_elt const * const end1 = p1 + size;
        T_elt const * const end2 = p2 + size;
        while (p1 != end1 && p2 != end2) {
            if (p1->v < p2->v) {
                *out++ = *p1++;
            } else {
                *out++ = *p2++;
            }
        }
        while (p1 != end1) {
            *out++ = *p1++;
        }
        while (p2 != end2) {
            *out++ = *p2++;
        }
    }

    /* Compute the sorted array T using the shift & merge algorithm.
     * Only consider index i in [start, end[.
     * scratch is used to hold temporary results.
     *
     * Assumptions:
     *  - T and scratch point to an allocated memory zone of size at least
     *    2^(end-start)
     *  - start < end
     *  - end <= crt_data_modp.size()
     */
    void set_Ti(
            T_elt * const T,
            size_t start,
            size_t end,
            T_elt * const scratch) const
    {
        ASSERT_EXPENSIVE(start < end);
        ASSERT_EXPENSIVE(end <= crt_data_modp.size());

        uint32_t mask = ((uint32_t) 1U << (start+1U)) - 1U;

        T_elt * const Tp =scratch;

        T[0] = { .v=0U, .l=0U };
        for (size_t k = start, s = 1U;
                k < end;
                ++k, mask=(mask << 1U) xor 1U, s<<=1U) {
            /* Invariants:
             *  - s = 2^(k-start)
             *  - mask = 0...01...1 with (k+1) 1s at the end
             */
            const uint32_t rk = crt_data_modp[k];
            const uint32_t pminusrk = pp-rk;
            T_elt * const Tm = scratch+s;

            shift(Tp, T, s, rk, pminusrk, mask);
            shift(Tm, T, s, pminusrk, rk);

            merge(T, Tp, Tm, s);
        }
    }

    std::size_t n1() const
    {
        return (n+1U)/2U;
    }

    std::size_t n2() const
    {
        return n/2U;
    }

    static std::span<std::byte> init_mem(
            std::span<std::byte> scratch,
            std::size_t memsize)
    {
        if (!scratch.empty()) {
            ASSERT_EXPENSIVE(scratch.size() >= memsize);
            return scratch;
        } else {
            return { new std::byte[memsize], memsize };
        }
    }

    std::span<T_elt> init_T1() const
    {
        std::byte * ptr = mem.data()+crt_data_modp.size_bytes();
        return { (T_elt *) ptr, 1U << n1() };
    }

    std::span<T_elt> init_T2() const
    {
        std::byte * ptr = mem.data()+crt_data_modp.size_bytes()+T1.size_bytes();
        return { (T_elt *) ptr, 1U << n2() };
    }

    T_elt * get_Tscratch() const
    {
        std::size_t delta = crt_data_modp.size_bytes()+T1.size_bytes()
                                                      +T2.size_bytes();
        return (T_elt *) (mem.data() + delta);
    }

  public:
    template<class FB_ENTRY_TYPE>
    siqs_largesieve(
            siqs_special_q_data const & Q,
            FB_ENTRY_TYPE const & e,
            size_t n,
            slice_offset_t const hint,
            std::span<std::byte> scratch = {})
        : pp(e.get_q())
        , mem(init_mem(scratch, memory_required(Q, n)))
        , mem_owner(scratch.empty())
        , crt_data_modp((uint32_t *) mem.data(), Q.crt_data_modq.size())
        , roots()
        , n(n)
        , T1(init_T1())
        , T2(init_T2())
        , hint(hint)
    {
        ASSERT_EXPENSIVE(2U <= n && n <= Q.nfactors());
        ASSERT_EXPENSIVE(Q.nfactors() <= 32);
        ASSERT_EXPENSIVE(mem.size() >= memory_required(Q, n));
        uint32_t invq; /* 1/q modulo pp */
        e.compute_crt_data_modp(invq, crt_data_modp, Q, false);

        /* copy roots and transform them: r/q mod pp */
        if constexpr (std::is_same_v<FB_ENTRY_TYPE, fb_entry_general>) {
            nroots = e.nr_roots;
        } else {
            nroots = e.roots.size();
        }
        ASSERT_EXPENSIVE(1 <= nroots && nroots <= 2);

        for(unsigned int i = 0U; i < nroots; ++i) {
            uint32_t r;
            if constexpr (std::is_same_v<FB_ENTRY_TYPE, fb_entry_general>) {
                r = e.roots[i].r;
            } else {
                r = e.roots[i];
            }

            if (UNLIKELY(!(pp & 1U))) { /* p power of 2 */
                roots[i] = (r * invq) & (pp - 1U);
            } else {
                roots[i] = mulmodredc_u32<true>(r, invq, pp, e.invq);
            }
        }

        T_elt * Tt = get_Tscratch();
        set_Ti(T1.data(), 0U, n1(), Tt);
        set_Ti(T2.data(), n1(), n, Tt);
    }

    siqs_largesieve(siqs_largesieve const &) = delete;
    siqs_largesieve & operator=(siqs_largesieve const &) = delete;

    /* The moved-from object loose ownership of the memory, i.e.,
     * other.mem_owner is set to false uncondionnally.
     */
    siqs_largesieve(siqs_largesieve && other) noexcept
        : pp(other.pp)
        , mem(other.mem)
        , mem_owner(std::exchange(other.mem_owner, false))
        , crt_data_modp(other.crt_data_modp)
        , roots(other.roots)
        , nroots(other.nroots)
        , n(other.n)
        , T1(other.T1)
        , T2(other.T2)
        , hint(other.hint)
    {
    }

    /* could be impleted (as above) but not needed for now */
    siqs_largesieve & operator=(siqs_largesieve &&) = delete;

    ~siqs_largesieve()
    {
        if (mem_owner) {
            delete[] mem.data();
        }
    }

    static std::size_t memory_required(
            siqs_special_q_data const & Q,
            std::size_t n)
    {
        /* 4 memory allocations:
         * - T1: size = 2^n1; type = T_elt;
         * - T2: size = 2^n2; type = T_elt;
         * - temporary vector: size = 2^n1 (max of T1 and T2); type = T_elt;
         * - crt_data_modp: size = #factors in sq; type = uint32_t.
         */
        return (2*(1U << (n+1U)/2U) + (1U << n/2U)) * sizeof(T_elt)
                + Q.crt_data_modq.size() * sizeof(uint32_t);
    }

    uint32_t prepare_for_root_mask(uint32_t const j_high) const
    {
        uint32_t mask = ((uint32_t) 1U << n);
        const bool is_nth_bit_set = j_high & mask;
        mask = mask - 1U; /* mask is 0...01..1 with n 1's */
        /* We do not need to put j_high in the mask because it would force us to
         * remove it later later to compute the bucket index. But we still need
         * to xor the mask if the nth bit of j_high is 1 in order to compute the
         * correct lowest bits of j later.
         */
        return is_nth_bit_set ? mask : 0U;
    }

    uint32_t prepare_for_root_precomp(
            int const logI,
            uint32_t const j_high) const
    {
        ASSERT_EXPENSIVE(!(j_high & (((uint32_t) 1U << n) - 1U)));

        uint32_t s = 1U << (logI-1); /* by hypothesis: 2^(I-1) < p */
        uint32_t g = siqs_special_q_data::gray_code_from_j(j_high >> n);
        for (size_t i = n; i < crt_data_modp.size(); ++i, g>>=1) {
            if (g & 1U) {
                s = addmod_u32(s, crt_data_modp[i], pp);
            } else {
                s = submod_u32(s, crt_data_modp[i], pp);
            }
        }
        ASSERT_EXPENSIVE(g == 0);
        return s;
    }

    std::span<T_elt> prepare_for_root(
            uint32_t const root,
            uint32_t s,
            uint32_t const mask) const
    {
        s = addmod_u32(root, s, pp);
        const std::span<T_elt> T2s { get_Tscratch(), T2.size() };
        if (mask) {
            shift(T2s.data(), T2.data(), T2.size(), s, pp-s, mask);
        } else {
            shift(T2s.data(), T2.data(), T2.size(), s, pp-s);
        }
        return T2s;
    }

    template <typename BA_t>
    friend void fill_in_buckets_siqs_compute_hits(
            siqs_largesieve & ple,
            BA_t & BA,
            int logI,
            int logB,
            uint32_t j_high,
            slice_index_t slice_index,
            where_am_I & w);
};

#endif /* CADO_SIQS_LARGESIEVE_HPP */
