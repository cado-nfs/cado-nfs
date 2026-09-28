#ifndef SIEVE_LAS_SUBLAT_HPP_
#define SIEVE_LAS_SUBLAT_HPP_

#include <vector>
#include <string>

#include "fmt/format.h"

#include "las-config.hpp"
#include "macros.h"

#include <cstdint>

// FIXME: could probably go somewhere else...
// Small struct for sublattice info:
// One sieves only positions congruent to (i0,j0) mod m
struct sublat_runtime_t;
template<uint32_t m>
struct sublat_t {
    static_assert(m == 1 || m == 2 || m == 3 || m == 6);
    static constexpr uint32_t modulus = m;
    uint32_t i0=0;
    uint32_t j0=0;

    void adjustIJ(int & i, unsigned int & j) const
    {
        i = i*m + i0;
        j = j*m + j0;
    }
    static bool not_coprime(uint32_t q);
    static bool coprime(uint32_t q) { return !not_coprime(q); }
    sublat_t() = default;
    inline explicit sublat_t(sublat_runtime_t const &);
    bool is_first() const { return m == 1 || (i0 == 0 && j0 == 1); }
};
template<>
inline void sublat_t<1>::adjustIJ(int &, unsigned int &) const {}
#ifndef BUCKET_SIEVE_POWERS
template<uint32_t M>
inline bool sublat_t<M>::not_coprime(uint32_t q) { ASSERT(q > M); return false; }
#else
template<>
inline bool sublat_t<1>::not_coprime(uint32_t) { return false; }
/* that is typically optimized away by serious compilers. mod 3 is
 * (q*0xAA...AAA) < 0x55..55 ; mod 6 is mod2 or mod3
 */
template<uint32_t M>
inline bool sublat_t<M>::not_coprime(uint32_t q) { return (M % 2 == 0 && q % 2 == 0) || (M % 3 == 0 && q % 3 == 0); }
#endif
struct sublat_runtime_t {
    uint32_t m=1;
    uint32_t i0=0;
    uint32_t j0=0;

    sublat_runtime_t() = default;

    constexpr sublat_runtime_t(uint32_t m, uint32_t i0, uint32_t j0)
        : m(m), i0(i0), j0(j0) {}

    /* The identity mapping, for the code paths that are only ever reached
     * without sublattices (the pattern-3 and pattern-5 variants). */
    static constexpr sublat_runtime_t no_sublat() { return { 1, 0, 0 }; }

    template<uint32_t M>
    explicit sublat_runtime_t(sublat_t<M> const & s)
    : m(M)
    , i0(s.i0)
    , j0(s.j0)
    {}
    void adjustIJ(int & i, unsigned int & j) const
    {
        if (m > 1) {
            i = i*static_cast<int>(m) + static_cast<int>(i0);
            j = j*m + j0;
        }
    }
    private:
    /* A class is useless when *both* i and j share a prime factor with m:
     * every (i,j) in it then has that factor in gcd(i,j), so it can never
     * be coprime. Testing i and j separately would be far too strong -- it
     * would throw away (0,1), and for m=6 it would keep 4 classes out of
     * the 24 that are worth sieving.
     *
     * This is slow, but it's only used in the non-critical function
     * sublattices()
     */
    bool degenerate(uint32_t i, uint32_t j) const {
        return (m % 2 == 0 && i % 2 == 0 && j % 2 == 0)
            || (m % 3 == 0 && i % 3 == 0 && j % 3 == 0);
    }
    public:
    /* the first class that sublattices() emits: (0,0) for m == 1, and
     * (0,1) otherwise. Getting this wrong for m == 1 means the toplevel
     * fill takes the "replay" branch and dereferences a null
     * precomputed-plattice vector.
     */
    bool is_first() const { return m == 1 || (i0 == 0 && j0 == 1); }

    std::vector<sublat_runtime_t> sublattices() const {
        ASSERT_ALWAYS(m == 1 || m == 2 || m == 3 || m == 6);
        if (m == 1)
            return { sublat_runtime_t(1, 0, 0) };
        std::vector<sublat_runtime_t> res;
        res.reserve(m*m-1);
        for(uint32_t i = 0 ; i < m ; i++)
            for(uint32_t j = 0 ; j < m ; j++)
                if ((i || j) && !degenerate(i, j)) res.emplace_back(m, i, j);
        return res;
    }

    /* the real row */
    unsigned int jj(unsigned int j) const { return m*j + j0; }

    /* |ii| for the sublattice abscissa x, given the region's i0 */
    unsigned int abs_ii(int region_i0, int x) const {
        return (unsigned int) std::abs((long) m * (region_i0 + x) + (long) i0);
    }

    /* Which x have ii even, on a row whose jj is even. The small sieve
     * leaves those alone because gcd(ii,jj) >= 2 there, and the survivor
     * search has to agree with it. 0 means "none of them": that is the
     * case for an even modulus, where the parity of ii is fixed at si0 and
     * the class with both parities even is never sieved. Otherwise i0 is
     * even, so ii has the parity of x + si0.
     *
     * Returns 0 for none, 1 for the even x, 2 for the odd x.
     */
    int parity_skip_class() const {
        if ((m & 1) == 0) return 0;
        return (i0 & 1) ? 2 : 1;
    }
};
template<uint32_t M>
inline sublat_t<M>::sublat_t(sublat_runtime_t const & sublat)
    : i0(sublat.i0)
    , j0(sublat.j0)
{
    ASSERT_ALWAYS(M == sublat.m);
}


inline std::string format_as(sublat_runtime_t const & S)
{
    return fmt::format("(i, j) == ({}, {}) mod {}", S.i0, S.j0, S.m);
}

template<uint32_t M>
inline std::string format_as(sublat_t<M> const & S)
{
    return format_as(sublat_runtime_t(S));
}

/* Call f with the compile-time-typed sublattice matching the runtime one.
 * Every site that needs the modulus as a constant goes through this, so the
 * list of supported moduli is written once. */
template<typename F> auto dispatch_sublat(sublat_runtime_t const & s, F && f) {
    switch (s.m) {
        case 1: return std::forward<F>(f)(sublat_t<1>(s));
        case 2: return std::forward<F>(f)(sublat_t<2>(s));
        case 3: return std::forward<F>(f)(sublat_t<3>(s));
        case 6: return std::forward<F>(f)(sublat_t<6>(s));
        default: ASSERT_ALWAYS(0);
    }
}

#endif	/* SIEVE_LAS_SUBLAT_HPP_ */
