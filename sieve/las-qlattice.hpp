#ifndef CADO_LAS_QLATTICE_HPP
#define CADO_LAS_QLATTICE_HPP

#include <exception>
#include <cstdint>
#include <iosfwd>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "cado_poly.h"
#include "special-q.hpp"
#include "fb-types.hpp"

struct special_q_data_base {
    special_q doing;
    // handy to have here.
    sublat_t sublat;

    special_q_data_base() = default;
    explicit special_q_data_base(special_q const & doing)
        : doing(doing)
    {}

    virtual ~special_q_data_base() = default;

    virtual void convert_ij_to_ab(
            int64_t & a,
            uint64_t & b,
            int i,
            unsigned int j) const = 0;

    virtual void convert_ij_to_ab(
            cxx_mpz & a,
            cxx_mpz & b,
            int i,
            unsigned int j) const = 0;

    virtual int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            int64_t a,
            uint64_t b) const = 0;

    virtual int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            cxx_mpz const & a,
            cxx_mpz const & b) const = 0;

    virtual std::ostream & print(std::ostream & os) const = 0;

    // Assumes ell is prime.
    bool is_coprime_to(unsigned long ell) const {
        return doing.is_coprime_to(ell);
    }
    size_t nfactors() const {
        return doing.prime_factors.size();
    }
};

std::ostream& operator<<(std::ostream& os, special_q_data_base const & Q);

struct qlattice_basis : public special_q_data_base {
    int64_t a0=0, b0=0, a1=0, b1=0;

    struct too_skewed : public std::exception { };

    qlattice_basis() = default;
    qlattice_basis(special_q const & doing, double skew);
    /* This ctor is used for tests */
    qlattice_basis(int64_t a0, int64_t b0, int64_t a1, int64_t b1)
        : a0(a0)
        , b0(b0)
        , a1(a1)
        , b1(b1)
    {}

    void convert_ij_to_ab(
            int64_t & a,
            uint64_t & b,
            int i,
            unsigned int j) const final;

    void convert_ij_to_ab(
            cxx_mpz & a,
            cxx_mpz & b,
            int i,
            unsigned int j) const final;

    int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            int64_t a,
            uint64_t b) const final;

    int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            cxx_mpz const & a,
            cxx_mpz const & b) const final;

    virtual std::ostream & print(std::ostream & os) const final;

    double skewed_norm0(double s) const { return a0*a0/s+b0*b0*s; }
    double skewed_norm1(double s) const { return a1*a1/s+b1*b1*s; }

    bool fits_31bits() const {
        constexpr int64_t t31 = int64_t(1) << 31;
        return a0 >= -t31 && a0 < t31 &&
               a1 >= -t31 && a1 < t31 &&
               b0 >= -t31 && b0 < t31 &&
               b1 >= -t31 && b1 < t31
               ;
    }
};

std::ostream& operator<<(std::ostream& os, qlattice_basis const & Q);

struct siqs_special_q_data : public special_q_data_base {
    std::vector<cxx_mpz> crt_data_modq;

    siqs_special_q_data() = default;
    siqs_special_q_data(special_q const & doing, cxx_cado_poly const & cpoly);

    void convert_ij_to_ab(
            int64_t & a,
            uint64_t & b,
            int i,
            unsigned int j) const final;

    void convert_ij_to_ab(
            cxx_mpz & a,
            cxx_mpz & b,
            int i,
            unsigned int j) const final;

    int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            int64_t a,
            uint64_t b) const final;

    int convert_ab_to_ij(
            int & i,
            unsigned int & j,
            cxx_mpz const & a,
            cxx_mpz const & b) const final;

    virtual std::ostream & print(std::ostream & os) const final;

    static unsigned int gray_code_from_j(unsigned int j)
    {
        return j xor (j >> 1u);
    }

    /* Given j, compute
     *    g the gray code corresponding to j
     *    rj = add((-1)^(k-th bit of g)*Rk for (k, Rk) in
     *                                              enumerate(crt_data_modq))
     * Note: r0 == doing.r
     */
    cxx_mpz root_from_j(unsigned int j) const
    {
        unsigned int g = gray_code_from_j(j);
        cxx_mpz rj = 0U;
        for (auto const & Rk: crt_data_modq) {
            if (g & 1u) {
                mpz_sub(rj, rj, Rk);
            } else {
                mpz_add(rj, rj, Rk);
            }
            g >>= 1u;
        }
        ASSERT_EXPENSIVE(g == 0u);
        return rj;
    }

    /* Given j, compute delta_j = rj-doing.r (see above for definitions)
     *    delta_j = -2*add(Rk for (k, Rk) in enumerate(crt_data_modq)
     *                                                  if k-th bit of g == 1)
     *    with g the gray code corresponding to j
     * Note: As r0 == doing.r, delta_0 = 0.
     */
    cxx_mpz delta_to_r(unsigned int j) const
    {
        unsigned int g = gray_code_from_j(j);
        cxx_mpz dj = 0U;
        for (auto const & Rk: crt_data_modq) {
            if (g & 1u) {
                mpz_add(dj, dj, Rk);
            }
            g >>= 1u;
        }
        ASSERT_EXPENSIVE(g == 0u);
        mpz_mul_2exp(dj, dj, 1);
        mpz_neg(dj, dj);
        return dj;
    }
};

std::ostream& operator<<(std::ostream& os, siqs_special_q_data const & Q);

namespace fmt {
    template <> struct formatter<special_q_data_base>: ostream_formatter {};
    template <> struct formatter<qlattice_basis>: ostream_formatter {};
    template <> struct formatter<siqs_special_q_data>: ostream_formatter {};
}


#endif	/* CADO_LAS_QLATTICE_HPP */
