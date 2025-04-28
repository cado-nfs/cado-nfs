#ifndef CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP
#define CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP

#include <string>
#include <vector>

#include "fmt/format.h"

#include "gmp_aux.h"
#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field.hpp"
#include "cxx_mpz.hpp"
#include "mpz_mat.h"
#include "numbertheory/fmt_helpers.hpp"

class number_field_order {
    friend class number_field;
    friend class number_field_fractional_ideal;
    friend class number_field_prime_ideal;
    friend class number_field_order_element;
    class number_field const & K;
    public:
    std::string name;
    private:
    cxx_mpq_mat basis_matrix;
    cxx_mpq_mat inv_basis_matrix;
    cxx_mpz_mat multiplication_table;

    /* return the order given by this basis (n*n matrix with respect to
     * the number field polynomial basis in alpha)
     */
    number_field_order(class number_field const &, cxx_mpq_mat);
    public:
    class number_field const & number_field() const { return K; }

    number_field_order_element operator()(cxx_mpz_mat const &) const;

    /* return the i-th basis element */
    number_field_element operator[](int i) const;

    std::vector<number_field_element> basis() const;

    void bless(std::string const & name) { this->name = name; }

    number_field_fractional_ideal fractional_ideal(std::vector<number_field_element> const & gens) const;

    number_field_fractional_ideal p_radical(cxx_mpz const& p) const;

    std::vector<std::pair<number_field_prime_ideal, int>> factor(cxx_mpz const &, cxx_gmp_randstate &) const;
    std::vector<std::pair<number_field_prime_ideal, int>> factor(cxx_mpz const & p) const;
    std::vector<number_field_prime_ideal> factor_radical(cxx_mpz const &, cxx_gmp_randstate &) const;
    std::vector<number_field_prime_ideal> factor_radical(cxx_mpz const & p) const;

    number_field_order(number_field_order const & a)
        : K(a.K)
        , basis_matrix(a.basis_matrix)
        , inv_basis_matrix(a.inv_basis_matrix)
        , multiplication_table(a.multiplication_table)
    {}

    number_field_order(number_field_order && a) noexcept
        : K(a.K)
        , basis_matrix(std::move(a.basis_matrix))
        , inv_basis_matrix(std::move(a.inv_basis_matrix))
        , multiplication_table(std::move(a.multiplication_table))
    {}
};

namespace fmt {
    template <>
    struct formatter<number_field_order>
        : formatter<string_view>
        , fmt_helper_sagemath<number_field_order>
    {
        using fmt_helper_sagemath::parse;
        static constexpr const decltype(custom_format) custom_format_default = TEXT;
        auto format(number_field_order const & O, format_context& ctx) const
            -> format_context::iterator;
    };
}
inline std::ostream& operator<<(std::ostream& os, number_field_order const & K) { return os << fmt::format("{}", K); }



#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP_ */
