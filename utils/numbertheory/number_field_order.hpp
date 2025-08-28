#ifndef CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP
#define CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP

#include <string>
#include <vector>
#include <ostream>
#include <utility>

#include "fmt/format.h"

#include "gmp_aux.h"
#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field.hpp"
#include "cxx_mpz.hpp"
#include "mpz_mat.h"
#include "fmt_helper_sagemath.hpp"

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

    public:
    /* return the order given by this basis (n*n matrix with respect to
     * the number field polynomial basis in alpha)
     */
    number_field_order(class number_field const &, cxx_mpq_mat);

    class number_field const & number_field() const { return K; }

    number_field_order_element operator()(cxx_mpz_mat const &) const;

    number_field_order p_maximal_order(cxx_mpz const & p) const;

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

    /* returns the smallest positive integer z such that z*a is in the
     * order */
    cxx_mpz index(number_field_element const & a) const;

    /*
    number_field_order(number_field_order const & a) = default;
    number_field_order& operator=(number_field_order const & a) = delete;
    number_field_order(number_field_order && a) noexcept = default;
    number_field_order& operator=(number_field_order && a) = default;
    ~number_field_order() = default;
    */

    /* True if the transformation matrix is in SLn(Zp) */
    bool equal_mod(number_field_order const & O, cxx_mpz const & p) const;
};

namespace fmt {
    template <>
    struct formatter<number_field_order>
        : fmt_helper_sagemath<number_field_order>
    {
        static constexpr const decltype(custom_format) custom_format_default = TEXT;
        auto format(number_field_order const & O, format_context& ctx) const
            -> format_context::iterator;
    };
}       /* namespace fmt */


inline std::ostream& operator<<(std::ostream& os, number_field_order const & K) { return os << fmt::format("{}", K); }



#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_HPP_ */
