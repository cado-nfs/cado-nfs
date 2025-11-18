#ifndef CADO_UTILS_MATRIX_HPP
#define CADO_UTILS_MATRIX_HPP

/* This defines matrices over arbitrary types, provided that these have
 * standard operator overloads defined. This should provide equivalent
 * interfaces to cxx_mpz_mat, cxx_mpq_mat, and also provide working
 * equivalent for other types as well
 */

#include <ostream>
#include <type_traits>
#include <vector>
#include <utility>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "macros.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "meta_pow.hpp"

namespace cado {

    namespace arithmetic_reductions {
        struct mod_mpz;
        struct mod_fx_mod_p;
    } /* namespace arithmetic_reductions */

template<typename T>
struct matrix {
    friend struct arithmetic_reductions::mod_mpz;
    friend struct arithmetic_reductions::mod_fx_mod_p;

    private:
    std::vector<T> coeffs;
    unsigned int m = 0;
    unsigned int n = 0;

    void resize(unsigned int m, unsigned int n)
    {
        this->m = m;
        this->n = n;
        set_zero();
    }

    public:
    unsigned int nrows() const { return m; }
    unsigned int ncols() const { return n; }

    matrix() = default;
    matrix(unsigned int m, unsigned int n)
        : coeffs(m * n)
        , m(m)
        , n(n)
    { }
    matrix(unsigned int m, unsigned int n, unsigned long z)
        : matrix(m, n)
    {
        for(unsigned int i = 0 ; i < m && i < n ; i++)
            (*this)(i, i) = z;
    }
    explicit matrix(cxx_mpz_mat const & o)
        requires std::is_same_v<T, cxx_mpz>
        : matrix(o.nrows(), o.ncols())
        {
            for(unsigned int i = 0 ; i < m ; i++)
                for(unsigned int j = 0 ; j < n ; j++)
                    mpz_set((*this)(i, j), o(i, j));
        }
    explicit matrix(cxx_mpz_mat && o)
        requires std::is_same_v<T, cxx_mpz>
        : matrix(o.nrows(), o.ncols())
        {
            for(unsigned int i = 0 ; i < m ; i++)
                for(unsigned int j = 0 ; j < n ; j++)
                    mpz_swap((*this)(i, j), o(i, j));
        }

    explicit operator cxx_mpz_mat() const
        requires std::is_same_v<T, cxx_mpz>
        {
            cxx_mpz_mat Z(m, n);
            for(unsigned int i = 0 ; i < m ; i++)
                for(unsigned int j = 0 ; j < n ; j++)
                    mpz_set(Z(i, j), (*this)(i, j));
            return Z;
        }

    ~matrix() = default;
    matrix(matrix const&) = default;
    matrix(matrix &&) = default;
    matrix& operator=(matrix const&) = default;
    matrix& operator=(matrix &&) = default;

    /* all instantations love each other */
    template<typename U> friend struct matrix;
    template<typename U>
        explicit matrix(matrix<U> const & a)
        requires (!std::is_convertible_v<U, T>)
        : coeffs { a.coeffs.begin(), a.coeffs.end() }
        , m(a.m)
        , n(a.n)
    {}

    void set_zero() { coeffs.assign(m * n, 0); }

    matrix& operator=(T v) {
        for(unsigned int i = 0 ; i < m && i < n ; i++)
            (*this)(i, i) = v;
        return *this;
    }

    /* TODO: creation from string? */

    public:

    T const & operator()(unsigned int i, unsigned int j) const
    {
        ASSERT(i < m);
        ASSERT(j < n);
        return coeffs[i * n + j];
    }
    T & operator()(unsigned int i, unsigned int j)
    {
        ASSERT(i < m);
        ASSERT(j < n);
        return coeffs[i * n + j];
    }

    matrix& operator*=(T a) { for(auto & x : coeffs) x *= a; return *this; }
    matrix operator*(T a) const { matrix h = *this; return h *= a; }
    matrix& operator/=(T a) { for(auto & x : coeffs) x /= a; return *this; }
    matrix operator/(T a) const { matrix h = *this; return h /= a; }

    matrix& operator-=(T const & a)
    {
        for(unsigned int i = 0 ; i < m && i < n ; i++)
            (*this)(i, i) -= a;
        return *this;
    }
    matrix operator-(T a) const { matrix h = *this; return h -= a; }

    matrix& operator+=(T const & a)
    {
        for(unsigned int i = 0 ; i < m && i < n ; i++)
            (*this)(i, i) += a;
        return *this;
    }
    matrix operator+(T a) const { matrix h = *this; return h += a; }

    matrix operator*(matrix const & b) const
    {
        matrix const & a = *this;
        matrix c(a.m, b.n);
        ASSERT_ALWAYS(a.n == b.m);
        for(unsigned int i = 0 ; i < a.m ; i++) {
            for(unsigned int j = 0 ; j < b.n ; j++) {
                for(unsigned int k = 0 ; k < a.n ; k++) {
                    c(i,j) += a(i,k) * b(k,j);
                }
            }
        }
        return c;
    }
    matrix& operator*=(matrix const & b) const
    {
        return (*this) = (*this) * b;
    }
    matrix& operator+=(matrix const & b) const
    {
        matrix const & a = *this;
        ASSERT_ALWAYS(a.m == b.m);
        ASSERT_ALWAYS(a.n == b.n);
        for(unsigned int i = 0 ; i < a.m * a.n ; i++)
            a.coeffs[i] += b.coeffs[i];
        return *this;
    }
    matrix operator+(matrix const & a) const {
        matrix h = *this;
        return h += a;
    }
    matrix& operator-=(matrix const & b) const
    {
        matrix const & a = *this;
        ASSERT_ALWAYS(a.m == b.m);
        ASSERT_ALWAYS(a.n == b.n);
        for(unsigned int i = 0 ; i < a.m * a.n ; i++)
            a.coeffs[i] -= b.coeffs[i];
        return *this;
    }
    matrix operator-(matrix const & a) const {
        matrix h = *this;
        return h -= a;
    }

    matrix transpose() const {
        matrix t(n, m);
        matrix const & a = *this;
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                t(j,i) = a(i,j);
            }
        }
        return t;
    }

    static void mul(matrix & c, matrix const & a, matrix const & b) {
        c = a * b;
    }

    template<typename R>
    static void pow(matrix & b, matrix const & a, unsigned long n, R const & r)
    {
        cado::meta_pow<matrix, R> {r} (b, a, n);
    }

    decltype(T() <=> T()) operator<=>(matrix const & a)
    {
        if (auto r = m <=> a.m; r != 0)
            return r;
        if (auto r = n <=> a.n; r != 0)
            return r;
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++)
                if (auto r = (*this)(i, j) <=> a(i, j) ; r != 0)
                    return r;
        return std::strong_ordering::equal;
    }

    template<typename U>
    decltype(T() <=> U()) operator<=>(U const & a) const
    requires std::is_convertible_v<U, T>
    {
        if (auto r = m <=> a.m; r != 0)
            return r;
        if (auto r = n <=> a.n; r != 0)
            return r;
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++)
                if (auto r = (*this)(i, j) <=> (i == j ? a : 0); r != 0)
                    return r;
        return std::strong_ordering::equal;
    }
    template<typename U> bool operator==(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) == 0; }
    /*
    template<typename U> bool operator!=(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) != 0; }
    template<typename U> bool operator<(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) < 0; }
    template<typename U> bool operator>(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) > 0; }
    template<typename U> bool operator<=(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) <= 0; }
    template<typename U> bool operator>=(U const & a) const requires std::is_convertible_v<U, T> { return operator<=>(a) >= 0; }
    */

    template<typename U>
    matrix& operator=(U const & a)
    requires std::is_convertible_v<U, T>
    {
        ASSERT_ALWAYS(m == n);
        set_zero();
        for(unsigned int i = 0 ; i < m && i < n ; i++)
            (*this)(i, i) = a;
        return *this;
    }
};

namespace matrix_reductions {
    struct noop {
    };
} /* namespace matrix_reductions */

} /* namespace cado */

/* we do have a default behaviour, though */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, cado::matrix<T> const & f);

namespace fmt {
    template<typename T>
    struct formatter<cado::matrix<T>>: ostream_formatter {};
}


#endif	/* CADO_UTILS_MATRIX_HPP */
