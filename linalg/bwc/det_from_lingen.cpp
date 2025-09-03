#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <iostream>
#include <sstream>

#include "fmt/base.h"
#include "fmt/format.h"

#include "arith-generic.hpp"
#include "bw-common.h"
#include "bwc_filenames.hpp"
#include "params.h"


/* Some utility classes and functions to compute the constant and leading terms
 * of the linear generator from the output of lingen.
 * The output of lingen is a polynomial matrix whose determinant is the linear
 * generator of the sequence computed by krylov.
 * Efficiency is not a goal as it will be used with very small matrix.
 */
struct Monomial
{
    cxx_mpz value;
    int exp;

    Monomial operator*(Monomial const & other) const
    {
        Monomial r;
        mpz_mul(r.value, value, other.value);
        r.exp = exp + other.exp;
        return r;
    }

    Monomial operator-(Monomial const & other) const
    {
        Monomial r;
        if (exp == other.exp) {
            mpz_sub(r.value, value, other.value);
            r.exp = exp;
        } else if (exp > other.exp) {
            r = *this;
        } else {
            r = other;
        }
        return r;
    }

    Monomial operator%(cxx_mpz const & m) const
    {
        Monomial r;
        mpz_mod(r.value, value, m);
        r.exp = exp;
        return r;
    }

    Monomial invmod(cxx_mpz const & p) const
    {
        Monomial r;
        r.value = value.invmod(p);
        r.exp = -exp;
        return r;
    }

    bool is_zero() const
    {
        return mpz_sgn(value) == 0;
    }

    friend std::ostream& operator<<(std::ostream& os, Monomial const & e)
    {
        return os << e.value << "*X^" << e.exp;
    }
};


bool is_zero(cxx_mpz const & v)
{
    return mpz_sgn(v) == 0;
}

bool is_zero(Monomial const & v)
{
    return v.is_zero();
}


template<typename T>
class SquareMatrixModp
{
    unsigned int n;
    cxx_mpz p;
    std::vector<T> coeffs;

    /* Set d to the determinant of the submatrix (M[j,k])_{i<=j,k<n} using
     * in-place Gauss elimination.
     */
    void determinant_inner(T & d, unsigned int i)
    {
        if (i+1U == n) {
            d = (*this)(i,i) % p;
        } else if (i+2U == n) {
            d = ((*this)(i,i)*(*this)(i+1,i+1)-(*this)(i,i+1)*(*this)(i+1,i)) % p;
        } else {
            unsigned int pivot_idx = i;
            for (; pivot_idx < n; ++pivot_idx) {
                if (!is_zero((*this)(pivot_idx, i))) {
                    break;
                }
            }
            if (pivot_idx == n) {
                d = (*this)(i, i); /* M(i, i) is zero in this case */
                return;
            }
            /* swap row i and pivot_idx */
            if (pivot_idx != i) {
                for (unsigned int j = i; j < n; ++j) {
                    std::swap((*this)(i, j), (*this)(pivot_idx, j));
                }
            }

            /* substract new ith row to all the other rows below */
            T inv = (*this)(i, i).invmod(p);
            for (unsigned int j = i+1; j < n; ++j) {
                T c = ((*this)(j, i) * inv) % p;
                for (unsigned int k = i; k < n; ++k) {
                    (*this)(j,k) = ((*this)(j,k) - c * (*this)(i, k)) % p;
                }
            }
            determinant_inner(d, i+1);
            d = (d * (*this)(i, i)) % p;
        }
    }

    public:

    using Elt = T;

    SquareMatrixModp(unsigned int n, cxx_mpz const & p)
        : n(n), p(p), coeffs(n*n)
    {
    }

    T const & operator()(unsigned int i, unsigned int j) const
    {
        ASSERT_ALWAYS(i < n && j < n);
        return coeffs[i*n+j];
    }

    T & operator()(unsigned int i, unsigned int j)
    {
        ASSERT_ALWAYS(i < n && j < n);
        return coeffs[i*n+j];
    }

    T determinant() const
    {
        SquareMatrixModp<T> tmp = *this;
        T d;
        tmp.determinant_inner(d, 0U);
        return d;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    SquareMatrixModp<T> const & M)
    {
        for (unsigned int i = 0; i < M.n; ++i) {
            os << (i == 0 ? "[ " : "  ");
            for (unsigned int j = 0; j < M.n; ++j) {
                os << (j == 0 ? "[" : ", ") << M(i, j);
            }
            os << (i+1 == M.n ? "] ]" : "],\n");
        }
        return os;
    }
};


namespace fmt {
    template <> struct formatter<Monomial>: ostream_formatter {};
    template <> struct formatter<SquareMatrixModp<cxx_mpz>>: ostream_formatter {};
    template <> struct formatter<SquareMatrixModp<Monomial>>: ostream_formatter {};
}

static int
det_from_lingen_prog(unsigned int n, cxx_mpz const & p,
                     std::string const & ffile, unsigned int charpoly_deg)
{
    ASSERT_ALWAYS(p != 2U);

    unsigned int const width = 1U;

    /* Define and initialize our arithmetic back-ends. */
    std::unique_ptr<arith_generic> A(arith_generic::instance(p, width));
    fmt::print("# Using '{}' for arith\n", A->impl_name());

    size_t const one_fcoeff = A->vec_elt_stride(A->simd_groupsize());

    SquareMatrixModp<cxx_mpz> M0(n, p);
    SquareMatrixModp<Monomial> Mlt(n, p);

    for (unsigned int i = 0; i < n; i+=width) {
        for (unsigned int j = 0; j < n; j+=width) {
            auto tmp = A->alloc_vector(one_fcoeff);
            auto lt = A->alloc_vector(one_fcoeff);
            size_t lt_deg = SIZE_MAX;

            std::string const f_name = fmt::format("{}.sols{}-{}.{}-{}",
                                                    ffile, i, (i+1)*width,
                                                           j, (j+1)*width);
            fmt::print("# reading from {}\n", f_name);
            auto f = fopen_helper(f_name, "rb");
            size_t nread = 0;
            while (true) {
                int rc = std::fread(tmp.get(), one_fcoeff, 1U, f.get());
                if (rc < 1) {
                    break;
                }
                if (!nread) {
                    std::stringstream s;
                    A->cxx_out(s, *tmp);
                    s >> M0(i, j);
                }
                if (!A->is_zero(*tmp)) {
                    A->set(*lt, *tmp);
                    lt_deg = nread;
                }
                nread++;
            }
            ASSERT_ALWAYS(nread > 0); /* File should not be empty */
            std::stringstream s;
            A->cxx_out(s, *lt);
            s >> Mlt(i, j).value;
            Mlt(i, j).exp = lt_deg;
        }
    }

    fmt::print("Matrix of degree-0 terms:\n{}\n", M0);
    fmt::print("Matrix of leading terms:\n{}\n", Mlt);

    auto d0 = M0.determinant();
    auto dlt = Mlt.determinant();

    fmt::print("Determinant of matrix of degree-0 terms: {}\n", d0);
    fmt::print("Determinant of matrix of leading terms: {}\n", dlt);

    if (dlt.exp < 0 || (unsigned int) dlt.exp > charpoly_deg) {
        fmt::print(stderr, "Error, the degree of the determinant of the matrix "
                           "of leading terms should be in [0, {}]\n",
                           charpoly_deg);
        return EXIT_FAILURE;
    } else if ((unsigned int) dlt.exp < charpoly_deg) {
        fmt::print("# Degree is smaller than {}, cannot compute the "
                   "determinant\n", charpoly_deg);
        fmt::print("determinant modulo p: inconclusive\n");
        return EXIT_FAILURE;
    } else {
        cxx_mpz inv = dlt.value.invmod(p);
        mpz_mul(d0, d0, inv);
        mpz_mod(d0, d0, p);
        fmt::print("determinant modulo p: {}\n", d0);
        return EXIT_SUCCESS;
    }
}


// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

    bw_common_decl_usage(pl);
    /* {{{ declare local parameters and switches */
    param_list_decl_usage(pl, "ffile", "generator file");
    param_list_decl_usage(pl, "charpoly-degree", "degree of the characteristic "
                                                 "polynomial");
    /* }}} */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);

    /* {{{ interpret our parameters */
    std::string ffile;
    param_list_parse(pl, "ffile", ffile);
    if (ffile.empty()) {
        fmt::print(stderr, "Error, -ffile must be a nonempty string\n");
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    unsigned int charpoly_deg = UINT_MAX;
    param_list_parse(pl, "charpoly-degree", charpoly_deg);
    if (charpoly_deg == UINT_MAX) {
        fmt::print(stderr, "Error, -charpoly-degree must be specified\n");
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }
    /* }}} */

    if (mpz_cmp_ui(bw->p, 2U) == 0) {
        fmt::print(stderr, "Error, -p cannot be 2 in this binary\n");
        exit(EXIT_FAILURE);
    }

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    int ret = det_from_lingen_prog(bw->n, bw->p, ffile, charpoly_deg);

    bw_common_clear(bw);

    return ret;
}
