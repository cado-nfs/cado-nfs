#ifndef CADO_MPZ_POLY_BIVARIATE_HPP
#define CADO_MPZ_POLY_BIVARIATE_HPP

/*
 * A bivariate poly structure is used for two things.
 *
 *  - univariate polynomials with coefficients in a finite field, which
 *    may be an extension of GF(p).  In some occasions, we need to do
 *    root finding over GF(p^k). This occurs very rarely in NFS contexts,
 *    and typically only for small p^k. The quantity of such computations
 *    is bound to be a lot smaller than what we do over "normal" finite
 *    fields!
 *
 *  - true bivariate stuff, just in case.
 *
 * This implementation is for dense bivariate polynomials. We build them
 * as a tower (polynomials of polynomials, which are represented as
 * mpz_polys).
 *
 * This is expected to be SLOW ! We do not care a lot about speed, here.
 * Because we're chiefly interested in getting this done quick, this is
 * C++ only.
 */

#include "mpz_poly.h"

#include <istream>
#include <ostream>
#include <type_traits>
#include <vector>
#include <algorithm>

#include <gmp.h>

#include "macros.h"
#include "named_proxy.hpp"
#include "runtime_numeric_cast.hpp"

class cxx_mpz_poly_bivariate : private std::vector<cxx_mpz_poly>
{
    using super = std::vector<cxx_mpz_poly>;
    using self = cxx_mpz_poly_bivariate;

  public:
    static constexpr int number_of_variables = 2;
    cxx_mpz_poly_bivariate() {}
#if 0
    cxx_mpz_poly_bivariate(int d)
        ATTRIBUTE_DEPRECATED /* it's too dangerous */
        : super(d+1) {}
#endif
    int degree() const { return ((int)size()) - 1; }
    int degree_y() const { return degree(); }
    int degree_x() const
    {
        int d = -1;
        for (auto const & c: *this) {
            d = std::max(d, c->deg);
        }
        return d;
    }
    int degree_in_yi(unsigned int i) const
    {
        if ((int)i >= degree())
            return -1;
        return mpz_poly_degree((*this)[i]);
    }
    int degree_in_xi(unsigned int i) const
    {
        for (int d = degree_y(); d >= 0; d--) {
            if (mpz_poly_degree((*this)[d]) >= (int)i &&
                mpz_cmp_ui(mpz_poly_coeff_const((*this)[d], i), 0) != 0)
                return d;
        }
        return -1;
    }
    self & operator=(self const & o) = default;
    cxx_mpz_poly_bivariate(self const &) = default;
    ~cxx_mpz_poly_bivariate() = default;
    cxx_mpz_poly_bivariate(self && o) = default;
    cxx_mpz_poly_bivariate & operator=(self && o) = default;
#if 0
    private:
    /* non-const accesses have the potential to ruin the consistency. So
     * this must be restricted to the implementation.
     */
    cxx_mpz_poly & operator[](unsigned int d) { return super::operator[](d); }
    public:
#endif
    cxx_mpz_poly const & operator[](unsigned int d) const
    {
        return super::operator[](d);
    }
    std::string print_poly(std::string const & var) const;

    void swap(self & a) noexcept { ((super &)*this).swap((super &)a); }

  private:
    void cleandeg(int deg)
    {
        ASSERT_ALWAYS(deg >= -1);
        if (runtime_numeric_cast<size_t>(deg + 1) >= size())
            deg = runtime_numeric_cast<int>(size() - 1);
        erase(begin() + (deg + 1), end());
        while (!empty() && super::back() == 0)
            pop_back();
    }

  public:
    /* setcoeffs from an array ? Do we need that ? */

    template <typename T>
    self & operator=(T c)
    requires std::is_integral_v<T>
    {
        assign(1, cxx_mpz_poly(c));
        cleandeg(0);
        return *this;
    }
    self & operator=(mpz_srcptr c)
    {
        assign(1, cxx_mpz_poly(c));
        cleandeg(0);
        return *this;
    }

    struct lifted_x {
        mpz_poly_srcptr c;
        lifted_x(mpz_poly_srcptr c)
            : c(c)
        {
        }
        // operator mpz_poly_srcptr() const { return c; }
        mpz_poly_srcptr operator->() const { return c; }
    };

    struct lifted_y {
        mpz_poly_srcptr c;
        lifted_y(mpz_poly_srcptr c)
            : c(c)
        {
        }
        // operator mpz_poly_srcptr() const { return c; }
        mpz_poly_srcptr operator->() const { return c; }
    };

    self & operator=(lifted_x const & c)
    {
        assign(1, c.c);
        cleandeg(0);
        return *this;
    }

    self & operator=(lifted_y const & c)
    {
        clear();
        reserve(c->deg + 1);
        for (int i = 0; i <= c->deg; i++) {
            emplace_back(mpz_poly_coeff_const(c.c, i));
        }
        return *this;
    }

    /* This is debatable. Do we need to pay the price of having a
     * potentially ambiguous interface here ?
     * Implicitly, we understand c as meaning c(x). But it's also
     * possible to think of it as c(y)...
     */
    self & operator=(mpz_poly_srcptr c) { return *this = lifted_x(c); }

    template <typename T>
    cxx_mpz_poly_bivariate(T c)
    requires std::is_integral_v<T>
    {
        *this = c;
    }
    cxx_mpz_poly_bivariate(mpz_srcptr c) { *this = c; }
    cxx_mpz_poly_bivariate(lifted_x const & c) { *this = c; }
    cxx_mpz_poly_bivariate(lifted_y const & c) { *this = c; }
    cxx_mpz_poly_bivariate(mpz_poly_srcptr c) { *this = c; }

    cado::named_proxy<self &> named(std::string const & x, std::string const & y)
    {
        return { *this, x, y };
    }
    cado::named_proxy<self const &> named(std::string const & x,
                                    std::string const & y) const
    {
        return { *this, x, y };
    }

    static self yi(unsigned int i)
    {
        self f;
        return f.set_xi(i);
    }

    self & set_yi(unsigned int i)
    {
        super::assign(i + 1, {});
        mpz_poly_set_xi(back(), 0);
        return *this;
    }

    static self xi(unsigned int i)
    {
        self f;
        return f.set_xi(i);
    }

    self & set_xi(unsigned int i)
    {
        super::assign(1, {});
        mpz_poly_set_xi(back(), i);
        return *this;
    }

    /* Do we want to keep analogues of the set_ab functions ? */
    /*
    void mpz_poly_bivariate_set_ab (mpz_poly_bivariate_ptr rel, int64_t a,
    uint64_t b); void mpz_poly_bivariate_init_set_mpz_ab (mpz_poly_bivariate_ptr
    rel, mpz_srcptr a, mpz_srcptr b);
    */

#if 0
    private:
    cxx_mpz_poly & lc() { return back(); }
    public:
#endif
    cxx_mpz_poly const & lc() const { return back(); }

    bool normalized_p() const
    {
        return empty() || (mpz_poly_normalized_p(lc()) && lc()->deg != -1);
    }

    bool monic_p() const
    {
        return empty() || (lc()->deg == 0 &&
                           mpz_cmp_ui(mpz_poly_coeff_const(lc(), 0), 1) == 0);
    }

    int cmp(self const & o) const
    {
        int r = (size() > o.size()) - (o.size() > size());
        if (r)
            return r;
        for (int d = size(); --d > 0;) {
            r = mpz_poly_cmp((*this)[d], o[d]);
            if (r)
                return r;
        }
        return 0;
    }

    bool operator<(self const & o) const { return cmp(o) < 0; }
    bool operator==(self const & o) const { return cmp(o) == 0; }

    /* We don't have any operator overloads, on purpose. No reason to
     * have this one specifically
    self operator-() const {
        self f(size());
        for(size_t i = 0 ; i < size() ; i++) {
            mpz_poly_neg(f[i], (*this)[i]);
        }
        return f;
    }
     */

    static void neg(self &, self const &);
    static void add(self &, self const &, self const &);
    static void sub(self &, self const &, self const &);
    /*
     * we might want to special-case these, but it does not seem super
     * important.
    static void add(self &, self const &, self::lifted_x const &);
    static void sub(self &, self const &, self::lifted_x const &);
    static void add(self &, self const &, self::lifted_y const &);
    static void sub(self &, self const &, self::lifted_y const &);
    */
    static void mod_mpz(self &, self const &, mpz_srcptr);
    static void mod_fx(self &, self const &, mpz_poly_srcptr);
    static void mod_fy(self &, self const &, self const &);
    static void mul(self &, self const &, self const &);
    static void mul(self &, self const &, mpz_poly_srcptr);
    static void mul(self &, self const &, mpz_srcptr);
    static void pow_ui(self &, self const &, unsigned long);
    static void mod_fy(self & a, self const & b, mpz_poly_srcptr fy)
    {
        self FY {lifted_y(fy)};
        mod_fy(a, b, FY);
    }
    static void div_qr(self & q, self & r, self const & f, self const & g);

    /* this substitutes the variable y with the given evaluation
     * polynomial, and returns a cxx_mpz_poly */
    static void eval_fy(cxx_mpz_poly & a, self const & f,
                        cxx_mpz_poly const & e);
    /* this substitutes the variable x with the given evaluation
     * integer, and returns a cxx_mpz_poly */
    static void eval_fx(cxx_mpz_poly & a, self const & f, mpz_srcptr e);

    /* computes f(y,x) from f(x,y) */
    static void transpose(self &, self &&);
    static void transpose(self &, self const &);

    static void resultant_y(cxx_mpz_poly & R, self const & f, self const & g);
    static void resultant_x(cxx_mpz_poly & R, self const & f, self const & g);

    static void set_rrandomb(self & f, int dx, int dy, int bits,
                             gmp_randstate_ptr rstate);
    static void set_rrandomb_cab(self & f, int dx, int dy, int bits,
                                 gmp_randstate_ptr rstate);
};

/* printing needs a way to specify the variables... */
namespace cado {
std::ostream & operator<<(
    std::ostream & o,
    cado::named_proxy<cxx_mpz_poly_bivariate const &> const & f);

std::istream &
operator>>(std::istream & in,
           cado::named_proxy<cxx_mpz_poly_bivariate &> f);
} /* namespace cado */

/* we do have a default behaviour, though */
inline std::ostream & operator<<(std::ostream & o,
                                 cxx_mpz_poly_bivariate const & f)
{
    return o << f.named("x", "y");
}

inline std::istream & operator>>(std::istream & in, cxx_mpz_poly_bivariate & f)
{
    return in >> f.named("x", "y");
}

/*
void mpz_poly_bivariate_div_xi(mpz_poly_bivariate_ptr g,
mpz_poly_bivariate_srcptr f, int i); void
mpz_poly_bivariate_mul_xi(mpz_poly_bivariate_ptr g, mpz_poly_bivariate_srcptr f,
int i); void mpz_poly_bivariate_mul_xplusa(mpz_poly_bivariate_ptr g,
mpz_poly_bivariate_srcptr f, mpz_srcptr a);


void mpz_poly_bivariate_eval(mpz_ptr res, mpz_poly_bivariate_srcptr f,
mpz_srcptr x); void mpz_poly_bivariate_eval_ui (mpz_ptr res,
mpz_poly_bivariate_srcptr f, unsigned long x); void
mpz_poly_bivariate_eval_diff_ui (mpz_ptr res, mpz_poly_bivariate_srcptr f,
unsigned long x); void mpz_poly_bivariate_eval_diff (mpz_ptr res,
mpz_poly_bivariate_srcptr f, mpz_srcptr x); void
mpz_poly_bivariate_eval_poly(mpz_poly_bivariate_ptr res,
mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr x); void
mpz_poly_bivariate_eval_diff_poly (mpz_poly_bivariate_ptr res,
mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr x); void
mpz_poly_bivariate_eval_mod_mpz(mpz_t res, mpz_poly_bivariate_srcptr f,
mpz_srcptr x, mpz_srcptr m); int
mpz_poly_bivariate_is_root(mpz_poly_bivariate_srcptr poly, mpz_srcptr root,
mpz_srcptr modulus); void mpz_poly_bivariate_eval_several_mod_mpz(mpz_ptr * res,
mpz_poly_bivariate_srcptr * f, int k, mpz_srcptr x, mpz_srcptr m); void
mpz_poly_bivariate_sqr_mod_f_mod_mpz(mpz_poly_bivariate_ptr Q,
mpz_poly_bivariate_srcptr P, mpz_poly_bivariate_srcptr f, mpz_srcptr m,
mpz_srcptr invf, mpz_srcptr invm); void
mpz_poly_bivariate_pow_ui(mpz_poly_bivariate_ptr B, mpz_poly_bivariate_srcptr A,
unsigned long n); void mpz_poly_bivariate_pow_ui_mod_f(mpz_poly_bivariate_ptr B,
mpz_poly_bivariate_srcptr A, unsigned long n, mpz_poly_bivariate_srcptr f); void
mpz_poly_bivariate_pow_mod_f_mod_ui(mpz_poly_bivariate_ptr Q,
mpz_poly_bivariate_srcptr P, mpz_poly_bivariate_srcptr f, mpz_srcptr a, unsigned
long p); void mpz_poly_bivariate_pow_mod_f_mod_mpz(mpz_poly_bivariate_ptr Q,
mpz_poly_bivariate_srcptr P, mpz_poly_bivariate_srcptr f, mpz_srcptr a,
mpz_srcptr p); void mpz_poly_bivariate_pow_ui_mod_f_mod_mpz
(mpz_poly_bivariate_ptr Q, mpz_poly_bivariate_srcptr P,
mpz_poly_bivariate_srcptr f, unsigned long a, mpz_srcptr p); void
mpz_poly_bivariate_derivative(mpz_poly_bivariate_ptr df,
mpz_poly_bivariate_srcptr f); mpz_poly_bivariate*
mpz_poly_bivariate_base_modp_init (mpz_poly_bivariate_srcptr P0, int p, unsigned
long *K, int l); void mpz_poly_bivariate_base_modp_clear (mpz_poly_bivariate *P,
int l); void mpz_poly_bivariate_base_modp_lift(mpz_poly_bivariate_ptr a,
mpz_poly_bivariate *P, int k, mpz_srcptr pk); size_t
mpz_poly_bivariate_sizeinbase (mpz_poly_bivariate_srcptr f, int base); size_t
mpz_poly_bivariate_size (mpz_poly_bivariate_srcptr f); void
mpz_poly_bivariate_infinity_norm(mpz_ptr in, mpz_poly_bivariate_srcptr f);
size_t mpz_poly_bivariate_totalsize (mpz_poly_bivariate_srcptr f);
void mpz_poly_bivariate_gcd_mpz (mpz_poly_bivariate_ptr h,
mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g, mpz_srcptr p);
// compute f = GCD(f,g) mod N. If this fails, put the factor in the last
// given argument.
int mpz_poly_bivariate_pseudogcd_mpz(mpz_poly_bivariate_ptr ,
mpz_poly_bivariate_ptr , mpz_srcptr , mpz_ptr); void
mpz_poly_bivariate_xgcd_mpz(mpz_poly_bivariate_ptr gcd,
mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g, mpz_poly_bivariate_ptr
u, mpz_poly_bivariate_ptr v, mpz_srcptr p); void mpz_poly_bivariate_homography
(mpz_poly_bivariate_ptr Fij, mpz_poly_bivariate_srcptr F, int64_t H[4]); void
mpz_poly_bivariate_homogeneous_eval_siui (mpz_ptr v, mpz_poly_bivariate_srcptr
f, const int64_t i, const uint64_t j); void mpz_poly_bivariate_content (mpz_ptr
c, mpz_poly_bivariate_srcptr F); int mpz_poly_bivariate_has_trivial_content
(mpz_poly_bivariate_srcptr F); void mpz_poly_bivariate_resultant(mpz_ptr res,
mpz_poly_bivariate_srcptr p, mpz_poly_bivariate_srcptr q); void
mpz_poly_bivariate_discriminant(mpz_ptr res, mpz_poly_bivariate_srcptr f); int
mpz_poly_bivariate_squarefree_p(mpz_poly_bivariate_srcptr f); int
mpz_poly_bivariate_is_irreducible_z(mpz_poly_bivariate_srcptr f);

int mpz_poly_bivariate_number_of_real_roots(mpz_poly_bivariate_srcptr f);

struct mpz_poly_bivariate_with_m_s {
    mpz_poly_bivariate f;
    int m;
};
typedef struct mpz_poly_bivariate_with_m_s mpz_poly_bivariate_with_m[1];
typedef struct mpz_poly_bivariate_with_m_s * mpz_poly_bivariate_with_m_ptr;
typedef const struct mpz_poly_bivariate_with_m_s *
mpz_poly_bivariate_with_m_srcptr;

struct mpz_poly_bivariate_factor_list_s {
    mpz_poly_bivariate_with_m * factors;
    int alloc;
    int size;
};
typedef struct mpz_poly_bivariate_factor_list_s
mpz_poly_bivariate_factor_list[1]; typedef struct
mpz_poly_bivariate_factor_list_s * mpz_poly_bivariate_factor_list_ptr; typedef
const struct mpz_poly_bivariate_factor_list_s *
mpz_poly_bivariate_factor_list_srcptr;

void mpz_poly_bivariate_factor_list_init(mpz_poly_bivariate_factor_list_ptr l);
void mpz_poly_bivariate_factor_list_clear(mpz_poly_bivariate_factor_list_ptr l);
void mpz_poly_bivariate_factor_list_flush(mpz_poly_bivariate_factor_list_ptr l);
void mpz_poly_bivariate_factor_list_push(mpz_poly_bivariate_factor_list_ptr l,
mpz_poly_bivariate_srcptr f, int m); void
mpz_poly_bivariate_factor_list_fprintf(FILE* fp,
mpz_poly_bivariate_factor_list_srcptr l); int
mpz_poly_bivariate_factor_sqf(mpz_poly_bivariate_factor_list_ptr lf,
mpz_poly_bivariate_srcptr f, mpz_srcptr p); int
mpz_poly_bivariate_factor_ddf(mpz_poly_bivariate_factor_list_ptr lf,
mpz_poly_bivariate_srcptr f0, mpz_srcptr p); int
mpz_poly_bivariate_factor_edf(mpz_poly_bivariate_factor_list_ptr lf,
mpz_poly_bivariate_srcptr f, int k, mpz_srcptr p, gmp_randstate_t rstate);

// output is sorted by degree and lexicographically
int mpz_poly_bivariate_factor(mpz_poly_bivariate_factor_list lf,
mpz_poly_bivariate_srcptr f, mpz_srcptr p, gmp_randstate_t rstate); int
mpz_poly_bivariate_is_irreducible(mpz_poly_bivariate_srcptr f, mpz_srcptr p);

*/

#endif /* CADO_MPZ_POLY_BIVARIATE_HPP */
