#ifndef MPZ_POLY_H_
#define MPZ_POLY_H_

// IWYU pragma: no_include "double_poly.h"
// (only the fwd-decl is needed)
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <istream>      // std::istream // IWYU pragma: keep
#include <ostream>      // std::ostream // IWYU pragma: keep
#include <type_traits>
#include <string>
#include <vector>
#include <istream>      // std::istream // IWYU pragma: keep
#include <ostream>      // std::ostream // IWYU pragma: keep
#endif

#ifdef __cplusplus
#include "fmt/ostream.h"
#include "fmt/format.h"
#endif
#include "gmp_aux.h"

#ifdef __cplusplus
#include "cxx_mpz.hpp"
#endif
#include "macros.h"

#define xxxMPZ_POLY_TIMINGS
// for timings of roots mod p (beware, this is not thread-safe)

#ifndef DOUBLE_POLY_H_
typedef struct double_poly_s * double_poly_ptr;
typedef const struct double_poly_s * double_poly_srcptr;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* maximum degree we can reconstruct using mpz_poly_mul_tc_interpolate */
#define MAX_TC_DEGREE 19


/* Note, deg = -1 means P=0; otherwise, one should have coeff[deg] != 0.
   Warning: a polynomial of degree d needs d+1 allocation. */

struct mpz_poly_s {
  unsigned int alloc;
  int deg;
  mpz_t *_coeff;
};
#ifndef DOUBLE_POLY_H_
/* double_poly.h forward-declares these. Don't do it twice */
typedef struct mpz_poly_s * mpz_poly_ptr;
typedef const struct mpz_poly_s * mpz_poly_srcptr;
#endif
typedef struct mpz_poly_s mpz_poly[1];

/* Note on parallelism.
 *
 * Some functions in this API can optionally use openmp. We want to make
 * sure that the code that calls these openmp-enabled functions does so
 * **willingly**. And we don't want to pollute the "normal" interface.
 *
 * The chosen solution is as follows.
 *
 * - the functions declared here, and use as plain (e.g.) mpz_poly_mul
 *   resolve to something that does **NOT** use openmp.
 *
 * - to use openmp, include mpz_poly_parallel.hpp instead, and call the
 *   member functions of an mpz_poly_parallel_info object. E.g.
 *
 *   mpz_poly_parallel_info inf;
 *   // (maybe add some configuration code for the inf object, if the
 *   // need for that ever appears)
 *   inf.mpz_poly_mul(....)
 *
 * More detail on can be found in mpz_poly_parallel.hpp and mpz_poly.cpp
 */

/* We should not access the polynomial coefficients directly. At least we
 * need specific accessors for them.
 * The accessors below behave in the exact same way as old direct
 * accesses that were within (allocation) bounds. Beyond bounds, they
 * force a realloc if we want a read-write access, and they return a
 * const static zero for const accesses.
 */
mpz_ptr mpz_poly_coeff(mpz_poly_ptr, int i);
mpz_srcptr mpz_poly_coeff_const(mpz_poly_srcptr, int i);

/* Management of the structure, set and print coefficients. */
void mpz_poly_init(mpz_poly_ptr, int d);
void mpz_poly_realloc (mpz_poly_ptr f, unsigned int nc);
void mpz_poly_set(mpz_poly_ptr g, mpz_poly_srcptr f);
void mpz_poly_swap (mpz_poly_ptr f, mpz_poly_ptr g);
void mpz_poly_clear(mpz_poly_ptr f);
static inline int mpz_poly_degree(mpz_poly_srcptr f) { return f->deg; }
int mpz_poly_valuation(mpz_poly_srcptr f);

void mpz_poly_cleandeg(mpz_poly_ptr f, int deg);
void mpz_poly_setcoeffs(mpz_poly_ptr f, mpz_t * coeffs, int d);
void mpz_poly_setcoeffs_si(mpz_poly_ptr f, const long int * h, int d);
void mpz_poly_setcoeffs_ui(mpz_poly_ptr f, const unsigned long int * h, int d);
void mpz_poly_set_zero(mpz_poly_ptr f);
void mpz_poly_set_ui(mpz_poly_ptr f, unsigned long z);
void mpz_poly_set_si(mpz_poly_ptr f, long z);
void mpz_poly_set_xi(mpz_poly_ptr f, int i);
void mpz_poly_set_mpz(mpz_poly_ptr f, mpz_srcptr z);
void mpz_poly_set_double_poly(mpz_poly_ptr g, double_poly_srcptr f);
/* returns 1 if parsing was successful */
int mpz_poly_set_from_expression(mpz_poly_ptr f, const char * value);

void mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b);
void mpz_poly_init_set_mpz_ab (mpz_poly_ptr rel, mpz_srcptr a, mpz_srcptr b);
void mpz_poly_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b);
void mpz_poly_set_mpz_ab (mpz_poly_ptr rel, mpz_srcptr a, mpz_srcptr b);

void mpz_poly_setcoeff(mpz_poly_ptr f, int i, mpz_srcptr z);
void mpz_poly_setcoeff_si(mpz_poly_ptr f, int i, long z);
void mpz_poly_setcoeff_ui(mpz_poly_ptr f, int i, unsigned long z);
void mpz_poly_setcoeff_int64(mpz_poly_ptr f, int i, int64_t z);
void mpz_poly_setcoeff_uint64(mpz_poly_ptr f, int i, uint64_t z);
void mpz_poly_setcoeff_double(mpz_poly_ptr f, int i, double z);

void mpz_poly_set_signed_urandomb (mpz_poly_ptr f, int d, gmp_randstate_ptr state, int k);
void mpz_poly_set_signed_rrandomb (mpz_poly_ptr f, int d, gmp_randstate_ptr state, int k);
void mpz_poly_set_signed_urandomm (mpz_poly_ptr f, int d, gmp_randstate_ptr state, mpz_srcptr N);
void mpz_poly_set_urandomb (mpz_poly_ptr f, int d, gmp_randstate_ptr state, int k);
void mpz_poly_set_rrandomb (mpz_poly_ptr f, int d, gmp_randstate_ptr state, int k);
void mpz_poly_set_urandomm (mpz_poly_ptr f, int d, gmp_randstate_ptr state, mpz_srcptr N);
void mpz_poly_set_urandomm_ui (mpz_poly_ptr f, int d, gmp_randstate_ptr state, unsigned long m);

/* functions for Joux-Lercier and generalized Joux-Lercier */
int mpz_poly_setcoeffs_counter(mpz_poly_ptr f, int* max_abs_coeffs, unsigned long *next_counter, int deg, unsigned long counter, unsigned int bound);
void  mpz_poly_setcoeffs_counter_print_error_code(int error_code);
unsigned long mpz_poly_getcounter(mpz_poly_ptr f, unsigned int bound);
unsigned long mpz_poly_cardinality(int deg, unsigned int bound);

/* return the leading coefficient of f */
static inline mpz_srcptr mpz_poly_lc (mpz_poly_srcptr f) {
    ASSERT(f->deg >= 0);
    return mpz_poly_coeff_const(f, f->deg);
}

static inline mpz_ptr mpz_poly_lc_w (mpz_poly_ptr f) {
    ASSERT(f->deg >= 0);
    return mpz_poly_coeff(f, f->deg);
}


/* Print functions */
int mpz_poly_asprintf(char ** res, mpz_poly_srcptr f);
/* Print coefficients of f.
 * endl = 1 if "\n" at the end of fprintf. */
void mpz_poly_fprintf_endl (FILE *fp, mpz_poly_srcptr f, int endl);
void mpz_poly_fprintf(FILE *fp, mpz_poly_srcptr f);
void mpz_poly_fprintf_coeffs (FILE *fp, mpz_poly_srcptr f, const char sep);
void mpz_poly_fscanf_coeffs (FILE *fp, mpz_poly_ptr f, const char sep);
void mpz_poly_fprintf_cado_format (FILE *fp, mpz_poly_srcptr f,
                                   const char letter, const char *pre);
void mpz_poly_asprintf_cado_format (char **pstr, mpz_poly_srcptr f, const char letter,
                              const char *prefix);
void mpz_poly_print_raw(mpz_poly_srcptr f);
#ifdef MPZ_POLY_TIMINGS
  void print_timings_pow_mod_f_mod_p();
#endif
/* Tests and comparison functions */
int mpz_poly_cmp (mpz_poly_srcptr, mpz_poly_srcptr);
int mpz_poly_normalized_p (mpz_poly_srcptr f);
int mpz_poly_is_monic (mpz_poly_srcptr f);
int mpz_poly_is_monomial_multiple(mpz_poly_srcptr f);
int mpz_poly_is_monomial(mpz_poly_srcptr f);

/* Polynomial arithmetic */
void mpz_poly_to_monic(mpz_poly_ptr g, mpz_poly_srcptr f);
void mpz_poly_neg(mpz_poly_ptr f, mpz_poly_srcptr g);
void mpz_poly_add(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_sub(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_add_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a);
void mpz_poly_sub_ui(mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a);
void mpz_poly_add_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a);
void mpz_poly_sub_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr a);
void mpz_poly_sub_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h, mpz_srcptr m);
void mpz_poly_mul(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h);
void mpz_poly_mul_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
void mpz_poly_divexact_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a);
int mpz_poly_divisible_mpz (mpz_poly_srcptr P, mpz_srcptr a);
void mpz_poly_translation (mpz_poly_ptr, mpz_poly_srcptr, mpz_srcptr);

void mpz_poly_rotation (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, mpz_srcptr, int);
void mpz_poly_rotation_si (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, long int, int);
void mpz_poly_rotation_ui (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, unsigned long int, int);
void mpz_poly_rotation_int64 (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, int64_t, int);
void mpz_poly_reverse_rotation (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, mpz_srcptr, int);
void mpz_poly_reverse_rotation_si (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, long int, int);
void mpz_poly_reverse_rotation_ui (mpz_poly_ptr, mpz_poly_srcptr, mpz_poly_srcptr, unsigned long int, int);

void mpz_poly_addmul_si (mpz_poly_ptr, mpz_poly_srcptr, long);
void mpz_poly_mul_si (mpz_poly_ptr, mpz_poly_srcptr, long);
void mpz_poly_divexact_ui (mpz_poly_ptr, mpz_poly_srcptr, unsigned long);
void mpz_poly_makemonic_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr m);
void barrett_precompute_inverse (mpz_ptr invm, mpz_srcptr m);
int mpz_poly_mod_f_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
int mpz_poly_mod_mpz(mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m, mpz_srcptr invm);
int mpz_poly_mod_mpz_lazy (mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m);
void mpz_poly_mul_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
void mpz_poly_mul_mod_f (mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2, mpz_poly_srcptr f);
void mpz_poly_reduce_frac_mod_f_mod_mpz (mpz_poly_ptr num, mpz_poly_ptr denom, mpz_poly_srcptr F, mpz_srcptr m);
int mpz_poly_div_qr_mod_mpz (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_srcptr p);
int mpz_poly_div_r_mod_mpz_clobber (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_srcptr p);
int mpz_poly_div_qr (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g);
int mpz_poly_div_r (mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g);
int mpz_poly_mod (mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g);
int mpz_poly_divexact (mpz_poly_ptr q, mpz_poly_srcptr h, mpz_poly_srcptr f, mpz_srcptr p);
void mpz_poly_div_2_mod_mpz(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr m);
void mpz_poly_div_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i);
void mpz_poly_mul_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i);
void mpz_poly_mul_xplusa(mpz_poly_ptr g, mpz_poly_srcptr f, mpz_srcptr a);

  
void mpz_poly_eval(mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x);
void mpz_poly_eval_ui (mpz_ptr res, mpz_poly_srcptr f, unsigned long x);
void mpz_poly_eval_diff_ui (mpz_ptr res, mpz_poly_srcptr f, unsigned long x);
void mpz_poly_eval_diff (mpz_ptr res, mpz_poly_srcptr f, mpz_srcptr x);
void mpz_poly_eval_poly(mpz_poly_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr x);
void mpz_poly_eval_diff_poly (mpz_poly_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr x);
void mpz_poly_eval_mod_mpz(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x, mpz_srcptr m);
int mpz_poly_is_root(mpz_poly_srcptr poly, mpz_srcptr root, mpz_srcptr modulus);
void mpz_poly_eval_several_mod_mpz(mpz_ptr * res, mpz_poly_srcptr * f, int k, mpz_srcptr x, mpz_srcptr m);
void mpz_poly_sqr_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invf, mpz_srcptr invm);
void mpz_poly_pow_ui(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n);
void mpz_poly_pow_ui_mod_f(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n, mpz_poly_srcptr f);
void mpz_poly_pow_mod_f_mod_ui(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, unsigned long p);
void mpz_poly_pow_mod_f_mod_mpz(mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, mpz_srcptr a, mpz_srcptr p);
void mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f, unsigned long a, mpz_srcptr p);
void mpz_poly_derivative(mpz_poly_ptr df, mpz_poly_srcptr f);
mpz_poly* mpz_poly_base_modp_init (mpz_poly_srcptr P0, unsigned long p, unsigned long *K, int l);
void mpz_poly_base_modp_clear (mpz_poly *P, int l);
void mpz_poly_base_modp_lift(mpz_poly_ptr a, mpz_poly *P, int k, mpz_srcptr pk);
size_t mpz_poly_sizeinbase (mpz_poly_srcptr f, int base);
size_t mpz_poly_size (mpz_poly_srcptr f);
void mpz_poly_infinity_norm(mpz_ptr in, mpz_poly_srcptr f);
size_t mpz_poly_totalsize (mpz_poly_srcptr f);
void mpz_poly_gcd_mpz (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_srcptr p);
// compute f = GCD(f,g) mod N. If this fails, put the factor in the last
// given argument.
int mpz_poly_pseudogcd_mpz(mpz_poly_ptr , mpz_poly_ptr , mpz_srcptr , mpz_ptr);
/* lc(b)*(deg(a)-deg(b)+1) = q*b+r */
void mpz_poly_pseudo_division(mpz_poly_ptr q, mpz_poly_ptr r,
    mpz_poly_srcptr a, mpz_poly_srcptr b);
/* lc(b)*(deg(a)-deg(b)+1) = q*b+r */
void mpz_poly_pseudo_remainder(mpz_poly_ptr r,
    mpz_poly_srcptr a, mpz_poly_srcptr b);
void mpz_poly_xgcd_mpz(mpz_poly_ptr gcd, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_poly_ptr u, mpz_poly_ptr v, mpz_srcptr p);
void mpz_poly_homography (mpz_poly_ptr Fij, mpz_poly_srcptr F, int64_t H[4]);
void mpz_poly_homogeneous_eval_siui (mpz_ptr v, mpz_poly_srcptr f, int64_t i, uint64_t j);
void mpz_poly_content (mpz_ptr c, mpz_poly_srcptr F);
int mpz_poly_has_trivial_content (mpz_poly_srcptr F);
int mpz_poly_divide_by_content (mpz_poly_ptr F);
void mpz_poly_resultant(mpz_ptr res, mpz_poly_srcptr p, mpz_poly_srcptr q);
void mpz_poly_discriminant(mpz_ptr res, mpz_poly_srcptr f);
int mpz_poly_squarefree_p(mpz_poly_srcptr f);
int mpz_poly_is_irreducible_z(mpz_poly_srcptr f);

int mpz_poly_number_of_real_roots(mpz_poly_srcptr f);
void mpz_poly_discriminant_of_linear_combination (mpz_poly_ptr D, mpz_poly_srcptr f0, mpz_poly_srcptr g);

struct mpz_poly_with_m_s {
    mpz_poly f;
    int m;
};
typedef struct mpz_poly_with_m_s mpz_poly_with_m[1];
typedef struct mpz_poly_with_m_s * mpz_poly_with_m_ptr;
typedef const struct mpz_poly_with_m_s * mpz_poly_with_m_srcptr;

struct mpz_poly_factor_list_s {
    mpz_poly_with_m * factors;
    int alloc;
    int size;
};
typedef struct mpz_poly_factor_list_s mpz_poly_factor_list[1];
typedef struct mpz_poly_factor_list_s * mpz_poly_factor_list_ptr;
typedef const struct mpz_poly_factor_list_s * mpz_poly_factor_list_srcptr;

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l);
void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly_srcptr f, int m);
void mpz_poly_factor_list_fprintf(FILE* fp, mpz_poly_factor_list_srcptr l);
void mpz_poly_factor_list_accumulate(mpz_poly_ptr f, mpz_poly_factor_list_srcptr l);
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, mpz_srcptr p);
int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0, mpz_srcptr p);
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int k, mpz_srcptr p, gmp_randstate_t rstate);

/* output is sorted by degree and lexicographically */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_t rstate);

int mpz_poly_is_irreducible(mpz_poly_srcptr f, mpz_srcptr p);

/* This computes the ell-adic lifts of the factors of f, assuming
 * we have no multiplicities, using Newton lifting.
 * This requires that f be monic 
 *
 * The output is sorted based on the order of the factors mod p (that is,
 * factors are the lifts of the factors returned by mpz_poly_factor mod
 * p, in the same order).
 */
int mpz_poly_factor_and_lift_padically(mpz_poly_factor_list_ptr fac, mpz_poly_srcptr f, mpz_srcptr ell, int prec, gmp_randstate_t rstate);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* This is sort of a generic way to write a c++ equivalent to the C type.
 * The first-class citizen in the cado-nfs code is (still) the C type, so
 * we're definitely bound to have a few infelicities here:
 *  - the type name can't be the same because of the size-1 array trick
 *    in C.
 *  - the C type is embedded as a member x for the same reason.
 *  - most operations on the C type should go through the member x
 *    (however, the conversions we have to _ptr and _srcptr can ease
 *    things a bit).
 */
struct cxx_mpz_poly {
    mpz_poly x;

    ATTRIBUTE_NODISCARD
    inline int degree() const { return x->deg; } /* handy */

    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    cxx_mpz_poly() { mpz_poly_init(x, -1); }

    cxx_mpz_poly(mpz_poly_srcptr f) { mpz_poly_init(x, -1); mpz_poly_set(x, f); }
    cxx_mpz_poly(cxx_mpz_poly const & o) {
        mpz_poly_init(x, -1);
        mpz_poly_set(x, o.x);
    }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    cxx_mpz_poly (const T & rhs) {
        mpz_poly_init(x, -1);
        *this = rhs;
    }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    cxx_mpz_poly & operator=(const T a) {
        cxx_mpz b;
        gmp_auxx::mpz_set(b, a);
        mpz_poly_set_mpz(x, b);
        return *this;
    }
    cxx_mpz_poly (const cxx_mpz & rhs) {
        mpz_poly_init(x, -1);
        *this = rhs;
    }
    cxx_mpz_poly & operator=(cxx_mpz const & a) {
        mpz_poly_set_mpz(x, a);
        return *this;
    }
#if __cplusplus >= 201103L
    cxx_mpz_poly(cxx_mpz_poly && o) {
        mpz_poly_init(x, -1);
        mpz_poly_swap(x, o.x);
    }
    cxx_mpz_poly& operator=(cxx_mpz_poly && o) {
        mpz_poly_swap(x, o.x);
        return *this;
    }
#endif
    cxx_mpz_poly(std::string const & e) : cxx_mpz_poly() {
        mpz_poly_set_from_expression(x, e.c_str());
    }
    cxx_mpz_poly& operator=(std::string const & e) {
        mpz_poly_set_from_expression(x, e.c_str());
        return *this;
    }
    cxx_mpz_poly(const char * e) : cxx_mpz_poly() {
        mpz_poly_set_from_expression(x, e);
    }
    cxx_mpz_poly& operator=(const char * e) {
        mpz_poly_set_from_expression(x, e);
        return *this;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

    ~cxx_mpz_poly() { mpz_poly_clear(x); }

    cxx_mpz_poly & operator=(cxx_mpz_poly const & o) {
        mpz_poly_set(x, o.x);
        return *this;
    }
    operator mpz_poly_ptr() { return x; }
    operator mpz_poly_srcptr() const { return x; }
    mpz_poly_ptr operator->() { return x; }
    mpz_poly_srcptr operator->() const { return x; }
    std::string print_poly(std::string const& var) const;

    template<typename T>
    class named_proxy {
        static_assert(std::is_reference<T>::value, "T must be a reference");
        typedef typename std::remove_reference<T>::type V;
        typedef typename std::remove_const<V>::type Vnc;
        typedef named_proxy<Vnc &> nc;
        static constexpr const bool is_c = std::is_const<V>::value;
        public:
        T c;
        std::string x;
        named_proxy(T c, std::string const & x)
            : c(c), x(x)
        {}
        template<
                typename U = T,
                typename = typename std::enable_if<
                    std::is_same<U, nc>::value
                >::type
            >
        named_proxy(U const & c) : c(c.c), x(c.x) {}
    };

    named_proxy<cxx_mpz_poly &> named(std::string const & x) {
        return named_proxy<cxx_mpz_poly &>(*this, x);
    }
    named_proxy<cxx_mpz_poly const &> named(std::string const & x) const {
        return named_proxy<cxx_mpz_poly const &>(*this, x);
    }

    /* A few initializers and convenience functions */
    cxx_mpz_poly& operator=(mpz_srcptr a)
    {
        mpz_poly_realloc(x, 1);
        mpz_poly_setcoeff(x, 0, a);
        mpz_poly_cleandeg(x, 0);
        return *this;
    }
    inline cxx_mpz_poly(mpz_srcptr c) : cxx_mpz_poly() { *this = c; }
    bool operator==(mpz_srcptr a) const {
        if (mpz_cmp_ui(a, 0) == 0) return x->deg == -1;
        return x->deg == 0 && mpz_cmp(mpz_poly_coeff_const(x, 0), a) == 0;
    }
    inline bool operator==(mpz_poly_srcptr a) const { return mpz_poly_cmp(x, a) == 0; }
    inline bool operator!=(mpz_poly_srcptr a) const { return mpz_poly_cmp(x, a) != 0; }
    inline bool operator<(mpz_poly_srcptr a) const { return mpz_poly_cmp(x, a) < 0; }
    /* we need to add these explicitly in order to resolve ambiguities */
    inline bool operator==(cxx_mpz_poly const & a) const { return mpz_poly_cmp(x, a) == 0; }
    inline bool operator!=(cxx_mpz_poly const & a) const { return mpz_poly_cmp(x, a) != 0; }
    inline bool operator<(cxx_mpz_poly const & a) const { return mpz_poly_cmp(x, a) < 0; }


    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    bool operator==(T a) const {
        if (a == 0) return x->deg == -1;
        return x->deg == 0 && gmp_auxx::mpz_cmp(mpz_poly_coeff_const(x, 0), a) == 0;
    }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0 >
    inline bool operator!=(T a) const { return !((*this) == a); }

    // mpz_ptr operator[](unsigned int i) { return mpz_poly_coeff_const(x, i); }
    mpz_srcptr
    coeff(unsigned int i) const { return mpz_poly_coeff_const(x, i); }
    
    /* We don't want to use it. I got "ISO C++ says that these are
     * ambiguous, even though the worst conversion for the first is
     * better than the worst conversion for the second:", because
     * operator[] can also be first a conversion to mpz_poly_srcptr, and
     * then []. Not taking sides here, but better avoid the issue. => use
     * .coeff() everywhere.
     */
    mpz_srcptr
    operator[](unsigned int i) const ATTRIBUTE_DEPRECATED { return this->coeff(i); }
};



/* printing needs a way to specify the variables... */
std::ostream& operator<<(std::ostream& o, cxx_mpz_poly::named_proxy<cxx_mpz_poly const &> f);
std::istream& operator>>(std::istream& in, cxx_mpz_poly::named_proxy<cxx_mpz_poly &> f);

/* we do have a default behaviour, though */
inline std::ostream& operator<<(std::ostream& o, cxx_mpz_poly const & f)
{
    return o << f.named("x");
}

inline std::istream& operator>>(std::istream& in, cxx_mpz_poly & f)
{
    return in >> f.named("x");
}

/* If there is an integer that takes the value evaluations[i] at
 * points[i] for all i, store it to resultant and return 1. Otherwise 0
 * is return, and the resultant value is clobbered.
 */
int mpz_poly_interpolate(mpz_poly_ptr resultant,
        std::vector<cxx_mpz> const & points,
        std::vector<cxx_mpz> const & evaluations);

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpz_poly_init(cxx_mpz_poly & pl, int) __attribute__((error("mpz_poly_init must not be called on a mpz_poly reference -- it is the caller's business (via a ctor)")));
extern void mpz_poly_clear(cxx_mpz_poly & pl) __attribute__((error("mpz_poly_clear must not be called on a mpz_poly reference -- it is the caller's business (via a dtor)")));
#endif

struct mpz_poly_coeff_list {
    cxx_mpz_poly const & P;
    std::string sep;
    mpz_poly_coeff_list(cxx_mpz_poly const & P, std::string const & sep = ", "): P(P), sep(sep) {}
};
std::ostream& operator<<(std::ostream& os, mpz_poly_coeff_list const & P);

namespace fmt {
    template <> struct formatter<cxx_mpz_poly>: ostream_formatter {};
    template <> struct formatter<mpz_poly_coeff_list>: ostream_formatter {};
}


cxx_mpz_poly prod(std::vector<std::pair<cxx_mpz_poly, int>> const &lf, mpz_srcptr modulus = NULL, mpz_srcptr invm = NULL);
std::vector<std::pair<cxx_mpz_poly, int>> mpz_poly_factor(mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_t rstate);
std::vector<std::pair<cxx_mpz_poly, int>> 
mpz_poly_factor_and_lift_padically(mpz_poly_srcptr f, mpz_srcptr ell, int prec, gmp_randstate_t rstate);

/* returns 1 if parsing was successful */
inline int mpz_poly_set_from_expression(mpz_poly_ptr f, std::string const & value)
{
    return mpz_poly_set_from_expression(f, value.c_str());
}
#endif

#endif	/* MPZ_POLY_H_ */
