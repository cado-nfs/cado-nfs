#ifndef CADO_MPZ_MAT_H
#define CADO_MPZ_MAT_H

#include <stdio.h>
#ifdef __cplusplus
#include <iosfwd>      // for ostream
#endif

#include <gmp.h>
#include "macros.h"
#include "mpz_poly.h"

#ifdef __cplusplus
#include "fmt_helper_sagemath.hpp"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* typedefs */

// mpz_t : multiple precision Z-elements ; same for mpq_t
// a matric on mpz_t type
struct mpz_mat_s {
    mpz_t *x;
    unsigned int m;     // number of rows
    unsigned int n;     // number of columns
};
typedef struct mpz_mat_s mpz_mat[1];
typedef struct mpz_mat_s * mpz_mat_ptr;
typedef const struct mpz_mat_s * mpz_mat_srcptr;

// same on mpq_t type
struct mpq_mat_s {
    mpq_t *x;
    unsigned int m;     // number of rows
    unsigned int n;     // number of columns
};
typedef struct mpq_mat_s mpq_mat[1];
typedef struct mpq_mat_s * mpq_mat_ptr;
typedef const struct mpq_mat_s * mpq_mat_srcptr;

/* entry access*/
mpz_ptr mpz_mat_entry(mpz_mat_ptr M, unsigned int i, unsigned int j);
mpz_srcptr mpz_mat_entry_const(mpz_mat_srcptr M, unsigned int i, unsigned int j);
mpq_ptr mpq_mat_entry(mpq_mat_ptr M, unsigned int i, unsigned int j);
mpq_srcptr mpq_mat_entry_const(mpq_mat_srcptr M, unsigned int i, unsigned int j);

/* init/clear/realloc*/
void mpz_mat_init(mpz_mat_ptr M, unsigned int m, unsigned int n);
void mpz_mat_clear(mpz_mat_ptr M);
void mpz_mat_realloc(mpz_mat_ptr M, unsigned int m, unsigned int n);
void mpq_mat_init(mpq_mat_ptr M, unsigned int m, unsigned int n);
void mpq_mat_clear(mpq_mat_ptr M);
void mpq_mat_realloc(mpq_mat_ptr M, unsigned int m, unsigned int n);

/* operations on submatrices, and swaps*/
void mpz_mat_swap(mpz_mat_ptr A, mpz_mat_ptr B);
void mpq_mat_swap(mpq_mat_ptr A, mpq_mat_ptr B);
void mpz_mat_submat_swap(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);
void mpq_mat_submat_swap(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);
void mpz_mat_submat_set(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_srcptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);
void mpq_mat_submat_set(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_srcptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);

/* set/set_ui*/
void mpz_mat_set(mpz_mat_ptr dst, mpz_mat_srcptr src);
void mpz_mat_set_ui(mpz_mat_ptr M, unsigned long a);
void mpz_mat_set_mpz(mpz_mat_ptr M, mpz_srcptr a);
void mpz_mat_add_ui(mpz_mat_ptr M, unsigned long a);
void mpz_mat_add_mpz(mpz_mat_ptr M, mpz_srcptr a);
void mpq_mat_set_ui(mpq_mat_ptr M, unsigned long a);
void mpq_mat_set(mpq_mat_ptr dst, mpq_mat_srcptr src);
void mpz_mat_urandomm(mpz_mat_ptr M, gmp_randstate_t state, mpz_srcptr p);
void mpq_mat_urandomm(mpq_mat_ptr M, gmp_randstate_t state, mpz_srcptr p);

/* Joining matrices */
void mpz_mat_vertical_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2);
void mpq_mat_vertical_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2);
void mpz_mat_horizontal_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2);
void mpq_mat_horizontal_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2);

/* determinant and trace */
// We assume that M is square
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M);
// We assume that M is triangular (and square)
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M);
/* this one does the HNF first */
void mpz_mat_determinant(mpz_ptr d, mpz_mat_srcptr M);
// We assume that M is square
void mpq_mat_trace(mpq_ptr t, mpq_mat_srcptr M);
// We assume that M is triangular (and square)
void mpq_mat_determinant_triangular(mpq_ptr d, mpq_mat_srcptr M);
/* uses HNF first */
void mpq_mat_determinant(mpq_ptr d, mpq_mat_srcptr M);
void mpz_mat_transpose(mpz_mat_ptr D, mpz_mat_srcptr M);
void mpq_mat_transpose(mpq_mat_ptr D, mpq_mat_srcptr M);
void mpz_mat_reverse_rows(mpz_mat_ptr B, mpz_mat_srcptr A);
void mpz_mat_reverse_columns(mpz_mat_ptr B, mpz_mat_srcptr A);

/* miscellaneous */
int mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M);
void mpq_mat_set_mpz_mat(mpq_mat_ptr N, mpz_mat_srcptr M);
void mpq_mat_set_mpz_mat_denom(mpq_mat_ptr N, mpz_mat_srcptr M, mpz_srcptr d);
int mpz_mat_p_valuation(mpz_mat_srcptr A, mpz_srcptr p);
int mpz_mat_p_valuation_ui(mpz_mat_srcptr A, unsigned long p);

void mpz_mat_mod_ui(mpz_mat_ptr dst, mpz_mat_srcptr src, unsigned long p);
void mpz_mat_mod_mpz(mpz_mat_ptr dst, mpz_mat_srcptr src, mpz_srcptr p);

/*  row-level operations */
// Return 1 if the k-th line of M is null, 0 otherwise
int mpz_mat_isnull_row(mpz_mat_srcptr M, unsigned int k);
void mpz_mat_swaprows(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
void mpq_mat_swaprows(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* add lambda times row i1 to row i0 */
void mpz_mat_addmulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda);
/* add lambda times row i1 to row i0 */
void mpq_mat_addmulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda);
/* subtract lambda times row i1 to row i0 */
void mpz_mat_submulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda);
/* subtract lambda times row i1 to row i0 */
void mpq_mat_submulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda);
/* add row i1 to row i0 */
void mpz_mat_addrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
/* add row i1 to row i0 */
void mpq_mat_addrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* subtract row i1 to row i0 */
void mpz_mat_subrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
/* subtract row i1 to row i0 */
void mpq_mat_subrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* multiply row i0 by lambda */
void mpz_mat_mulrow(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda);
/* multiply row i0 by lambda */
void mpq_mat_mulrow(mpq_mat_ptr M, unsigned int i0, mpq_srcptr lambda);


/* multiplication */
void mpz_mat_mul(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B);
void mpz_mat_mul_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B, unsigned long p);
void mpz_mat_mul_mod_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B, mpz_srcptr p);
void mpz_mat_mul_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_ptr k);
void mpz_mat_divexact_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, mpz_srcptr k);
void mpz_mat_divexact_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long k);
void mpz_mat_mul_si(mpz_mat_ptr C, mpz_mat_srcptr A, long k);
void mpz_mat_mul_ui(mpz_mat_ptr C, mpz_mat_srcptr A, unsigned long k);
// Returns A^n for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_pow_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n);
// Returns A^n mod p for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_pow_ui_mod_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n, mpz_srcptr p);
void mpz_mat_pow_ui_mod_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n, unsigned long p);
void mpq_mat_mul(mpq_mat_ptr D, mpq_mat_srcptr A, mpq_mat_srcptr B);
void mpq_mat_mul_mpq(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_srcptr k);
void mpq_mat_mul_si(mpq_mat_ptr C, mpq_mat_srcptr A, long k);
void mpq_mat_mul_ui(mpq_mat_ptr C, mpq_mat_srcptr A, unsigned long k);
void mpq_mat_mul_mpz(mpq_mat_ptr C, mpq_mat_srcptr A, mpz_srcptr k);
void mpq_mat_div_mpq(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_srcptr k);
void mpq_mat_div_si(mpq_mat_ptr C, mpq_mat_srcptr A, long k);
void mpq_mat_div_ui(mpq_mat_ptr C, mpq_mat_srcptr A, unsigned long k);
void mpq_mat_div_mpz(mpq_mat_ptr C, mpq_mat_srcptr A, mpz_srcptr k);
void mpz_mat_add(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B);
void mpq_mat_add(mpq_mat_ptr D, mpq_mat_srcptr A, mpq_mat_srcptr B);
void mpz_mat_sub(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B);
void mpq_mat_sub(mpq_mat_ptr D, mpq_mat_srcptr A, mpq_mat_srcptr B);
void mpz_poly_eval_mpz_mat(mpz_mat_ptr D, mpz_poly_srcptr f, mpz_mat_srcptr M);
void mpz_poly_eval_mpz_mat_mod_ui(mpz_mat_ptr D, mpz_poly_srcptr f, mpz_mat_srcptr M, unsigned long p);
void mpz_poly_eval_mpz_mat_mod_mpz(mpz_mat_ptr D, mpz_poly_srcptr f, mpz_mat_srcptr M, mpz_srcptr p);

/*  gaussian reduction over the rationals
 * this is a backend for row gaussian reduction. T receives the
 * transformation matrix. M is modified in place.
 * this includes reordering of the rows.
 *
 * this works reasonably well over Fp, but for rationals we suffer from
 * coefficient blowup.
 *
 * T may be NULL in case we don't care.
 */
void mpq_mat_gauss_backend(mpq_mat_ptr M, mpq_mat_ptr T);


/* Gaussian reduction over Z/pZ
 *
 * We sort of return something even in the case where p is not prime.
 *
 * M is the reduced matrix (modification in place of the input), T the
 * transformation matrix (so that T * input-M = output-M)
 */
void mpz_mat_gauss_backend_mod_mpz(mpz_mat_ptr M, mpz_mat_ptr T, mpz_srcptr p);
void mpz_mat_gauss_backend_mod_ui(mpz_mat_ptr M, mpz_mat_ptr T, unsigned long p);

/* this computes the row hnf on a column C.
 * The result is always of the form C'=(gcd(C),0,0,...,0). The transform
 * matrix T such that C'=T*C is most interesting.  a is modified in
 * place.
 */
void mpz_gcd_many(mpz_mat_ptr dT, mpz_ptr * a, unsigned int n);

/* return +1 or -1, which is the determinant of the transformation matrix */
int mpz_mat_hermite_form(mpz_mat_ptr M);
int mpz_mat_hermite_form_rev(mpz_mat_ptr M);

/* kernel*/
// This is supposed to compute the Kernel of M mod p and to store it in the matrix K. If r is the rank of M, and M is a square matrix n*n, K is a n*(n-r) matrix
void mpz_mat_kernel_mod_ui(mpz_mat_ptr K, mpz_mat_srcptr M, unsigned long p);
void mpz_mat_kernel_mod_mpz(mpz_mat_ptr K, mpz_mat_srcptr M, mpz_srcptr p);

void mpq_mat_inv(mpq_mat_ptr dst, mpq_mat_srcptr src);

int mpz_mat_cmp(mpz_mat_srcptr M, mpz_mat_srcptr N);
int mpz_mat_is_zero(mpz_mat_srcptr M);
int mpq_mat_cmp(mpq_mat_srcptr M, mpq_mat_srcptr N);
void mpz_mat_fprint(FILE * stream, mpz_mat_srcptr M);
void mpq_mat_fprint(FILE * stream, mpq_mat_srcptr M);

void mpz_mat_LLL(mpz_ptr det, mpz_mat_ptr M, mpz_mat_ptr U, mpz_srcptr a,
    mpz_srcptr b);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

struct cxx_mpz_mat {
    using cxx_coeff_type = cxx_mpz;
    mpz_mat x;
    cxx_mpz_mat() { mpz_mat_init(x, 0, 0); }
    cxx_mpz_mat(unsigned int m, unsigned int n) { mpz_mat_init(x, m, n); }
    cxx_mpz_mat(unsigned int m, unsigned int n, unsigned long z) {
        mpz_mat_init(x, m, n);
        mpz_mat_set_ui(x, z);
    }
    cxx_mpz_mat(unsigned int m, unsigned int n, mpz_srcptr z) {
        mpz_mat_init(x, m, n);
        mpz_mat_set_mpz(x, z);
    }
    ~cxx_mpz_mat() { mpz_mat_clear(x); }
    cxx_mpz_mat(cxx_mpz_mat const & o) {
        mpz_mat_init(x, 0, 0);
        mpz_mat_set(x, o.x);
    }
    cxx_mpz_mat & operator=(cxx_mpz_mat const & o) {
        if (this != &o)
            mpz_mat_set(x, o.x);
        return *this;
    }
    cxx_mpz_mat(cxx_mpz_mat && o) {
        mpz_mat_init(x, 0, 0);
        mpz_mat_swap(x, o.x);
    }
    cxx_mpz_mat& operator=(cxx_mpz_mat && o) noexcept {
        mpz_mat_swap(x, o.x);
        return *this;
    }
    bool operator==(cxx_mpz_mat const & o) const {
        return mpz_mat_cmp(x, o.x) == 0;
    }
    bool operator!=(cxx_mpz_mat const & o) const {
        return mpz_mat_cmp(x, o.x) != 0;
    }
    unsigned int nrows() const { return x->m; }
    unsigned int ncols() const { return x->n; }
    operator mpz_mat_ptr() { return x; }
    operator mpz_mat_srcptr() const { return x; }
    mpz_mat_ptr operator->() { return x; }
    mpz_mat_srcptr operator->() const { return x; }
    mpz_ptr operator ()(unsigned int i, unsigned int j) { return mpz_mat_entry(x, i, j); }
    mpz_srcptr operator ()(unsigned int i, unsigned int j) const { return mpz_mat_entry_const(x, i, j); }
};
#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpz_mat_init(cxx_mpz_mat & pl, unsigned int, unsigned int) __attribute__((error("mpz_mat_init must not be called on a mpz_mat reference -- it is the caller's business (via a ctor)")));
extern void mpz_mat_clear(cxx_mpz_mat & pl) __attribute__((error("mpz_mat_clear must not be called on a mpz_mat reference -- it is the caller's business (via a dtor)")));
#endif
struct cxx_mpq_mat {
    using cxx_coeff_type = cxx_mpq;
    mpq_mat x;
    cxx_mpq_mat() { mpq_mat_init(x, 0, 0); }
    cxx_mpq_mat(unsigned int m, unsigned int n) { mpq_mat_init(x, m, n); }
    cxx_mpq_mat(unsigned int m, unsigned int n, unsigned long z) {
        mpq_mat_init(x, m, n);
        mpq_mat_set_ui(x, z);
    }
    ~cxx_mpq_mat() { mpq_mat_clear(x); }
    cxx_mpq_mat(cxx_mpq_mat const & o) {
        mpq_mat_init(x, 0, 0);
        mpq_mat_set(x, o.x);
    }
    cxx_mpq_mat & operator=(cxx_mpq_mat const & o) {
        mpq_mat_set(x, o.x);
        return *this;
    }
    cxx_mpq_mat(cxx_mpq_mat && o) {
        mpq_mat_init(x, 0, 0);
        mpq_mat_swap(x, o.x);
    }
    cxx_mpq_mat& operator=(cxx_mpq_mat && o) {
        mpq_mat_swap(x, o.x);
        return *this;
    }
    /* As a bonus, allow conversion and assignment from mpz matrices,
     * too. Move ctor and move assignment do not make sense in that case,
     * of course. */
    cxx_mpq_mat(cxx_mpz_mat const & o) {
        mpq_mat_init(x, 0, 0);
        mpq_mat_set_mpz_mat(x, o.x);
    }
    cxx_mpq_mat & operator=(cxx_mpz_mat const & o) {
        mpq_mat_set_mpz_mat(x, o.x);
        return *this;
    }
    bool operator==(cxx_mpq_mat const & o) const {
        return mpq_mat_cmp(x, o.x) == 0;
    }
    bool operator!=(cxx_mpq_mat const & o) const {
        return mpq_mat_cmp(x, o.x) != 0;
    }
    unsigned int nrows() const { return x->m; }
    unsigned int ncols() const { return x->n; }
    operator mpq_mat_ptr() { return x; }
    operator mpq_mat_srcptr() const { return x; }
    mpq_mat_ptr operator->() { return x; }
    mpq_mat_srcptr operator->() const { return x; }
    mpq_ptr operator ()(unsigned int i, unsigned int j) { return mpq_mat_entry(x, i, j); }
    mpq_srcptr operator ()(unsigned int i, unsigned int j) const { return mpq_mat_entry_const(x, i, j); }
};

extern std::ostream& operator<<(std::ostream& os, cxx_mpz_mat const& M);
extern std::ostream& operator<<(std::ostream& os, cxx_mpq_mat const& M);
#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpq_mat_init(cxx_mpq_mat & pl, unsigned int, unsigned int) __attribute__((error("mpq_mat_init must not be called on a mpq_mat reference -- it is the caller's business (via a ctor)")));
extern void mpq_mat_clear(cxx_mpq_mat & pl) __attribute__((error("mpq_mat_clear must not be called on a mpq_mat reference -- it is the caller's business (via a dtor)")));
#endif

namespace fmt {
    template <>
    struct formatter<cxx_mpz_mat> : fmt_helper_sagemath<cxx_mpz_mat>
    {
        static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        auto format(cxx_mpz_mat const & e, format_context& ctx) const -> format_context::iterator;
    };
    template<>
    struct formatter<cxx_mpq_mat> : fmt_helper_sagemath<cxx_mpq_mat>
    {
        static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        auto format(cxx_mpq_mat const & e, format_context& ctx) const -> format_context::iterator;
    };
} /* namespace fmt */
#endif

#endif	/* CADO_MPZ_MAT_H */
