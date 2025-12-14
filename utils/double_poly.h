#ifndef CADO_DOUBLE_POLY_H
#define CADO_DOUBLE_POLY_H

// IWYU pragma: no_include "mpz_poly.h"
// the fwd-decl is enough

#include <stdint.h>    // for uint32_t
#include <stdio.h>     // for FILE

#ifdef __cplusplus
#include <string>      // for string
#include <ostream>
#include <type_traits>
#include <utility>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "named_proxy.hpp"
#endif

#include "macros.h"    // for GNUC_VERSION_ATLEAST

#ifndef CADO_MPZ_POLY_H
typedef struct mpz_poly_s * mpz_poly_ptr;
typedef const struct mpz_poly_s * mpz_poly_srcptr;
#endif

/* floating point polynomials */

#ifdef __cplusplus
extern "C" {
#endif

struct double_poly_s {
    unsigned int alloc;
    int deg;
    double *coeff;
};

typedef struct double_poly_s double_poly[1];
#ifndef CADO_MPZ_POLY_H
/* mpz_poly.h forward-declares these. Don't do it twice */
typedef struct double_poly_s * double_poly_ptr;
typedef const struct double_poly_s * double_poly_srcptr;
#endif

/* double_poly.c */
void double_poly_init (double_poly_ptr, int deg);
void double_poly_realloc (double_poly_ptr, unsigned int nc);
void double_poly_clear (double_poly_ptr);
void double_poly_swap(double_poly_ptr p, double_poly_ptr q);
void double_poly_set (double_poly_ptr, double_poly_srcptr);

void double_poly_set_zero (double_poly_ptr r);
void double_poly_set_xi(double_poly_ptr s, int i);

void double_poly_cleandeg(double_poly_ptr f, int deg);

int double_poly_cmp(double_poly_srcptr a, double_poly_srcptr b);
double double_poly_lc(double_poly_srcptr f);

double double_poly_eval (double_poly_srcptr, double);
double double_poly_eval_homogeneous (double_poly_srcptr p, double x, double y);
double double_poly_eval_safe (double_poly_srcptr, double);
double double_poly_dichotomy (double_poly_srcptr, double, double, double,
                              unsigned int);
void double_poly_derivative(double_poly_ptr, double_poly_srcptr);
void double_poly_mul(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_add(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_sub(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_addmul(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_submul(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_addmul_double(double_poly_ptr, double_poly_srcptr, double);
void double_poly_submul_double(double_poly_ptr, double_poly_srcptr, double);
void double_poly_neg(double_poly_ptr, double_poly_srcptr);
void double_poly_revert(double_poly_ptr, double_poly_srcptr);
double double_poly_bound_roots (double_poly_srcptr p);
unsigned int double_poly_compute_roots(double *, double_poly_srcptr, double);
unsigned int double_poly_compute_all_roots_with_bound(double *,
                                                      double_poly_srcptr,
                                                      double);
unsigned int double_poly_compute_all_roots(double *, double_poly_srcptr);
void double_poly_print (FILE *, double_poly_srcptr, const char *variable_name);
int double_poly_asprint (char **t, double_poly_srcptr p, const char *variable_name);
void double_poly_set_mpz_poly (double_poly_ptr p, mpz_poly_srcptr q);

void double_poly_mul_double(double_poly_ptr f, double_poly_srcptr g,
    double mul);
double double_poly_div_linear(double_poly_ptr q, double_poly_srcptr p, double r);
void double_poly_set_string(double_poly_ptr poly, const char *str);

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
struct cxx_double_poly {
    double_poly x;
    static constexpr int number_of_variables = 1;
    cxx_double_poly(int deg = -1) { double_poly_init(x, deg); }
    cxx_double_poly(double_poly_srcptr f) { double_poly_init(x, -1); double_poly_set(x, f); }
    ~cxx_double_poly() { double_poly_clear(x); }
    cxx_double_poly(cxx_double_poly const & o) {
        double_poly_init(x, -1);
        double_poly_set(x, o.x);
    }
    cxx_double_poly & operator=(cxx_double_poly const & o) {
        double_poly_set(x, o.x);
        return *this;
    }
    cxx_double_poly(cxx_double_poly && o) {
        double_poly_init(x, -1);
        double_poly_swap(x, o.x);
    }
    cxx_double_poly& operator=(cxx_double_poly && o) {
        double_poly_swap(x, o.x);
        return *this;
    }
    operator double_poly_ptr() { return x; }
    operator double_poly_srcptr() const { return x; }
    double_poly_ptr operator->() { return x; }
    double_poly_srcptr operator->() const { return x; }
    std::string print_poly(std::string const& var) const;

    cado::named_proxy<cxx_double_poly &> named(std::string const & x) {
        return { *this, x };
    }
    cado::named_proxy<cxx_double_poly const &> named(std::string const & x) const {
        return { *this, x };
    }

};

/* printing needs a way to specify the variables... */
inline std::ostream& operator<<(std::ostream& o, cado::named_proxy<cxx_double_poly const &> const & f)
{
    return o << f.c.print_poly(f.x());
}

/* we do have a default behaviour, though */
inline std::ostream& operator<<(std::ostream& o, cxx_double_poly const & f)
{
    return o << f.named("x");
}

namespace fmt {
    template <> struct formatter<cxx_double_poly>: ostream_formatter {};
}

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void double_poly_init(cxx_double_poly & pl, int) __attribute__((error("double_poly_init must not be called on a double_poly reference -- it is the caller's business (via a ctor)")));
extern void double_poly_clear(cxx_double_poly & pl) __attribute__((error("double_poly_clear must not be called on a double_poly reference -- it is the caller's business (via a dtor)")));
#endif

#endif


#endif	/* CADO_DOUBLE_POLY_H */
