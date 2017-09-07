#ifndef DOUBLE_POLY_H_
#define DOUBLE_POLY_H_

#include <stdio.h>
#include <limits.h>

#ifdef DOUBLE_POLY_EXPOSE_COMPLEX_FUNCTIONS
#include <complex.h>
#endif

/* forward-declare our type before inclusion by mpz_poly.h, since we
 * include eachother
 */
struct double_poly_s;
typedef struct double_poly_s * double_poly_ptr;
typedef const struct double_poly_s * double_poly_srcptr;

#include "mpz_poly.h"

/* floating point polynomials */

#ifdef __cplusplus
extern "C" {
#endif

struct double_poly_s {
    int alloc;
    int deg;
    double *coeff;
};

typedef struct double_poly_s double_poly[1];

/* double_poly.c */
void double_poly_init (double_poly_ptr, int deg);
void double_poly_realloc (double_poly_ptr, int nc);
void double_poly_clear (double_poly_ptr);
void double_poly_swap(double_poly_ptr p, double_poly_ptr q);
void double_poly_set (double_poly_ptr, double_poly_srcptr);

void double_poly_set_zero (double_poly_ptr r);
void double_poly_set_xi(double_poly_ptr s, int i);

void double_poly_cleandeg(double_poly_ptr f, int deg);

int double_poly_cmp(double_poly_srcptr a, double_poly_srcptr b);

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
void double_poly_print (FILE *, double_poly_srcptr, char *variable_name);
int double_poly_asprint (char **t, double_poly_srcptr p, char *variable_name);
void double_poly_set_mpz_poly (double_poly_ptr p, mpz_poly_srcptr q);

double double_poly_resultant(double_poly_srcptr p, double_poly_srcptr q);
void double_poly_mul_double(double_poly_ptr f, double_poly_srcptr g,
    double mul);
double double_poly_div_linear(double_poly_ptr q, double_poly_srcptr p, const double r);
void double_poly_set_string(double_poly_ptr poly, const char *str);

#ifdef DOUBLE_POLY_EXPOSE_COMPLEX_FUNCTIONS
/* we use _Complex and not complex below, because it's mandatory in order
 * to be nice to C++ code */
void double_poly_complex_roots(double _Complex *roots, double_poly_srcptr);
void double_poly_complex_roots_long(long double _Complex *roots, double_poly_srcptr);
/* The implementation of these in in polyroots.c -- we expose them here,
 * but new implementation should prefer using the functions above.
 */
uint32_t poly_roots_double(double *poly, uint32_t degree, double _Complex *roots);
uint32_t poly_roots_longdouble(double *poly, uint32_t degree, long double _Complex *roots);
#endif

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
#if __cplusplus >= 201103L
    cxx_double_poly(cxx_double_poly && o) {
        double_poly_init(x, -1);
        double_poly_swap(x, o.x);
    }
    cxx_double_poly& operator=(cxx_double_poly && o) {
        double_poly_swap(x, o.x);
        return *this;
    }
#endif
    operator double_poly_ptr() { return x; }
    operator double_poly_srcptr() const { return x; }
    double_poly_ptr operator->() { return x; }
    double_poly_srcptr operator->() const { return x; }
    std::string print_poly(std::string const& var) const;
};

#endif


#endif	/* DOUBLE_POLY_H_ */
