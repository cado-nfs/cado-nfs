#ifndef MPFQ_2_64_H_
#define MPFQ_2_64_H_

/* MPFQ generated file -- do not edit */

#include "gf2x.h"
#include "gf2x/gf2x-small.h"

#include "mpfq.h"
#include "mpfq_gf2n_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>

#include "assert.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      2_64

/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::gf2n::field */
/* Automatically generated code for GF(2^64) */
/* Definition polynomial P = X^64 + X^4 + X^3 + X + 1 */
/* Active handler: Mpfq::gf2n::trivialities */
/* Active handler: Mpfq::gf2n::io */
/* Active handler: Mpfq::gf2n::linearops */
/* Active handler: Mpfq::gf2n::inversion */
/* Active handler: Mpfq::gf2n::reduction */
/* Active handler: Mpfq::gf2n::mul */
/* Active handler: Mpfq::defaults::poly */
/* Options used:{
   coeffs=[ 64, 4, 3, 1, 0, ],
   helper=/tmp/mpfq-cado/gf2n/helper/helper,
   n=64,
   no_gmp=1,
   output_path=i386,
   slice=4,
   table=/tmp/mpfq-cado/gf2x/wizard.table,
   tag=2_64,
   w=32,
   } */

typedef mpfq_2_field mpfq_2_64_field;
typedef mpfq_2_dst_field mpfq_2_64_dst_field;

typedef unsigned long mpfq_2_64_elt[2];
typedef unsigned long * mpfq_2_64_dst_elt;
typedef const unsigned long * mpfq_2_64_src_elt;

typedef unsigned long mpfq_2_64_elt_ur[4];
typedef unsigned long * mpfq_2_64_dst_elt_ur;
typedef const unsigned long * mpfq_2_64_src_elt_ur;

typedef mpfq_2_64_elt * mpfq_2_64_vec;
typedef mpfq_2_64_elt * mpfq_2_64_dst_vec;
typedef mpfq_2_64_elt * mpfq_2_64_src_vec;

typedef mpfq_2_64_elt_ur * mpfq_2_64_vec_ur;
typedef mpfq_2_64_elt_ur * mpfq_2_64_dst_vec_ur;
typedef mpfq_2_64_elt_ur * mpfq_2_64_src_vec_ur;

typedef struct {
  mpfq_2_64_vec c;
  unsigned int alloc;
  unsigned int size;
} mpfq_2_64_poly_struct;
typedef mpfq_2_64_poly_struct mpfq_2_64_poly [1];
typedef mpfq_2_64_poly_struct * mpfq_2_64_dst_poly;
typedef mpfq_2_64_poly_struct * mpfq_2_64_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name */
#define mpfq_2_64_impl_name()	"2_64"
/* *Mpfq::gf2n::field::code_for_impl_max_characteristic_bits */
#define mpfq_2_64_impl_max_characteristic_bits()	2
/* *Mpfq::gf2n::field::code_for_impl_max_degree */
#define mpfq_2_64_impl_max_degree()	1

/* Functions operating on the field structure */
/* *Mpfq::gf2n::field::code_for_field_characteristic_bits */
#define mpfq_2_64_field_characteristic_bits(f)	1
/* *Mpfq::gf2n::field::code_for_field_degree */
#define mpfq_2_64_field_degree(f)	64
static inline
void mpfq_2_64_field_init(mpfq_2_64_dst_field);
/* *Mpfq::gf2n::field::code_for_field_clear */
#define mpfq_2_64_field_clear(K)	/**/
/* *Mpfq::gf2n::field::code_for_field_specify */
#define mpfq_2_64_field_specify(k, dummy, vp)	/**/
static inline
void mpfq_2_64_field_setopt(mpfq_2_64_dst_field, unsigned long, void *);

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, Mpfq::gf2n::trivialities */
#define mpfq_2_64_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, Mpfq::gf2n::trivialities */
#define mpfq_2_64_clear(f, px)	/**/

/* Elementary assignment functions */
static inline
void mpfq_2_64_set(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_set_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, unsigned long);
static inline
void mpfq_2_64_set_zero(mpfq_2_64_dst_field, mpfq_2_64_dst_elt);
static inline
unsigned long mpfq_2_64_get_ui(mpfq_2_64_dst_field, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_set_uipoly(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, unsigned long);
static inline
void mpfq_2_64_set_uipoly_wide(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, const unsigned long *, unsigned int);
static inline
unsigned long mpfq_2_64_get_uipoly(mpfq_2_64_dst_field, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_get_uipoly_wide(mpfq_2_64_dst_field, unsigned long *, mpfq_2_64_src_elt);

/* Assignment of random values */

/* Arithmetic operations on elements */
static inline
void mpfq_2_64_add(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, mpfq_2_64_src_elt);
/* *Mpfq::gf2n::trivialities::code_for_sub */
#define mpfq_2_64_sub(K, r, s1, s2)	mpfq_2_64_add(K,r,s1,s2)
/* *Mpfq::gf2n::trivialities::code_for_neg */
#define mpfq_2_64_neg(K, r, s)	mpfq_2_64_set(K,r,s)
static inline
void mpfq_2_64_mul(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_sqr(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt);
/* *Mpfq::gf2n::trivialities::code_for_is_sqr */
#define mpfq_2_64_is_sqr(f, p)	1
static inline
int mpfq_2_64_sqrt(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_pow(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, unsigned long *, size_t);
/* missing powz */
/* *Mpfq::gf2n::trivialities::code_for_frobenius */
#define mpfq_2_64_frobenius(K, r, s)	mpfq_2_64_sqr(K,r,s)
static inline
void mpfq_2_64_add_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, unsigned long);
/* *Mpfq::gf2n::trivialities::code_for_sub_ui */
#define mpfq_2_64_sub_ui(K, r, s1, s2)	mpfq_2_64_add_ui(K,r,s1,s2)
static inline
void mpfq_2_64_mul_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, unsigned long);
static inline
void mpfq_2_64_add_uipoly(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, unsigned long);
/* *Mpfq::gf2n::trivialities::code_for_sub_uipoly */
#define mpfq_2_64_sub_uipoly(K, r, s1, s2)	mpfq_2_64_add_uipoly(K,r,s1,s2)
static inline
void mpfq_2_64_mul_uipoly(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt, unsigned long);
static inline
void mpfq_2_64_longshift_left(unsigned long *, const unsigned long *, int, int);
static inline
void mpfq_2_64_longaddshift_left(unsigned long *, const unsigned long *, int, int);
static inline
int mpfq_2_64_inv(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_as_solve(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_elt);
static inline
unsigned long mpfq_2_64_trace(mpfq_2_64_dst_field, mpfq_2_64_src_elt);

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, Mpfq::gf2n::trivialities */
#define mpfq_2_64_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, Mpfq::gf2n::trivialities */
#define mpfq_2_64_elt_ur_clear(f, px)	/**/
static inline
void mpfq_2_64_elt_ur_set(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_elt_ur);
static inline
void mpfq_2_64_elt_ur_set_elt(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_elt_ur_set_zero(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur);
static inline
void mpfq_2_64_elt_ur_set_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, unsigned long);
static inline
void mpfq_2_64_elt_ur_add(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_elt_ur, mpfq_2_64_src_elt_ur);
/* *Mpfq::gf2n::trivialities::code_for_elt_ur_neg */
#define mpfq_2_64_elt_ur_neg(K, r, s)	mpfq_2_64_elt_ur_set(K,r,s)
/* *Mpfq::gf2n::trivialities::code_for_elt_ur_sub */
#define mpfq_2_64_elt_ur_sub(K, r, s1, s2)	mpfq_2_64_elt_ur_add(K,r,s1,s2)
static inline
void mpfq_2_64_mul_ur(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_elt, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_sqr_ur(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_reduce(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_dst_elt_ur);

/* Comparison functions */
static inline
int mpfq_2_64_cmp(mpfq_2_64_dst_field, mpfq_2_64_src_elt, mpfq_2_64_src_elt);
static inline
int mpfq_2_64_cmp_ui(mpfq_2_64_dst_field, mpfq_2_64_src_elt, unsigned long);
static inline
int mpfq_2_64_is_zero(mpfq_2_64_dst_field, mpfq_2_64_src_elt);

/* Input/output functions */
/* missing fprint */
/* missing print */
/* missing scan */

/* Vector functions */
void mpfq_2_64_vec_init(mpfq_2_64_dst_field, mpfq_2_64_vec *, unsigned int);
void mpfq_2_64_vec_reinit(mpfq_2_64_dst_field, mpfq_2_64_vec *, unsigned int, unsigned int);
void mpfq_2_64_vec_clear(mpfq_2_64_dst_field, mpfq_2_64_vec *, unsigned int);
static inline
void mpfq_2_64_vec_set(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_set_zero(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, unsigned int);
static inline
void mpfq_2_64_vec_setcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_elt, unsigned int);
static inline
void mpfq_2_64_vec_setcoeff_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, unsigned long, unsigned int);
static inline
void mpfq_2_64_vec_getcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_add(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_neg(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_rev(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_sub(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_scal_mul(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, mpfq_2_64_src_elt, unsigned int);
static inline
void mpfq_2_64_vec_conv(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_src_vec, unsigned int, mpfq_2_64_src_vec, unsigned int);
/* missing vec_random */
/* missing vec_random2 */
static inline
int mpfq_2_64_vec_cmp(mpfq_2_64_dst_field, mpfq_2_64_src_vec, mpfq_2_64_src_vec, unsigned int);
static inline
int mpfq_2_64_vec_is_zero(mpfq_2_64_dst_field, mpfq_2_64_src_vec, unsigned int);
static inline
mpfq_2_64_dst_vec mpfq_2_64_vec_subvec(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, int);
static inline
mpfq_2_64_src_vec mpfq_2_64_vec_subvec_const(mpfq_2_64_dst_field, mpfq_2_64_src_vec, int);
static inline
mpfq_2_64_dst_elt mpfq_2_64_vec_coeff_ptr(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, int);
static inline
mpfq_2_64_src_elt mpfq_2_64_vec_coeff_ptr_const(mpfq_2_64_dst_field, mpfq_2_64_src_vec, int);
/* missing vec_asprint */
/* missing vec_fprint */
/* missing vec_print */
/* missing vec_sscan */
/* missing vec_fscan */
/* missing vec_scan */
void mpfq_2_64_vec_ur_init(mpfq_2_64_dst_field, mpfq_2_64_vec_ur *, unsigned int);
static inline
void mpfq_2_64_vec_ur_set_zero(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_set_vec(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec, unsigned int);
void mpfq_2_64_vec_ur_reinit(mpfq_2_64_dst_field, mpfq_2_64_vec_ur *, unsigned int, unsigned int);
void mpfq_2_64_vec_ur_clear(mpfq_2_64_dst_field, mpfq_2_64_vec_ur *, unsigned int);
static inline
void mpfq_2_64_vec_ur_set(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_setcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_elt_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_getcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_elt_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_add(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_sub(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_neg(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_ur_rev(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec_ur, unsigned int);
static inline
void mpfq_2_64_vec_scal_mul_ur(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec, mpfq_2_64_src_elt, unsigned int);
static inline
void mpfq_2_64_vec_conv_ur_n(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_conv_ur(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, mpfq_2_64_src_vec, unsigned int, mpfq_2_64_src_vec, unsigned int);
static inline
void mpfq_2_64_vec_reduce(mpfq_2_64_dst_field, mpfq_2_64_dst_vec, mpfq_2_64_dst_vec_ur, unsigned int);
static inline
mpfq_2_64_dst_vec_ur mpfq_2_64_vec_ur_subvec(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, int);
static inline
mpfq_2_64_src_vec_ur mpfq_2_64_vec_ur_subvec_const(mpfq_2_64_dst_field, mpfq_2_64_src_vec_ur, int);
static inline
mpfq_2_64_dst_elt mpfq_2_64_vec_ur_coeff_ptr(mpfq_2_64_dst_field, mpfq_2_64_dst_vec_ur, int);
static inline
mpfq_2_64_src_elt mpfq_2_64_vec_ur_coeff_ptr_const(mpfq_2_64_dst_field, mpfq_2_64_src_vec_ur, int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, Mpfq::gf2n::trivialities */
#define mpfq_2_64_vec_elt_stride(K, n)	((n)*sizeof(mpfq_2_64_elt))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, Mpfq::gf2n::trivialities */
#define mpfq_2_64_vec_ur_elt_stride(K, n)	((n)*sizeof(mpfq_2_64_elt_ur))

/* Polynomial functions */
static inline
void mpfq_2_64_poly_init(mpfq_2_64_dst_field, mpfq_2_64_poly, unsigned int);
static inline
void mpfq_2_64_poly_clear(mpfq_2_64_dst_field, mpfq_2_64_poly);
static inline
void mpfq_2_64_poly_set(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly);
void mpfq_2_64_poly_setmonic(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_setcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_elt, unsigned int);
static inline
void mpfq_2_64_poly_setcoeff_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, unsigned long, unsigned int);
static inline
void mpfq_2_64_poly_getcoeff(mpfq_2_64_dst_field, mpfq_2_64_dst_elt, mpfq_2_64_src_poly, unsigned int);
static inline
int mpfq_2_64_poly_deg(mpfq_2_64_dst_field, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_add(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_sub(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_set_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, unsigned long);
static inline
void mpfq_2_64_poly_add_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, unsigned long);
static inline
void mpfq_2_64_poly_sub_ui(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, unsigned long);
static inline
void mpfq_2_64_poly_neg(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_scal_mul(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_elt);
static inline
void mpfq_2_64_poly_mul(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
int mpfq_2_64_poly_divmod(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
void mpfq_2_64_poly_precomp_mod(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly);
void mpfq_2_64_poly_mod_pre(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_gcd(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
static inline
void mpfq_2_64_poly_xgcd(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_dst_poly, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
/* missing poly_random */
/* missing poly_random2 */
static inline
int mpfq_2_64_poly_cmp(mpfq_2_64_dst_field, mpfq_2_64_src_poly, mpfq_2_64_src_poly);
/* missing poly_asprint */
/* missing poly_fprint */
/* missing poly_print */
/* missing poly_sscan */
/* missing poly_fscan */
/* missing poly_scan */
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::gf2n::field::code_for_field_init */
static inline
void mpfq_2_64_field_init(mpfq_2_64_dst_field f)
{
    f->io_type=16;
}

/* *Mpfq::gf2n::field::code_for_field_setopt */
static inline
void mpfq_2_64_field_setopt(mpfq_2_64_dst_field f, unsigned long x MAYBE_UNUSED, void * y)
{
    assert(x == MPFQ_IO_TYPE);
    f->io_type=((unsigned long*)y)[0];
}

/* *Mpfq::defaults::flatdata::code_for_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_set(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_64_elt));
}

/* *Mpfq::gf2n::trivialities::code_for_set_ui */
static inline
void mpfq_2_64_set_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, unsigned long x)
{
    r[0] = x & 1UL;
    memset(r + 1, 0, sizeof(mpfq_2_64_elt) - sizeof(unsigned long));
}

/* *Mpfq::defaults::flatdata::code_for_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_set_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r)
{
    mpfq_2_64_vec_set_zero(K,(mpfq_2_64_dst_vec)r,1);
}

/* *Mpfq::gf2n::trivialities::code_for_get_ui */
static inline
unsigned long mpfq_2_64_get_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt r)
{
    return r[0] & 1UL;
}

/* *Mpfq::gf2n::trivialities::code_for_set_uipoly */
static inline
void mpfq_2_64_set_uipoly(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, unsigned long x)
{
    r[0] = x;
    memset(r + 1, 0, sizeof(mpfq_2_64_elt) - sizeof(unsigned long));
}

/* *Mpfq::gf2n::trivialities::code_for_set_uipoly_wide */
static inline
void mpfq_2_64_set_uipoly_wide(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, const unsigned long * x, unsigned int n)
{
    unsigned int i;
    for (i = 0 ; i < n && i < 2 ; i++)
        r[i] = x[i];
}

/* *Mpfq::gf2n::trivialities::code_for_get_uipoly */
static inline
unsigned long mpfq_2_64_get_uipoly(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt r)
{
    return r[0];
}

/* *Mpfq::gf2n::trivialities::code_for_get_uipoly_wide */
static inline
void mpfq_2_64_get_uipoly_wide(mpfq_2_64_dst_field K MAYBE_UNUSED, unsigned long * r, mpfq_2_64_src_elt x)
{
    unsigned int i;
    for(i = 0 ; i < 2 ; i++) r[i] = x[i];
}

/* *Mpfq::gf2n::trivialities::code_for_add */
static inline
void mpfq_2_64_add(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s1, mpfq_2_64_src_elt s2)
{
    int i;
    for(i = 0 ; i < 2 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::trivialities::code_for_mul */
static inline
void mpfq_2_64_mul(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s1, mpfq_2_64_src_elt s2)
{
    mpfq_2_64_elt_ur t;
    mpfq_2_64_mul_ur(K, t, s1, s2);
    mpfq_2_64_reduce(K, r, t);
}

/* *Mpfq::gf2n::trivialities::code_for_sqr */
static inline
void mpfq_2_64_sqr(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s)
{
    mpfq_2_64_elt_ur t;
    mpfq_2_64_sqr_ur(K, t, s);
    mpfq_2_64_reduce(K, r, t);
}

/* *Mpfq::gf2n::linearops::code_for_sqrt */
static inline
int mpfq_2_64_sqrt(mpfq_2_64_dst_field K, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s)
{
    static const unsigned long shuffle_table[256] = {
        0, 1, 16, 17, 2, 3, 18, 19,
        32, 33, 48, 49, 34, 35, 50, 51,
        4, 5, 20, 21, 6, 7, 22, 23,
        36, 37, 52, 53, 38, 39, 54, 55,
        64, 65, 80, 81, 66, 67, 82, 83,
        96, 97, 112, 113, 98, 99, 114, 115,
        68, 69, 84, 85, 70, 71, 86, 87,
        100, 101, 116, 117, 102, 103, 118, 119,
        8, 9, 24, 25, 10, 11, 26, 27,
        40, 41, 56, 57, 42, 43, 58, 59,
        12, 13, 28, 29, 14, 15, 30, 31,
        44, 45, 60, 61, 46, 47, 62, 63,
        72, 73, 88, 89, 74, 75, 90, 91,
        104, 105, 120, 121, 106, 107, 122, 123,
        76, 77, 92, 93, 78, 79, 94, 95,
        108, 109, 124, 125, 110, 111, 126, 127,
        128, 129, 144, 145, 130, 131, 146, 147,
        160, 161, 176, 177, 162, 163, 178, 179,
        132, 133, 148, 149, 134, 135, 150, 151,
        164, 165, 180, 181, 166, 167, 182, 183,
        192, 193, 208, 209, 194, 195, 210, 211,
        224, 225, 240, 241, 226, 227, 242, 243,
        196, 197, 212, 213, 198, 199, 214, 215,
        228, 229, 244, 245, 230, 231, 246, 247,
        136, 137, 152, 153, 138, 139, 154, 155,
        168, 169, 184, 185, 170, 171, 186, 187,
        140, 141, 156, 157, 142, 143, 158, 159,
        172, 173, 188, 189, 174, 175, 190, 191,
        200, 201, 216, 217, 202, 203, 218, 219,
        232, 233, 248, 249, 234, 235, 250, 251,
        204, 205, 220, 221, 206, 207, 222, 223,
        236, 237, 252, 253, 238, 239, 254, 255,
    };
    
    mpfq_2_64_elt sqrt_t ={ 0xaUL, 0xffffffffUL, }
    ;
    
    mpfq_2_64_elt odd, even;
    mpfq_2_64_elt_ur odd_t;
    unsigned long t;
#define	EVEN_MASK	(((unsigned long)-1)/3UL)
#define	ODD_MASK	((EVEN_MASK)<<1)
    unsigned int i;
    for(i = 0 ; i < 2 ; i++) {
        even[i] = s[i] & EVEN_MASK;
    }
    for(i = 0 ; i < 1 ; i++) {
        t = even[2*i];   t |= t >> 7;
              even[i]  = shuffle_table[t & 255];
        t >>= 16; even[i] |= shuffle_table[t & 255] << 8;
        t = even[2*i+1]; t |= t >> 7;
              even[i] |= shuffle_table[t & 255] << 16;
        t >>= 16; even[i] |= shuffle_table[t & 255] << 24;
    }
    memset(even + 1, 0, 1 * sizeof(unsigned long));
    for(i = 0 ; i < 2 ; i++) {
        odd[i] = (s[i] & ODD_MASK) >> 1;
    }
    for(i = 0 ; i < 1 ; i++) {
        t = odd[2*i];   t |= t >> 7;
              odd[i]  = shuffle_table[t & 255];
        t >>= 16; odd[i] |= shuffle_table[t & 255] << 8;
        t = odd[2*i+1]; t |= t >> 7;
              odd[i] |= shuffle_table[t & 255] << 16;
        t >>= 16; odd[i] |= shuffle_table[t & 255] << 24;
    }
    memset(odd + 1, 0, 1 * sizeof(unsigned long));
    mpfq_2_64_mul_ur(K, odd_t, odd, sqrt_t);
    for(i = 0 ; i < (2+1)/2 ; i++) {
        odd_t[i] ^= even[i];
    }
    mpfq_2_64_reduce(K, r, odd_t);
    return 1;
}

/* *Mpfq::defaults::pow::code_for_pow, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_pow(mpfq_2_64_dst_field k, mpfq_2_64_dst_elt res, mpfq_2_64_src_elt r, unsigned long * x, size_t n)
{
    mpfq_2_64_elt u, a;
    long i, j, lead;     /* it is a signed type */
    unsigned long mask;
    
    /* get the correct (i,j) position of the most significant bit in x */
    for(i = ((long)n)-1; i>=0 && x[i]==0; i--)
        ;
    if (i < 0) {
        /* power zero gets 1 */
        mpfq_2_64_set_ui(k, res, 1);
        return;
    }
    j = 32 - 1;
    mask = (1UL<<j);
    for( ; (x[i]&mask)==0 ;j--, mask>>=1)
        ;
    lead = i*32+j;      /* Ensured. */
    
    mpfq_2_64_init(k, &u);
    mpfq_2_64_init(k, &a);
    mpfq_2_64_set(k, a, r);
    for( ; lead > 0; lead--) {
        if (j-- == 0) {
            i--;
            j = 32-1;
            mask = (1UL<<j);
        } else {
            mask >>= 1;
        }
        if (x[i]&mask) {
            mpfq_2_64_sqr(k, u, a);
            mpfq_2_64_mul(k, a, u, r);
        } else {
            mpfq_2_64_sqr(k, a,a);
        }
    }
    mpfq_2_64_set(k, res, a);
    mpfq_2_64_clear(k, &u);
    mpfq_2_64_clear(k, &a);
}

/* *Mpfq::gf2n::trivialities::code_for_add_ui */
static inline
void mpfq_2_64_add_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s, unsigned long x)
{
    mpfq_2_64_set(K, r, s);
    r[0] ^= x & 1UL;
}

/* *Mpfq::gf2n::trivialities::code_for_mul_ui */
static inline
void mpfq_2_64_mul_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s, unsigned long x)
{
    if (x & 1UL) {
        mpfq_2_64_set(K, r, s);
    } else {
        memset(r, 0, sizeof(mpfq_2_64_elt));
    }
}

/* *Mpfq::gf2n::trivialities::code_for_add_uipoly */
static inline
void mpfq_2_64_add_uipoly(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s, unsigned long x)
{
    mpfq_2_64_set(K, r, s);
    r[0] ^= x;
}

/* *Mpfq::gf2n::trivialities::code_for_mul_uipoly */
static inline
void mpfq_2_64_mul_uipoly(mpfq_2_64_dst_field k, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s, unsigned long x)
{
    mpfq_2_64_elt xx;
    mpfq_2_64_init(k, &xx);
    mpfq_2_64_set_uipoly(k, xx, x);
    mpfq_2_64_mul(k, r, s, xx);
    mpfq_2_64_clear(k, &xx);
}

/* *Mpfq::gf2n::inversion::code_for_inv */
/* Triggered by: inv */
static inline
void mpfq_2_64_longshift_left(unsigned long * dst, const unsigned long * src, int n, int s)
{
    int m = s / 32;
    int i;
    s %= 32;
    if (s > 0) {
        for(i = n-m-1 ; i > 0 ; i--) {
            dst[m+i] = src[i] << s ^ src[i-1] >> (32-s);
        }
        dst[m] = src[0] << s;
    } else {
        for(i = n-m-1 ; i > 0 ; i--) {
            dst[m+i] = src[i];
        }
        dst[m] = src[0];
    }
    for(i = m-1 ; i>= 0 ; i--) {
        dst[i] = 0UL;
    }
    
}

/* *Mpfq::gf2n::inversion::code_for_inv */
/* Triggered by: inv */
static inline
void mpfq_2_64_longaddshift_left(unsigned long * dst, const unsigned long * src, int n, int s)
{
    int m = s / 32;
    int i;
    s %= 32;
    if (s>0) {
        for(i = n-m-1 ; i > 0 ; i--) {
            dst[m+i] ^= src[i] << s ^ src[i-1] >> (32-s);
        }
        dst[m] ^= src[0] << s;
    } else {
        for(i = n-m-1 ; i > 0 ; i--) {
            dst[m+i] ^= src[i];
        }
        dst[m] ^= src[0];
    }
}

/* *Mpfq::gf2n::inversion::code_for_inv */
static inline
int mpfq_2_64_inv(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s)
{
    unsigned long a[3] = { 0x80000000UL, 0xdUL, 0x80000000UL, };
    unsigned long b[3];
    unsigned long u[3] = { 0, };
    unsigned long v[3] = { 1, 0, };
    unsigned long x;
    int ia, ib;
    int i,d;
    
    if (mpfq_2_64_cmp_ui(K, s, 0UL) == 0)
        return 0;
    {
        unsigned long z;
        z = s[0] << 31;
        b[0] = z;
        z = s[0] >>  1 ^ s[1] << 31;
        b[1] = z;
        z = s[1] >>  1;
        b[2] = z;
    }
    ib = mpfq_clzlx(b, 3);
    ia = 0;
    
    mpfq_2_64_longshift_left(b,b,3,ib);
    
    for(d = ib - ia ; ; ) {
            if (d == 0) {
                    for(i = 0 ; i < 3 ; i++) v[i] ^= u[i];
            b[0] ^= a[0]; x = b[0];
            b[1] ^= a[1]; x |= b[1];
            b[2] ^= a[2]; x |= b[2];
                    if (!x) { memcpy(r,u,2 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(b,3);
                    ib += t;
                    d += t;
                    mpfq_2_64_longshift_left(b,b,3,t);
            }
            for(;d > 0;) {
                    mpfq_2_64_longaddshift_left(u,v,3,d);
            a[0] ^= b[0]; x = a[0];
            a[1] ^= b[1]; x |= a[1];
            a[2] ^= b[2]; x |= a[2];
                    if (!x) { memcpy(r,v,2 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(a,3);
                    ia += t;
                    d -= t;
                    mpfq_2_64_longshift_left(a,a,3,t);
            } 
            if (d == 0) {
                    for(i = 0 ; i < 3 ; i++) u[i] ^= v[i];
            a[0] ^= b[0]; x = a[0];
            a[1] ^= b[1]; x |= a[1];
            a[2] ^= b[2]; x |= a[2];
                    if (!x) { memcpy(r,v,2 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(a,3);
                    ia += t;
                    d -= t;
                    mpfq_2_64_longshift_left(a,a,3,t);
            }
            for(;d < 0;) {
                    mpfq_2_64_longaddshift_left(v,u,3,-d);
            b[0] ^= a[0]; x = b[0];
            b[1] ^= a[1]; x |= b[1];
            b[2] ^= a[2]; x |= b[2];
                    if (!x) { memcpy(r,u,2 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(b,3);
                    ib += t;
                    d += t;
                    mpfq_2_64_longshift_left(b,b,3,t);
            }
    }
}

/* *Mpfq::gf2n::linearops::code_for_as_solve */
static inline
void mpfq_2_64_as_solve(mpfq_2_64_dst_field K, mpfq_2_64_dst_elt r, mpfq_2_64_src_elt s)
{
    static const mpfq_2_64_elt t[64] = {
        { 0x278adc02UL, 0x19c9369fUL, },
        { 0xb2383ee4UL, 0x84b2b22aUL, },
        { 0xb2383ee6UL, 0x84b2b22aUL, },
        { 0x95b3e3f6UL, 0x9d7b84b4UL, },
        { 0xb2383ee2UL, 0x84b2b22aUL, },
        { 0x9213f790UL, 0x37c470b4UL, },
        { 0x95b3e3feUL, 0x9d7b84b4UL, },
        { 0x105137cUL, 0x1000aUL, },
        { 0xb2383ef2UL, 0x84b2b22aUL, },
        { 0x8edce1fcUL, 0x368e964aUL, },
        { 0x9213f7b0UL, 0x37c470b4UL, },
        { 0x278fdf4cUL, 0x19c9368eUL, },
        { 0x95b3e3beUL, 0x9d7b84b4UL, },
        { 0xbc7d4570UL, 0x2e4da23cUL, },
        { 0x10513fcUL, 0x1000aUL, },
        { 0xbac24232UL, 0x84f35772UL, },
        { 0xb2383ff2UL, 0x84b2b22aUL, },
        { 0x9314e4fcUL, 0x37c570baUL, },
        { 0x8edce3fcUL, 0x368e964aUL, },
        { 0x213cdb0eUL, 0xb377c390UL, },
        { 0x9213f3b0UL, 0x37c470b4UL, },
        { 0xa99c24f2UL, 0x85ed5a3aUL, },
        { 0x278fd74cUL, 0x19c9368eUL, },
        { 0x780000eUL, 0xaabff000UL, },
        { 0x95b3f3beUL, 0x9d7b84b4UL, },
        { 0xb03038f2UL, 0x84b6b3daUL, },
        { 0xbc7d6570UL, 0x2e4da23cUL, },
        { 0x3494ffcUL, 0x511eaUL, },
        { 0x10553fcUL, 0x1000aUL, },
        { 0x343c6c0eUL, 0xae0c0220UL, },
        { 0xbac2c232UL, 0x84f35772UL, },
        { 0x8000000eUL, 0x80000000UL, },
        { 0xb2393ff2UL, 0x84b2b22aUL, },
        { 0x202bc97eUL, 0xb376c29cUL, },
        { 0x9316e4fcUL, 0x37c570baUL, },
        { 0x8879e6ceUL, 0x9c306248UL, },
        { 0x8ed8e3fcUL, 0x368e964aUL, },
        { 0x8e47e70UL, 0x41e42cUL, },
        { 0x2134db0eUL, 0xb377c390UL, },
        { 0xa60f56ceUL, 0x85b9b108UL, },
        { 0x9203f3b0UL, 0x37c470b4UL, },
        { 0x21f3cb4cUL, 0x19dd3b6eUL, },
        { 0xa9bc24f2UL, 0x85ed5a3aUL, },
        { 0x2c428ac0UL, 0x198ddf68UL, },
        { 0x27cfd74cUL, 0x19c9368eUL, },
        { 0x31ca84b0UL, 0x4b7c684UL, },
        { 0x700000eUL, 0xaabff000UL, },
        { 0x89ffefbeUL, 0x80406554UL, },
        { 0x94b3f3beUL, 0x9d7b84b4UL, },
        { 0x32bfa74cUL, 0x18c1354eUL, },
        { 0xb23038f2UL, 0x84b6b3daUL, },
        { 0xf74627eUL, 0xaaf613ccUL, },
        { 0xb87d6570UL, 0x2e4da23cUL, },
        { 0xb3342a8cUL, 0x3248b3d6UL, },
        { 0xb494ffcUL, 0x511eaUL, },
        { 0xe70700eUL, 0xb60813c0UL, },
        { 0x110553fcUL, 0x1000aUL, },
        { 0x5393ffcUL, 0x1e0d022aUL, },
        { 0x143c6c0eUL, 0xae0c0220UL, },
        { 0x143c6c00UL, 0xe0c0220UL, },
        { 0xfac2c232UL, 0x84f35772UL, },
        { 0x48fbfdceUL, 0xc041e559UL, },
        { 0xeUL, 0x80000000UL, },
        { 0x0UL, 0x0UL, },
    };
    const mpfq_2_64_elt * ptr = t;
    unsigned int i,j;
    memset(r, 0, sizeof(mpfq_2_64_elt));
    for(i = 0 ; i < 2 ; i++) {
        unsigned long a = s[i];
        for(j = 0 ; j < 32 && ptr != &t[64]; j++, ptr++) {
            if (a & 1UL) {
                mpfq_2_64_add(K, r, r, *ptr);
            }
            a >>= 1;
        }
    }
}

/* *Mpfq::gf2n::linearops::code_for_trace */
static inline
unsigned long mpfq_2_64_trace(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt s)
{
    return ((s[1]>>31) ^ (s[1]>>29)) & 1;
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_elt_ur_set(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur r, mpfq_2_64_src_elt_ur s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_64_elt_ur));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_elt_ur_set_elt(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur r, mpfq_2_64_src_elt s)
{
    memset(r, 0, sizeof(mpfq_2_64_elt_ur)); memcpy(r,s,sizeof(mpfq_2_64_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_elt_ur_set_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur r)
{
    memset(r, 0, sizeof(mpfq_2_64_elt_ur));
}

/* *Mpfq::gf2n::trivialities::code_for_elt_ur_set_ui */
static inline
void mpfq_2_64_elt_ur_set_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur r, unsigned long x)
{
    r[0] = x & 1UL;
    memset(r + 1, 0, sizeof(mpfq_2_64_elt_ur) - sizeof(unsigned long));
}

/* *Mpfq::gf2n::trivialities::code_for_elt_ur_add */
static inline
void mpfq_2_64_elt_ur_add(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur r, mpfq_2_64_src_elt_ur s1, mpfq_2_64_src_elt_ur s2)
{
    int i;
    for(i = 0 ; i < 4 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::mul::code_for_mul_ur */
static inline
void mpfq_2_64_mul_ur(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur t, mpfq_2_64_src_elt s1, mpfq_2_64_src_elt s2)
{
    gf2x_mul2(t, s1, s2);
}

/* *Mpfq::gf2n::squaring::code_for_sqr_ur */
static inline
void mpfq_2_64_sqr_ur(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur t, mpfq_2_64_src_elt s)
{
    static const unsigned long g[16] = {
        0, 1, 4, 5, 16, 17, 20, 21,
        64, 65, 68, 69, 80, 81, 84, 85,
    };
    {
        unsigned long u;
        u = g[s[0]       & 15];
    t[0]  = u;
        u = g[s[0] >>  4 & 15];
    t[0] ^= u <<  8;
        u = g[s[0] >>  8 & 15];
    t[0] ^= u << 16;
        u = g[s[0] >> 12 & 15];
    t[0] ^= u << 24;
        u = g[s[0] >> 16 & 15];
    t[1]  = u;
        u = g[s[0] >> 20 & 15];
    t[1] ^= u <<  8;
        u = g[s[0] >> 24 & 15];
    t[1] ^= u << 16;
        u = g[s[0] >> 28 & 15];
    t[1] ^= u << 24;
        u = g[s[1]       & 15];
    t[2]  = u;
        u = g[s[1] >>  4 & 15];
    t[2] ^= u <<  8;
        u = g[s[1] >>  8 & 15];
    t[2] ^= u << 16;
        u = g[s[1] >> 12 & 15];
    t[2] ^= u << 24;
        u = g[s[1] >> 16 & 15];
    t[3]  = u;
        u = g[s[1] >> 20 & 15];
    t[3] ^= u <<  8;
        u = g[s[1] >> 24 & 15];
    t[3] ^= u << 16;
        u = g[s[1] >> 28 & 15];
    t[3] ^= u << 24;
    }
}

/* *Mpfq::gf2n::reduction::code_for_reduce */
static inline
void mpfq_2_64_reduce(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt r, mpfq_2_64_dst_elt_ur t)
{
    {
        unsigned long s[3];
        /* 63 excess bits */
        {
            unsigned long z;
            z = t[0];
            s[0] = z;
            z = t[1];
            s[1] = z;
        }
        memset(s + 2, 0, 1 * sizeof(unsigned long));
        {
            unsigned long z;
            z = t[2];
            s[0]^= z <<  4;
            s[0]^= z <<  3;
            s[0]^= z <<  1;
            s[0]^= z;
            z >>= 28;
            z^= t[3] <<  4;
            s[1]^= z;
            z >>= 1;
            z^= t[3] >> 28 << 31;
            s[1]^= z;
            z >>= 2;
            z^= t[3] >> 29 << 30;
            s[1]^= z;
            z >>= 1;
            s[1]^= z;
            z >>= 28;
            s[2]^= z;
            z >>= 1;
            s[2]^= z;
        }
        /* 3 excess bits */
        {
            unsigned long z;
            z = s[0];
            r[0] = z;
            z = s[1];
            r[1] = z;
        }
        {
            unsigned long z;
            z = s[2];
            r[0]^= z <<  4;
            r[0]^= z <<  3;
            r[0]^= z <<  1;
            r[0]^= z;
        }
    }
}

/* *Mpfq::gf2n::trivialities::code_for_cmp */
static inline
int mpfq_2_64_cmp(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt a, mpfq_2_64_src_elt b)
{
    return memcmp(a, b, sizeof(mpfq_2_64_elt));
}

/* *Mpfq::gf2n::trivialities::code_for_cmp_ui */
static inline
int mpfq_2_64_cmp_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt r, unsigned long x)
{
    int i;
    if (r[0] < (x & 1UL)) return -1;
    if (r[0] > (x & 1UL)) return 1;
    for(i = 1 ; i < 2 ; i++) {
        if (r[i]) return 1;
    }
    return 0;
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gf2n::trivialities */
static inline
int mpfq_2_64_is_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(mpfq_2_64_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_vec_set(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec r, mpfq_2_64_src_vec s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_2_64_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_vec_set_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(mpfq_2_64_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_setcoeff(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_elt x, unsigned int i)
{
    mpfq_2_64_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoeff_ui, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_setcoeff_ui(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, unsigned long x, unsigned int i)
{
    mpfq_2_64_set_ui(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_getcoeff(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt x, mpfq_2_64_src_vec w, unsigned int i)
{
    mpfq_2_64_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_add(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, mpfq_2_64_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_neg(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_64_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_rev(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    mpfq_2_64_elt tmp[1];
    mpfq_2_64_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_2_64_set(K, tmp[0], u[i]);
        mpfq_2_64_set(K, w[i], u[n-1-i]);
        mpfq_2_64_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_2_64_set(K, w[nn], u[nn]);
    mpfq_2_64_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_sub(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, mpfq_2_64_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_64_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_scal_mul(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, mpfq_2_64_src_elt x, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_mul(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv */
static inline
void mpfq_2_64_vec_conv(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_src_vec u, unsigned int n, mpfq_2_64_src_vec v, unsigned int m)
{
    mpfq_2_64_vec_ur tmp;
    mpfq_2_64_vec_ur_init(K, &tmp, m+n-1);
    mpfq_2_64_vec_conv_ur(K, tmp, u, n, v, m);
    mpfq_2_64_vec_reduce(K, w, tmp, m+n-1);
    mpfq_2_64_vec_ur_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec */
static inline
int mpfq_2_64_vec_cmp(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec u, mpfq_2_64_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = mpfq_2_64_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_is_zero, Mpfq::defaults::vec */
static inline
int mpfq_2_64_vec_is_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!mpfq_2_64_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec */
static inline
mpfq_2_64_dst_vec mpfq_2_64_vec_subvec(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_2_64_src_vec mpfq_2_64_vec_subvec_const(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_2_64_dst_elt mpfq_2_64_vec_coeff_ptr(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_2_64_src_elt mpfq_2_64_vec_coeff_ptr_const(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_vec_ur_set_zero(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur r, unsigned int n)
{
    memset(r, 0, n*sizeof(mpfq_2_64_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_set_vec, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_set_vec(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_elt_ur_set_elt(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_64_vec_ur_set(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur r, mpfq_2_64_src_vec_ur s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_2_64_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_setcoeff(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_elt_ur x, unsigned int i)
{
    mpfq_2_64_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_getcoeff(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_elt_ur x, mpfq_2_64_src_vec_ur w, unsigned int i)
{
    mpfq_2_64_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_add(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec_ur u, mpfq_2_64_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_sub(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec_ur u, mpfq_2_64_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_neg(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_64_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_ur_rev(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec_ur u, unsigned int n)
{
    unsigned int nn = n >> 1;
    mpfq_2_64_elt_ur tmp[1];
    mpfq_2_64_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_2_64_elt_ur_set(K, tmp[0], u[i]);
        mpfq_2_64_elt_ur_set(K, w[i], u[n-1-i]);
        mpfq_2_64_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_2_64_elt_ur_set(K, w[nn], u[nn]);
    mpfq_2_64_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_scal_mul_ur(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec u, mpfq_2_64_src_elt x, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_mul_ur(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur */
/* Triggered by: vec_conv_ur */
static inline
void mpfq_2_64_vec_conv_ur_n(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec u, mpfq_2_64_src_vec v, unsigned int n)
{
    if (n == 0)
        return;
    if (n == 1) {
        mpfq_2_64_mul_ur(K, w[0], u[0], v[0]);
        return;
    }
    if (n == 2) {  // Kara 2
        mpfq_2_64_elt t1, t2;
        mpfq_2_64_init(K, &t1);
        mpfq_2_64_init(K, &t2);
        mpfq_2_64_mul_ur(K, w[0], u[0], v[0]);
        mpfq_2_64_mul_ur(K, w[2], u[1], v[1]);
        mpfq_2_64_add(K, t1, u[0], u[1]);
        mpfq_2_64_add(K, t2, v[0], v[1]);
        mpfq_2_64_mul_ur(K, w[1], t1, t2);
        mpfq_2_64_elt_ur_sub(K, w[1], w[1], w[0]);
        mpfq_2_64_elt_ur_sub(K, w[1], w[1], w[2]);
        mpfq_2_64_clear(K, &t1);
        mpfq_2_64_clear(K, &t2);
        return;
    }
    if (n == 3) {  // do it in 6
        mpfq_2_64_elt t1, t2;
        mpfq_2_64_elt_ur s;
        mpfq_2_64_init(K, &t1);
        mpfq_2_64_init(K, &t2);
        mpfq_2_64_elt_ur_init(K, &s);
        // a0*b0*(1 - X)
        mpfq_2_64_mul_ur(K, w[0], u[0], v[0]);
        mpfq_2_64_elt_ur_neg(K, w[1], w[0]);
        // a1*b1*(-X + 2*X^2 - X^3)
        mpfq_2_64_mul_ur(K, w[2], u[1], v[1]);
        mpfq_2_64_elt_ur_neg(K, w[3], w[2]);
        mpfq_2_64_elt_ur_add(K, w[2], w[2], w[2]);
        mpfq_2_64_elt_ur_add(K, w[1], w[1], w[3]);
        // a2*b2*(-X^3+X^4)
        mpfq_2_64_mul_ur(K, w[4], u[2], v[2]);
        mpfq_2_64_elt_ur_sub(K, w[3], w[3], w[4]);
        // (a0+a1)*(b0+b1)*(X - X^2)
        mpfq_2_64_add(K, t1, u[0], u[1]);
        mpfq_2_64_add(K, t2, v[0], v[1]);
        mpfq_2_64_mul_ur(K, s, t1, t2);
        mpfq_2_64_elt_ur_add(K, w[1], w[1], s);
        mpfq_2_64_elt_ur_sub(K, w[2], w[2], s);
        // (a1+a2)*(b1+b2)*(X^3 - X^2)
        mpfq_2_64_add(K, t1, u[1], u[2]);
        mpfq_2_64_add(K, t2, v[1], v[2]);
        mpfq_2_64_mul_ur(K, s, t1, t2);
        mpfq_2_64_elt_ur_add(K, w[3], w[3], s);
        mpfq_2_64_elt_ur_sub(K, w[2], w[2], s);
        // (a0+a1+a2)*(b0+b1+b2)* X^2
        mpfq_2_64_add(K, t1, u[0], t1);
        mpfq_2_64_add(K, t2, v[0], t2);
        mpfq_2_64_mul_ur(K, s, t1, t2);
        mpfq_2_64_elt_ur_add(K, w[2], w[2], s);
        return;
    }
    unsigned int n0, n1;
    n0 = n / 2;
    n1 = n - n0;
    mpfq_2_64_vec_conv_ur_n(K, w, u, v, n0);
    mpfq_2_64_vec_conv_ur_n(K, w + 2*n0, u + n0, v + n0, n1);
    mpfq_2_64_elt_ur_set_ui(K, w[2*n0-1], 0);
    
    mpfq_2_64_vec tmpu, tmpv;
    mpfq_2_64_vec_ur tmpw;
    mpfq_2_64_vec_init(K, &tmpu, n1);
    mpfq_2_64_vec_init(K, &tmpv, n1);
    mpfq_2_64_vec_ur_init(K, &tmpw, 2*n1-1);
    
    mpfq_2_64_vec_set(K, tmpu, u, n0);
    if (n1 != n0) 
        mpfq_2_64_set_ui(K, tmpu[n0], 0);
    mpfq_2_64_vec_add(K, tmpu, tmpu, u+n0, n1);
    mpfq_2_64_vec_set(K, tmpv, v, n0);
    if (n1 != n0) 
        mpfq_2_64_set_ui(K, tmpv[n0], 0);
    mpfq_2_64_vec_add(K, tmpv, tmpv, v+n0, n1);
    mpfq_2_64_vec_conv_ur_n(K, tmpw, tmpu, tmpv, n1);
    mpfq_2_64_vec_ur_sub(K, tmpw, tmpw, w, 2*n0-1);
    mpfq_2_64_vec_ur_sub(K, tmpw, tmpw, w + 2*n0, 2*n1-1);
    mpfq_2_64_vec_ur_add(K, w + n0, w + n0, tmpw, 2*n1-1);
    
    mpfq_2_64_vec_clear(K, &tmpu, n1);
    mpfq_2_64_vec_clear(K, &tmpv, n1);
    mpfq_2_64_vec_ur_clear(K, &tmpw, 2*n1-1);
    return;
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur */
static inline
void mpfq_2_64_vec_conv_ur(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur w, mpfq_2_64_src_vec u, unsigned int n, mpfq_2_64_src_vec v, unsigned int m)
{
    unsigned int i, j MAYBE_UNUSED, k;
    mpfq_2_64_elt_ur acc, z;
    mpfq_2_64_elt_ur_init(K, &acc);
    mpfq_2_64_elt_ur_init(K, &z);
    // swap pointers to have n <= m
    mpfq_2_64_src_vec uu, vv;
    if (n <= m) {
        uu = u; vv = v;
    } else {
        uu = v; vv = u;
        unsigned int tmp = n;
        n = m; m = tmp;
    }
    for(k = 0; k < n; ++k) {
        mpfq_2_64_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i <= k; ++i) {
            mpfq_2_64_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_64_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_64_elt_ur_set(K, w[k], acc);
    }
    for(k = n; k < m; ++k) {
        mpfq_2_64_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i < n; ++i) {
            mpfq_2_64_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_64_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_64_elt_ur_set(K, w[k], acc);
    }
    for(k = m; k < n+m-1; ++k) {
        mpfq_2_64_mul_ur(K, acc, uu[k-m+1], vv[m-1]);
        for(i = k-m+2; i < n; ++i) {
            mpfq_2_64_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_64_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_64_elt_ur_set(K, w[k], acc);
    }
    mpfq_2_64_elt_ur_clear(K, &acc);
    mpfq_2_64_elt_ur_clear(K, &z);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec */
static inline
void mpfq_2_64_vec_reduce(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec w, mpfq_2_64_dst_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_reduce(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec */
static inline
mpfq_2_64_dst_vec_ur mpfq_2_64_vec_ur_subvec(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_2_64_src_vec_ur mpfq_2_64_vec_ur_subvec_const(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_2_64_dst_elt mpfq_2_64_vec_ur_coeff_ptr(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_2_64_src_elt mpfq_2_64_vec_ur_coeff_ptr_const(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::poly::code_for_poly_init */
static inline
void mpfq_2_64_poly_init(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_poly p, unsigned int n)
{
    mpfq_2_64_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear */
static inline
void mpfq_2_64_poly_clear(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_poly p)
{
    mpfq_2_64_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set */
static inline
void mpfq_2_64_poly_set(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_2_64_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff */
static inline
void mpfq_2_64_poly_setcoeff(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_elt x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_2_64_vec_set_zero(k, mpfq_2_64_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_2_64_vec_setcoeff(k, w->c, x, i);
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff_ui */
static inline
void mpfq_2_64_poly_setcoeff_ui(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, unsigned long x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_2_64_vec_set_zero(k, mpfq_2_64_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_2_64_vec_setcoeff_ui(k, w->c, x, i);
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_getcoeff */
static inline
void mpfq_2_64_poly_getcoeff(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_elt x, mpfq_2_64_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       mpfq_2_64_set_ui(k,x,0);
    } else {
       mpfq_2_64_vec_getcoeff(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg */
static inline
int mpfq_2_64_poly_deg(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    mpfq_2_64_elt temp;
    mpfq_2_64_init(K, &temp);
    mpfq_2_64_vec_getcoeff(K, temp, w->c, deg);
    int comp=mpfq_2_64_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           mpfq_2_64_vec_getcoeff(K, temp, w->c, deg);
           comp=mpfq_2_64_cmp_ui(K, temp, 0);
        }
    }
    mpfq_2_64_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add */
static inline
void mpfq_2_64_poly_add(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, mpfq_2_64_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_2_64_vec_add(k, w->c, u->c, v->c, su);
        mpfq_2_64_vec_set(k, mpfq_2_64_vec_subvec(k, w->c, su), mpfq_2_64_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_2_64_vec_add(k, w->c, u->c, v->c, sv);
        mpfq_2_64_vec_set(k, mpfq_2_64_vec_subvec(k, w->c, sv), mpfq_2_64_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_sub */
static inline
void mpfq_2_64_poly_sub(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, mpfq_2_64_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_2_64_vec_sub(k, w->c, u->c, v->c, su);
        mpfq_2_64_vec_neg(k, mpfq_2_64_vec_subvec(k, w->c, su), mpfq_2_64_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_2_64_vec_sub(k, w->c, u->c, v->c, sv);
        mpfq_2_64_vec_set(k, mpfq_2_64_vec_subvec(k, w->c, sv), mpfq_2_64_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_set_ui */
static inline
void mpfq_2_64_poly_set_ui(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, unsigned long x)
{
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_64_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui */
static inline
void mpfq_2_64_poly_add_ui(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_64_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_2_64_poly_deg(k, w);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_2_64_vec_set(k, mpfq_2_64_vec_subvec(k, w->c, 1), mpfq_2_64_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_2_64_add_ui(k, mpfq_2_64_vec_coeff_ptr(k, w->c, 0), mpfq_2_64_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui */
static inline
void mpfq_2_64_poly_sub_ui(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_64_elt temp;
        mpfq_2_64_init(k, &temp);
        mpfq_2_64_set_ui(k, temp, x);
        mpfq_2_64_neg(k, mpfq_2_64_vec_coeff_ptr(k, w->c, 0), temp);
        w->size = mpfq_2_64_cmp_ui(k, temp, 0);
        mpfq_2_64_clear(k, &temp);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_2_64_vec_set(k, mpfq_2_64_vec_subvec(k, w->c, 1), mpfq_2_64_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_2_64_sub_ui(k, mpfq_2_64_vec_coeff_ptr(k, w->c, 0), mpfq_2_64_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_neg */
static inline
void mpfq_2_64_poly_neg(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_2_64_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul */
static inline
void mpfq_2_64_poly_scal_mul(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, mpfq_2_64_src_elt x)
{
    if (mpfq_2_64_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_2_64_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_mul */
static inline
void mpfq_2_64_poly_mul(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly w, mpfq_2_64_src_poly u, mpfq_2_64_src_poly v)
{
    unsigned int usize = mpfq_2_64_poly_deg(k, u)+1;
    unsigned int vsize = mpfq_2_64_poly_deg(k, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        mpfq_2_64_vec_reinit(k, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    mpfq_2_64_vec_conv(k, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
    w->size = 1 + mpfq_2_64_poly_deg(k, w);
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly */
static inline
void mpfq_2_64_poly_gcd(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly g, mpfq_2_64_src_poly a0, mpfq_2_64_src_poly b0)
{
    mpfq_2_64_poly a,b,q,r;
    int da0=mpfq_2_64_poly_deg(k,a0), db0=mpfq_2_64_poly_deg(k,b0);
    if (db0==-1)
     mpfq_2_64_poly_set(k,g,a0);
    else {
     mpfq_2_64_poly_init(k,a,da0+1);
     mpfq_2_64_poly_init(k,b,db0+1);
     mpfq_2_64_poly_init(k,q,1);
     mpfq_2_64_poly_init(k,r,db0);
     mpfq_2_64_poly_set(k,a,a0);
     mpfq_2_64_poly_set(k,b,b0);
     while (mpfq_2_64_poly_deg(k,b)>=0) {
      mpfq_2_64_poly_divmod(k,q,r,a,b);
      mpfq_2_64_poly_set(k,a,b);
      mpfq_2_64_poly_set(k,b,r); 
     }
     mpfq_2_64_poly_setmonic(k,g,a);
    mpfq_2_64_poly_clear(k,a);
    mpfq_2_64_poly_clear(k,b);
    mpfq_2_64_poly_clear(k,q);
    mpfq_2_64_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly */
static inline
void mpfq_2_64_poly_xgcd(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_dst_poly g, mpfq_2_64_dst_poly u0, mpfq_2_64_dst_poly v0, mpfq_2_64_src_poly a0, mpfq_2_64_src_poly b0)
{
    mpfq_2_64_poly a,b,u,v,w,x,q,r;
    mpfq_2_64_elt c;
    mpfq_2_64_init(k,&c);
    mpfq_2_64_set_ui(k,c,0);        /* placate gcc */
    int da0=mpfq_2_64_poly_deg(k,a0), db0=mpfq_2_64_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      mpfq_2_64_poly_set(k,u0,a0);
      mpfq_2_64_poly_set(k,v0,b0);
      mpfq_2_64_poly_set(k,g,a0);
     } else {
      mpfq_2_64_poly_getcoeff(k,c,a0,da0);
      mpfq_2_64_inv(k,c,c);
      mpfq_2_64_poly_scal_mul(k,g,a0,c);
      mpfq_2_64_poly_set(k,v0,b0);
      mpfq_2_64_poly_set(k,u0,b0);
      mpfq_2_64_poly_setcoeff(k,u0,c,0);
     }
    }
    else {
     mpfq_2_64_poly_init(k,a,da0+1);
     mpfq_2_64_poly_init(k,b,db0+1);
     mpfq_2_64_poly_init(k,q,1);
     mpfq_2_64_poly_init(k,r,db0);
     mpfq_2_64_poly_set(k,a,a0);
     mpfq_2_64_poly_set(k,b,b0);
     mpfq_2_64_poly_init(k,u,1);
     mpfq_2_64_poly_init(k,v,1);
     mpfq_2_64_poly_init(k,w,1);
     mpfq_2_64_poly_init(k,x,1);
     mpfq_2_64_poly_setcoeff_ui(k,u,1,0);
     mpfq_2_64_poly_setcoeff_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (mpfq_2_64_poly_deg(k,b)>=0) {
      mpfq_2_64_poly_divmod(k,q,r,a,b);
      mpfq_2_64_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      mpfq_2_64_poly_set(k,b,r);
      mpfq_2_64_poly_mul(k,r,q,w);
      mpfq_2_64_poly_sub(k,r,u,r);
      mpfq_2_64_poly_set(k,u,w);   /* u,w <- w,u-qw */
      mpfq_2_64_poly_set(k,w,r);
      mpfq_2_64_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      mpfq_2_64_poly_sub(k,r,v,r);
      mpfq_2_64_poly_set(k,v,x);
      mpfq_2_64_poly_set(k,x,r);
     }
     dega=mpfq_2_64_poly_deg(k,a);
     mpfq_2_64_poly_getcoeff(k,c,a,dega);
     mpfq_2_64_inv(k,c,c);
     mpfq_2_64_poly_scal_mul(k,g,a,c);
     mpfq_2_64_poly_scal_mul(k,u0,u,c);
     mpfq_2_64_poly_scal_mul(k,v0,v,c);
     mpfq_2_64_poly_clear(k,a);
     mpfq_2_64_poly_clear(k,b);
     mpfq_2_64_poly_clear(k,u);
     mpfq_2_64_poly_clear(k,v);
     mpfq_2_64_poly_clear(k,w);
     mpfq_2_64_poly_clear(k,x);
     mpfq_2_64_poly_clear(k,q);
     mpfq_2_64_poly_clear(k,r);
    }
    mpfq_2_64_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_cmp */
static inline
int mpfq_2_64_poly_cmp(mpfq_2_64_dst_field k MAYBE_UNUSED, mpfq_2_64_src_poly u, mpfq_2_64_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return mpfq_2_64_vec_cmp(k, u->c, v->c, u->size);
}


#endif  /* MPFQ_2_64_H_ */

/* vim:set ft=cpp: */
