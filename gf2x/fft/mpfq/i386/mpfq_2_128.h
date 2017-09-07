#ifndef MPFQ_2_128_H_
#define MPFQ_2_128_H_

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
#define MPFQ_LAST_GENERATED_TAG      2_128

/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::gf2n::field */
/* Automatically generated code for GF(2^128) */
/* Definition polynomial P = X^128 + X^7 + X^2 + X + 1 */
/* Active handler: Mpfq::gf2n::trivialities */
/* Active handler: Mpfq::gf2n::io */
/* Active handler: Mpfq::gf2n::linearops */
/* Active handler: Mpfq::gf2n::inversion */
/* Active handler: Mpfq::gf2n::reduction */
/* Active handler: Mpfq::gf2n::mul */
/* Active handler: Mpfq::defaults::poly */
/* Options used:{
   coeffs=[ 128, 7, 2, 1, 0, ],
   helper=/tmp/mpfq-cado/gf2n/helper/helper,
   n=128,
   no_gmp=1,
   output_path=i386,
   slice=4,
   table=/tmp/mpfq-cado/gf2x/wizard.table,
   tag=2_128,
   w=32,
   } */

typedef mpfq_2_field mpfq_2_128_field;
typedef mpfq_2_dst_field mpfq_2_128_dst_field;

typedef unsigned long mpfq_2_128_elt[4];
typedef unsigned long * mpfq_2_128_dst_elt;
typedef const unsigned long * mpfq_2_128_src_elt;

typedef unsigned long mpfq_2_128_elt_ur[8];
typedef unsigned long * mpfq_2_128_dst_elt_ur;
typedef const unsigned long * mpfq_2_128_src_elt_ur;

typedef mpfq_2_128_elt * mpfq_2_128_vec;
typedef mpfq_2_128_elt * mpfq_2_128_dst_vec;
typedef mpfq_2_128_elt * mpfq_2_128_src_vec;

typedef mpfq_2_128_elt_ur * mpfq_2_128_vec_ur;
typedef mpfq_2_128_elt_ur * mpfq_2_128_dst_vec_ur;
typedef mpfq_2_128_elt_ur * mpfq_2_128_src_vec_ur;

typedef struct {
  mpfq_2_128_vec c;
  unsigned int alloc;
  unsigned int size;
} mpfq_2_128_poly_struct;
typedef mpfq_2_128_poly_struct mpfq_2_128_poly [1];
typedef mpfq_2_128_poly_struct * mpfq_2_128_dst_poly;
typedef mpfq_2_128_poly_struct * mpfq_2_128_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name */
#define mpfq_2_128_impl_name()	"2_128"
/* *Mpfq::gf2n::field::code_for_impl_max_characteristic_bits */
#define mpfq_2_128_impl_max_characteristic_bits()	2
/* *Mpfq::gf2n::field::code_for_impl_max_degree */
#define mpfq_2_128_impl_max_degree()	1

/* Functions operating on the field structure */
/* *Mpfq::gf2n::field::code_for_field_characteristic_bits */
#define mpfq_2_128_field_characteristic_bits(f)	1
/* *Mpfq::gf2n::field::code_for_field_degree */
#define mpfq_2_128_field_degree(f)	128
static inline
void mpfq_2_128_field_init(mpfq_2_128_dst_field);
/* *Mpfq::gf2n::field::code_for_field_clear */
#define mpfq_2_128_field_clear(K)	/**/
/* *Mpfq::gf2n::field::code_for_field_specify */
#define mpfq_2_128_field_specify(k, dummy, vp)	/**/
static inline
void mpfq_2_128_field_setopt(mpfq_2_128_dst_field, unsigned long, void *);

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, Mpfq::gf2n::trivialities */
#define mpfq_2_128_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, Mpfq::gf2n::trivialities */
#define mpfq_2_128_clear(f, px)	/**/

/* Elementary assignment functions */
static inline
void mpfq_2_128_set(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_set_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, unsigned long);
static inline
void mpfq_2_128_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_elt);
static inline
unsigned long mpfq_2_128_get_ui(mpfq_2_128_dst_field, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_set_uipoly(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, unsigned long);
static inline
void mpfq_2_128_set_uipoly_wide(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, const unsigned long *, unsigned int);
static inline
unsigned long mpfq_2_128_get_uipoly(mpfq_2_128_dst_field, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_get_uipoly_wide(mpfq_2_128_dst_field, unsigned long *, mpfq_2_128_src_elt);

/* Assignment of random values */

/* Arithmetic operations on elements */
static inline
void mpfq_2_128_add(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
/* *Mpfq::gf2n::trivialities::code_for_sub */
#define mpfq_2_128_sub(K, r, s1, s2)	mpfq_2_128_add(K,r,s1,s2)
/* *Mpfq::gf2n::trivialities::code_for_neg */
#define mpfq_2_128_neg(K, r, s)	mpfq_2_128_set(K,r,s)
static inline
void mpfq_2_128_mul(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_sqr(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
/* *Mpfq::gf2n::trivialities::code_for_is_sqr */
#define mpfq_2_128_is_sqr(f, p)	1
static inline
int mpfq_2_128_sqrt(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_pow(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, unsigned long *, size_t);
/* missing powz */
/* *Mpfq::gf2n::trivialities::code_for_frobenius */
#define mpfq_2_128_frobenius(K, r, s)	mpfq_2_128_sqr(K,r,s)
static inline
void mpfq_2_128_add_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, unsigned long);
/* *Mpfq::gf2n::trivialities::code_for_sub_ui */
#define mpfq_2_128_sub_ui(K, r, s1, s2)	mpfq_2_128_add_ui(K,r,s1,s2)
static inline
void mpfq_2_128_mul_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, unsigned long);
static inline
void mpfq_2_128_add_uipoly(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, unsigned long);
/* *Mpfq::gf2n::trivialities::code_for_sub_uipoly */
#define mpfq_2_128_sub_uipoly(K, r, s1, s2)	mpfq_2_128_add_uipoly(K,r,s1,s2)
static inline
void mpfq_2_128_mul_uipoly(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt, unsigned long);
static inline
void mpfq_2_128_longshift_left(unsigned long *, const unsigned long *, int, int);
static inline
void mpfq_2_128_longaddshift_left(unsigned long *, const unsigned long *, int, int);
static inline
int mpfq_2_128_inv(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_as_solve(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_elt);
static inline
unsigned long mpfq_2_128_trace(mpfq_2_128_dst_field, mpfq_2_128_src_elt);

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, Mpfq::gf2n::trivialities */
#define mpfq_2_128_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, Mpfq::gf2n::trivialities */
#define mpfq_2_128_elt_ur_clear(f, px)	/**/
static inline
void mpfq_2_128_elt_ur_set(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt_ur);
static inline
void mpfq_2_128_elt_ur_set_elt(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_elt_ur_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur);
static inline
void mpfq_2_128_elt_ur_set_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, unsigned long);
static inline
void mpfq_2_128_elt_ur_add(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt_ur, mpfq_2_128_src_elt_ur);
/* *Mpfq::gf2n::trivialities::code_for_elt_ur_neg */
#define mpfq_2_128_elt_ur_neg(K, r, s)	mpfq_2_128_elt_ur_set(K,r,s)
/* *Mpfq::gf2n::trivialities::code_for_elt_ur_sub */
#define mpfq_2_128_elt_ur_sub(K, r, s1, s2)	mpfq_2_128_elt_ur_add(K,r,s1,s2)
static inline
void mpfq_2_128_mul_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_sqr_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_reduce(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_dst_elt_ur);

/* Comparison functions */
static inline
int mpfq_2_128_cmp(mpfq_2_128_dst_field, mpfq_2_128_src_elt, mpfq_2_128_src_elt);
static inline
int mpfq_2_128_cmp_ui(mpfq_2_128_dst_field, mpfq_2_128_src_elt, unsigned long);
static inline
int mpfq_2_128_is_zero(mpfq_2_128_dst_field, mpfq_2_128_src_elt);

/* Input/output functions */
/* missing fprint */
/* missing print */
/* missing scan */

/* Vector functions */
void mpfq_2_128_vec_init(mpfq_2_128_dst_field, mpfq_2_128_vec *, unsigned int);
void mpfq_2_128_vec_reinit(mpfq_2_128_dst_field, mpfq_2_128_vec *, unsigned int, unsigned int);
void mpfq_2_128_vec_clear(mpfq_2_128_dst_field, mpfq_2_128_vec *, unsigned int);
static inline
void mpfq_2_128_vec_set(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, unsigned int);
static inline
void mpfq_2_128_vec_setcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_elt, unsigned int);
static inline
void mpfq_2_128_vec_setcoeff_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, unsigned long, unsigned int);
static inline
void mpfq_2_128_vec_getcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_add(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_neg(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_rev(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_sub(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_scal_mul(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, mpfq_2_128_src_elt, unsigned int);
static inline
void mpfq_2_128_vec_conv(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_src_vec, unsigned int, mpfq_2_128_src_vec, unsigned int);
/* missing vec_random */
/* missing vec_random2 */
static inline
int mpfq_2_128_vec_cmp(mpfq_2_128_dst_field, mpfq_2_128_src_vec, mpfq_2_128_src_vec, unsigned int);
static inline
int mpfq_2_128_vec_is_zero(mpfq_2_128_dst_field, mpfq_2_128_src_vec, unsigned int);
static inline
mpfq_2_128_dst_vec mpfq_2_128_vec_subvec(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, int);
static inline
mpfq_2_128_src_vec mpfq_2_128_vec_subvec_const(mpfq_2_128_dst_field, mpfq_2_128_src_vec, int);
static inline
mpfq_2_128_dst_elt mpfq_2_128_vec_coeff_ptr(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, int);
static inline
mpfq_2_128_src_elt mpfq_2_128_vec_coeff_ptr_const(mpfq_2_128_dst_field, mpfq_2_128_src_vec, int);
/* missing vec_asprint */
/* missing vec_fprint */
/* missing vec_print */
/* missing vec_sscan */
/* missing vec_fscan */
/* missing vec_scan */
void mpfq_2_128_vec_ur_init(mpfq_2_128_dst_field, mpfq_2_128_vec_ur *, unsigned int);
static inline
void mpfq_2_128_vec_ur_set_zero(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_set_vec(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec, unsigned int);
void mpfq_2_128_vec_ur_reinit(mpfq_2_128_dst_field, mpfq_2_128_vec_ur *, unsigned int, unsigned int);
void mpfq_2_128_vec_ur_clear(mpfq_2_128_dst_field, mpfq_2_128_vec_ur *, unsigned int);
static inline
void mpfq_2_128_vec_ur_set(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_setcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_elt_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_getcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_elt_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_add(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_sub(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_neg(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_ur_rev(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec_ur, unsigned int);
static inline
void mpfq_2_128_vec_scal_mul_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec, mpfq_2_128_src_elt, unsigned int);
static inline
void mpfq_2_128_vec_conv_ur_n(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_conv_ur(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, mpfq_2_128_src_vec, unsigned int, mpfq_2_128_src_vec, unsigned int);
static inline
void mpfq_2_128_vec_reduce(mpfq_2_128_dst_field, mpfq_2_128_dst_vec, mpfq_2_128_dst_vec_ur, unsigned int);
static inline
mpfq_2_128_dst_vec_ur mpfq_2_128_vec_ur_subvec(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, int);
static inline
mpfq_2_128_src_vec_ur mpfq_2_128_vec_ur_subvec_const(mpfq_2_128_dst_field, mpfq_2_128_src_vec_ur, int);
static inline
mpfq_2_128_dst_elt mpfq_2_128_vec_ur_coeff_ptr(mpfq_2_128_dst_field, mpfq_2_128_dst_vec_ur, int);
static inline
mpfq_2_128_src_elt mpfq_2_128_vec_ur_coeff_ptr_const(mpfq_2_128_dst_field, mpfq_2_128_src_vec_ur, int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, Mpfq::gf2n::trivialities */
#define mpfq_2_128_vec_elt_stride(K, n)	((n)*sizeof(mpfq_2_128_elt))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, Mpfq::gf2n::trivialities */
#define mpfq_2_128_vec_ur_elt_stride(K, n)	((n)*sizeof(mpfq_2_128_elt_ur))

/* Polynomial functions */
static inline
void mpfq_2_128_poly_init(mpfq_2_128_dst_field, mpfq_2_128_poly, unsigned int);
static inline
void mpfq_2_128_poly_clear(mpfq_2_128_dst_field, mpfq_2_128_poly);
static inline
void mpfq_2_128_poly_set(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly);
void mpfq_2_128_poly_setmonic(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_setcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_elt, unsigned int);
static inline
void mpfq_2_128_poly_setcoeff_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, unsigned long, unsigned int);
static inline
void mpfq_2_128_poly_getcoeff(mpfq_2_128_dst_field, mpfq_2_128_dst_elt, mpfq_2_128_src_poly, unsigned int);
static inline
int mpfq_2_128_poly_deg(mpfq_2_128_dst_field, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_add(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_sub(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_set_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, unsigned long);
static inline
void mpfq_2_128_poly_add_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, unsigned long);
static inline
void mpfq_2_128_poly_sub_ui(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, unsigned long);
static inline
void mpfq_2_128_poly_neg(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_scal_mul(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_elt);
static inline
void mpfq_2_128_poly_mul(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
int mpfq_2_128_poly_divmod(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
void mpfq_2_128_poly_precomp_mod(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly);
void mpfq_2_128_poly_mod_pre(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_gcd(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
static inline
void mpfq_2_128_poly_xgcd(mpfq_2_128_dst_field, mpfq_2_128_dst_poly, mpfq_2_128_dst_poly, mpfq_2_128_dst_poly, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
/* missing poly_random */
/* missing poly_random2 */
static inline
int mpfq_2_128_poly_cmp(mpfq_2_128_dst_field, mpfq_2_128_src_poly, mpfq_2_128_src_poly);
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
void mpfq_2_128_field_init(mpfq_2_128_dst_field f)
{
    f->io_type=16;
}

/* *Mpfq::gf2n::field::code_for_field_setopt */
static inline
void mpfq_2_128_field_setopt(mpfq_2_128_dst_field f, unsigned long x MAYBE_UNUSED, void * y)
{
    assert(x == MPFQ_IO_TYPE);
    f->io_type=((unsigned long*)y)[0];
}

/* *Mpfq::defaults::flatdata::code_for_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_128_elt));
}

/* *Mpfq::gf2n::trivialities::code_for_set_ui */
static inline
void mpfq_2_128_set_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, unsigned long x)
{
    r[0] = x & 1UL;
    memset(r + 1, 0, sizeof(mpfq_2_128_elt) - sizeof(unsigned long));
}

/* *Mpfq::defaults::flatdata::code_for_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r)
{
    mpfq_2_128_vec_set_zero(K,(mpfq_2_128_dst_vec)r,1);
}

/* *Mpfq::gf2n::trivialities::code_for_get_ui */
static inline
unsigned long mpfq_2_128_get_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt r)
{
    return r[0] & 1UL;
}

/* *Mpfq::gf2n::trivialities::code_for_set_uipoly */
static inline
void mpfq_2_128_set_uipoly(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, unsigned long x)
{
    r[0] = x;
    memset(r + 1, 0, sizeof(mpfq_2_128_elt) - sizeof(unsigned long));
}

/* *Mpfq::gf2n::trivialities::code_for_set_uipoly_wide */
static inline
void mpfq_2_128_set_uipoly_wide(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, const unsigned long * x, unsigned int n)
{
    unsigned int i;
    for (i = 0 ; i < n && i < 4 ; i++)
        r[i] = x[i];
}

/* *Mpfq::gf2n::trivialities::code_for_get_uipoly */
static inline
unsigned long mpfq_2_128_get_uipoly(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt r)
{
    return r[0];
}

/* *Mpfq::gf2n::trivialities::code_for_get_uipoly_wide */
static inline
void mpfq_2_128_get_uipoly_wide(mpfq_2_128_dst_field K MAYBE_UNUSED, unsigned long * r, mpfq_2_128_src_elt x)
{
    unsigned int i;
    for(i = 0 ; i < 4 ; i++) r[i] = x[i];
}

/* *Mpfq::gf2n::trivialities::code_for_add */
static inline
void mpfq_2_128_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    int i;
    for(i = 0 ; i < 4 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::trivialities::code_for_mul */
static inline
void mpfq_2_128_mul(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    mpfq_2_128_elt_ur t;
    mpfq_2_128_mul_ur(K, t, s1, s2);
    mpfq_2_128_reduce(K, r, t);
}

/* *Mpfq::gf2n::trivialities::code_for_sqr */
static inline
void mpfq_2_128_sqr(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    mpfq_2_128_elt_ur t;
    mpfq_2_128_sqr_ur(K, t, s);
    mpfq_2_128_reduce(K, r, t);
}

/* *Mpfq::gf2n::linearops::code_for_sqrt */
static inline
int mpfq_2_128_sqrt(mpfq_2_128_dst_field K, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
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
    
    mpfq_2_128_elt sqrt_t ={ 0xb6db6da4UL, 0x6db6db6dUL, 0x92492492UL, 0x24924924UL, }
    ;
    
    mpfq_2_128_elt odd, even;
    mpfq_2_128_elt_ur odd_t;
    unsigned long t;
#define	EVEN_MASK	(((unsigned long)-1)/3UL)
#define	ODD_MASK	((EVEN_MASK)<<1)
    unsigned int i;
    for(i = 0 ; i < 4 ; i++) {
        even[i] = s[i] & EVEN_MASK;
    }
    for(i = 0 ; i < 2 ; i++) {
        t = even[2*i];   t |= t >> 7;
              even[i]  = shuffle_table[t & 255];
        t >>= 16; even[i] |= shuffle_table[t & 255] << 8;
        t = even[2*i+1]; t |= t >> 7;
              even[i] |= shuffle_table[t & 255] << 16;
        t >>= 16; even[i] |= shuffle_table[t & 255] << 24;
    }
    memset(even + 2, 0, 2 * sizeof(unsigned long));
    for(i = 0 ; i < 4 ; i++) {
        odd[i] = (s[i] & ODD_MASK) >> 1;
    }
    for(i = 0 ; i < 2 ; i++) {
        t = odd[2*i];   t |= t >> 7;
              odd[i]  = shuffle_table[t & 255];
        t >>= 16; odd[i] |= shuffle_table[t & 255] << 8;
        t = odd[2*i+1]; t |= t >> 7;
              odd[i] |= shuffle_table[t & 255] << 16;
        t >>= 16; odd[i] |= shuffle_table[t & 255] << 24;
    }
    memset(odd + 2, 0, 2 * sizeof(unsigned long));
    mpfq_2_128_mul_ur(K, odd_t, odd, sqrt_t);
    for(i = 0 ; i < (4+1)/2 ; i++) {
        odd_t[i] ^= even[i];
    }
    mpfq_2_128_reduce(K, r, odd_t);
    return 1;
}

/* *Mpfq::defaults::pow::code_for_pow, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_pow(mpfq_2_128_dst_field k, mpfq_2_128_dst_elt res, mpfq_2_128_src_elt r, unsigned long * x, size_t n)
{
    mpfq_2_128_elt u, a;
    long i, j, lead;     /* it is a signed type */
    unsigned long mask;
    
    /* get the correct (i,j) position of the most significant bit in x */
    for(i = ((long)n)-1; i>=0 && x[i]==0; i--)
        ;
    if (i < 0) {
        /* power zero gets 1 */
        mpfq_2_128_set_ui(k, res, 1);
        return;
    }
    j = 32 - 1;
    mask = (1UL<<j);
    for( ; (x[i]&mask)==0 ;j--, mask>>=1)
        ;
    lead = i*32+j;      /* Ensured. */
    
    mpfq_2_128_init(k, &u);
    mpfq_2_128_init(k, &a);
    mpfq_2_128_set(k, a, r);
    for( ; lead > 0; lead--) {
        if (j-- == 0) {
            i--;
            j = 32-1;
            mask = (1UL<<j);
        } else {
            mask >>= 1;
        }
        if (x[i]&mask) {
            mpfq_2_128_sqr(k, u, a);
            mpfq_2_128_mul(k, a, u, r);
        } else {
            mpfq_2_128_sqr(k, a,a);
        }
    }
    mpfq_2_128_set(k, res, a);
    mpfq_2_128_clear(k, &u);
    mpfq_2_128_clear(k, &a);
}

/* *Mpfq::gf2n::trivialities::code_for_add_ui */
static inline
void mpfq_2_128_add_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s, unsigned long x)
{
    mpfq_2_128_set(K, r, s);
    r[0] ^= x & 1UL;
}

/* *Mpfq::gf2n::trivialities::code_for_mul_ui */
static inline
void mpfq_2_128_mul_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s, unsigned long x)
{
    if (x & 1UL) {
        mpfq_2_128_set(K, r, s);
    } else {
        memset(r, 0, sizeof(mpfq_2_128_elt));
    }
}

/* *Mpfq::gf2n::trivialities::code_for_add_uipoly */
static inline
void mpfq_2_128_add_uipoly(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s, unsigned long x)
{
    mpfq_2_128_set(K, r, s);
    r[0] ^= x;
}

/* *Mpfq::gf2n::trivialities::code_for_mul_uipoly */
static inline
void mpfq_2_128_mul_uipoly(mpfq_2_128_dst_field k, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s, unsigned long x)
{
    mpfq_2_128_elt xx;
    mpfq_2_128_init(k, &xx);
    mpfq_2_128_set_uipoly(k, xx, x);
    mpfq_2_128_mul(k, r, s, xx);
    mpfq_2_128_clear(k, &xx);
}

/* *Mpfq::gf2n::inversion::code_for_inv */
/* Triggered by: inv */
static inline
void mpfq_2_128_longshift_left(unsigned long * dst, const unsigned long * src, int n, int s)
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
void mpfq_2_128_longaddshift_left(unsigned long * dst, const unsigned long * src, int n, int s)
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
int mpfq_2_128_inv(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    unsigned long a[5] = { 0x80000000UL, 0x43UL, 0x0UL, 0x0UL, 0x80000000UL, };
    unsigned long b[5];
    unsigned long u[5] = { 0, };
    unsigned long v[5] = { 1, 0, };
    unsigned long x;
    int ia, ib;
    int i,d;
    
    if (mpfq_2_128_cmp_ui(K, s, 0UL) == 0)
        return 0;
    {
        unsigned long z;
        z = s[0] << 31;
        b[0] = z;
        z = s[0] >>  1 ^ s[1] << 31;
        b[1] = z;
        z = s[1] >>  1 ^ s[2] << 31;
        b[2] = z;
        z = s[2] >>  1 ^ s[3] << 31;
        b[3] = z;
        z = s[3] >>  1;
        b[4] = z;
    }
    ib = mpfq_clzlx(b, 5);
    ia = 0;
    
    mpfq_2_128_longshift_left(b,b,5,ib);
    
    for(d = ib - ia ; ; ) {
            if (d == 0) {
                    for(i = 0 ; i < 5 ; i++) v[i] ^= u[i];
            b[0] ^= a[0]; x = b[0];
            b[1] ^= a[1]; x |= b[1];
            b[2] ^= a[2]; x |= b[2];
            b[3] ^= a[3]; x |= b[3];
            b[4] ^= a[4]; x |= b[4];
                    if (!x) { memcpy(r,u,4 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(b,5);
                    ib += t;
                    d += t;
                    mpfq_2_128_longshift_left(b,b,5,t);
            }
            for(;d > 0;) {
                    mpfq_2_128_longaddshift_left(u,v,5,d);
            a[0] ^= b[0]; x = a[0];
            a[1] ^= b[1]; x |= a[1];
            a[2] ^= b[2]; x |= a[2];
            a[3] ^= b[3]; x |= a[3];
            a[4] ^= b[4]; x |= a[4];
                    if (!x) { memcpy(r,v,4 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(a,5);
                    ia += t;
                    d -= t;
                    mpfq_2_128_longshift_left(a,a,5,t);
            } 
            if (d == 0) {
                    for(i = 0 ; i < 5 ; i++) u[i] ^= v[i];
            a[0] ^= b[0]; x = a[0];
            a[1] ^= b[1]; x |= a[1];
            a[2] ^= b[2]; x |= a[2];
            a[3] ^= b[3]; x |= a[3];
            a[4] ^= b[4]; x |= a[4];
                    if (!x) { memcpy(r,v,4 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(a,5);
                    ia += t;
                    d -= t;
                    mpfq_2_128_longshift_left(a,a,5,t);
            }
            for(;d < 0;) {
                    mpfq_2_128_longaddshift_left(v,u,5,-d);
            b[0] ^= a[0]; x = b[0];
            b[1] ^= a[1]; x |= b[1];
            b[2] ^= a[2]; x |= b[2];
            b[3] ^= a[3]; x |= b[3];
            b[4] ^= a[4]; x |= b[4];
                    if (!x) { memcpy(r,u,4 * sizeof(unsigned long)); return 1; }
                    unsigned long t = mpfq_clzlx(b,5);
                    ib += t;
                    d += t;
                    mpfq_2_128_longshift_left(b,b,5,t);
            }
    }
}

/* *Mpfq::gf2n::linearops::code_for_as_solve */
static inline
void mpfq_2_128_as_solve(mpfq_2_128_dst_field K, mpfq_2_128_dst_elt r, mpfq_2_128_src_elt s)
{
    static const mpfq_2_128_elt t[128] = {
        { 0xa4b20b08UL, 0x676aac9fUL, 0xf4731af9UL, 0x295ac0b1UL, },
        { 0x766af404UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0x766af406UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0x3404a20cUL, 0x7783f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0x766af402UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0x83074e6aUL, 0x7dc714cbUL, 0xd7da34b0UL, 0x2aaea79aUL, },
        { 0x3404a204UL, 0x7783f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0xd2d9fe18UL, 0xf5f3f531UL, 0x1ae41171UL, 0xcc809524UL, },
        { 0x766af412UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0x8ac85b44UL, 0x77cb862eUL, 0xd159ce96UL, 0x4e983687UL, },
        { 0x83074e4aUL, 0x7dc714cbUL, 0xd7da34b0UL, 0x2aaea79aUL, },
        { 0x6aaff56UL, 0xf73fcb3cUL, 0xaf144aafUL, 0xaae2ab6aUL, },
        { 0x3404a244UL, 0x7783f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0x40603f62UL, 0x893fe9caUL, 0xcd5792ddUL, 0xe62f327cUL, },
        { 0xd2d9fe98UL, 0xf5f3f531UL, 0x1ae41171UL, 0xcc809524UL, },
        { 0xd453b626UL, 0xe1336bc7UL, 0x3f47c90bUL, 0xab52f032UL, },
        { 0x766af512UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0xabf4573aUL, 0x79a6a089UL, 0xd75338b4UL, 0x2abe34baUL, },
        { 0x8ac85944UL, 0x77cb862eUL, 0xd159ce96UL, 0x4e983687UL, },
        { 0x232562fcUL, 0x1a9cbef8UL, 0xdcc72f2eUL, 0x1e8ccb9UL, },
        { 0x83074a4aUL, 0x7dc714cbUL, 0xd7da34b0UL, 0x2aaea79aUL, },
        { 0x4a80e48cUL, 0xce869b2UL, 0x21f3dc1aUL, 0x2442c32UL, },
        { 0x6aaf756UL, 0xf73fcb3cUL, 0xaf144aafUL, 0xaae2ab6aUL, },
        { 0xdd9aa228UL, 0xeb3ff936UL, 0x39c4323dUL, 0xcf64612fUL, },
        { 0x3404b244UL, 0x7783f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0x57572458UL, 0xf0d5382aUL, 0x880d4909UL, 0xcc911696UL, },
        { 0x40601f62UL, 0x893fe9caUL, 0xcd5792ddUL, 0xe62f327cUL, },
        { 0xe5db2834UL, 0x697f6478UL, 0x29d28cf2UL, 0x4d7cd7ccUL, },
        { 0xd2d9be98UL, 0xf5f3f531UL, 0x1ae41171UL, 0xcc809524UL, },
        { 0x1a5b86a4UL, 0xfc0956caUL, 0x569146a7UL, 0xccc99127UL, },
        { 0xd4533626UL, 0xe1336bc7UL, 0x3f47c90bUL, 0xab52f032UL, },
        { 0x1eeb82UL, 0x14f5986bUL, 0xdacd5b1cUL, 0x65cece94UL, },
        { 0x766bf512UL, 0x929959afUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0xc7c3e73aUL, 0x6fe27169UL, 0xd509ca74UL, 0x2b0e7fa2UL, },
        { 0xabf6573aUL, 0x79a6a089UL, 0xd75338b4UL, 0x2abe34baUL, },
        { 0x85980142UL, 0x11d35538UL, 0x48241324UL, 0x65df4d26UL, },
        { 0x8acc5944UL, 0x77cb862eUL, 0xd159ce96UL, 0x4e983687UL, },
        { 0x2af671a4UL, 0xef6fd21dUL, 0x70ef2af7UL, 0xcf74f70eUL, },
        { 0x232d62fcUL, 0x1a9cbef8UL, 0xdcc72f2eUL, 0x1e8ccb9UL, },
        { 0x99490f58UL, 0xf91a9d84UL, 0x3b16cdf9UL, 0xcec4b917UL, },
        { 0x83174a4aUL, 0x7dc714cbUL, 0xd7da34b0UL, 0x2aaea79aUL, },
        { 0x6b513bdcUL, 0x99a6cf20UL, 0x175f6247UL, 0x83f16f99UL, },
        { 0x4aa0e48cUL, 0xce869b2UL, 0x21f3dc1aUL, 0x2442c32UL, },
        { 0xe84671c0UL, 0x1b9e27ebUL, 0x6a2aa908UL, 0x3f562a9UL, },
        { 0x6eaf756UL, 0xf73fcb3cUL, 0xaf144aafUL, 0xaae2ab6aUL, },
        { 0xcddeebd4UL, 0xfcd6d526UL, 0xe11bcce3UL, 0xccc9c07aUL, },
        { 0xdd1aa228UL, 0xeb3ff936UL, 0x39c4323dUL, 0xcf64612fUL, },
        { 0xe8467196UL, 0xe461dfebUL, 0xc09556f7UL, 0xa95fc803UL, },
        { 0x3504b244UL, 0x7783f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0xc15c7384UL, 0x7b20ee03UL, 0xf0ab5b2eUL, 0x4cdc1ab0UL, },
        { 0x55572458UL, 0xf0d5382aUL, 0x880d4909UL, 0xcc911696UL, },
        { 0xfe68cf6aUL, 0xfb48c8c7UL, 0x1da11d79UL, 0xaaba6896UL, },
        { 0x44601f62UL, 0x893fe9caUL, 0xcd5792ddUL, 0xe62f327cUL, },
        { 0x638d3beeUL, 0x938b5c70UL, 0xee209847UL, 0xe5db5416UL, },
        { 0xeddb2834UL, 0x697f6478UL, 0x29d28cf2UL, 0x4d7cd7ccUL, },
        { 0xa7e512ceUL, 0xb49fa26UL, 0x95e2e76eUL, 0x6437818bUL, },
        { 0xc2d9be98UL, 0xf5f3f531UL, 0x1ae41171UL, 0xcc809524UL, },
        { 0x22ddf34UL, 0x7868a299UL, 0x46184ba2UL, 0x4cd90207UL, },
        { 0x3a5b86a4UL, 0xfc0956caUL, 0x569146a7UL, 0xccc99127UL, },
        { 0x7e8e76cUL, 0x8f1d24e4UL, 0xefbe4a8bUL, 0x824171dcUL, },
        { 0x94533626UL, 0xe1336bc7UL, 0x3f47c90bUL, 0xab52f032UL, },
        { 0x2cdc3b5eUL, 0x88711e20UL, 0x1f6bf32bUL, 0xe62e22afUL, },
        { 0x801eeb82UL, 0x14f5986bUL, 0xdacd5b1cUL, 0x65cece94UL, },
        { 0x27b54520UL, 0x9aadb854UL, 0x23a92e49UL, 0x83f4672bUL, },
        { 0x766bf512UL, 0x929959aeUL, 0xee970b89UL, 0xe5da5595UL, },
        { 0xbeccf94cUL, 0x48711dUL, 0x92097172UL, 0x1040UL, },
        { 0xc7c3e73aUL, 0x6fe2716bUL, 0xd509ca74UL, 0x2b0e7fa2UL, },
        { 0xd0d7976cUL, 0x99d6b267UL, 0x7a7437cbUL, 0x81edd40aUL, },
        { 0xabf6573aUL, 0x79a6a08dUL, 0xd75338b4UL, 0x2abe34baUL, },
        { 0x1a5d4462UL, 0x84dd884dUL, 0xeccdf96dUL, 0xe46a1e8dUL, },
        { 0x85980142UL, 0x11d35530UL, 0x48241324UL, 0x65df4d26UL, },
        { 0x3a605920UL, 0x8ab1aeebUL, 0xeac70f69UL, 0x804c1cb1UL, },
        { 0x8acc5944UL, 0x77cb863eUL, 0xd159ce96UL, 0x4e983687UL, },
        { 0x7dbd55aaUL, 0xe0451437UL, 0x52499e01UL, 0xa94f4b32UL, },
        { 0x2af671a4UL, 0xef6fd23dUL, 0x70ef2af7UL, 0xcf74f70eUL, },
        { 0xf56f619eUL, 0x978a93f3UL, 0x831189d3UL, 0xe7d77da5UL, },
        { 0x232d62fcUL, 0x1a9cbeb8UL, 0xdcc72f2eUL, 0x1e8ccb9UL, },
        { 0xb0a16ef4UL, 0x7a5109eaUL, 0x2b68e52aUL, 0x4cdd0f36UL, },
        { 0x99490f58UL, 0xf91a9d04UL, 0x3b16cdf9UL, 0xcec4b917UL, },
        { 0x88d05912UL, 0x8830781eUL, 0x7bf23369UL, 0xe4329c3dUL, },
        { 0x83174a4aUL, 0x7dc715cbUL, 0xd7da34b0UL, 0x2aaea79aUL, },
        { 0xcff334e0UL, 0x8618bcbfUL, 0x59738741UL, 0x82082082UL, },
        { 0x6b513bdcUL, 0x99a6cd20UL, 0x175f6247UL, 0x83f16f99UL, },
        { 0xf4adcc3aUL, 0x6297e291UL, 0xd1181a84UL, 0x294b124eUL, },
        { 0x4aa0e48cUL, 0xce86db2UL, 0x21f3dc1aUL, 0x2442c32UL, },
        { 0x6550119eUL, 0x9b696b66UL, 0x14466493UL, 0xe78fff41UL, },
        { 0xe84671c0UL, 0x1b9e2febUL, 0x6a2aa908UL, 0x3f562a9UL, },
        { 0x63d55a14UL, 0xeb568d9eUL, 0x54a1436bUL, 0xcd78dbfdUL, },
        { 0x6eaf756UL, 0xf73fdb3cUL, 0xaf144aafUL, 0xaae2ab6aUL, },
        { 0xed722b64UL, 0xfacabbbeUL, 0x1ce38d0fUL, 0xccdd4e7bUL, },
        { 0xcddeebd4UL, 0xfcd6f526UL, 0xe11bcce3UL, 0xccc9c07aUL, },
        { 0xa43bbf3cUL, 0x1fb4754fUL, 0xb1bf89c6UL, 0x3e5e1ccUL, },
        { 0xdd1aa228UL, 0xeb3fb936UL, 0x39c4323dUL, 0xcf64612fUL, },
        { 0x172caea2UL, 0x926f68f4UL, 0x3461b455UL, 0xe5da0089UL, },
        { 0xe8467196UL, 0xe4615febUL, 0xc09556f7UL, 0xa95fc803UL, },
        { 0xcb9d1fd2UL, 0x981cf1dcUL, 0x7b3e87b1UL, 0xe78ae2a7UL, },
        { 0x3504b244UL, 0x7782f733UL, 0x4350bfe6UL, 0x4e9826c7UL, },
        { 0xc15f72c6UL, 0x6fd777b6UL, 0x2a675b0aUL, 0x2912d422UL, },
        { 0xc15c7384UL, 0x7b22ee03UL, 0xf0ab5b2eUL, 0x4cdc1ab0UL, },
        { 0x85956fd2UL, 0x96f8740fUL, 0x58d237f1UL, 0xe7d6682eUL, },
        { 0x55572458UL, 0xf0d1382aUL, 0x880d4909UL, 0xcc911696UL, },
        { 0x33474626UL, 0xe46143c9UL, 0x764a0b4bUL, 0xab4366b0UL, },
        { 0xfe68cf6aUL, 0xfb40c8c7UL, 0x1da11d79UL, 0xaaba6896UL, },
        { 0xf7de0d9eUL, 0x97849bb7UL, 0x8378aeb3UL, 0xe7d77d71UL, },
        { 0x44601f62UL, 0x892fe9caUL, 0xcd5792ddUL, 0xe62f327cUL, },
        { 0x382e380eUL, 0xd13f9bdUL, 0xfb5729d6UL, 0x64231d9eUL, },
        { 0x638d3beeUL, 0x93ab5c70UL, 0xee209847UL, 0xe5db5416UL, },
        { 0x9a557220UL, 0x959811bdUL, 0xa5fe1499UL, 0x81b51ab7UL, },
        { 0xeddb2834UL, 0x693f6478UL, 0x29d28cf2UL, 0x4d7cd7ccUL, },
        { 0x75a08406UL, 0x678a4ee9UL, 0x43072c12UL, 0x295ab4bbUL, },
        { 0xa7e512ceUL, 0xbc9fa26UL, 0x95e2e76eUL, 0x6437818bUL, },
        { 0x25989a26UL, 0xef9e1acdUL, 0x8eb165ebUL, 0xab0eead3UL, },
        { 0xc2d9be98UL, 0xf4f3f531UL, 0x1ae41171UL, 0xcc809524UL, },
        { 0x29b00000UL, 0x10185b9bUL, 0xdb230c00UL, 0x3ba7996UL, },
        { 0x22ddf34UL, 0x7a68a299UL, 0x46184ba2UL, 0x4cd90207UL, },
        { 0x38765990UL, 0x8261f453UL, 0x10890d05UL, 0x82109320UL, },
        { 0x3a5b86a4UL, 0xf80956caUL, 0x569146a7UL, 0xccc99127UL, },
        { 0x3db361faUL, 0x7f14722eUL, 0xb92f0c2cUL, 0x2a88e0fbUL, },
        { 0x7e8e76cUL, 0x871d24e4UL, 0xefbe4a8bUL, 0x824171dcUL, },
        { 0x7d17578aUL, 0x6fb067f4UL, 0xb8989fd8UL, 0x29b11e87UL, },
        { 0x94533626UL, 0xf1336bc7UL, 0x3f47c90bUL, 0xab52f032UL, },
        { 0x3891e6faUL, 0xdb7ed8cUL, 0xfae1613cUL, 0x38b21c09UL, },
        { 0x2cdc3b5eUL, 0xa8711e20UL, 0x1f6bf32bUL, 0xe62e22afUL, },
        { 0x2cdc3b6cUL, 0xa8711e20UL, 0x1f6bf32bUL, 0x862e22afUL, },
        { 0x801eeb82UL, 0x54f5986bUL, 0xdacd5b1cUL, 0x65cece94UL, },
        { 0xc270bdacUL, 0xb1ef36f7UL, 0x770aef73UL, 0x8e8cbdc6UL, },
        { 0x27b54520UL, 0x1aadb854UL, 0x23a92e49UL, 0x83f4672bUL, },
        { 0x0UL, 0x0UL, 0x0UL, 0x0UL, },
    };
    const mpfq_2_128_elt * ptr = t;
    unsigned int i,j;
    memset(r, 0, sizeof(mpfq_2_128_elt));
    for(i = 0 ; i < 4 ; i++) {
        unsigned long a = s[i];
        for(j = 0 ; j < 32 && ptr != &t[128]; j++, ptr++) {
            if (a & 1UL) {
                mpfq_2_128_add(K, r, r, *ptr);
            }
            a >>= 1;
        }
    }
}

/* *Mpfq::gf2n::linearops::code_for_trace */
static inline
unsigned long mpfq_2_128_trace(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt s)
{
    return ((s[3]>>31) ^ (s[3]>>25)) & 1;
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt_ur s)
{
    if (r != s) memcpy(r,s,sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set_elt(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt s)
{
    memset(r, 0, sizeof(mpfq_2_128_elt_ur)); memcpy(r,s,sizeof(mpfq_2_128_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_elt_ur_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r)
{
    memset(r, 0, sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::gf2n::trivialities::code_for_elt_ur_set_ui */
static inline
void mpfq_2_128_elt_ur_set_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, unsigned long x)
{
    r[0] = x & 1UL;
    memset(r + 1, 0, sizeof(mpfq_2_128_elt_ur) - sizeof(unsigned long));
}

/* *Mpfq::gf2n::trivialities::code_for_elt_ur_add */
static inline
void mpfq_2_128_elt_ur_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur r, mpfq_2_128_src_elt_ur s1, mpfq_2_128_src_elt_ur s2)
{
    int i;
    for(i = 0 ; i < 8 ; i++)
        r[i] = s1[i] ^ s2[i];
}

/* *Mpfq::gf2n::mul::code_for_mul_ur */
static inline
void mpfq_2_128_mul_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur t, mpfq_2_128_src_elt s1, mpfq_2_128_src_elt s2)
{
    gf2x_mul4(t, s1, s2);
}

/* *Mpfq::gf2n::squaring::code_for_sqr_ur */
static inline
void mpfq_2_128_sqr_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur t, mpfq_2_128_src_elt s)
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
        u = g[s[2]       & 15];
    t[4]  = u;
        u = g[s[2] >>  4 & 15];
    t[4] ^= u <<  8;
        u = g[s[2] >>  8 & 15];
    t[4] ^= u << 16;
        u = g[s[2] >> 12 & 15];
    t[4] ^= u << 24;
        u = g[s[2] >> 16 & 15];
    t[5]  = u;
        u = g[s[2] >> 20 & 15];
    t[5] ^= u <<  8;
        u = g[s[2] >> 24 & 15];
    t[5] ^= u << 16;
        u = g[s[2] >> 28 & 15];
    t[5] ^= u << 24;
        u = g[s[3]       & 15];
    t[6]  = u;
        u = g[s[3] >>  4 & 15];
    t[6] ^= u <<  8;
        u = g[s[3] >>  8 & 15];
    t[6] ^= u << 16;
        u = g[s[3] >> 12 & 15];
    t[6] ^= u << 24;
        u = g[s[3] >> 16 & 15];
    t[7]  = u;
        u = g[s[3] >> 20 & 15];
    t[7] ^= u <<  8;
        u = g[s[3] >> 24 & 15];
    t[7] ^= u << 16;
        u = g[s[3] >> 28 & 15];
    t[7] ^= u << 24;
    }
}

/* *Mpfq::gf2n::reduction::code_for_reduce */
static inline
void mpfq_2_128_reduce(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt r, mpfq_2_128_dst_elt_ur t)
{
    {
        unsigned long s[5];
        /* 127 excess bits */
        {
            unsigned long z;
            z = t[0];
            s[0] = z;
            z = t[1];
            s[1] = z;
            z = t[2];
            s[2] = z;
            z = t[3];
            s[3] = z;
        }
        memset(s + 4, 0, 1 * sizeof(unsigned long));
        {
            unsigned long z;
            z = t[4];
            s[0]^= z <<  7;
            s[0]^= z <<  2;
            s[0]^= z <<  1;
            s[0]^= z;
            z >>= 25;
            z^= t[5] <<  7;
            s[1]^= z;
            z >>= 5;
            z^= t[5] >> 25 << 27;
            s[1]^= z;
            z >>= 1;
            z^= t[5] >> 30 << 31;
            s[1]^= z;
            z >>= 1;
            z^= (t[5] & ~0x7fffffffUL);
            s[1]^= z;
            z >>= 25;
            z^= t[6] <<  7;
            s[2]^= z;
            z >>= 5;
            z^= t[6] >> 25 << 27;
            s[2]^= z;
            z >>= 1;
            z^= t[6] >> 30 << 31;
            s[2]^= z;
            z >>= 1;
            z^= (t[6] & ~0x7fffffffUL);
            s[2]^= z;
            z >>= 25;
            z^= t[7] <<  7;
            s[3]^= z;
            z >>= 5;
            z^= t[7] >> 25 << 27;
            s[3]^= z;
            z >>= 1;
            z^= t[7] >> 30 << 31;
            s[3]^= z;
            z >>= 1;
            s[3]^= z;
            z >>= 25;
            s[4]^= z;
            z >>= 5;
            s[4]^= z;
        }
        /* 6 excess bits */
        {
            unsigned long z;
            z = s[0];
            r[0] = z;
            z = s[1];
            r[1] = z;
            z = s[2];
            r[2] = z;
            z = s[3];
            r[3] = z;
        }
        {
            unsigned long z;
            z = s[4];
            r[0]^= z <<  7;
            r[0]^= z <<  2;
            r[0]^= z <<  1;
            r[0]^= z;
        }
    }
}

/* *Mpfq::gf2n::trivialities::code_for_cmp */
static inline
int mpfq_2_128_cmp(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt a, mpfq_2_128_src_elt b)
{
    return memcmp(a, b, sizeof(mpfq_2_128_elt));
}

/* *Mpfq::gf2n::trivialities::code_for_cmp_ui */
static inline
int mpfq_2_128_cmp_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt r, unsigned long x)
{
    int i;
    if (r[0] < (x & 1UL)) return -1;
    if (r[0] > (x & 1UL)) return 1;
    for(i = 1 ; i < 4 ; i++) {
        if (r[i]) return 1;
    }
    return 0;
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gf2n::trivialities */
static inline
int mpfq_2_128_is_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(mpfq_2_128_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_vec_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec r, mpfq_2_128_src_vec s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_2_128_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_vec_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(mpfq_2_128_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_setcoeff(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_elt x, unsigned int i)
{
    mpfq_2_128_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoeff_ui, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_setcoeff_ui(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, unsigned long x, unsigned int i)
{
    mpfq_2_128_set_ui(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_getcoeff(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt x, mpfq_2_128_src_vec w, unsigned int i)
{
    mpfq_2_128_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, mpfq_2_128_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_neg(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_128_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_rev(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    mpfq_2_128_elt tmp[1];
    mpfq_2_128_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_2_128_set(K, tmp[0], u[i]);
        mpfq_2_128_set(K, w[i], u[n-1-i]);
        mpfq_2_128_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_2_128_set(K, w[nn], u[nn]);
    mpfq_2_128_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_sub(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, mpfq_2_128_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_128_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_scal_mul(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, mpfq_2_128_src_elt x, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_mul(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv */
static inline
void mpfq_2_128_vec_conv(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_src_vec u, unsigned int n, mpfq_2_128_src_vec v, unsigned int m)
{
    mpfq_2_128_vec_ur tmp;
    mpfq_2_128_vec_ur_init(K, &tmp, m+n-1);
    mpfq_2_128_vec_conv_ur(K, tmp, u, n, v, m);
    mpfq_2_128_vec_reduce(K, w, tmp, m+n-1);
    mpfq_2_128_vec_ur_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec */
static inline
int mpfq_2_128_vec_cmp(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec u, mpfq_2_128_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = mpfq_2_128_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_is_zero, Mpfq::defaults::vec */
static inline
int mpfq_2_128_vec_is_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!mpfq_2_128_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec */
static inline
mpfq_2_128_dst_vec mpfq_2_128_vec_subvec(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_2_128_src_vec mpfq_2_128_vec_subvec_const(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_2_128_dst_elt mpfq_2_128_vec_coeff_ptr(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_2_128_src_elt mpfq_2_128_vec_coeff_ptr_const(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_vec_ur_set_zero(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur r, unsigned int n)
{
    memset(r, 0, n*sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_set_vec, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_set_vec(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_elt_ur_set_elt(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, Mpfq::gf2n::trivialities */
static inline
void mpfq_2_128_vec_ur_set(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur r, mpfq_2_128_src_vec_ur s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_2_128_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_setcoeff(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_elt_ur x, unsigned int i)
{
    mpfq_2_128_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_getcoeff(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_elt_ur x, mpfq_2_128_src_vec_ur w, unsigned int i)
{
    mpfq_2_128_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_add(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec_ur u, mpfq_2_128_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_sub(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec_ur u, mpfq_2_128_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_neg(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_2_128_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_ur_rev(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec_ur u, unsigned int n)
{
    unsigned int nn = n >> 1;
    mpfq_2_128_elt_ur tmp[1];
    mpfq_2_128_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_2_128_elt_ur_set(K, tmp[0], u[i]);
        mpfq_2_128_elt_ur_set(K, w[i], u[n-1-i]);
        mpfq_2_128_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_2_128_elt_ur_set(K, w[nn], u[nn]);
    mpfq_2_128_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_scal_mul_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec u, mpfq_2_128_src_elt x, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_mul_ur(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur */
/* Triggered by: vec_conv_ur */
static inline
void mpfq_2_128_vec_conv_ur_n(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec u, mpfq_2_128_src_vec v, unsigned int n)
{
    if (n == 0)
        return;
    if (n == 1) {
        mpfq_2_128_mul_ur(K, w[0], u[0], v[0]);
        return;
    }
    if (n == 2) {  // Kara 2
        mpfq_2_128_elt t1, t2;
        mpfq_2_128_init(K, &t1);
        mpfq_2_128_init(K, &t2);
        mpfq_2_128_mul_ur(K, w[0], u[0], v[0]);
        mpfq_2_128_mul_ur(K, w[2], u[1], v[1]);
        mpfq_2_128_add(K, t1, u[0], u[1]);
        mpfq_2_128_add(K, t2, v[0], v[1]);
        mpfq_2_128_mul_ur(K, w[1], t1, t2);
        mpfq_2_128_elt_ur_sub(K, w[1], w[1], w[0]);
        mpfq_2_128_elt_ur_sub(K, w[1], w[1], w[2]);
        mpfq_2_128_clear(K, &t1);
        mpfq_2_128_clear(K, &t2);
        return;
    }
    if (n == 3) {  // do it in 6
        mpfq_2_128_elt t1, t2;
        mpfq_2_128_elt_ur s;
        mpfq_2_128_init(K, &t1);
        mpfq_2_128_init(K, &t2);
        mpfq_2_128_elt_ur_init(K, &s);
        // a0*b0*(1 - X)
        mpfq_2_128_mul_ur(K, w[0], u[0], v[0]);
        mpfq_2_128_elt_ur_neg(K, w[1], w[0]);
        // a1*b1*(-X + 2*X^2 - X^3)
        mpfq_2_128_mul_ur(K, w[2], u[1], v[1]);
        mpfq_2_128_elt_ur_neg(K, w[3], w[2]);
        mpfq_2_128_elt_ur_add(K, w[2], w[2], w[2]);
        mpfq_2_128_elt_ur_add(K, w[1], w[1], w[3]);
        // a2*b2*(-X^3+X^4)
        mpfq_2_128_mul_ur(K, w[4], u[2], v[2]);
        mpfq_2_128_elt_ur_sub(K, w[3], w[3], w[4]);
        // (a0+a1)*(b0+b1)*(X - X^2)
        mpfq_2_128_add(K, t1, u[0], u[1]);
        mpfq_2_128_add(K, t2, v[0], v[1]);
        mpfq_2_128_mul_ur(K, s, t1, t2);
        mpfq_2_128_elt_ur_add(K, w[1], w[1], s);
        mpfq_2_128_elt_ur_sub(K, w[2], w[2], s);
        // (a1+a2)*(b1+b2)*(X^3 - X^2)
        mpfq_2_128_add(K, t1, u[1], u[2]);
        mpfq_2_128_add(K, t2, v[1], v[2]);
        mpfq_2_128_mul_ur(K, s, t1, t2);
        mpfq_2_128_elt_ur_add(K, w[3], w[3], s);
        mpfq_2_128_elt_ur_sub(K, w[2], w[2], s);
        // (a0+a1+a2)*(b0+b1+b2)* X^2
        mpfq_2_128_add(K, t1, u[0], t1);
        mpfq_2_128_add(K, t2, v[0], t2);
        mpfq_2_128_mul_ur(K, s, t1, t2);
        mpfq_2_128_elt_ur_add(K, w[2], w[2], s);
        return;
    }
    unsigned int n0, n1;
    n0 = n / 2;
    n1 = n - n0;
    mpfq_2_128_vec_conv_ur_n(K, w, u, v, n0);
    mpfq_2_128_vec_conv_ur_n(K, w + 2*n0, u + n0, v + n0, n1);
    mpfq_2_128_elt_ur_set_ui(K, w[2*n0-1], 0);
    
    mpfq_2_128_vec tmpu, tmpv;
    mpfq_2_128_vec_ur tmpw;
    mpfq_2_128_vec_init(K, &tmpu, n1);
    mpfq_2_128_vec_init(K, &tmpv, n1);
    mpfq_2_128_vec_ur_init(K, &tmpw, 2*n1-1);
    
    mpfq_2_128_vec_set(K, tmpu, u, n0);
    if (n1 != n0) 
        mpfq_2_128_set_ui(K, tmpu[n0], 0);
    mpfq_2_128_vec_add(K, tmpu, tmpu, u+n0, n1);
    mpfq_2_128_vec_set(K, tmpv, v, n0);
    if (n1 != n0) 
        mpfq_2_128_set_ui(K, tmpv[n0], 0);
    mpfq_2_128_vec_add(K, tmpv, tmpv, v+n0, n1);
    mpfq_2_128_vec_conv_ur_n(K, tmpw, tmpu, tmpv, n1);
    mpfq_2_128_vec_ur_sub(K, tmpw, tmpw, w, 2*n0-1);
    mpfq_2_128_vec_ur_sub(K, tmpw, tmpw, w + 2*n0, 2*n1-1);
    mpfq_2_128_vec_ur_add(K, w + n0, w + n0, tmpw, 2*n1-1);
    
    mpfq_2_128_vec_clear(K, &tmpu, n1);
    mpfq_2_128_vec_clear(K, &tmpv, n1);
    mpfq_2_128_vec_ur_clear(K, &tmpw, 2*n1-1);
    return;
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur */
static inline
void mpfq_2_128_vec_conv_ur(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur w, mpfq_2_128_src_vec u, unsigned int n, mpfq_2_128_src_vec v, unsigned int m)
{
    unsigned int i, j MAYBE_UNUSED, k;
    mpfq_2_128_elt_ur acc, z;
    mpfq_2_128_elt_ur_init(K, &acc);
    mpfq_2_128_elt_ur_init(K, &z);
    // swap pointers to have n <= m
    mpfq_2_128_src_vec uu, vv;
    if (n <= m) {
        uu = u; vv = v;
    } else {
        uu = v; vv = u;
        unsigned int tmp = n;
        n = m; m = tmp;
    }
    for(k = 0; k < n; ++k) {
        mpfq_2_128_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i <= k; ++i) {
            mpfq_2_128_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_128_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_128_elt_ur_set(K, w[k], acc);
    }
    for(k = n; k < m; ++k) {
        mpfq_2_128_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i < n; ++i) {
            mpfq_2_128_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_128_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_128_elt_ur_set(K, w[k], acc);
    }
    for(k = m; k < n+m-1; ++k) {
        mpfq_2_128_mul_ur(K, acc, uu[k-m+1], vv[m-1]);
        for(i = k-m+2; i < n; ++i) {
            mpfq_2_128_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_2_128_elt_ur_add(K, acc, acc, z);
        }
        mpfq_2_128_elt_ur_set(K, w[k], acc);
    }
    mpfq_2_128_elt_ur_clear(K, &acc);
    mpfq_2_128_elt_ur_clear(K, &z);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec */
static inline
void mpfq_2_128_vec_reduce(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec w, mpfq_2_128_dst_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_128_reduce(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec */
static inline
mpfq_2_128_dst_vec_ur mpfq_2_128_vec_ur_subvec(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_2_128_src_vec_ur mpfq_2_128_vec_ur_subvec_const(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_2_128_dst_elt mpfq_2_128_vec_ur_coeff_ptr(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_dst_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_2_128_src_elt mpfq_2_128_vec_ur_coeff_ptr_const(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::poly::code_for_poly_init */
static inline
void mpfq_2_128_poly_init(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_poly p, unsigned int n)
{
    mpfq_2_128_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear */
static inline
void mpfq_2_128_poly_clear(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_poly p)
{
    mpfq_2_128_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set */
static inline
void mpfq_2_128_poly_set(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_2_128_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff */
static inline
void mpfq_2_128_poly_setcoeff(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_elt x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_2_128_vec_set_zero(k, mpfq_2_128_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_2_128_vec_setcoeff(k, w->c, x, i);
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff_ui */
static inline
void mpfq_2_128_poly_setcoeff_ui(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, unsigned long x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_2_128_vec_set_zero(k, mpfq_2_128_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_2_128_vec_setcoeff_ui(k, w->c, x, i);
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_getcoeff */
static inline
void mpfq_2_128_poly_getcoeff(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_elt x, mpfq_2_128_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       mpfq_2_128_set_ui(k,x,0);
    } else {
       mpfq_2_128_vec_getcoeff(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg */
static inline
int mpfq_2_128_poly_deg(mpfq_2_128_dst_field K MAYBE_UNUSED, mpfq_2_128_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    mpfq_2_128_elt temp;
    mpfq_2_128_init(K, &temp);
    mpfq_2_128_vec_getcoeff(K, temp, w->c, deg);
    int comp=mpfq_2_128_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           mpfq_2_128_vec_getcoeff(K, temp, w->c, deg);
           comp=mpfq_2_128_cmp_ui(K, temp, 0);
        }
    }
    mpfq_2_128_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add */
static inline
void mpfq_2_128_poly_add(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, mpfq_2_128_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_2_128_vec_add(k, w->c, u->c, v->c, su);
        mpfq_2_128_vec_set(k, mpfq_2_128_vec_subvec(k, w->c, su), mpfq_2_128_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_2_128_vec_add(k, w->c, u->c, v->c, sv);
        mpfq_2_128_vec_set(k, mpfq_2_128_vec_subvec(k, w->c, sv), mpfq_2_128_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_sub */
static inline
void mpfq_2_128_poly_sub(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, mpfq_2_128_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_2_128_vec_sub(k, w->c, u->c, v->c, su);
        mpfq_2_128_vec_neg(k, mpfq_2_128_vec_subvec(k, w->c, su), mpfq_2_128_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_2_128_vec_sub(k, w->c, u->c, v->c, sv);
        mpfq_2_128_vec_set(k, mpfq_2_128_vec_subvec(k, w->c, sv), mpfq_2_128_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_set_ui */
static inline
void mpfq_2_128_poly_set_ui(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, unsigned long x)
{
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_128_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui */
static inline
void mpfq_2_128_poly_add_ui(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_128_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_2_128_poly_deg(k, w);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_2_128_vec_set(k, mpfq_2_128_vec_subvec(k, w->c, 1), mpfq_2_128_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_2_128_add_ui(k, mpfq_2_128_vec_coeff_ptr(k, w->c, 0), mpfq_2_128_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui */
static inline
void mpfq_2_128_poly_sub_ui(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_2_128_elt temp;
        mpfq_2_128_init(k, &temp);
        mpfq_2_128_set_ui(k, temp, x);
        mpfq_2_128_neg(k, mpfq_2_128_vec_coeff_ptr(k, w->c, 0), temp);
        w->size = mpfq_2_128_cmp_ui(k, temp, 0);
        mpfq_2_128_clear(k, &temp);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_2_128_vec_set(k, mpfq_2_128_vec_subvec(k, w->c, 1), mpfq_2_128_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_2_128_sub_ui(k, mpfq_2_128_vec_coeff_ptr(k, w->c, 0), mpfq_2_128_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_neg */
static inline
void mpfq_2_128_poly_neg(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_2_128_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul */
static inline
void mpfq_2_128_poly_scal_mul(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, mpfq_2_128_src_elt x)
{
    if (mpfq_2_128_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_2_128_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_mul */
static inline
void mpfq_2_128_poly_mul(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly w, mpfq_2_128_src_poly u, mpfq_2_128_src_poly v)
{
    unsigned int usize = mpfq_2_128_poly_deg(k, u)+1;
    unsigned int vsize = mpfq_2_128_poly_deg(k, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        mpfq_2_128_vec_reinit(k, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    mpfq_2_128_vec_conv(k, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
    w->size = 1 + mpfq_2_128_poly_deg(k, w);
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly */
static inline
void mpfq_2_128_poly_gcd(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly g, mpfq_2_128_src_poly a0, mpfq_2_128_src_poly b0)
{
    mpfq_2_128_poly a,b,q,r;
    int da0=mpfq_2_128_poly_deg(k,a0), db0=mpfq_2_128_poly_deg(k,b0);
    if (db0==-1)
     mpfq_2_128_poly_set(k,g,a0);
    else {
     mpfq_2_128_poly_init(k,a,da0+1);
     mpfq_2_128_poly_init(k,b,db0+1);
     mpfq_2_128_poly_init(k,q,1);
     mpfq_2_128_poly_init(k,r,db0);
     mpfq_2_128_poly_set(k,a,a0);
     mpfq_2_128_poly_set(k,b,b0);
     while (mpfq_2_128_poly_deg(k,b)>=0) {
      mpfq_2_128_poly_divmod(k,q,r,a,b);
      mpfq_2_128_poly_set(k,a,b);
      mpfq_2_128_poly_set(k,b,r); 
     }
     mpfq_2_128_poly_setmonic(k,g,a);
    mpfq_2_128_poly_clear(k,a);
    mpfq_2_128_poly_clear(k,b);
    mpfq_2_128_poly_clear(k,q);
    mpfq_2_128_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly */
static inline
void mpfq_2_128_poly_xgcd(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_dst_poly g, mpfq_2_128_dst_poly u0, mpfq_2_128_dst_poly v0, mpfq_2_128_src_poly a0, mpfq_2_128_src_poly b0)
{
    mpfq_2_128_poly a,b,u,v,w,x,q,r;
    mpfq_2_128_elt c;
    mpfq_2_128_init(k,&c);
    mpfq_2_128_set_ui(k,c,0);        /* placate gcc */
    int da0=mpfq_2_128_poly_deg(k,a0), db0=mpfq_2_128_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      mpfq_2_128_poly_set(k,u0,a0);
      mpfq_2_128_poly_set(k,v0,b0);
      mpfq_2_128_poly_set(k,g,a0);
     } else {
      mpfq_2_128_poly_getcoeff(k,c,a0,da0);
      mpfq_2_128_inv(k,c,c);
      mpfq_2_128_poly_scal_mul(k,g,a0,c);
      mpfq_2_128_poly_set(k,v0,b0);
      mpfq_2_128_poly_set(k,u0,b0);
      mpfq_2_128_poly_setcoeff(k,u0,c,0);
     }
    }
    else {
     mpfq_2_128_poly_init(k,a,da0+1);
     mpfq_2_128_poly_init(k,b,db0+1);
     mpfq_2_128_poly_init(k,q,1);
     mpfq_2_128_poly_init(k,r,db0);
     mpfq_2_128_poly_set(k,a,a0);
     mpfq_2_128_poly_set(k,b,b0);
     mpfq_2_128_poly_init(k,u,1);
     mpfq_2_128_poly_init(k,v,1);
     mpfq_2_128_poly_init(k,w,1);
     mpfq_2_128_poly_init(k,x,1);
     mpfq_2_128_poly_setcoeff_ui(k,u,1,0);
     mpfq_2_128_poly_setcoeff_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (mpfq_2_128_poly_deg(k,b)>=0) {
      mpfq_2_128_poly_divmod(k,q,r,a,b);
      mpfq_2_128_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      mpfq_2_128_poly_set(k,b,r);
      mpfq_2_128_poly_mul(k,r,q,w);
      mpfq_2_128_poly_sub(k,r,u,r);
      mpfq_2_128_poly_set(k,u,w);   /* u,w <- w,u-qw */
      mpfq_2_128_poly_set(k,w,r);
      mpfq_2_128_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      mpfq_2_128_poly_sub(k,r,v,r);
      mpfq_2_128_poly_set(k,v,x);
      mpfq_2_128_poly_set(k,x,r);
     }
     dega=mpfq_2_128_poly_deg(k,a);
     mpfq_2_128_poly_getcoeff(k,c,a,dega);
     mpfq_2_128_inv(k,c,c);
     mpfq_2_128_poly_scal_mul(k,g,a,c);
     mpfq_2_128_poly_scal_mul(k,u0,u,c);
     mpfq_2_128_poly_scal_mul(k,v0,v,c);
     mpfq_2_128_poly_clear(k,a);
     mpfq_2_128_poly_clear(k,b);
     mpfq_2_128_poly_clear(k,u);
     mpfq_2_128_poly_clear(k,v);
     mpfq_2_128_poly_clear(k,w);
     mpfq_2_128_poly_clear(k,x);
     mpfq_2_128_poly_clear(k,q);
     mpfq_2_128_poly_clear(k,r);
    }
    mpfq_2_128_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_cmp */
static inline
int mpfq_2_128_poly_cmp(mpfq_2_128_dst_field k MAYBE_UNUSED, mpfq_2_128_src_poly u, mpfq_2_128_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return mpfq_2_128_vec_cmp(k, u->c, v->c, u->size);
}


#endif  /* MPFQ_2_128_H_ */

/* vim:set ft=cpp: */
