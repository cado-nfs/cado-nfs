#ifndef MPFQ_P_3_H_
#define MPFQ_P_3_H_

/* MPFQ generated file -- do not edit */

#include "mpfq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include "assert.h"
#include <limits.h>
#include "mpfq_fixmp.h"
#include "mpfq_gfp_common.h"
#include "mpfq_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      p_3

/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
   fieldtype=prime,
   n=3,
   nn=7,
   opthw=,
   tag=p_3,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_10, tag=p_10, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_11, tag=p_11, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_12, tag=p_12, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_13, tag=p_13, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_14, tag=p_14, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_15, tag=p_15, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_5, tag=p_5, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_6, tag=p_6, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_7, tag=p_7, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_9, tag=p_9, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     ],
    member_templates_restrict={
     m128=[
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     p_10=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_10, tag=p_10, }, ],
     p_11=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_11, tag=p_11, }, ],
     p_12=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_12, tag=p_12, }, ],
     p_13=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_13, tag=p_13, }, ],
     p_14=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_14, tag=p_14, }, ],
     p_15=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_15, tag=p_15, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     p_5=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_5, tag=p_5, }, ],
     p_6=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_6, tag=p_6, }, ],
     p_7=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_7, tag=p_7, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     p_9=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_9, tag=p_9, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     u64k1=[
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     u64k2=[
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     u64k3=[
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     u64k4=[
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
      { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
      ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=mpfq_vbase,
    global_prefix=mpfq_,
    name=mpfq_vbase,
    substitutions=[
     [ (?^:mpfq_p_3_elt \*), void *, ],
     [ (?^:mpfq_p_3_src_elt\b), const void *, ],
     [ (?^:mpfq_p_3_elt\b), void *, ],
     [ (?^:mpfq_p_3_dst_elt\b), void *, ],
     [ (?^:mpfq_p_3_elt_ur \*), void *, ],
     [ (?^:mpfq_p_3_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_p_3_elt_ur\b), void *, ],
     [ (?^:mpfq_p_3_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_p_3_vec \*), void *, ],
     [ (?^:mpfq_p_3_src_vec\b), const void *, ],
     [ (?^:mpfq_p_3_vec\b), void *, ],
     [ (?^:mpfq_p_3_dst_vec\b), void *, ],
     [ (?^:mpfq_p_3_vec_ur \*), void *, ],
     [ (?^:mpfq_p_3_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_p_3_vec_ur\b), void *, ],
     [ (?^:mpfq_p_3_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_p_3_poly \*), void *, ],
     [ (?^:mpfq_p_3_src_poly\b), const void *, ],
     [ (?^:mpfq_p_3_poly\b), void *, ],
     [ (?^:mpfq_p_3_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */

typedef mpfq_p_field mpfq_p_3_field;
typedef mpfq_p_dst_field mpfq_p_3_dst_field;
typedef mpfq_p_src_field mpfq_p_3_src_field;

typedef unsigned long mpfq_p_3_elt[3];
typedef unsigned long * mpfq_p_3_dst_elt;
typedef const unsigned long * mpfq_p_3_src_elt;

typedef unsigned long mpfq_p_3_elt_ur[7];
typedef unsigned long * mpfq_p_3_dst_elt_ur;
typedef const unsigned long * mpfq_p_3_src_elt_ur;

typedef mpfq_p_3_elt * mpfq_p_3_vec;
typedef mpfq_p_3_elt * mpfq_p_3_dst_vec;
typedef mpfq_p_3_elt * mpfq_p_3_src_vec;

typedef mpfq_p_3_elt_ur * mpfq_p_3_vec_ur;
typedef mpfq_p_3_elt_ur * mpfq_p_3_dst_vec_ur;
typedef mpfq_p_3_elt_ur * mpfq_p_3_src_vec_ur;

typedef struct {
  mpfq_p_3_vec c;
  unsigned long alloc;
  unsigned long size;
} mpfq_p_3_poly_struct;
typedef mpfq_p_3_poly_struct mpfq_p_3_poly [1];
typedef mpfq_p_3_poly_struct * mpfq_p_3_dst_poly;
typedef mpfq_p_3_poly_struct * mpfq_p_3_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name, Mpfq::gfp */
#define mpfq_p_3_impl_name()	"p_3"
/* *Mpfq::gfp::field::code_for_impl_max_characteristic_bits, Mpfq::gfp */
#define mpfq_p_3_impl_max_characteristic_bits()	192
/* *Mpfq::gfp::field::code_for_impl_max_degree, Mpfq::gfp */
#define mpfq_p_3_impl_max_degree()	1

/* Functions operating on the field structure */
static inline
void mpfq_p_3_field_characteristic(mpfq_p_3_src_field, mpz_ptr);
static inline
mpz_srcptr mpfq_p_3_field_characteristic_srcptr(mpfq_p_3_src_field);
static inline
unsigned long mpfq_p_3_field_characteristic_bits(mpfq_p_3_src_field);
/* *Mpfq::gfp::field::code_for_field_degree, Mpfq::gfp */
#define mpfq_p_3_field_degree(K)	1
static inline
void mpfq_p_3_field_init(mpfq_p_3_dst_field);
void mpfq_p_3_field_clear(mpfq_p_3_dst_field);
void mpfq_p_3_field_specify(mpfq_p_3_dst_field, unsigned long, const void *);
/* *Mpfq::gfp::field::code_for_field_setopt, Mpfq::gfp */
#define mpfq_p_3_field_setopt(f, x, y)	/**/

/* Element allocation functions */
static inline
void mpfq_p_3_init(mpfq_p_3_dst_field, mpfq_p_3_elt *);
static inline
void mpfq_p_3_clear(mpfq_p_3_dst_field, mpfq_p_3_elt *);
/* *Mpfq::defaults::flatdata::code_for_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define mpfq_p_3_elt_stride(k)	sizeof(mpfq_p_3_elt)

/* Elementary assignment functions */
static inline
void mpfq_p_3_set(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_set_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, unsigned long);
static inline
void mpfq_p_3_set_zero(mpfq_p_3_dst_field, mpfq_p_3_dst_elt);
static inline
unsigned long mpfq_p_3_get_ui(mpfq_p_3_dst_field, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_set_mpn(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mp_limb_t *, size_t);
static inline
void mpfq_p_3_set_mpz(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpz_srcptr);
static inline
void mpfq_p_3_get_mpn(mpfq_p_3_dst_field, mp_limb_t *, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_get_mpz(mpfq_p_3_dst_field, mpz_ptr, mpfq_p_3_src_elt);

/* Assignment of random values */
static inline
void mpfq_p_3_random(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, gmp_randstate_t);
static inline
void mpfq_p_3_random2(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
static inline
void mpfq_p_3_add(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_sub(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_neg(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_mul(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_sqr(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt);
static inline
int mpfq_p_3_is_sqr(mpfq_p_3_dst_field, mpfq_p_3_src_elt);
int mpfq_p_3_sqrt(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_pow(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, unsigned long *, size_t);
void mpfq_p_3_powz(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, mpz_srcptr);
/* *Mpfq::gfp::elt::code_for_frobenius, Mpfq::gfp */
#define mpfq_p_3_frobenius(k, x, y)	mpfq_p_3_set(k, x, y)
static inline
void mpfq_p_3_add_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, unsigned long);
static inline
void mpfq_p_3_sub_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, unsigned long);
static inline
void mpfq_p_3_mul_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, unsigned long);
#define HAVE_mpfq_p_3_normalize
static inline
void mpfq_p_3_normalize(mpfq_p_3_dst_field, mpfq_p_3_dst_elt);
static inline
int mpfq_p_3_inv(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt);
#define HAVE_mpfq_p_3_hadamard
static inline
void mpfq_p_3_hadamard(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_dst_elt, mpfq_p_3_dst_elt, mpfq_p_3_dst_elt);

/* Operations involving unreduced elements */
static inline
void mpfq_p_3_elt_ur_init(mpfq_p_3_dst_field, mpfq_p_3_elt_ur *);
static inline
void mpfq_p_3_elt_ur_clear(mpfq_p_3_dst_field, mpfq_p_3_elt_ur *);
/* *Mpfq::defaults::flatdata::code_for_elt_ur_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define mpfq_p_3_elt_ur_stride(k)	sizeof(mpfq_p_3_elt_ur)
static inline
void mpfq_p_3_elt_ur_set(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt_ur);
static inline
void mpfq_p_3_elt_ur_set_elt(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_elt_ur_set_zero(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur);
static inline
void mpfq_p_3_elt_ur_set_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, unsigned long);
static inline
void mpfq_p_3_elt_ur_add(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt_ur, mpfq_p_3_src_elt_ur);
static inline
void mpfq_p_3_elt_ur_neg(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt_ur);
static inline
void mpfq_p_3_elt_ur_sub(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt_ur, mpfq_p_3_src_elt_ur);
static inline
void mpfq_p_3_mul_ur(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_sqr_ur(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_reduce(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_dst_elt_ur);
#define HAVE_mpfq_p_3_addmul_si_ur
static inline
void mpfq_p_3_addmul_si_ur(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_elt, long);

/* Comparison functions */
static inline
int mpfq_p_3_cmp(mpfq_p_3_dst_field, mpfq_p_3_src_elt, mpfq_p_3_src_elt);
static inline
int mpfq_p_3_cmp_ui(mpfq_p_3_dst_field, mpfq_p_3_src_elt, unsigned long);
static inline
int mpfq_p_3_is_zero(mpfq_p_3_dst_field, mpfq_p_3_src_elt);

/* Input/output functions */
int mpfq_p_3_asprint(mpfq_p_3_dst_field, char * *, mpfq_p_3_src_elt);
int mpfq_p_3_fprint(mpfq_p_3_dst_field, FILE *, mpfq_p_3_src_elt);
/* *Mpfq::defaults::code_for_print, Mpfq::gfp */
#define mpfq_p_3_print(k, x)	mpfq_p_3_fprint(k,stdout,x)
int mpfq_p_3_sscan(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, const char *);
int mpfq_p_3_fscan(mpfq_p_3_dst_field, FILE *, mpfq_p_3_dst_elt);
/* *Mpfq::defaults::code_for_scan, Mpfq::gfp */
#define mpfq_p_3_scan(k, x)	mpfq_p_3_fscan(k,stdin,x)

/* Vector functions */
void mpfq_p_3_vec_init(mpfq_p_3_dst_field, mpfq_p_3_vec *, unsigned long);
void mpfq_p_3_vec_reinit(mpfq_p_3_dst_field, mpfq_p_3_vec *, unsigned long, unsigned long);
void mpfq_p_3_vec_clear(mpfq_p_3_dst_field, mpfq_p_3_vec *, unsigned long);
static inline
void mpfq_p_3_vec_set(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_set_zero(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, unsigned long);
static inline
void mpfq_p_3_vec_setcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_elt, unsigned long);
static inline
void mpfq_p_3_vec_setcoeff_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, unsigned long, unsigned long);
static inline
void mpfq_p_3_vec_getcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_add(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_neg(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_rev(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_sub(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_scal_mul(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, mpfq_p_3_src_elt, unsigned long);
static inline
void mpfq_p_3_vec_conv(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, unsigned long, mpfq_p_3_src_vec, unsigned long);
void mpfq_p_3_vec_random(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, unsigned long, gmp_randstate_t);
void mpfq_p_3_vec_random2(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, unsigned long, gmp_randstate_t);
long mpfq_p_3_vec_cmp(mpfq_p_3_dst_field, mpfq_p_3_src_vec, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_is_zero(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
static inline
mpfq_p_3_dst_vec mpfq_p_3_vec_subvec(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, long);
static inline
mpfq_p_3_src_vec mpfq_p_3_vec_subvec_const(mpfq_p_3_dst_field, mpfq_p_3_src_vec, long);
static inline
mpfq_p_3_dst_elt mpfq_p_3_vec_coeff_ptr(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, long);
static inline
mpfq_p_3_src_elt mpfq_p_3_vec_coeff_ptr_const(mpfq_p_3_dst_field, mpfq_p_3_src_vec, long);
long mpfq_p_3_vec_asprint(mpfq_p_3_dst_field, char * *, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_fprint(mpfq_p_3_dst_field, FILE *, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_print(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_sscan(mpfq_p_3_dst_field, mpfq_p_3_vec *, unsigned long *, const char *);
long mpfq_p_3_vec_fscan(mpfq_p_3_dst_field, FILE *, mpfq_p_3_vec *, unsigned long *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec, Mpfq::gfp */
#define mpfq_p_3_vec_scan(K, w, n)	mpfq_p_3_vec_fscan(K,stdin,w,n)
long mpfq_p_3_vec_hamming_weight(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_find_first_set(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_simd_hamming_weight(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
long mpfq_p_3_vec_simd_find_first_set(mpfq_p_3_dst_field, mpfq_p_3_src_vec, unsigned long);
void mpfq_p_3_vec_ur_init(mpfq_p_3_dst_field, mpfq_p_3_vec_ur *, unsigned long);
static inline
void mpfq_p_3_vec_ur_set_zero(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_set_vec(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec, unsigned long);
void mpfq_p_3_vec_ur_reinit(mpfq_p_3_dst_field, mpfq_p_3_vec_ur *, unsigned long, unsigned long);
void mpfq_p_3_vec_ur_clear(mpfq_p_3_dst_field, mpfq_p_3_vec_ur *, unsigned long);
static inline
void mpfq_p_3_vec_ur_set(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_setcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_elt_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_getcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_elt_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_add(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_sub(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_neg(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_ur_rev(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec_ur, unsigned long);
static inline
void mpfq_p_3_vec_scal_mul_ur(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec, mpfq_p_3_src_elt, unsigned long);
void mpfq_p_3_vec_conv_ur_ks(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec, unsigned int, mpfq_p_3_src_vec, unsigned int);
static inline
void mpfq_p_3_vec_conv_ur_n(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec, mpfq_p_3_src_vec, unsigned int);
static inline
void mpfq_p_3_vec_conv_ur(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, mpfq_p_3_src_vec, unsigned long, mpfq_p_3_src_vec, unsigned long);
static inline
void mpfq_p_3_vec_reduce(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_dst_vec_ur, unsigned long);
static inline
mpfq_p_3_dst_vec_ur mpfq_p_3_vec_ur_subvec(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, long);
static inline
mpfq_p_3_src_vec_ur mpfq_p_3_vec_ur_subvec_const(mpfq_p_3_dst_field, mpfq_p_3_src_vec_ur, long);
static inline
mpfq_p_3_dst_elt mpfq_p_3_vec_ur_coeff_ptr(mpfq_p_3_dst_field, mpfq_p_3_dst_vec_ur, long);
static inline
mpfq_p_3_src_elt mpfq_p_3_vec_ur_coeff_ptr_const(mpfq_p_3_dst_field, mpfq_p_3_src_vec_ur, long);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define mpfq_p_3_vec_elt_stride(k, n)	((n) * mpfq_p_3_elt_stride((k)))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define mpfq_p_3_vec_ur_elt_stride(k, n)	((n) * mpfq_p_3_elt_ur_stride((k)))

/* Polynomial functions */
static inline
void mpfq_p_3_poly_init(mpfq_p_3_dst_field, mpfq_p_3_poly, unsigned int);
static inline
void mpfq_p_3_poly_clear(mpfq_p_3_dst_field, mpfq_p_3_poly);
static inline
void mpfq_p_3_poly_set(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly);
void mpfq_p_3_poly_setmonic(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_setcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_elt, unsigned int);
static inline
void mpfq_p_3_poly_setcoeff_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, unsigned long, unsigned int);
static inline
void mpfq_p_3_poly_getcoeff(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_poly, unsigned int);
static inline
int mpfq_p_3_poly_deg(mpfq_p_3_dst_field, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_add(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_sub(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_set_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, unsigned long);
static inline
void mpfq_p_3_poly_add_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, unsigned long);
static inline
void mpfq_p_3_poly_sub_ui(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, unsigned long);
static inline
void mpfq_p_3_poly_neg(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_scal_mul(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_elt);
static inline
void mpfq_p_3_poly_mul(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
int mpfq_p_3_poly_divmod(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
void mpfq_p_3_poly_precomp_mod(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly);
void mpfq_p_3_poly_mod_pre(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_gcd(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_xgcd(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, mpfq_p_3_dst_poly, mpfq_p_3_dst_poly, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
void mpfq_p_3_poly_random(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, unsigned int, gmp_randstate_t);
static inline
void mpfq_p_3_poly_random2(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, unsigned int, gmp_randstate_t);
static inline
int mpfq_p_3_poly_cmp(mpfq_p_3_dst_field, mpfq_p_3_src_poly, mpfq_p_3_src_poly);
static inline
int mpfq_p_3_poly_asprint(mpfq_p_3_dst_field, char * *, mpfq_p_3_src_poly);
static inline
int mpfq_p_3_poly_fprint(mpfq_p_3_dst_field, FILE *, mpfq_p_3_src_poly);
static inline
int mpfq_p_3_poly_print(mpfq_p_3_dst_field, mpfq_p_3_src_poly);
static inline
int mpfq_p_3_poly_sscan(mpfq_p_3_dst_field, mpfq_p_3_dst_poly, const char *);
static inline
int mpfq_p_3_poly_fscan(mpfq_p_3_dst_field, FILE *, mpfq_p_3_dst_poly);
static inline
int mpfq_p_3_poly_scan(mpfq_p_3_dst_field, mpfq_p_3_dst_poly);

/* Functions related to SIMD operation */
/* *simd_gfp::code_for_simd_groupsize */
#define mpfq_p_3_simd_groupsize(K)	1
static inline
int mpfq_p_3_simd_hamming_weight(mpfq_p_3_dst_field, mpfq_p_3_src_elt);
static inline
int mpfq_p_3_simd_find_first_set(mpfq_p_3_dst_field, mpfq_p_3_src_elt);
static inline
unsigned long mpfq_p_3_simd_get_ui_at(mpfq_p_3_dst_field, mpfq_p_3_src_elt, int);
static inline
void mpfq_p_3_simd_set_ui_at(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, int, unsigned long);
static inline
void mpfq_p_3_simd_add_ui_at(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, mpfq_p_3_src_elt, int, unsigned long);
static inline
void mpfq_p_3_simd_set_ui_all(mpfq_p_3_dst_field, mpfq_p_3_dst_elt, unsigned long);
void mpfq_p_3_add_dotprod(mpfq_p_3_dst_field, mpfq_p_3_dst_vec, mpfq_p_3_src_vec, mpfq_p_3_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* Object-oriented interface */
static inline
void mpfq_p_3_oo_field_clear(mpfq_vbase_ptr);
void mpfq_p_3_oo_field_init(mpfq_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::gfp::field::code_for_field_characteristic, Mpfq::gfp */
static inline
void mpfq_p_3_field_characteristic(mpfq_p_3_src_field k, mpz_ptr z)
{
        mpz_set(z, k->p);
}

/* *Mpfq::gfp::field::code_for_field_characteristic_srcptr, Mpfq::gfp */
static inline
mpz_srcptr mpfq_p_3_field_characteristic_srcptr(mpfq_p_3_src_field k)
{
    return k->p;
}

/* *Mpfq::gfp::field::code_for_field_characteristic_bits, Mpfq::gfp */
static inline
unsigned long mpfq_p_3_field_characteristic_bits(mpfq_p_3_src_field k)
{
        return mpz_sizeinbase(k->p, 2);
}

/* *Mpfq::gfp::field::code_for_field_init, Mpfq::gfp */
static inline
void mpfq_p_3_field_init(mpfq_p_3_dst_field k)
{
    mpz_init(k->p);
    mpz_init(k->bigmul_p);
    k->io_base = 10;
    mpz_init(k->factor);
    k->ts_info.e=0;
}

/* *Mpfq::gfp::elt::code_for_init, Mpfq::gfp */
static inline
void mpfq_p_3_init(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_clear, Mpfq::gfp */
static inline
void mpfq_p_3_clear(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_set, Mpfq::gfp */
static inline
void mpfq_p_3_set(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt r, mpfq_p_3_src_elt s)
{
    mpfq_copy(r,s,3);
}

/* *Mpfq::gfp::elt::code_for_set_ui, Mpfq::gfp */
static inline
void mpfq_p_3_set_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt r, unsigned long x)
{
    assert (r);
    r[0] = x;
    mpfq_zero(r + 1, 3 - 1);
}

/* *Mpfq::gfp::elt::code_for_set_zero, Mpfq::gfp */
static inline
void mpfq_p_3_set_zero(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt r)
{
    mpfq_zero(r, 3);
}

/* *Mpfq::gfp::elt::code_for_get_ui, Mpfq::gfp */
static inline
unsigned long mpfq_p_3_get_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_src_elt x)
{
    return x[0];
}

/* *Mpfq::gfp::elt::code_for_set_mpn, Mpfq::gfp */
static inline
void mpfq_p_3_set_mpn(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt r, mp_limb_t * x, size_t n)
{
    if (n < 3) {
        mpfq_copy(r, x, n);
        mpfq_zero(r + n, 3 - n);
    } else {
        mp_limb_t tmp[n-3+1];
        mpn_tdiv_qr(tmp, r, 0, x, n, k->p->_mp_d, 3);
    }
}

/* *Mpfq::gfp::elt::code_for_set_mpz, Mpfq::gfp */
static inline
void mpfq_p_3_set_mpz(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt r, mpz_srcptr z)
{
    if (z->_mp_size < 0) {
        mpfq_p_3_set_mpn(k, r, z->_mp_d, -z->_mp_size);
        mpfq_p_3_neg(k, r, r);
    } else {
        mpfq_p_3_set_mpn(k, r, z->_mp_d, z->_mp_size);
    }
}

/* *Mpfq::gfp::elt::code_for_get_mpn, Mpfq::gfp */
static inline
void mpfq_p_3_get_mpn(mpfq_p_3_dst_field k MAYBE_UNUSED, mp_limb_t * r, mpfq_p_3_src_elt x)
{
    mpfq_copy(r, x, 3);
}

/* *Mpfq::gfp::elt::code_for_get_mpz, Mpfq::gfp */
static inline
void mpfq_p_3_get_mpz(mpfq_p_3_dst_field k MAYBE_UNUSED, mpz_ptr z, mpfq_p_3_src_elt y)
{
    int i; 
    mpz_realloc2(z, 3*64);
    for (i = 0; i < 3; ++i)
        z->_mp_d[i] = y[i];
    i = 3;
    while (i>=1 && z->_mp_d[i-1] == 0)
        i--;
    z->_mp_size = i;
}

/* *Mpfq::gfp::elt::code_for_random, Mpfq::gfp */
static inline
void mpfq_p_3_random(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt x, gmp_randstate_t state)
{
      mpz_t z;
      mpz_init(z);
      mpz_urandomb(z, state, 3 * GMP_LIMB_BITS);
      mpfq_copy(x, z->_mp_d, 3);
      mpz_clear(z);
    mpfq_p_3_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_random2, Mpfq::gfp */
static inline
void mpfq_p_3_random2(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt x, gmp_randstate_t state)
{
      mpz_t z;
      mpz_init(z);
      mpz_rrandomb(z, state, 3 * GMP_LIMB_BITS);
      mpfq_copy(x, z->_mp_d, 3);
      mpz_clear(z);
    mpfq_p_3_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_add, Mpfq::gfp */
static inline
void mpfq_p_3_add(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, mpfq_p_3_src_elt y)
{
    mp_limb_t cy;
    cy = mpfq_fixmp_3_add(z, x, y);
    if (cy || (mpfq_fixmp_3_cmp(z, k->p->_mp_d) >= 0))
        mpfq_fixmp_3_sub(z, z, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_sub, Mpfq::gfp */
static inline
void mpfq_p_3_sub(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, mpfq_p_3_src_elt y)
{
    mp_limb_t cy;
    cy = mpfq_fixmp_3_sub(z, x, y);
    if (cy) // negative result
        mpfq_fixmp_3_add(z, z, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_neg, Mpfq::gfp */
static inline
void mpfq_p_3_neg(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x)
{
    if (mpfq_fixmp_3_cmp_ui(x, 0))
        mpfq_fixmp_3_sub(z, k->p->_mp_d, x);
    else {
        int i;
        for (i = 0; i < 3; ++i)
            z[i] = 0;
        }
}

/* *Mpfq::gfp::elt::code_for_mul, Mpfq::gfp */
static inline
void mpfq_p_3_mul(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, mpfq_p_3_src_elt y)
{
    mp_limb_t tmp[6];
    mpfq_fixmp_3_mul(tmp, x, y);
    mpfq_fixmp_3_mod(z, tmp, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_sqr, Mpfq::gfp */
static inline
void mpfq_p_3_sqr(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x)
{
    mp_limb_t tmp[6];
    mpfq_fixmp_3_sqr(tmp, x);
    mpfq_fixmp_3_mod(z, tmp, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_is_sqr, Mpfq::gfp */
static inline
int mpfq_p_3_is_sqr(mpfq_p_3_dst_field k, mpfq_p_3_src_elt x)
{
    mp_limb_t pp[3];
    mpfq_p_3_elt y;
    mpfq_fixmp_3_sub_ui_nc(pp, k->p->_mp_d, 1);
    mpfq_fixmp_3_rshift(pp, 1);
    mpfq_p_3_init(k, &y);
    mpfq_p_3_pow(k, y, x, pp, 3);
    int res = mpfq_p_3_cmp_ui(k, y, 1);
    mpfq_p_3_clear(k, &y);
    if (res == 0)
        return 1;
    else 
        return 0;
}

/* *Mpfq::defaults::pow::code_for_pow, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void mpfq_p_3_pow(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt res, mpfq_p_3_src_elt r, unsigned long * x, size_t n)
{
    mpfq_p_3_elt u, a;
    long i, j, lead;     /* it is a signed type */
    unsigned long mask;
    
    /* get the correct (i,j) position of the most significant bit in x */
    for(i = ((long)n)-1; i>=0 && x[i]==0; i--)
        ;
    if (i < 0) {
        /* power zero gets 1 */
        mpfq_p_3_set_ui(k, res, 1);
        return;
    }
    j = 64 - 1;
    mask = (1UL<<j);
    for( ; (x[i]&mask)==0 ;j--, mask>>=1)
        ;
    lead = i*64+j;      /* Ensured. */
    
    mpfq_p_3_init(k, &u);
    mpfq_p_3_init(k, &a);
    mpfq_p_3_set(k, a, r);
    for( ; lead > 0; lead--) {
        if (j-- == 0) {
            i--;
            j = 64-1;
            mask = (1UL<<j);
        } else {
            mask >>= 1;
        }
        if (x[i]&mask) {
            mpfq_p_3_sqr(k, u, a);
            mpfq_p_3_mul(k, a, u, r);
        } else {
            mpfq_p_3_sqr(k, a,a);
        }
    }
    mpfq_p_3_set(k, res, a);
    mpfq_p_3_clear(k, &u);
    mpfq_p_3_clear(k, &a);
}

/* *Mpfq::gfp::elt::code_for_add_ui, Mpfq::gfp */
static inline
void mpfq_p_3_add_ui(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = mpfq_fixmp_3_add_ui(z, x, y);
    if (cy || (mpfq_fixmp_3_cmp(z, k->p->_mp_d) >= 0))
        mpfq_fixmp_3_sub(z, z, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_sub_ui, Mpfq::gfp */
static inline
void mpfq_p_3_sub_ui(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = mpfq_fixmp_3_sub_ui(z, x, y);
    if (cy) // negative result
        mpfq_fixmp_3_add(z, z, k->p->_mp_d);
}

/* *Mpfq::gfp::elt::code_for_mul_ui, Mpfq::gfp */
static inline
void mpfq_p_3_mul_ui(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x, unsigned long y)
{
    mp_limb_t tmp[3+1], q[2];
    mpfq_fixmp_3_mul1(tmp,x,y);
    mpn_tdiv_qr(q, z, 0, tmp, 3+1, k->p->_mp_d, 3);
}

/* *Mpfq::gfp::elt::code_for_normalize, Mpfq::gfp */
static inline
void mpfq_p_3_normalize(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt x)
{
    if (mpfq_fixmp_3_cmp(x,k->p->_mp_d)>=0) {
      mp_limb_t q[3+1];
      mpfq_p_3_elt r;
      mpn_tdiv_qr(q, r, 0, x, 3, k->p->_mp_d, 3);
      mpfq_p_3_set(k, x, r);
    }
}

/* *Mpfq::gfp::elt::code_for_inv, Mpfq::gfp */
static inline
int mpfq_p_3_inv(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_src_elt x)
{
    int ret=mpfq_fixmp_3_invmod(z, x, k->p->_mp_d);
    if (!ret)
        mpfq_p_3_get_mpz(k, k->factor, z);
    return ret;
}

/* *Mpfq::gfp::elt::code_for_hadamard, Mpfq::gfp */
static inline
void mpfq_p_3_hadamard(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt x, mpfq_p_3_dst_elt y, mpfq_p_3_dst_elt z, mpfq_p_3_dst_elt t)
{
    mpfq_p_3_elt tmp;
    mpfq_p_3_init(k, &tmp);
    mpfq_p_3_add(k, tmp, x, y);
    mpfq_p_3_sub(k, y, x, y);
    mpfq_p_3_set(k, x, tmp);
    mpfq_p_3_add(k, tmp, z, t);
    mpfq_p_3_sub(k, t, z, t);
    mpfq_p_3_set(k, z, tmp);
    mpfq_p_3_sub(k, tmp, x, z);
    mpfq_p_3_add(k, x, x, z);
    mpfq_p_3_add(k, z, y, t);
    mpfq_p_3_sub(k, t, y, t);
    mpfq_p_3_set(k, y, tmp);
    mpfq_p_3_clear(k, &tmp); 
}

/* *Mpfq::gfp::elt::code_for_elt_ur_init, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_init(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_clear, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_clear(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_set(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt_ur x)
{
    mpfq_copy(z, x, 7);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set_elt, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_set_elt(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt_ur r, mpfq_p_3_src_elt s)
{
    mpfq_copy(r, s, 3); mpfq_zero(r+3, 7-3);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set_zero, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_set_zero(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt_ur r)
{
    mpfq_zero(r, 7);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set_ui, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_set_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur r, unsigned long x)
{
    assert (r); 
    r[0] = x;
    mpfq_zero(r + 1, 7 - 1);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_add, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_add(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt_ur x, mpfq_p_3_src_elt_ur y)
{
    mpn_add_n(z, x, y, 7);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_neg, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_neg(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt_ur x)
{
    mpfq_p_3_elt_ur tmp;
    mpfq_p_3_elt_ur_init(k, &tmp);
    mpfq_zero(tmp, 7);
    mpn_sub_n(z, tmp, x, 7);
    mpfq_p_3_elt_ur_clear(k, &tmp);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_sub, Mpfq::gfp */
static inline
void mpfq_p_3_elt_ur_sub(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt_ur x, mpfq_p_3_src_elt_ur y)
{
    mpn_sub_n(z, x, y, 7);
}

/* *Mpfq::gfp::elt::code_for_mul_ur, Mpfq::gfp */
static inline
void mpfq_p_3_mul_ur(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt x, mpfq_p_3_src_elt y)
{
    mpfq_fixmp_3_mul(z, x, y);
    mpfq_zero(z + 6, 7 - 6);
}

/* *Mpfq::gfp::elt::code_for_sqr_ur, Mpfq::gfp */
static inline
void mpfq_p_3_sqr_ur(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt_ur z, mpfq_p_3_src_elt x)
{
    mpfq_fixmp_3_sqr(z, x);
    mpfq_zero(z + 6, 7 - 6);
}

/* *Mpfq::gfp::elt::code_for_reduce, Mpfq::gfp */
static inline
void mpfq_p_3_reduce(mpfq_p_3_dst_field k, mpfq_p_3_dst_elt z, mpfq_p_3_dst_elt_ur x)
{
    mp_limb_t q[7+1];
    if (x[7-1]>>(64-1)) {
        // negative number, add bigmul_p to make it positive before reduction
        mpn_add_n(x, x, k->bigmul_p->_mp_d, 7);
    }
    mpn_tdiv_qr(q, z, 0, x, 7, k->p->_mp_d, 3);
}

/* *simd_gfp::code_for_addmul_si_ur */
static inline
void mpfq_p_3_addmul_si_ur(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt_ur w, mpfq_p_3_src_elt u, long v)
{
        mpfq_p_3_elt_ur s;
        mpfq_p_3_elt vx;
        mpfq_p_3_elt_ur_init(K, &s);
        mpfq_p_3_init(K, &vx);
        if (v>0) {
            mpfq_p_3_set_ui(K, vx, v);
            mpfq_p_3_mul_ur(K, s, u, vx);
            mpfq_p_3_elt_ur_add(K, w, w, s);
        } else {
            mpfq_p_3_set_ui(K, vx, -v);
            mpfq_p_3_mul_ur(K, s, u, vx);
            mpfq_p_3_elt_ur_sub(K, w, w, s);
        }
        mpfq_p_3_clear(K, &vx);
        mpfq_p_3_elt_ur_clear(K, &s);
}

/* *Mpfq::gfp::elt::code_for_cmp, Mpfq::gfp */
static inline
int mpfq_p_3_cmp(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_src_elt x, mpfq_p_3_src_elt y)
{
    return mpfq_fixmp_3_cmp(x,y);
}

/* *Mpfq::gfp::elt::code_for_cmp_ui, Mpfq::gfp */
static inline
int mpfq_p_3_cmp_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_src_elt x, unsigned long y)
{
    return mpfq_fixmp_3_cmp_ui(x,y);
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
int mpfq_p_3_is_zero(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(mpfq_p_3_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void mpfq_p_3_vec_set(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec r, mpfq_p_3_src_vec s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_p_3_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void mpfq_p_3_vec_set_zero(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_p_3_elt));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_setcoeff, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_setcoeff(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_elt x, unsigned long i)
{
            mpfq_p_3_dst_elt y = mpfq_p_3_vec_coeff_ptr(K, w, i);
            mpfq_p_3_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_setcoeff_ui, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_setcoeff_ui(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, unsigned long x, unsigned long i)
{
            mpfq_p_3_dst_elt y = mpfq_p_3_vec_coeff_ptr(K, w, i);
            mpfq_p_3_set_ui(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_getcoeff, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_getcoeff(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt x, mpfq_p_3_src_vec w, unsigned long i)
{
            mpfq_p_3_src_elt y = mpfq_p_3_vec_coeff_ptr_const(K, w, i);
            mpfq_p_3_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_add(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, mpfq_p_3_src_vec v, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_p_3_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_neg(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_p_3_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_rev(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_p_3_elt tmp[1];
    mpfq_p_3_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_p_3_set(K, tmp[0], u[i]);
        mpfq_p_3_set(K, w[i], u[n-1-i]);
        mpfq_p_3_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_p_3_set(K, w[nn], u[nn]);
    mpfq_p_3_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_sub(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, mpfq_p_3_src_vec v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_p_3_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_scal_mul(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, mpfq_p_3_src_elt c, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_p_3_src_elt x = mpfq_p_3_vec_coeff_ptr_const(K, u, i);
        mpfq_p_3_dst_elt y = mpfq_p_3_vec_coeff_ptr(K, w, i);
        mpfq_p_3_mul(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv, Mpfq::gfp */
static inline
void mpfq_p_3_vec_conv(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_src_vec u, unsigned long n, mpfq_p_3_src_vec v, unsigned long m)
{
    mpfq_p_3_vec_ur tmp;
    mpfq_p_3_vec_ur_init(K, &tmp, m+n-1);
    mpfq_p_3_vec_conv_ur(K, tmp, u, n, v, m);
    mpfq_p_3_vec_reduce(K, w, tmp, m+n-1);
    mpfq_p_3_vec_ur_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_dst_vec mpfq_p_3_vec_subvec(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_src_vec mpfq_p_3_vec_subvec_const(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_dst_elt mpfq_p_3_vec_coeff_ptr(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_src_elt mpfq_p_3_vec_coeff_ptr_const(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_set_zero(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_p_3_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_set_vec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_set_vec(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec u, unsigned long n)
{
            unsigned int i;
            for(i = 0; i < n; ++i) {
                mpfq_p_3_src_elt x = mpfq_p_3_vec_coeff_ptr_const(K, u, i);
                mpfq_p_3_dst_elt_ur y = mpfq_p_3_vec_ur_coeff_ptr(K, w, i);
                mpfq_p_3_elt_ur_set_elt(K, y, x);
            }
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_set(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur r, mpfq_p_3_src_vec_ur s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_p_3_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_setcoeff, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_setcoeff(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_elt_ur x, unsigned long i)
{
            mpfq_p_3_dst_elt_ur y = mpfq_p_3_vec_ur_coeff_ptr(K, w, i);
            mpfq_p_3_elt_ur_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_getcoeff, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_getcoeff(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt_ur x, mpfq_p_3_src_vec_ur w, unsigned long i)
{
            mpfq_p_3_src_elt_ur y = mpfq_p_3_vec_ur_coeff_ptr_const(K, w, i);
            mpfq_p_3_elt_ur_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_add(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec_ur u, mpfq_p_3_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_p_3_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_sub(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec_ur u, mpfq_p_3_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_p_3_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_neg(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_p_3_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_ur_rev(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec_ur u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_p_3_elt_ur tmp[1];
    mpfq_p_3_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_p_3_elt_ur_set(K, tmp[0], u[i]);
        mpfq_p_3_elt_ur_set(K, w[i], u[n-1-i]);
        mpfq_p_3_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_p_3_elt_ur_set(K, w[nn], u[nn]);
    mpfq_p_3_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_scal_mul_ur(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec u, mpfq_p_3_src_elt c, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_p_3_src_elt x = mpfq_p_3_vec_coeff_ptr_const(K, u, i);
        mpfq_p_3_dst_elt_ur y = mpfq_p_3_vec_ur_coeff_ptr(K, w, i);
        mpfq_p_3_mul_ur(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
/* Triggered by: vec_conv_ur */
static inline
void mpfq_p_3_vec_conv_ur_n(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec u, mpfq_p_3_src_vec v, unsigned int n)
{
    if (n == 0)
        return;
    if (n == 1) {
        mpfq_p_3_mul_ur(K, w[0], u[0], v[0]);
        return;
    }
    if (n == 2) {  // Kara 2
        mpfq_p_3_elt t1, t2;
        mpfq_p_3_init(K, &t1);
        mpfq_p_3_init(K, &t2);
        mpfq_p_3_mul_ur(K, w[0], u[0], v[0]);
        mpfq_p_3_mul_ur(K, w[2], u[1], v[1]);
        mpfq_p_3_add(K, t1, u[0], u[1]);
        mpfq_p_3_add(K, t2, v[0], v[1]);
        mpfq_p_3_mul_ur(K, w[1], t1, t2);
        mpfq_p_3_elt_ur_sub(K, w[1], w[1], w[0]);
        mpfq_p_3_elt_ur_sub(K, w[1], w[1], w[2]);
        mpfq_p_3_clear(K, &t1);
        mpfq_p_3_clear(K, &t2);
        return;
    }
    if (n == 3) {  // do it in 6
        mpfq_p_3_elt t1, t2;
        mpfq_p_3_elt_ur s;
        mpfq_p_3_init(K, &t1);
        mpfq_p_3_init(K, &t2);
        mpfq_p_3_elt_ur_init(K, &s);
        // a0*b0*(1 - X)
        mpfq_p_3_mul_ur(K, w[0], u[0], v[0]);
        mpfq_p_3_elt_ur_neg(K, w[1], w[0]);
        // a1*b1*(-X + 2*X^2 - X^3)
        mpfq_p_3_mul_ur(K, w[2], u[1], v[1]);
        mpfq_p_3_elt_ur_neg(K, w[3], w[2]);
        mpfq_p_3_elt_ur_add(K, w[2], w[2], w[2]);
        mpfq_p_3_elt_ur_add(K, w[1], w[1], w[3]);
        // a2*b2*(-X^3+X^4)
        mpfq_p_3_mul_ur(K, w[4], u[2], v[2]);
        mpfq_p_3_elt_ur_sub(K, w[3], w[3], w[4]);
        // (a0+a1)*(b0+b1)*(X - X^2)
        mpfq_p_3_add(K, t1, u[0], u[1]);
        mpfq_p_3_add(K, t2, v[0], v[1]);
        mpfq_p_3_mul_ur(K, s, t1, t2);
        mpfq_p_3_elt_ur_add(K, w[1], w[1], s);
        mpfq_p_3_elt_ur_sub(K, w[2], w[2], s);
        // (a1+a2)*(b1+b2)*(X^3 - X^2)
        mpfq_p_3_add(K, t1, u[1], u[2]);
        mpfq_p_3_add(K, t2, v[1], v[2]);
        mpfq_p_3_mul_ur(K, s, t1, t2);
        mpfq_p_3_elt_ur_add(K, w[3], w[3], s);
        mpfq_p_3_elt_ur_sub(K, w[2], w[2], s);
        // (a0+a1+a2)*(b0+b1+b2)* X^2
        mpfq_p_3_add(K, t1, u[0], t1);
        mpfq_p_3_add(K, t2, v[0], t2);
        mpfq_p_3_mul_ur(K, s, t1, t2);
        mpfq_p_3_elt_ur_add(K, w[2], w[2], s);
        return;
    }
    unsigned int n0, n1;
    n0 = n / 2;
    n1 = n - n0;
    mpfq_p_3_vec_conv_ur_n(K, w, u, v, n0);
    mpfq_p_3_vec_conv_ur_n(K, w + 2*n0, u + n0, v + n0, n1);
    mpfq_p_3_elt_ur_set_ui(K, w[2*n0-1], 0);
    
    mpfq_p_3_vec tmpu, tmpv;
    mpfq_p_3_vec_ur tmpw;
    mpfq_p_3_vec_init(K, &tmpu, n1);
    mpfq_p_3_vec_init(K, &tmpv, n1);
    mpfq_p_3_vec_ur_init(K, &tmpw, 2*n1-1);
    
    mpfq_p_3_vec_set(K, tmpu, u, n0);
    if (n1 != n0) 
        mpfq_p_3_set_ui(K, tmpu[n0], 0);
    mpfq_p_3_vec_add(K, tmpu, tmpu, u+n0, n1);
    mpfq_p_3_vec_set(K, tmpv, v, n0);
    if (n1 != n0) 
        mpfq_p_3_set_ui(K, tmpv[n0], 0);
    mpfq_p_3_vec_add(K, tmpv, tmpv, v+n0, n1);
    mpfq_p_3_vec_conv_ur_n(K, tmpw, tmpu, tmpv, n1);
    mpfq_p_3_vec_ur_sub(K, tmpw, tmpw, w, 2*n0-1);
    mpfq_p_3_vec_ur_sub(K, tmpw, tmpw, w + 2*n0, 2*n1-1);
    mpfq_p_3_vec_ur_add(K, w + n0, w + n0, tmpw, 2*n1-1);
    
    mpfq_p_3_vec_clear(K, &tmpu, n1);
    mpfq_p_3_vec_clear(K, &tmpv, n1);
    mpfq_p_3_vec_ur_clear(K, &tmpw, 2*n1-1);
    return;
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
static inline
void mpfq_p_3_vec_conv_ur(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur w, mpfq_p_3_src_vec u, unsigned long n, mpfq_p_3_src_vec v, unsigned long m)
{
    if ((n > 1) && (m > 1) && (n+m > 15)) {
        mpfq_p_3_vec_conv_ur_ks(K, w, u, n, v, m);
        return;
    }
    if (n == m) {
        mpfq_p_3_vec_conv_ur_n(K, w, u, v, n);
        return;
    }
    unsigned int i, j MAYBE_UNUSED, k;
    mpfq_p_3_elt_ur acc, z;
    mpfq_p_3_elt_ur_init(K, &acc);
    mpfq_p_3_elt_ur_init(K, &z);
    // swap pointers to have n <= m
    mpfq_p_3_src_vec uu, vv;
    if (n <= m) {
        uu = u; vv = v;
    } else {
        uu = v; vv = u;
        unsigned int tmp = n;
        n = m; m = tmp;
    }
    for(k = 0; k < n; ++k) {
        mpfq_p_3_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i <= k; ++i) {
            mpfq_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_p_3_elt_ur_add(K, acc, acc, z);
        }
        mpfq_p_3_elt_ur_set(K, w[k], acc);
    }
    for(k = n; k < m; ++k) {
        mpfq_p_3_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i < n; ++i) {
            mpfq_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_p_3_elt_ur_add(K, acc, acc, z);
        }
        mpfq_p_3_elt_ur_set(K, w[k], acc);
    }
    for(k = m; k < n+m-1; ++k) {
        mpfq_p_3_mul_ur(K, acc, uu[k-m+1], vv[m-1]);
        for(i = k-m+2; i < n; ++i) {
            mpfq_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            mpfq_p_3_elt_ur_add(K, acc, acc, z);
        }
        mpfq_p_3_elt_ur_set(K, w[k], acc);
    }
    mpfq_p_3_elt_ur_clear(K, &acc);
    mpfq_p_3_elt_ur_clear(K, &z);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void mpfq_p_3_vec_reduce(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec w, mpfq_p_3_dst_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_p_3_dst_elt_ur x = mpfq_p_3_vec_ur_coeff_ptr(K, u, i);
        mpfq_p_3_dst_elt y = mpfq_p_3_vec_coeff_ptr(K, w, i);
        mpfq_p_3_reduce(K, y, x);
    }
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_dst_vec_ur mpfq_p_3_vec_ur_subvec(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_src_vec_ur mpfq_p_3_vec_ur_subvec_const(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_dst_elt mpfq_p_3_vec_ur_coeff_ptr(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_vec_ur v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
mpfq_p_3_src_elt mpfq_p_3_vec_ur_coeff_ptr_const(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_vec_ur v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::poly::code_for_poly_init, Mpfq::gfp */
static inline
void mpfq_p_3_poly_init(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_poly p, unsigned int n)
{
    mpfq_p_3_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear, Mpfq::gfp */
static inline
void mpfq_p_3_poly_clear(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_poly p)
{
    mpfq_p_3_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set, Mpfq::gfp */
static inline
void mpfq_p_3_poly_set(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_p_3_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff, Mpfq::gfp */
static inline
void mpfq_p_3_poly_setcoeff(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_elt x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_p_3_vec_set_zero(k, mpfq_p_3_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_p_3_vec_setcoeff(k, w->c, x, i);
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff_ui, Mpfq::gfp */
static inline
void mpfq_p_3_poly_setcoeff_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, unsigned long x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_p_3_vec_set_zero(k, mpfq_p_3_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_p_3_vec_setcoeff_ui(k, w->c, x, i);
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_getcoeff, Mpfq::gfp */
static inline
void mpfq_p_3_poly_getcoeff(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_elt x, mpfq_p_3_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       mpfq_p_3_set_ui(k,x,0);
    } else {
       mpfq_p_3_vec_getcoeff(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg, Mpfq::gfp */
static inline
int mpfq_p_3_poly_deg(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    mpfq_p_3_elt temp;
    mpfq_p_3_init(K, &temp);
    mpfq_p_3_vec_getcoeff(K, temp, w->c, deg);
    int comp=mpfq_p_3_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           mpfq_p_3_vec_getcoeff(K, temp, w->c, deg);
           comp=mpfq_p_3_cmp_ui(K, temp, 0);
        }
    }
    mpfq_p_3_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add, Mpfq::gfp */
static inline
void mpfq_p_3_poly_add(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, mpfq_p_3_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_p_3_vec_add(k, w->c, u->c, v->c, su);
        mpfq_p_3_vec_set(k, mpfq_p_3_vec_subvec(k, w->c, su), mpfq_p_3_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_p_3_vec_add(k, w->c, u->c, v->c, sv);
        mpfq_p_3_vec_set(k, mpfq_p_3_vec_subvec(k, w->c, sv), mpfq_p_3_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_sub, Mpfq::gfp */
static inline
void mpfq_p_3_poly_sub(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, mpfq_p_3_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_p_3_vec_sub(k, w->c, u->c, v->c, su);
        mpfq_p_3_vec_neg(k, mpfq_p_3_vec_subvec(k, w->c, su), mpfq_p_3_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_p_3_vec_sub(k, w->c, u->c, v->c, sv);
        mpfq_p_3_vec_set(k, mpfq_p_3_vec_subvec(k, w->c, sv), mpfq_p_3_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_set_ui, Mpfq::gfp */
static inline
void mpfq_p_3_poly_set_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, unsigned long x)
{
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_p_3_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui, Mpfq::gfp */
static inline
void mpfq_p_3_poly_add_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_p_3_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_p_3_poly_deg(k, w);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_p_3_vec_set(k, mpfq_p_3_vec_subvec(k, w->c, 1), mpfq_p_3_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_p_3_add_ui(k, mpfq_p_3_vec_coeff_ptr(k, w->c, 0), mpfq_p_3_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui, Mpfq::gfp */
static inline
void mpfq_p_3_poly_sub_ui(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_p_3_elt temp;
        mpfq_p_3_init(k, &temp);
        mpfq_p_3_set_ui(k, temp, x);
        mpfq_p_3_neg(k, mpfq_p_3_vec_coeff_ptr(k, w->c, 0), temp);
        w->size = mpfq_p_3_cmp_ui(k, temp, 0);
        mpfq_p_3_clear(k, &temp);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_p_3_vec_set(k, mpfq_p_3_vec_subvec(k, w->c, 1), mpfq_p_3_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_p_3_sub_ui(k, mpfq_p_3_vec_coeff_ptr(k, w->c, 0), mpfq_p_3_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_neg, Mpfq::gfp */
static inline
void mpfq_p_3_poly_neg(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_p_3_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul, Mpfq::gfp */
static inline
void mpfq_p_3_poly_scal_mul(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, mpfq_p_3_src_elt x)
{
    if (mpfq_p_3_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_p_3_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_mul, Mpfq::gfp */
static inline
void mpfq_p_3_poly_mul(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, mpfq_p_3_src_poly u, mpfq_p_3_src_poly v)
{
    unsigned int usize = mpfq_p_3_poly_deg(k, u)+1;
    unsigned int vsize = mpfq_p_3_poly_deg(k, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    mpfq_p_3_vec_conv(k, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
    w->size = 1 + mpfq_p_3_poly_deg(k, w);
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly, Mpfq::gfp */
static inline
void mpfq_p_3_poly_gcd(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly g, mpfq_p_3_src_poly a0, mpfq_p_3_src_poly b0)
{
    mpfq_p_3_poly a,b,q,r;
    int da0=mpfq_p_3_poly_deg(k,a0), db0=mpfq_p_3_poly_deg(k,b0);
    if (db0==-1)
     mpfq_p_3_poly_set(k,g,a0);
    else {
     mpfq_p_3_poly_init(k,a,da0+1);
     mpfq_p_3_poly_init(k,b,db0+1);
     mpfq_p_3_poly_init(k,q,1);
     mpfq_p_3_poly_init(k,r,db0);
     mpfq_p_3_poly_set(k,a,a0);
     mpfq_p_3_poly_set(k,b,b0);
     while (mpfq_p_3_poly_deg(k,b)>=0) {
      mpfq_p_3_poly_divmod(k,q,r,a,b);
      mpfq_p_3_poly_set(k,a,b);
      mpfq_p_3_poly_set(k,b,r); 
     }
     mpfq_p_3_poly_setmonic(k,g,a);
    mpfq_p_3_poly_clear(k,a);
    mpfq_p_3_poly_clear(k,b);
    mpfq_p_3_poly_clear(k,q);
    mpfq_p_3_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly, Mpfq::gfp */
static inline
void mpfq_p_3_poly_xgcd(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly g, mpfq_p_3_dst_poly u0, mpfq_p_3_dst_poly v0, mpfq_p_3_src_poly a0, mpfq_p_3_src_poly b0)
{
    mpfq_p_3_poly a,b,u,v,w,x,q,r;
    mpfq_p_3_elt c;
    mpfq_p_3_init(k,&c);
    mpfq_p_3_set_ui(k,c,0);        /* placate gcc */
    int da0=mpfq_p_3_poly_deg(k,a0), db0=mpfq_p_3_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      mpfq_p_3_poly_set(k,u0,a0);
      mpfq_p_3_poly_set(k,v0,b0);
      mpfq_p_3_poly_set(k,g,a0);
     } else {
      mpfq_p_3_poly_getcoeff(k,c,a0,da0);
      mpfq_p_3_inv(k,c,c);
      mpfq_p_3_poly_scal_mul(k,g,a0,c);
      mpfq_p_3_poly_set(k,v0,b0);
      mpfq_p_3_poly_set(k,u0,b0);
      mpfq_p_3_poly_setcoeff(k,u0,c,0);
     }
    }
    else {
     mpfq_p_3_poly_init(k,a,da0+1);
     mpfq_p_3_poly_init(k,b,db0+1);
     mpfq_p_3_poly_init(k,q,1);
     mpfq_p_3_poly_init(k,r,db0);
     mpfq_p_3_poly_set(k,a,a0);
     mpfq_p_3_poly_set(k,b,b0);
     mpfq_p_3_poly_init(k,u,1);
     mpfq_p_3_poly_init(k,v,1);
     mpfq_p_3_poly_init(k,w,1);
     mpfq_p_3_poly_init(k,x,1);
     mpfq_p_3_poly_setcoeff_ui(k,u,1,0);
     mpfq_p_3_poly_setcoeff_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (mpfq_p_3_poly_deg(k,b)>=0) {
      mpfq_p_3_poly_divmod(k,q,r,a,b);
      mpfq_p_3_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      mpfq_p_3_poly_set(k,b,r);
      mpfq_p_3_poly_mul(k,r,q,w);
      mpfq_p_3_poly_sub(k,r,u,r);
      mpfq_p_3_poly_set(k,u,w);   /* u,w <- w,u-qw */
      mpfq_p_3_poly_set(k,w,r);
      mpfq_p_3_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      mpfq_p_3_poly_sub(k,r,v,r);
      mpfq_p_3_poly_set(k,v,x);
      mpfq_p_3_poly_set(k,x,r);
     }
     dega=mpfq_p_3_poly_deg(k,a);
     ASSERT_FOR_STATIC_ANALYZER(dega >= 0);
     mpfq_p_3_poly_getcoeff(k,c,a,dega);
     mpfq_p_3_inv(k,c,c);
     mpfq_p_3_poly_scal_mul(k,g,a,c);
     mpfq_p_3_poly_scal_mul(k,u0,u,c);
     mpfq_p_3_poly_scal_mul(k,v0,v,c);
     mpfq_p_3_poly_clear(k,a);
     mpfq_p_3_poly_clear(k,b);
     mpfq_p_3_poly_clear(k,u);
     mpfq_p_3_poly_clear(k,v);
     mpfq_p_3_poly_clear(k,w);
     mpfq_p_3_poly_clear(k,x);
     mpfq_p_3_poly_clear(k,q);
     mpfq_p_3_poly_clear(k,r);
    }
    mpfq_p_3_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_random, Mpfq::gfp */
static inline
void mpfq_p_3_poly_random(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_p_3_vec_random(k, w->c, n,state);
    w->size=n;
    int wdeg = mpfq_p_3_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_random2, Mpfq::gfp */
static inline
void mpfq_p_3_poly_random2(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        mpfq_p_3_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_p_3_vec_random2(k, w->c, n,state);
    w->size=n;
    int wdeg = mpfq_p_3_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_cmp, Mpfq::gfp */
static inline
int mpfq_p_3_poly_cmp(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_src_poly u, mpfq_p_3_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return mpfq_p_3_vec_cmp(k, u->c, v->c, u->size);
}

/* *Mpfq::defaults::poly::code_for_poly_asprint, Mpfq::gfp */
static inline
int mpfq_p_3_poly_asprint(mpfq_p_3_dst_field k MAYBE_UNUSED, char * * pstr, mpfq_p_3_src_poly w)
{
    return mpfq_p_3_vec_asprint(k, pstr, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_fprint, Mpfq::gfp */
static inline
int mpfq_p_3_poly_fprint(mpfq_p_3_dst_field k MAYBE_UNUSED, FILE * file, mpfq_p_3_src_poly w)
{
    return mpfq_p_3_vec_fprint(k, file, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_print, Mpfq::gfp */
static inline
int mpfq_p_3_poly_print(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_src_poly w)
{
    return mpfq_p_3_vec_print(k, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_sscan, Mpfq::gfp */
static inline
int mpfq_p_3_poly_sscan(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w, const char * str)
{
    int ret;
    ret = mpfq_p_3_vec_sscan(k, &(w->c), &(w->alloc), str);
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_fscan, Mpfq::gfp */
static inline
int mpfq_p_3_poly_fscan(mpfq_p_3_dst_field k MAYBE_UNUSED, FILE * file, mpfq_p_3_dst_poly w)
{
    int ret;
    ret = mpfq_p_3_vec_fscan(k, file, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_scan, Mpfq::gfp */
static inline
int mpfq_p_3_poly_scan(mpfq_p_3_dst_field k MAYBE_UNUSED, mpfq_p_3_dst_poly w)
{
    int ret;
    ret = mpfq_p_3_vec_scan(k, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *simd_gfp::code_for_simd_hamming_weight */
static inline
int mpfq_p_3_simd_hamming_weight(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_elt p)
{
    return !mpfq_p_3_is_zero(K,p);
}

/* *simd_gfp::code_for_simd_find_first_set */
static inline
int mpfq_p_3_simd_find_first_set(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_elt p)
{
    return mpfq_p_3_simd_hamming_weight(K,p) - 1;
}

/* *simd_gfp::code_for_simd_get_ui_at */
static inline
unsigned long mpfq_p_3_simd_get_ui_at(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_src_elt p, int k MAYBE_UNUSED)
{
    return mpfq_p_3_get_ui(K,p);
}

/* *simd_gfp::code_for_simd_set_ui_at */
static inline
void mpfq_p_3_simd_set_ui_at(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    mpfq_p_3_set_ui(K,p,v);
}

/* *simd_gfp::code_for_simd_add_ui_at */
static inline
void mpfq_p_3_simd_add_ui_at(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt p, mpfq_p_3_src_elt p0, int k MAYBE_UNUSED, unsigned long v)
{
    mpfq_p_3_add_ui(K,p,p0,v);
}

/* *simd_gfp::code_for_simd_set_ui_all */
static inline
void mpfq_p_3_simd_set_ui_all(mpfq_p_3_dst_field K MAYBE_UNUSED, mpfq_p_3_dst_elt p, unsigned long v)
{
    mpfq_p_3_set_ui(K,p,v);
}

/* Mpfq::engine::oo::oo_field_clear */
/* Triggered by: oo */
static inline
void mpfq_p_3_oo_field_clear(mpfq_vbase_ptr f)
{
    mpfq_p_3_field_clear((mpfq_p_3_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* MPFQ_P_3_H_ */

/* vim:set ft=cpp: */
