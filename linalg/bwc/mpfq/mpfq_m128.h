#ifndef MPFQ_M128_H_
#define MPFQ_M128_H_

/* MPFQ generated file -- do not edit */

#include "mpfq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <x86intrin.h>
#include <stddef.h>
#include <stdio.h>
#include "assert.h"
#include "mpfq_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      m128

/* Active handler: simd_m128 */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: simd_char2 */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Options used:{
   family=[
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
    ],
   k=2,
   tag=m128,
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
     [ (?^:mpfq_m128_elt \*), void *, ],
     [ (?^:mpfq_m128_src_elt\b), const void *, ],
     [ (?^:mpfq_m128_elt\b), void *, ],
     [ (?^:mpfq_m128_dst_elt\b), void *, ],
     [ (?^:mpfq_m128_elt_ur \*), void *, ],
     [ (?^:mpfq_m128_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_m128_elt_ur\b), void *, ],
     [ (?^:mpfq_m128_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_m128_vec \*), void *, ],
     [ (?^:mpfq_m128_src_vec\b), const void *, ],
     [ (?^:mpfq_m128_vec\b), void *, ],
     [ (?^:mpfq_m128_dst_vec\b), void *, ],
     [ (?^:mpfq_m128_vec_ur \*), void *, ],
     [ (?^:mpfq_m128_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_m128_vec_ur\b), void *, ],
     [ (?^:mpfq_m128_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_m128_poly \*), void *, ],
     [ (?^:mpfq_m128_src_poly\b), const void *, ],
     [ (?^:mpfq_m128_poly\b), void *, ],
     [ (?^:mpfq_m128_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */

typedef struct {
    /* empty struct is not allowed in C. No, really it's not allowed. */
    char c;
} mpfq_m128_field_struct;
typedef mpfq_m128_field_struct mpfq_m128_field [1];
typedef mpfq_m128_field_struct * mpfq_m128_dst_field;
typedef const mpfq_m128_field_struct * mpfq_m128_src_field;

typedef __m128i mpfq_m128_elt[1];
typedef __m128i * mpfq_m128_dst_elt;
typedef const __m128i * mpfq_m128_src_elt;

typedef __m128i mpfq_m128_elt_ur[1];
typedef __m128i * mpfq_m128_dst_elt_ur;
typedef const __m128i * mpfq_m128_src_elt_ur;

typedef mpfq_m128_elt * mpfq_m128_vec;
typedef mpfq_m128_elt * mpfq_m128_dst_vec;
typedef mpfq_m128_elt * mpfq_m128_src_vec;

typedef mpfq_m128_elt_ur * mpfq_m128_vec_ur;
typedef mpfq_m128_elt_ur * mpfq_m128_dst_vec_ur;
typedef mpfq_m128_elt_ur * mpfq_m128_src_vec_ur;

typedef struct {
  mpfq_m128_vec c;
  unsigned long alloc;
  unsigned long size;
} mpfq_m128_poly_struct;
typedef mpfq_m128_poly_struct mpfq_m128_poly [1];
typedef mpfq_m128_poly_struct * mpfq_m128_dst_poly;
typedef const mpfq_m128_poly_struct * mpfq_m128_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name */
#define mpfq_m128_impl_name()	"m128"
/* *simd_char2::code_for_impl_max_characteristic_bits */
#define mpfq_m128_impl_max_characteristic_bits()	2
/* *simd_char2::code_for_impl_max_degree */
#define mpfq_m128_impl_max_degree()	1

/* Functions operating on the field structure */
/* *simd_char2::code_for_field_characteristic */
#define mpfq_m128_field_characteristic(K, z)	mpz_set_ui(z,2)
mpz_srcptr mpfq_m128_field_characteristic_srcptr(mpfq_m128_src_field);
/* *simd_char2::code_for_field_degree */
#define mpfq_m128_field_degree(K)	1
/* *simd_char2::code_for_field_init */
#define mpfq_m128_field_init(f)	((mpfq_m128_dst_field) (f))->c=0
/* *simd_char2::code_for_field_clear */
#define mpfq_m128_field_clear(K)	/**/
void mpfq_m128_field_specify(mpfq_m128_dst_field, unsigned long, const void *);
/* *simd_char2::code_for_field_setopt */
#define mpfq_m128_field_setopt(f, x, y)	/**/

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, simd_flat, simd_char2 */
#define mpfq_m128_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, simd_flat, simd_char2 */
#define mpfq_m128_clear(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_stride, simd_flat, simd_char2 */
#define mpfq_m128_elt_stride(k)	sizeof(mpfq_m128_elt)

/* Elementary assignment functions */
static inline
void mpfq_m128_set(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt);
static inline
void mpfq_m128_set_zero(mpfq_m128_dst_field, mpfq_m128_dst_elt);

/* Assignment of random values */
static inline
void mpfq_m128_random(mpfq_m128_dst_field, mpfq_m128_dst_elt, gmp_randstate_t);
static inline
void mpfq_m128_random2(mpfq_m128_dst_field, mpfq_m128_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
static inline
void mpfq_m128_add(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt, mpfq_m128_src_elt);
static inline
void mpfq_m128_sub(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt, mpfq_m128_src_elt);
static inline
void mpfq_m128_neg(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt);
static inline
void mpfq_m128_mul(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt, mpfq_m128_src_elt);
static inline
int mpfq_m128_inv(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt);

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, simd_flat, simd_char2 */
#define mpfq_m128_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, simd_flat, simd_char2 */
#define mpfq_m128_elt_ur_clear(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_stride, simd_flat, simd_char2 */
#define mpfq_m128_elt_ur_stride(k)	sizeof(mpfq_m128_elt_ur)
static inline
void mpfq_m128_elt_ur_set(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt_ur);
static inline
void mpfq_m128_elt_ur_set_elt(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt);
static inline
void mpfq_m128_elt_ur_set_zero(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur);
static inline
void mpfq_m128_elt_ur_add(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt_ur, mpfq_m128_src_elt_ur);
static inline
void mpfq_m128_elt_ur_neg(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt_ur);
static inline
void mpfq_m128_elt_ur_sub(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt_ur, mpfq_m128_src_elt_ur);
static inline
void mpfq_m128_mul_ur(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_elt, mpfq_m128_src_elt);
static inline
void mpfq_m128_reduce(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_dst_elt_ur);

/* Comparison functions */
static inline
int mpfq_m128_cmp(mpfq_m128_dst_field, mpfq_m128_src_elt, mpfq_m128_src_elt);
static inline
int mpfq_m128_is_zero(mpfq_m128_dst_field, mpfq_m128_src_elt);

/* Input/output functions */
int mpfq_m128_asprint(mpfq_m128_dst_field, char * *, mpfq_m128_src_elt);
int mpfq_m128_fprint(mpfq_m128_dst_field, FILE *, mpfq_m128_src_elt);
/* *io::code_for_print */
#define mpfq_m128_print(k, x)	mpfq_m128_fprint(k,stdout,x)
int mpfq_m128_sscan(mpfq_m128_dst_field, mpfq_m128_dst_elt, const char *);
int mpfq_m128_fscan(mpfq_m128_dst_field, FILE *, mpfq_m128_dst_elt);
/* *Mpfq::defaults::code_for_scan */
#define mpfq_m128_scan(k, x)	mpfq_m128_fscan(k,stdin,x)

/* Vector functions */
void mpfq_m128_vec_init(mpfq_m128_dst_field, mpfq_m128_vec *, unsigned long);
void mpfq_m128_vec_reinit(mpfq_m128_dst_field, mpfq_m128_vec *, unsigned long, unsigned long);
void mpfq_m128_vec_clear(mpfq_m128_dst_field, mpfq_m128_vec *, unsigned long);
static inline
void mpfq_m128_vec_set(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_set_zero(mpfq_m128_dst_field, mpfq_m128_dst_vec, unsigned long);
static inline
void mpfq_m128_vec_setcoeff(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_elt, unsigned long);
/* missing vec_setcoeff_ui */
static inline
void mpfq_m128_vec_getcoeff(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_add(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_neg(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_rev(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_sub(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, mpfq_m128_src_vec, unsigned long);
static inline
void mpfq_m128_vec_scal_mul(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, mpfq_m128_src_elt, unsigned long);
void mpfq_m128_vec_random(mpfq_m128_dst_field, mpfq_m128_dst_vec, unsigned long, gmp_randstate_t);
void mpfq_m128_vec_random2(mpfq_m128_dst_field, mpfq_m128_dst_vec, unsigned long, gmp_randstate_t);
long mpfq_m128_vec_cmp(mpfq_m128_dst_field, mpfq_m128_src_vec, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_is_zero(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
static inline
mpfq_m128_dst_vec mpfq_m128_vec_subvec(mpfq_m128_dst_field, mpfq_m128_dst_vec, long);
static inline
mpfq_m128_src_vec mpfq_m128_vec_subvec_const(mpfq_m128_dst_field, mpfq_m128_src_vec, long);
static inline
mpfq_m128_dst_elt mpfq_m128_vec_coeff_ptr(mpfq_m128_dst_field, mpfq_m128_dst_vec, long);
static inline
mpfq_m128_src_elt mpfq_m128_vec_coeff_ptr_const(mpfq_m128_dst_field, mpfq_m128_src_vec, long);
long mpfq_m128_vec_asprint(mpfq_m128_dst_field, char * *, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_fprint(mpfq_m128_dst_field, FILE *, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_print(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_sscan(mpfq_m128_dst_field, mpfq_m128_vec *, unsigned long *, const char *);
long mpfq_m128_vec_fscan(mpfq_m128_dst_field, FILE *, mpfq_m128_vec *, unsigned long *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec */
#define mpfq_m128_vec_scan(K, w, n)	mpfq_m128_vec_fscan(K,stdin,w,n)
long mpfq_m128_vec_hamming_weight(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_find_first_set(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_simd_hamming_weight(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
long mpfq_m128_vec_simd_find_first_set(mpfq_m128_dst_field, mpfq_m128_src_vec, unsigned long);
void mpfq_m128_vec_ur_init(mpfq_m128_dst_field, mpfq_m128_vec_ur *, unsigned long);
static inline
void mpfq_m128_vec_ur_set_zero(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_set_vec(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec, unsigned long);
void mpfq_m128_vec_ur_reinit(mpfq_m128_dst_field, mpfq_m128_vec_ur *, unsigned long, unsigned long);
void mpfq_m128_vec_ur_clear(mpfq_m128_dst_field, mpfq_m128_vec_ur *, unsigned long);
static inline
void mpfq_m128_vec_ur_set(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_setcoeff(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_elt_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_getcoeff(mpfq_m128_dst_field, mpfq_m128_dst_elt_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_add(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_sub(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_neg(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_ur_rev(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec_ur, unsigned long);
static inline
void mpfq_m128_vec_scal_mul_ur(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, mpfq_m128_src_vec, mpfq_m128_src_elt, unsigned long);
static inline
void mpfq_m128_vec_reduce(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_dst_vec_ur, unsigned long);
static inline
mpfq_m128_dst_vec_ur mpfq_m128_vec_ur_subvec(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, long);
static inline
mpfq_m128_src_vec_ur mpfq_m128_vec_ur_subvec_const(mpfq_m128_dst_field, mpfq_m128_src_vec_ur, long);
static inline
mpfq_m128_dst_elt mpfq_m128_vec_ur_coeff_ptr(mpfq_m128_dst_field, mpfq_m128_dst_vec_ur, long);
static inline
mpfq_m128_src_elt mpfq_m128_vec_ur_coeff_ptr_const(mpfq_m128_dst_field, mpfq_m128_src_vec_ur, long);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, simd_flat, simd_char2 */
#define mpfq_m128_vec_elt_stride(k, n)	((n) * mpfq_m128_elt_stride((k)))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, simd_flat, simd_char2 */
#define mpfq_m128_vec_ur_elt_stride(k, n)	((n) * mpfq_m128_elt_ur_stride((k)))

/* Polynomial functions */

/* Functions related to SIMD operation */
/* *simd_m128::code_for_simd_groupsize */
#define mpfq_m128_simd_groupsize(K)	128
static inline
int mpfq_m128_simd_hamming_weight(mpfq_m128_dst_field, mpfq_m128_src_elt);
static inline
int mpfq_m128_simd_find_first_set(mpfq_m128_dst_field, mpfq_m128_src_elt);
static inline
unsigned long mpfq_m128_simd_get_ui_at(mpfq_m128_dst_field, mpfq_m128_src_elt, int);
static inline
void mpfq_m128_simd_set_ui_at(mpfq_m128_dst_field, mpfq_m128_dst_elt, int, unsigned long);
static inline
void mpfq_m128_simd_add_ui_at(mpfq_m128_dst_field, mpfq_m128_dst_elt, mpfq_m128_src_elt, int, unsigned long);
static inline
void mpfq_m128_simd_set_ui_all(mpfq_m128_dst_field, mpfq_m128_dst_elt, unsigned long);
void mpfq_m128_add_dotprod(mpfq_m128_dst_field, mpfq_m128_dst_vec, mpfq_m128_src_vec, mpfq_m128_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* Object-oriented interface */
static inline
void mpfq_m128_oo_field_clear(mpfq_vbase_ptr);
void mpfq_m128_oo_field_init(mpfq_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *simd_m128::code_for_set */
static inline
void mpfq_m128_set(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s)
{
    *r=*s;
}

/* *simd_m128::code_for_set_zero */
static inline
void mpfq_m128_set_zero(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r)
{
    *r=_mm_setzero_si128();
}

/* *simd_flat::code_for_random, simd_char2 */
static inline
void mpfq_m128_random(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, gmp_randstate_t state)
{
        mpz_t ugly;
        ugly->_mp_d = (mp_limb_t *) r;
        ugly->_mp_alloc = sizeof(mpfq_m128_elt) / (sizeof(mp_limb_t));
        ugly->_mp_size = sizeof(mpfq_m128_elt) / (sizeof(mp_limb_t));
        mpz_urandomb(ugly, state, mpfq_m128_simd_groupsize(K));
}

/* *simd_flat::code_for_random2, simd_char2 */
static inline
void mpfq_m128_random2(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, gmp_randstate_t state)
{
        mpz_t ugly;
        ugly->_mp_d = (mp_limb_t *) r;
        ugly->_mp_alloc = sizeof(mpfq_m128_elt) / (sizeof(mp_limb_t));
        ugly->_mp_size = sizeof(mpfq_m128_elt) / (sizeof(mp_limb_t));
        mpz_rrandomb(ugly, state, mpfq_m128_simd_groupsize(K));
}

/* *simd_m128::code_for_add */
static inline
void mpfq_m128_add(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s1, mpfq_m128_src_elt s2)
{
    *r = _mm_xor_si128(*s1, *s2);
}

/* *simd_m128::code_for_sub */
static inline
void mpfq_m128_sub(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s1, mpfq_m128_src_elt s2)
{
    *r = _mm_xor_si128(*s1, *s2);
}

/* *simd_m128::code_for_neg */
static inline
void mpfq_m128_neg(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s)
{
    *r=*s;
}

/* *simd_m128::code_for_mul */
static inline
void mpfq_m128_mul(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s1, mpfq_m128_src_elt s2)
{
    *r = _mm_and_si128(*s1, *s2);
}

/* *simd_m128::code_for_inv */
static inline
int mpfq_m128_inv(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_src_elt s)
{
        *r = *s;
        return 1;
}

/* *simd_m128::code_for_elt_ur_set */
static inline
void mpfq_m128_elt_ur_set(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt_ur s)
{
    *r=*s;
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, simd_flat, simd_char2 */
static inline
void mpfq_m128_elt_ur_set_elt(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt s)
{
    memset(r, 0, sizeof(mpfq_m128_elt_ur)); memcpy(r,s,sizeof(mpfq_m128_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, simd_flat, simd_char2 */
static inline
void mpfq_m128_elt_ur_set_zero(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r)
{
    memset(r, 0, sizeof(mpfq_m128_elt_ur));
}

/* *simd_m128::code_for_elt_ur_add */
static inline
void mpfq_m128_elt_ur_add(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt_ur s1, mpfq_m128_src_elt_ur s2)
{
    *r = _mm_xor_si128(*s1, *s2);
}

/* *simd_m128::code_for_elt_ur_neg */
static inline
void mpfq_m128_elt_ur_neg(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt_ur s)
{
    *r=*s;
}

/* *simd_m128::code_for_elt_ur_sub */
static inline
void mpfq_m128_elt_ur_sub(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt_ur s1, mpfq_m128_src_elt_ur s2)
{
    *r = _mm_xor_si128(*s1, *s2);
}

/* *simd_m128::code_for_mul_ur */
static inline
void mpfq_m128_mul_ur(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur r, mpfq_m128_src_elt s1, mpfq_m128_src_elt s2)
{
    *r = _mm_and_si128(*s1, *s2);
}

/* *simd_m128::code_for_reduce */
static inline
void mpfq_m128_reduce(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, mpfq_m128_dst_elt_ur s)
{
    *r=*s;
}

/* *Mpfq::defaults::flatdata::code_for_cmp, simd_flat, simd_char2 */
static inline
int mpfq_m128_cmp(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_elt r, mpfq_m128_src_elt s)
{
    return memcmp(r,s,sizeof(mpfq_m128_elt));
}

/* *simd_m128::code_for_is_zero */
static inline
int mpfq_m128_is_zero(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_elt r)
{
            return _mm_extract_epi64(*r, 0) == 0 && 
                   _mm_extract_epi64(*r, 1) == 0;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_m128_vec_set(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec r, mpfq_m128_src_vec s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_m128_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_m128_vec_set_zero(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_m128_elt));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_setcoeff(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_elt x, unsigned long i)
{
            mpfq_m128_dst_elt y = mpfq_m128_vec_coeff_ptr(K, w, i);
            mpfq_m128_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_getcoeff(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt x, mpfq_m128_src_vec w, unsigned long i)
{
            mpfq_m128_src_elt y = mpfq_m128_vec_coeff_ptr_const(K, w, i);
            mpfq_m128_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_add(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_vec u, mpfq_m128_src_vec v, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_m128_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_neg(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_vec u, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_m128_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_rev(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_vec u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_m128_elt tmp[1];
    mpfq_m128_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_m128_set(K, tmp[0], u[i]);
        mpfq_m128_set(K, w[i], u[n-1-i]);
        mpfq_m128_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_m128_set(K, w[nn], u[nn]);
    mpfq_m128_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_sub(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_vec u, mpfq_m128_src_vec v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_m128_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_scal_mul(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_src_vec u, mpfq_m128_src_elt c, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_m128_src_elt x = mpfq_m128_vec_coeff_ptr_const(K, u, i);
        mpfq_m128_dst_elt y = mpfq_m128_vec_coeff_ptr(K, w, i);
        mpfq_m128_mul(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec */
static inline
mpfq_m128_dst_vec mpfq_m128_vec_subvec(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_m128_src_vec mpfq_m128_vec_subvec_const(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_m128_dst_elt mpfq_m128_vec_coeff_ptr(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_m128_src_elt mpfq_m128_vec_coeff_ptr_const(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_m128_vec_ur_set_zero(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_m128_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_set_vec, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_set_vec(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec u, unsigned long n)
{
            unsigned int i;
            for(i = 0; i < n; ++i) {
                mpfq_m128_src_elt x = mpfq_m128_vec_coeff_ptr_const(K, u, i);
                mpfq_m128_dst_elt_ur y = mpfq_m128_vec_ur_coeff_ptr(K, w, i);
                mpfq_m128_elt_ur_set_elt(K, y, x);
            }
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_m128_vec_ur_set(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur r, mpfq_m128_src_vec_ur s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_m128_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_setcoeff(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_elt_ur x, unsigned long i)
{
            mpfq_m128_dst_elt_ur y = mpfq_m128_vec_ur_coeff_ptr(K, w, i);
            mpfq_m128_elt_ur_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_getcoeff(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt_ur x, mpfq_m128_src_vec_ur w, unsigned long i)
{
            mpfq_m128_src_elt_ur y = mpfq_m128_vec_ur_coeff_ptr_const(K, w, i);
            mpfq_m128_elt_ur_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_add(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec_ur u, mpfq_m128_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_m128_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_sub(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec_ur u, mpfq_m128_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_m128_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_neg(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_m128_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_ur_rev(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec_ur u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_m128_elt_ur tmp[1];
    mpfq_m128_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_m128_elt_ur_set(K, tmp[0], u[i]);
        mpfq_m128_elt_ur_set(K, w[i], u[n-1-i]);
        mpfq_m128_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_m128_elt_ur_set(K, w[nn], u[nn]);
    mpfq_m128_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_scal_mul_ur(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur w, mpfq_m128_src_vec u, mpfq_m128_src_elt c, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_m128_src_elt x = mpfq_m128_vec_coeff_ptr_const(K, u, i);
        mpfq_m128_dst_elt_ur y = mpfq_m128_vec_ur_coeff_ptr(K, w, i);
        mpfq_m128_mul_ur(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec */
static inline
void mpfq_m128_vec_reduce(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_m128_dst_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_m128_dst_elt_ur x = mpfq_m128_vec_ur_coeff_ptr(K, u, i);
        mpfq_m128_dst_elt y = mpfq_m128_vec_coeff_ptr(K, w, i);
        mpfq_m128_reduce(K, y, x);
    }
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec */
static inline
mpfq_m128_dst_vec_ur mpfq_m128_vec_ur_subvec(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_m128_src_vec_ur mpfq_m128_vec_ur_subvec_const(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_m128_dst_elt mpfq_m128_vec_ur_coeff_ptr(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_vec_ur v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_m128_src_elt mpfq_m128_vec_ur_coeff_ptr_const(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_vec_ur v, long i)
{
    return v[i];
}

/* *simd_m128::code_for_simd_hamming_weight */
static inline
int mpfq_m128_simd_hamming_weight(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_elt r)
{
            return _mm_popcnt_u64(_mm_extract_epi64(*r, 0)) +
                   _mm_popcnt_u64(_mm_extract_epi64(*r, 1));
}

/* *simd_m128::code_for_simd_find_first_set */
static inline
int mpfq_m128_simd_find_first_set(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_elt p)
{
        int f = 0;
#if GNUC_VERSION_ATLEAST(3,4,0)
            unsigned long * xp = (unsigned long *) p;
            for(size_t c = 0 ; c < 2 ; c++, f += 64) {
                if (!xp[c]) continue;
                return f + __builtin_ctzl(xp[c]);
            }
#else
            uint8_t * xp = (uint8_t *) p;
            int tab[16] = { -1,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0 };
            for(size_t c = 0 ; c < 16 ; c++, f += 8) {
                if (!xp[c]) continue;
                if (xp[c] & 15) return f + tab[xp[c] & 15];
                return f + 4 + tab[xp[c] >> 4];
            }
#endif
        return -1;
}

/* *simd_flat::code_for_simd_get_ui_at, simd_char2 */
static inline
unsigned long mpfq_m128_simd_get_ui_at(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_src_elt p, int k)
{
        assert(k < mpfq_m128_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        return (xp[k/64] & mask) != 0;
}

/* *simd_flat::code_for_simd_set_ui_at, simd_char2 */
static inline
void mpfq_m128_simd_set_ui_at(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt p, int k, unsigned long v)
{
        assert(k < mpfq_m128_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
}

/* *simd_flat::code_for_simd_add_ui_at, simd_char2 */
static inline
void mpfq_m128_simd_add_ui_at(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt p, mpfq_m128_src_elt p0, int k, unsigned long v)
{
        mpfq_m128_set(K, p, p0);
        assert(k < mpfq_m128_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)(v&1)) << (k%64);
        xp[k/64] ^= mask;
}

/* *simd_m128::code_for_simd_set_ui_all */
static inline
void mpfq_m128_simd_set_ui_all(mpfq_m128_dst_field K MAYBE_UNUSED, mpfq_m128_dst_elt r, unsigned long v)
{
        *r = _mm_set1_epi64(_mm_cvtsi64_m64(-(v != 0)));
}

/* Mpfq::engine::oo::oo_field_clear */
/* Triggered by: oo */
static inline
void mpfq_m128_oo_field_clear(mpfq_vbase_ptr f)
{
    mpfq_m128_field_clear((mpfq_m128_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* MPFQ_M128_H_ */

/* vim:set ft=cpp: */
