#ifndef MPFQ_U64K2_H_
#define MPFQ_U64K2_H_

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
#include "mpfq_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      u64k2

/* Active handler: simd_u64k */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_char2 */
/* Options used:{
   family=[
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_m128, tag=m128, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k1, tag=u64k1, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k2, tag=u64k2, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k3, tag=u64k3, },
    { cpp_ifdef=COMPILE_MPFQ_BINARY_FIELD_u64k4, tag=u64k4, },
    ],
   k=2,
   tag=u64k2,
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
     [ (?^:mpfq_u64k2_elt \*), void *, ],
     [ (?^:mpfq_u64k2_src_elt\b), const void *, ],
     [ (?^:mpfq_u64k2_elt\b), void *, ],
     [ (?^:mpfq_u64k2_dst_elt\b), void *, ],
     [ (?^:mpfq_u64k2_elt_ur \*), void *, ],
     [ (?^:mpfq_u64k2_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_u64k2_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k2_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k2_vec \*), void *, ],
     [ (?^:mpfq_u64k2_src_vec\b), const void *, ],
     [ (?^:mpfq_u64k2_vec\b), void *, ],
     [ (?^:mpfq_u64k2_dst_vec\b), void *, ],
     [ (?^:mpfq_u64k2_vec_ur \*), void *, ],
     [ (?^:mpfq_u64k2_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_u64k2_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k2_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k2_poly \*), void *, ],
     [ (?^:mpfq_u64k2_src_poly\b), const void *, ],
     [ (?^:mpfq_u64k2_poly\b), void *, ],
     [ (?^:mpfq_u64k2_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */

typedef struct {
    /* empty struct is not allowed in C. No, really it's not allowed. */
    char c;
} mpfq_u64k2_field_struct;
typedef mpfq_u64k2_field_struct mpfq_u64k2_field [1];
typedef mpfq_u64k2_field_struct * mpfq_u64k2_dst_field;
typedef const mpfq_u64k2_field_struct * mpfq_u64k2_src_field;

typedef uint64_t mpfq_u64k2_elt[2];
typedef uint64_t * mpfq_u64k2_dst_elt;
typedef const uint64_t * mpfq_u64k2_src_elt;

typedef uint64_t mpfq_u64k2_elt_ur[2];
typedef uint64_t * mpfq_u64k2_dst_elt_ur;
typedef const uint64_t * mpfq_u64k2_src_elt_ur;

typedef mpfq_u64k2_elt * mpfq_u64k2_vec;
typedef mpfq_u64k2_elt * mpfq_u64k2_dst_vec;
typedef mpfq_u64k2_elt * mpfq_u64k2_src_vec;

typedef mpfq_u64k2_elt_ur * mpfq_u64k2_vec_ur;
typedef mpfq_u64k2_elt_ur * mpfq_u64k2_dst_vec_ur;
typedef mpfq_u64k2_elt_ur * mpfq_u64k2_src_vec_ur;

typedef struct {
  mpfq_u64k2_vec c;
  unsigned long alloc;
  unsigned long size;
} mpfq_u64k2_poly_struct;
typedef mpfq_u64k2_poly_struct mpfq_u64k2_poly [1];
typedef mpfq_u64k2_poly_struct * mpfq_u64k2_dst_poly;
typedef const mpfq_u64k2_poly_struct * mpfq_u64k2_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name */
#define mpfq_u64k2_impl_name()	"u64k2"
/* *simd_char2::code_for_impl_max_characteristic_bits */
#define mpfq_u64k2_impl_max_characteristic_bits()	2
/* *simd_char2::code_for_impl_max_degree */
#define mpfq_u64k2_impl_max_degree()	1

/* Functions operating on the field structure */
/* *simd_char2::code_for_field_characteristic */
#define mpfq_u64k2_field_characteristic(K, z)	mpz_set_ui(z,2)
mpz_srcptr mpfq_u64k2_field_characteristic_srcptr(mpfq_u64k2_src_field);
/* *simd_char2::code_for_field_degree */
#define mpfq_u64k2_field_degree(K)	1
/* *simd_char2::code_for_field_init */
#define mpfq_u64k2_field_init(f)	((mpfq_u64k2_dst_field) (f))->c=0
/* *simd_char2::code_for_field_clear */
#define mpfq_u64k2_field_clear(K)	/**/
void mpfq_u64k2_field_specify(mpfq_u64k2_dst_field, unsigned long, const void *);
/* *simd_char2::code_for_field_setopt */
#define mpfq_u64k2_field_setopt(f, x, y)	/**/

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, simd_flat, simd_char2 */
#define mpfq_u64k2_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, simd_flat, simd_char2 */
#define mpfq_u64k2_clear(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_stride, simd_flat, simd_char2 */
#define mpfq_u64k2_elt_stride(k)	sizeof(mpfq_u64k2_elt)

/* Elementary assignment functions */
static inline
void mpfq_u64k2_set(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_elt);
static inline
void mpfq_u64k2_set_zero(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt);

/* Assignment of random values */
static inline
void mpfq_u64k2_random(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, gmp_randstate_t);
static inline
void mpfq_u64k2_random2(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
static inline
void mpfq_u64k2_add(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_elt, mpfq_u64k2_src_elt);
/* *simd_char2::code_for_sub */
#define mpfq_u64k2_sub(K, r, s1, s2)	mpfq_u64k2_add(K,r,s1,s2)
/* *simd_char2::code_for_neg */
#define mpfq_u64k2_neg(K, r, s)	mpfq_u64k2_set(K,r,s)
static inline
void mpfq_u64k2_mul(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_elt, mpfq_u64k2_src_elt);
static inline
int mpfq_u64k2_inv(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_elt);

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, simd_flat, simd_char2 */
#define mpfq_u64k2_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, simd_flat, simd_char2 */
#define mpfq_u64k2_elt_ur_clear(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_stride, simd_flat, simd_char2 */
#define mpfq_u64k2_elt_ur_stride(k)	sizeof(mpfq_u64k2_elt_ur)
static inline
void mpfq_u64k2_elt_ur_set(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur, mpfq_u64k2_src_elt_ur);
static inline
void mpfq_u64k2_elt_ur_set_elt(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur, mpfq_u64k2_src_elt);
static inline
void mpfq_u64k2_elt_ur_set_zero(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur);
static inline
void mpfq_u64k2_elt_ur_add(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur, mpfq_u64k2_src_elt_ur, mpfq_u64k2_src_elt_ur);
/* *simd_char2::code_for_elt_ur_neg */
#define mpfq_u64k2_elt_ur_neg(K, r, s)	mpfq_u64k2_elt_ur_set(K,r,s)
/* *simd_char2::code_for_elt_ur_sub */
#define mpfq_u64k2_elt_ur_sub(K, r, s1, s2)	mpfq_u64k2_elt_ur_add(K,r,s1,s2)
static inline
void mpfq_u64k2_mul_ur(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur, mpfq_u64k2_src_elt, mpfq_u64k2_src_elt);
static inline
void mpfq_u64k2_reduce(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_dst_elt_ur);

/* Comparison functions */
static inline
int mpfq_u64k2_cmp(mpfq_u64k2_dst_field, mpfq_u64k2_src_elt, mpfq_u64k2_src_elt);
static inline
int mpfq_u64k2_is_zero(mpfq_u64k2_dst_field, mpfq_u64k2_src_elt);

/* Input/output functions */
int mpfq_u64k2_asprint(mpfq_u64k2_dst_field, char * *, mpfq_u64k2_src_elt);
int mpfq_u64k2_fprint(mpfq_u64k2_dst_field, FILE *, mpfq_u64k2_src_elt);
/* *io::code_for_print */
#define mpfq_u64k2_print(k, x)	mpfq_u64k2_fprint(k,stdout,x)
int mpfq_u64k2_sscan(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, const char *);
int mpfq_u64k2_fscan(mpfq_u64k2_dst_field, FILE *, mpfq_u64k2_dst_elt);
/* *Mpfq::defaults::code_for_scan */
#define mpfq_u64k2_scan(k, x)	mpfq_u64k2_fscan(k,stdin,x)

/* Vector functions */
void mpfq_u64k2_vec_init(mpfq_u64k2_dst_field, mpfq_u64k2_vec *, unsigned long);
void mpfq_u64k2_vec_reinit(mpfq_u64k2_dst_field, mpfq_u64k2_vec *, unsigned long, unsigned long);
void mpfq_u64k2_vec_clear(mpfq_u64k2_dst_field, mpfq_u64k2_vec *, unsigned long);
static inline
void mpfq_u64k2_vec_set(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_set_zero(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, unsigned long);
static inline
void mpfq_u64k2_vec_setcoeff(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_elt, unsigned long);
/* missing vec_setcoeff_ui */
static inline
void mpfq_u64k2_vec_getcoeff(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_add(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_neg(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_rev(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_sub(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, mpfq_u64k2_src_vec, unsigned long);
static inline
void mpfq_u64k2_vec_scal_mul(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, mpfq_u64k2_src_elt, unsigned long);
void mpfq_u64k2_vec_random(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, unsigned long, gmp_randstate_t);
void mpfq_u64k2_vec_random2(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, unsigned long, gmp_randstate_t);
long mpfq_u64k2_vec_cmp(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_is_zero(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
static inline
mpfq_u64k2_dst_vec mpfq_u64k2_vec_subvec(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, long);
static inline
mpfq_u64k2_src_vec mpfq_u64k2_vec_subvec_const(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, long);
static inline
mpfq_u64k2_dst_elt mpfq_u64k2_vec_coeff_ptr(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, long);
static inline
mpfq_u64k2_src_elt mpfq_u64k2_vec_coeff_ptr_const(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, long);
long mpfq_u64k2_vec_asprint(mpfq_u64k2_dst_field, char * *, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_fprint(mpfq_u64k2_dst_field, FILE *, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_print(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_sscan(mpfq_u64k2_dst_field, mpfq_u64k2_vec *, unsigned long *, const char *);
long mpfq_u64k2_vec_fscan(mpfq_u64k2_dst_field, FILE *, mpfq_u64k2_vec *, unsigned long *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec */
#define mpfq_u64k2_vec_scan(K, w, n)	mpfq_u64k2_vec_fscan(K,stdin,w,n)
long mpfq_u64k2_vec_hamming_weight(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_find_first_set(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_simd_hamming_weight(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
long mpfq_u64k2_vec_simd_find_first_set(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec, unsigned long);
void mpfq_u64k2_vec_ur_init(mpfq_u64k2_dst_field, mpfq_u64k2_vec_ur *, unsigned long);
static inline
void mpfq_u64k2_vec_ur_set_zero(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_set_vec(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec, unsigned long);
void mpfq_u64k2_vec_ur_reinit(mpfq_u64k2_dst_field, mpfq_u64k2_vec_ur *, unsigned long, unsigned long);
void mpfq_u64k2_vec_ur_clear(mpfq_u64k2_dst_field, mpfq_u64k2_vec_ur *, unsigned long);
static inline
void mpfq_u64k2_vec_ur_set(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_setcoeff(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_elt_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_getcoeff(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_add(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_sub(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_neg(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_ur_rev(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec_ur, unsigned long);
static inline
void mpfq_u64k2_vec_scal_mul_ur(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, mpfq_u64k2_src_vec, mpfq_u64k2_src_elt, unsigned long);
static inline
void mpfq_u64k2_vec_reduce(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_dst_vec_ur, unsigned long);
static inline
mpfq_u64k2_dst_vec_ur mpfq_u64k2_vec_ur_subvec(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, long);
static inline
mpfq_u64k2_src_vec_ur mpfq_u64k2_vec_ur_subvec_const(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec_ur, long);
static inline
mpfq_u64k2_dst_elt mpfq_u64k2_vec_ur_coeff_ptr(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec_ur, long);
static inline
mpfq_u64k2_src_elt mpfq_u64k2_vec_ur_coeff_ptr_const(mpfq_u64k2_dst_field, mpfq_u64k2_src_vec_ur, long);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, simd_flat, simd_char2 */
#define mpfq_u64k2_vec_elt_stride(k, n)	((n) * mpfq_u64k2_elt_stride((k)))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, simd_flat, simd_char2 */
#define mpfq_u64k2_vec_ur_elt_stride(k, n)	((n) * mpfq_u64k2_elt_ur_stride((k)))

/* Polynomial functions */

/* Functions related to SIMD operation */
/* *simd_u64k::code_for_simd_groupsize */
#define mpfq_u64k2_simd_groupsize(K)	128
static inline
int mpfq_u64k2_simd_hamming_weight(mpfq_u64k2_dst_field, mpfq_u64k2_src_elt);
static inline
int mpfq_u64k2_simd_find_first_set(mpfq_u64k2_dst_field, mpfq_u64k2_src_elt);
static inline
unsigned long mpfq_u64k2_simd_get_ui_at(mpfq_u64k2_dst_field, mpfq_u64k2_src_elt, int);
static inline
void mpfq_u64k2_simd_set_ui_at(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, int, unsigned long);
static inline
void mpfq_u64k2_simd_add_ui_at(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, mpfq_u64k2_src_elt, int, unsigned long);
static inline
void mpfq_u64k2_simd_set_ui_all(mpfq_u64k2_dst_field, mpfq_u64k2_dst_elt, unsigned long);
void mpfq_u64k2_add_dotprod(mpfq_u64k2_dst_field, mpfq_u64k2_dst_vec, mpfq_u64k2_src_vec, mpfq_u64k2_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* Object-oriented interface */
static inline
void mpfq_u64k2_oo_field_clear(mpfq_vbase_ptr);
void mpfq_u64k2_oo_field_init(mpfq_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *simd_flat::code_for_set, simd_char2 */
static inline
void mpfq_u64k2_set(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, mpfq_u64k2_src_elt s)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*s) ; i++) {
            r[i] = s[i];
        }
}

/* *simd_flat::code_for_set_zero, simd_char2 */
static inline
void mpfq_u64k2_set_zero(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r)
{
    memset(r, 0, sizeof(mpfq_u64k2_elt));
}

/* *simd_flat::code_for_random, simd_char2 */
static inline
void mpfq_u64k2_random(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, gmp_randstate_t state)
{
        mpz_t ugly;
        ugly->_mp_d = (mp_limb_t *) r;
        ugly->_mp_alloc = sizeof(mpfq_u64k2_elt) / (sizeof(mp_limb_t));
        ugly->_mp_size = sizeof(mpfq_u64k2_elt) / (sizeof(mp_limb_t));
        mpz_urandomb(ugly, state, mpfq_u64k2_simd_groupsize(K));
}

/* *simd_flat::code_for_random2, simd_char2 */
static inline
void mpfq_u64k2_random2(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, gmp_randstate_t state)
{
        mpz_t ugly;
        ugly->_mp_d = (mp_limb_t *) r;
        ugly->_mp_alloc = sizeof(mpfq_u64k2_elt) / (sizeof(mp_limb_t));
        ugly->_mp_size = sizeof(mpfq_u64k2_elt) / (sizeof(mp_limb_t));
        mpz_rrandomb(ugly, state, mpfq_u64k2_simd_groupsize(K));
}

/* *simd_char2::code_for_add */
static inline
void mpfq_u64k2_add(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, mpfq_u64k2_src_elt s1, mpfq_u64k2_src_elt s2)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] ^ s2[i];
        }
}

/* *simd_char2::code_for_mul */
static inline
void mpfq_u64k2_mul(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, mpfq_u64k2_src_elt s1, mpfq_u64k2_src_elt s2)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] & s2[i];
        }
}

/* *simd_char2::code_for_inv */
static inline
int mpfq_u64k2_inv(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, mpfq_u64k2_src_elt s)
{
        int rc = 0;
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = s[i];
            if (r[i]) rc = 1;
        }
        return rc;
}

/* *simd_flat::code_for_elt_ur_set, simd_char2 */
static inline
void mpfq_u64k2_elt_ur_set(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur r, mpfq_u64k2_src_elt_ur s)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt_ur)/sizeof(*s) ; i++) {
            r[i] = s[i];
        }
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_elt_ur_set_elt(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur r, mpfq_u64k2_src_elt s)
{
    memset(r, 0, sizeof(mpfq_u64k2_elt_ur)); memcpy(r,s,sizeof(mpfq_u64k2_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_elt_ur_set_zero(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur r)
{
    memset(r, 0, sizeof(mpfq_u64k2_elt_ur));
}

/* *simd_char2::code_for_elt_ur_add */
static inline
void mpfq_u64k2_elt_ur_add(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur r, mpfq_u64k2_src_elt_ur s1, mpfq_u64k2_src_elt_ur s2)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] ^ s2[i];
        }
}

/* *simd_char2::code_for_mul_ur */
static inline
void mpfq_u64k2_mul_ur(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur r, mpfq_u64k2_src_elt s1, mpfq_u64k2_src_elt s2)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] & s2[i];
        }
}

/* *simd_char2::code_for_reduce */
static inline
void mpfq_u64k2_reduce(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, mpfq_u64k2_dst_elt_ur s)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*s) ; i++) {
            r[i] = s[i];
        }
}

/* *Mpfq::defaults::flatdata::code_for_cmp, simd_flat, simd_char2 */
static inline
int mpfq_u64k2_cmp(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_elt r, mpfq_u64k2_src_elt s)
{
    return memcmp(r,s,sizeof(mpfq_u64k2_elt));
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, simd_flat, simd_char2 */
static inline
int mpfq_u64k2_is_zero(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_vec_set(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec r, mpfq_u64k2_src_vec s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_u64k2_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_vec_set_zero(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_u64k2_elt));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_setcoeff(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_elt x, unsigned long i)
{
            mpfq_u64k2_dst_elt y = mpfq_u64k2_vec_coeff_ptr(K, w, i);
            mpfq_u64k2_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_getcoeff(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt x, mpfq_u64k2_src_vec w, unsigned long i)
{
            mpfq_u64k2_src_elt y = mpfq_u64k2_vec_coeff_ptr_const(K, w, i);
            mpfq_u64k2_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_add(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_vec u, mpfq_u64k2_src_vec v, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_u64k2_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_neg(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_vec u, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_u64k2_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_rev(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_vec u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_u64k2_elt tmp[1];
    mpfq_u64k2_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_u64k2_set(K, tmp[0], u[i]);
        mpfq_u64k2_set(K, w[i], u[n-1-i]);
        mpfq_u64k2_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_u64k2_set(K, w[nn], u[nn]);
    mpfq_u64k2_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_sub(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_vec u, mpfq_u64k2_src_vec v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_u64k2_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_scal_mul(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_src_vec u, mpfq_u64k2_src_elt c, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_u64k2_src_elt x = mpfq_u64k2_vec_coeff_ptr_const(K, u, i);
        mpfq_u64k2_dst_elt y = mpfq_u64k2_vec_coeff_ptr(K, w, i);
        mpfq_u64k2_mul(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec */
static inline
mpfq_u64k2_dst_vec mpfq_u64k2_vec_subvec(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_u64k2_src_vec mpfq_u64k2_vec_subvec_const(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_vec v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_u64k2_dst_elt mpfq_u64k2_vec_coeff_ptr(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_u64k2_src_elt mpfq_u64k2_vec_coeff_ptr_const(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_vec v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_vec_ur_set_zero(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur r, unsigned long n)
{
    memset(r, 0, n*sizeof(mpfq_u64k2_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_set_vec, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_set_vec(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec u, unsigned long n)
{
            unsigned int i;
            for(i = 0; i < n; ++i) {
                mpfq_u64k2_src_elt x = mpfq_u64k2_vec_coeff_ptr_const(K, u, i);
                mpfq_u64k2_dst_elt_ur y = mpfq_u64k2_vec_ur_coeff_ptr(K, w, i);
                mpfq_u64k2_elt_ur_set_elt(K, y, x);
            }
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, simd_flat, simd_char2 */
static inline
void mpfq_u64k2_vec_ur_set(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur r, mpfq_u64k2_src_vec_ur s, unsigned long n)
{
    if (r != s) memmove(r, s, n*sizeof(mpfq_u64k2_elt_ur));
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_setcoeff, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_setcoeff(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_elt_ur x, unsigned long i)
{
            mpfq_u64k2_dst_elt_ur y = mpfq_u64k2_vec_ur_coeff_ptr(K, w, i);
            mpfq_u64k2_elt_ur_set(K, y, x);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_ur_getcoeff, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_getcoeff(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt_ur x, mpfq_u64k2_src_vec_ur w, unsigned long i)
{
            mpfq_u64k2_src_elt_ur y = mpfq_u64k2_vec_ur_coeff_ptr_const(K, w, i);
            mpfq_u64k2_elt_ur_set(K, x, y);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_add(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec_ur u, mpfq_u64k2_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_u64k2_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_sub(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec_ur u, mpfq_u64k2_src_vec_ur v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_u64k2_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_neg(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        mpfq_u64k2_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_ur_rev(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec_ur u, unsigned long n)
{
    unsigned int nn = n >> 1;
    mpfq_u64k2_elt_ur tmp[1];
    mpfq_u64k2_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        mpfq_u64k2_elt_ur_set(K, tmp[0], u[i]);
        mpfq_u64k2_elt_ur_set(K, w[i], u[n-1-i]);
        mpfq_u64k2_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        mpfq_u64k2_elt_ur_set(K, w[nn], u[nn]);
    mpfq_u64k2_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_scal_mul_ur(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur w, mpfq_u64k2_src_vec u, mpfq_u64k2_src_elt c, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_u64k2_src_elt x = mpfq_u64k2_vec_coeff_ptr_const(K, u, i);
        mpfq_u64k2_dst_elt_ur y = mpfq_u64k2_vec_ur_coeff_ptr(K, w, i);
        mpfq_u64k2_mul_ur(K, y, x, c);
    }
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec */
static inline
void mpfq_u64k2_vec_reduce(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k2_dst_vec_ur u, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_u64k2_dst_elt_ur x = mpfq_u64k2_vec_ur_coeff_ptr(K, u, i);
        mpfq_u64k2_dst_elt y = mpfq_u64k2_vec_coeff_ptr(K, w, i);
        mpfq_u64k2_reduce(K, y, x);
    }
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec */
static inline
mpfq_u64k2_dst_vec_ur mpfq_u64k2_vec_ur_subvec(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec */
static inline
mpfq_u64k2_src_vec_ur mpfq_u64k2_vec_ur_subvec_const(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_vec_ur v, long i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec */
static inline
mpfq_u64k2_dst_elt mpfq_u64k2_vec_ur_coeff_ptr(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_vec_ur v, long i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec */
static inline
mpfq_u64k2_src_elt mpfq_u64k2_vec_ur_coeff_ptr_const(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_vec_ur v, long i)
{
    return v[i];
}

/* *simd_flat::code_for_simd_hamming_weight, simd_char2 */
static inline
int mpfq_u64k2_simd_hamming_weight(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_elt p)
{
        int w = 0;
#if GNUC_VERSION_ATLEAST(3,4,0)
        unsigned long * xp = (unsigned long *) p;
        assert(mpfq_u64k2_elt_stride(K) % sizeof(unsigned long) == 0);
        // parentheses in the divisor are mandatory here, see
        // https://gcc.gnu.org/pipermail/gcc-patches/2020-September/553888.html
        for(size_t b = 0 ; b < mpfq_u64k2_elt_stride(K) / (sizeof(unsigned long)) ; b++) {
            w += __builtin_popcountl(xp[b]);
        }
#else
        int tab[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
        uint8_t * xp = (uint8_t *) p;
        for(size_t b = 0 ; b < mpfq_u64k2_elt_stride(K) ; b++) {
            w += tab[xp[b]&15] + tab[xp[b]>>4];
        }
#endif
        return w;
}

/* *simd_flat::code_for_simd_find_first_set, simd_char2 */
static inline
int mpfq_u64k2_simd_find_first_set(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_elt p)
{
        size_t bmax = mpfq_u64k2_elt_stride(K) / (sizeof(mpfq_u64k2_elt));
        assert(mpfq_u64k2_simd_groupsize(K) % bmax == 0);
        int g = mpfq_u64k2_simd_groupsize(K) / bmax;
        int f = 0;
        for(size_t b = 0 ; b < bmax ; b++, f+=g) {
            if (!p[b]) continue;
#if GNUC_VERSION_ATLEAST(3,4,0)
            unsigned long * xp = (unsigned long *) (p + b);
            // parentheses in the divisor are mandatory here, see
            // https://gcc.gnu.org/pipermail/gcc-patches/2020-September/553888.html
            for(size_t c = 0 ; c < sizeof(mpfq_u64k2_elt) / (sizeof(unsigned long)) ; c++, f += 64) {
                if (!xp[c]) continue;
                return f + __builtin_ctzl(xp[c]);
            }
#else
            uint8_t * xp = (uint8_t *) (p + b);
            int tab[16] = { -1,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0 };
            for(size_t c = 0 ; c < sizeof(mpfq_u64k2_elt) ; c++, f += 8) {
                if (!xp[c]) continue;
                if (xp[c] & 15) return f + tab[xp[c] & 15];
                return f + 4 + tab[xp[c] >> 4];
            }
#endif
            abort();
        }
        return -1;
}

/* *simd_flat::code_for_simd_get_ui_at, simd_char2 */
static inline
unsigned long mpfq_u64k2_simd_get_ui_at(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_src_elt p, int k)
{
        assert(k < mpfq_u64k2_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        return (xp[k/64] & mask) != 0;
}

/* *simd_flat::code_for_simd_set_ui_at, simd_char2 */
static inline
void mpfq_u64k2_simd_set_ui_at(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt p, int k, unsigned long v)
{
        assert(k < mpfq_u64k2_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
}

/* *simd_flat::code_for_simd_add_ui_at, simd_char2 */
static inline
void mpfq_u64k2_simd_add_ui_at(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt p, mpfq_u64k2_src_elt p0, int k, unsigned long v)
{
        mpfq_u64k2_set(K, p, p0);
        assert(k < mpfq_u64k2_simd_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)(v&1)) << (k%64);
        xp[k/64] ^= mask;
}

/* *simd_flat::code_for_simd_set_ui_all, simd_char2 */
static inline
void mpfq_u64k2_simd_set_ui_all(mpfq_u64k2_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_elt r, unsigned long v)
{
        for(unsigned int i = 0 ; i < sizeof(mpfq_u64k2_elt)/sizeof(*r) ; i++) {
            r[i] = -v;
        }
}

/* Mpfq::engine::oo::oo_field_clear */
/* Triggered by: oo */
static inline
void mpfq_u64k2_oo_field_clear(mpfq_vbase_ptr f)
{
    mpfq_u64k2_field_clear((mpfq_u64k2_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* MPFQ_U64K2_H_ */

/* vim:set ft=cpp: */
