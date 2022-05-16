#ifndef MPFQ_PZ_H_
#define MPFQ_PZ_H_

/* MPFQ generated file -- do not edit */

#include "mpfq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include "mpfq_gfp_common.h"
#include <limits.h>
#include "mpfq_fixmp.h"
#include "mpfq_gfp_common.h"
#include <stddef.h>
#include <stdio.h>
#include "assert.h"
#include "mpfq_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      pz

/* Active handler: simd_pz */
/* Automatically generated code  */
/* Active handler: pz */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::poly */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
   fieldtype=prime,
   n=mpz_size(k->p),
   nn=(2*mpz_size(k->p) + 1),
   tag=pz,
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
     [ (?^:mpfq_pz_elt \*), void *, ],
     [ (?^:mpfq_pz_src_elt\b), const void *, ],
     [ (?^:mpfq_pz_elt\b), void *, ],
     [ (?^:mpfq_pz_dst_elt\b), void *, ],
     [ (?^:mpfq_pz_elt_ur \*), void *, ],
     [ (?^:mpfq_pz_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_pz_elt_ur\b), void *, ],
     [ (?^:mpfq_pz_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_pz_vec \*), void *, ],
     [ (?^:mpfq_pz_src_vec\b), const void *, ],
     [ (?^:mpfq_pz_vec\b), void *, ],
     [ (?^:mpfq_pz_dst_vec\b), void *, ],
     [ (?^:mpfq_pz_vec_ur \*), void *, ],
     [ (?^:mpfq_pz_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_pz_vec_ur\b), void *, ],
     [ (?^:mpfq_pz_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_pz_poly \*), void *, ],
     [ (?^:mpfq_pz_src_poly\b), const void *, ],
     [ (?^:mpfq_pz_poly\b), void *, ],
     [ (?^:mpfq_pz_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */

typedef mpfq_p_field mpfq_pz_field;
typedef mpfq_p_dst_field mpfq_pz_dst_field;
typedef mpfq_p_src_field mpfq_pz_src_field;

typedef mp_limb_t * mpfq_pz_elt;
typedef mp_limb_t * mpfq_pz_dst_elt;
typedef const mp_limb_t * mpfq_pz_src_elt;

typedef mp_limb_t * mpfq_pz_elt_ur;
typedef mp_limb_t * mpfq_pz_dst_elt_ur;
typedef const mp_limb_t * mpfq_pz_src_elt_ur;

typedef mp_limb_t * mpfq_pz_vec;
typedef mp_limb_t * mpfq_pz_dst_vec;
typedef const mp_limb_t * mpfq_pz_src_vec;

typedef mp_limb_t * mpfq_pz_vec_ur;
typedef mp_limb_t * mpfq_pz_dst_vec_ur;
typedef const mp_limb_t * mpfq_pz_src_vec_ur;

typedef struct {
  mpfq_pz_vec c;
  unsigned long alloc;
  unsigned long size;
} mpfq_pz_poly_struct;
typedef mpfq_pz_poly_struct mpfq_pz_poly [1];
typedef mpfq_pz_poly_struct * mpfq_pz_dst_poly;
typedef mpfq_pz_poly_struct * mpfq_pz_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name, pz */
#define mpfq_pz_impl_name()	"pz"
/* *pz::code_for_impl_max_characteristic_bits */
#define mpfq_pz_impl_max_characteristic_bits()	UINT_MAX
/* *pz::code_for_impl_max_degree */
#define mpfq_pz_impl_max_degree()	1

/* Functions operating on the field structure */
static inline
void mpfq_pz_field_characteristic(mpfq_pz_src_field, mpz_ptr);
static inline
mpz_srcptr mpfq_pz_field_characteristic_srcptr(mpfq_pz_src_field);
static inline
unsigned long mpfq_pz_field_characteristic_bits(mpfq_pz_src_field);
/* *pz::code_for_field_degree */
#define mpfq_pz_field_degree(k)	1
static inline
void mpfq_pz_field_init(mpfq_pz_dst_field);
void mpfq_pz_field_clear(mpfq_pz_dst_field);
void mpfq_pz_field_specify(mpfq_pz_dst_field, unsigned long, const void *);
/* *pz::code_for_field_setopt */
#define mpfq_pz_field_setopt(k, x, y)	/**/

/* Element allocation functions */
void mpfq_pz_init(mpfq_pz_dst_field, mpfq_pz_elt *);
void mpfq_pz_clear(mpfq_pz_dst_field, mpfq_pz_elt *);
static inline
ptrdiff_t mpfq_pz_elt_stride(mpfq_pz_dst_field);

/* Elementary assignment functions */
static inline
void mpfq_pz_set(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static inline
void mpfq_pz_set_ui(mpfq_pz_dst_field, mpfq_pz_dst_elt, unsigned long);
static inline
void mpfq_pz_set_zero(mpfq_pz_dst_field, mpfq_pz_dst_elt);
static inline
unsigned long mpfq_pz_get_ui(mpfq_pz_dst_field, mpfq_pz_src_elt);
void mpfq_pz_set_mpn(mpfq_pz_dst_field, mpfq_pz_dst_elt, mp_limb_t *, size_t);
void mpfq_pz_set_mpz(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpz_srcptr);
static inline
void mpfq_pz_get_mpn(mpfq_pz_dst_field, mp_limb_t *, mpfq_pz_src_elt);
static inline
void mpfq_pz_get_mpz(mpfq_pz_dst_field, mpz_ptr, mpfq_pz_src_elt);

/* Assignment of random values */
void mpfq_pz_random(mpfq_pz_dst_field, mpfq_pz_dst_elt, gmp_randstate_t);
void mpfq_pz_random2(mpfq_pz_dst_field, mpfq_pz_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
void mpfq_pz_add(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
void mpfq_pz_sub(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
void mpfq_pz_neg(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt);
void mpfq_pz_mul(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
void mpfq_pz_sqr(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt);
int mpfq_pz_is_sqr(mpfq_pz_dst_field, mpfq_pz_src_elt);
void mpfq_pz_pow(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long *, size_t);
void mpfq_pz_powz(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpz_srcptr);
/* *pz::code_for_frobenius */
#define mpfq_pz_frobenius(k, x, y)	mpfq_pz_set(k,x,y)
void mpfq_pz_add_ui(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
void mpfq_pz_sub_ui(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
void mpfq_pz_mul_ui(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
#define HAVE_mpfq_pz_normalize
void mpfq_pz_normalize(mpfq_pz_dst_field, mpfq_pz_dst_elt);
int mpfq_pz_inv(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt);

/* Operations involving unreduced elements */
void mpfq_pz_elt_ur_init(mpfq_pz_dst_field, mpfq_pz_elt_ur *);
void mpfq_pz_elt_ur_clear(mpfq_pz_dst_field, mpfq_pz_elt_ur *);
static inline
ptrdiff_t mpfq_pz_elt_ur_stride(mpfq_pz_dst_field);
static inline
void mpfq_pz_elt_ur_set(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur);
static inline
void mpfq_pz_elt_ur_set_elt(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt);
static inline
void mpfq_pz_elt_ur_set_zero(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur);
static inline
void mpfq_pz_elt_ur_set_ui(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, unsigned long);
void mpfq_pz_elt_ur_add(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur, mpfq_pz_src_elt_ur);
void mpfq_pz_elt_ur_neg(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur);
void mpfq_pz_elt_ur_sub(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur, mpfq_pz_src_elt_ur);
void mpfq_pz_mul_ur(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt, mpfq_pz_src_elt);
void mpfq_pz_sqr_ur(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt);
void mpfq_pz_reduce(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_dst_elt_ur);
#define HAVE_mpfq_pz_addmul_si_ur
static inline
void mpfq_pz_addmul_si_ur(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt, long);

/* Comparison functions */
static inline
int mpfq_pz_cmp(mpfq_pz_dst_field, mpfq_pz_src_elt, mpfq_pz_src_elt);
static inline
int mpfq_pz_cmp_ui(mpfq_pz_dst_field, mpfq_pz_src_elt, unsigned long);
static inline
int mpfq_pz_is_zero(mpfq_pz_dst_field, mpfq_pz_src_elt);

/* Input/output functions */
int mpfq_pz_asprint(mpfq_pz_dst_field, char * *, mpfq_pz_src_elt);
int mpfq_pz_fprint(mpfq_pz_dst_field, FILE *, mpfq_pz_src_elt);
/* *Mpfq::defaults::code_for_print, pz */
#define mpfq_pz_print(k, x)	mpfq_pz_fprint(k,stdout,x)
int mpfq_pz_sscan(mpfq_pz_dst_field, mpfq_pz_dst_elt, const char *);
int mpfq_pz_fscan(mpfq_pz_dst_field, FILE *, mpfq_pz_dst_elt);
/* *Mpfq::defaults::code_for_scan, pz */
#define mpfq_pz_scan(k, x)	mpfq_pz_fscan(k,stdin,x)

/* Vector functions */
void mpfq_pz_vec_init(mpfq_pz_dst_field, mpfq_pz_vec *, unsigned long);
void mpfq_pz_vec_reinit(mpfq_pz_dst_field, mpfq_pz_vec *, unsigned long, unsigned long);
void mpfq_pz_vec_clear(mpfq_pz_dst_field, mpfq_pz_vec *, unsigned long);
void mpfq_pz_vec_set(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_set_zero(mpfq_pz_dst_field, mpfq_pz_dst_vec, unsigned long);
void mpfq_pz_vec_setcoeff(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_elt, unsigned long);
void mpfq_pz_vec_setcoeff_ui(mpfq_pz_dst_field, mpfq_pz_dst_vec, unsigned long, unsigned long);
void mpfq_pz_vec_getcoeff(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_add(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_neg(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_rev(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_sub(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
static inline
void mpfq_pz_vec_scal_mul(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_elt, unsigned long);
void mpfq_pz_vec_conv(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_random(mpfq_pz_dst_field, mpfq_pz_dst_vec, unsigned long, gmp_randstate_t);
void mpfq_pz_vec_random2(mpfq_pz_dst_field, mpfq_pz_dst_vec, unsigned long, gmp_randstate_t);
long mpfq_pz_vec_cmp(mpfq_pz_dst_field, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_is_zero(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
static inline
mpfq_pz_dst_vec mpfq_pz_vec_subvec(mpfq_pz_dst_field, mpfq_pz_dst_vec, long);
static inline
mpfq_pz_src_vec mpfq_pz_vec_subvec_const(mpfq_pz_dst_field, mpfq_pz_src_vec, long);
static inline
mpfq_pz_dst_elt mpfq_pz_vec_coeff_ptr(mpfq_pz_dst_field, mpfq_pz_dst_vec, long);
static inline
mpfq_pz_src_elt mpfq_pz_vec_coeff_ptr_const(mpfq_pz_dst_field, mpfq_pz_src_vec, long);
long mpfq_pz_vec_asprint(mpfq_pz_dst_field, char * *, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_fprint(mpfq_pz_dst_field, FILE *, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_print(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_sscan(mpfq_pz_dst_field, mpfq_pz_vec *, unsigned long *, const char *);
long mpfq_pz_vec_fscan(mpfq_pz_dst_field, FILE *, mpfq_pz_vec *, unsigned long *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, pz */
#define mpfq_pz_vec_scan(K, w, n)	mpfq_pz_vec_fscan(K,stdin,w,n)
long mpfq_pz_vec_hamming_weight(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_find_first_set(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_simd_hamming_weight(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
long mpfq_pz_vec_simd_find_first_set(mpfq_pz_dst_field, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_ur_init(mpfq_pz_dst_field, mpfq_pz_vec_ur *, unsigned long);
void mpfq_pz_vec_ur_set_zero(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, unsigned long);
void mpfq_pz_vec_ur_set_vec(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_ur_reinit(mpfq_pz_dst_field, mpfq_pz_vec_ur *, unsigned long, unsigned long);
void mpfq_pz_vec_ur_clear(mpfq_pz_dst_field, mpfq_pz_vec_ur *, unsigned long);
void mpfq_pz_vec_ur_set(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
void mpfq_pz_vec_ur_setcoeff(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_elt_ur, unsigned long);
void mpfq_pz_vec_ur_getcoeff(mpfq_pz_dst_field, mpfq_pz_dst_elt_ur, mpfq_pz_src_vec_ur, unsigned long);
void mpfq_pz_vec_ur_add(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
void mpfq_pz_vec_ur_sub(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
void mpfq_pz_vec_ur_neg(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
void mpfq_pz_vec_ur_rev(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static inline
void mpfq_pz_vec_scal_mul_ur(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, mpfq_pz_src_elt, unsigned long);
void mpfq_pz_vec_conv_ur(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, unsigned long, mpfq_pz_src_vec, unsigned long);
void mpfq_pz_vec_reduce(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_dst_vec_ur, unsigned long);
static inline
mpfq_pz_dst_vec_ur mpfq_pz_vec_ur_subvec(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, long);
static inline
mpfq_pz_src_vec_ur mpfq_pz_vec_ur_subvec_const(mpfq_pz_dst_field, mpfq_pz_src_vec_ur, long);
static inline
mpfq_pz_dst_elt mpfq_pz_vec_ur_coeff_ptr(mpfq_pz_dst_field, mpfq_pz_dst_vec_ur, long);
static inline
mpfq_pz_src_elt mpfq_pz_vec_ur_coeff_ptr_const(mpfq_pz_dst_field, mpfq_pz_src_vec_ur, long);
/* *pz::code_for_vec_elt_stride */
#define mpfq_pz_vec_elt_stride(k, n)	((n) * mpfq_pz_elt_stride((k)))
/* *pz::code_for_vec_ur_elt_stride */
#define mpfq_pz_vec_ur_elt_stride(k, n)	((n) * mpfq_pz_elt_ur_stride((k)))

/* Polynomial functions */
static inline
void mpfq_pz_poly_init(mpfq_pz_dst_field, mpfq_pz_poly, unsigned int);
static inline
void mpfq_pz_poly_clear(mpfq_pz_dst_field, mpfq_pz_poly);
static inline
void mpfq_pz_poly_set(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly);
void mpfq_pz_poly_setmonic(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_setcoeff(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_elt, unsigned int);
static inline
void mpfq_pz_poly_setcoeff_ui(mpfq_pz_dst_field, mpfq_pz_dst_poly, unsigned long, unsigned int);
static inline
void mpfq_pz_poly_getcoeff(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_poly, unsigned int);
static inline
int mpfq_pz_poly_deg(mpfq_pz_dst_field, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_add(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_sub(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_set_ui(mpfq_pz_dst_field, mpfq_pz_dst_poly, unsigned long);
static inline
void mpfq_pz_poly_add_ui(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, unsigned long);
static inline
void mpfq_pz_poly_sub_ui(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, unsigned long);
static inline
void mpfq_pz_poly_neg(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_scal_mul(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_elt);
static inline
void mpfq_pz_poly_mul(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
int mpfq_pz_poly_divmod(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
void mpfq_pz_poly_precomp_mod(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly);
void mpfq_pz_poly_mod_pre(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_gcd(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_xgcd(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
void mpfq_pz_poly_random(mpfq_pz_dst_field, mpfq_pz_dst_poly, unsigned int, gmp_randstate_t);
static inline
void mpfq_pz_poly_random2(mpfq_pz_dst_field, mpfq_pz_dst_poly, unsigned int, gmp_randstate_t);
static inline
int mpfq_pz_poly_cmp(mpfq_pz_dst_field, mpfq_pz_src_poly, mpfq_pz_src_poly);
static inline
int mpfq_pz_poly_asprint(mpfq_pz_dst_field, char * *, mpfq_pz_src_poly);
static inline
int mpfq_pz_poly_fprint(mpfq_pz_dst_field, FILE *, mpfq_pz_src_poly);
static inline
int mpfq_pz_poly_print(mpfq_pz_dst_field, mpfq_pz_src_poly);
static inline
int mpfq_pz_poly_sscan(mpfq_pz_dst_field, mpfq_pz_dst_poly, const char *);
static inline
int mpfq_pz_poly_fscan(mpfq_pz_dst_field, FILE *, mpfq_pz_dst_poly);
static inline
int mpfq_pz_poly_scan(mpfq_pz_dst_field, mpfq_pz_dst_poly);

/* Functions related to SIMD operation */
/* *simd_pz::code_for_simd_groupsize */
#define mpfq_pz_simd_groupsize(K)	1
static inline
int mpfq_pz_simd_hamming_weight(mpfq_pz_dst_field, mpfq_pz_src_elt);
static inline
int mpfq_pz_simd_find_first_set(mpfq_pz_dst_field, mpfq_pz_src_elt);
static inline
unsigned long mpfq_pz_simd_get_ui_at(mpfq_pz_dst_field, mpfq_pz_src_elt, int);
static inline
void mpfq_pz_simd_set_ui_at(mpfq_pz_dst_field, mpfq_pz_dst_elt, int, unsigned long);
static inline
void mpfq_pz_simd_add_ui_at(mpfq_pz_dst_field, mpfq_pz_dst_elt, mpfq_pz_src_elt, int, unsigned long);
static inline
void mpfq_pz_simd_set_ui_all(mpfq_pz_dst_field, mpfq_pz_dst_elt, unsigned long);
void mpfq_pz_add_dotprod(mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* Object-oriented interface */
static inline
void mpfq_pz_oo_field_clear(mpfq_vbase_ptr);
void mpfq_pz_oo_field_init(mpfq_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::gfp::field::code_for_field_characteristic, pz */
static inline
void mpfq_pz_field_characteristic(mpfq_pz_src_field k, mpz_ptr z)
{
        mpz_set(z, k->p);
}

/* *Mpfq::gfp::field::code_for_field_characteristic_srcptr, pz */
static inline
mpz_srcptr mpfq_pz_field_characteristic_srcptr(mpfq_pz_src_field k)
{
    return k->p;
}

/* *Mpfq::gfp::field::code_for_field_characteristic_bits, pz */
static inline
unsigned long mpfq_pz_field_characteristic_bits(mpfq_pz_src_field k)
{
        return mpz_sizeinbase(k->p, 2);
}

/* *Mpfq::gfp::field::code_for_field_init, pz */
static inline
void mpfq_pz_field_init(mpfq_pz_dst_field k)
{
    mpz_init(k->p);
    mpz_init(k->bigmul_p);
    k->io_base = 10;
    mpz_init(k->factor);
    k->ts_info.e=0;
}

/* *pz::code_for_elt_stride */
static inline
ptrdiff_t mpfq_pz_elt_stride(mpfq_pz_dst_field k MAYBE_UNUSED)
{
    return mpz_size(k->p) * sizeof(mp_limb_t);
}

/* *pz::code_for_set */
static inline
void mpfq_pz_set(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x)
{
        mpfq_copy(z, x, mpz_size(k->p));
}

/* *pz::code_for_set_ui */
static inline
void mpfq_pz_set_ui(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, unsigned long x0)
{
        ASSERT_FOR_STATIC_ANALYZER(mpz_size(k->p) > 1 || mpz_getlimbn(k->p, 0) > 0);
        z[0] = mpz_size(k->p) == 1 ? x0 % mpz_getlimbn(k->p, 0) : x0;
        if (mpz_size(k->p) >= 1) {
            /* Sigh. mpz_size(k->p) is a size_t. As such, it is bounded by
             * SIZE_MAX, and gcc 12+ sees it as such. It is apparently not
             * possible to use this bound in a memset directly.
             */
#if GNUC_VERSION_ATLEAST(7,1,0)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
            mpfq_zero(z + 1, (mpz_size(k->p) - 1));
#if GNUC_VERSION_ATLEAST(7,1,0)
#pragma GCC diagnostic pop
#endif
        }
}

/* *pz::code_for_set_zero */
static inline
void mpfq_pz_set_zero(mpfq_pz_dst_field k, mpfq_pz_dst_elt z)
{
        mpfq_zero(z, (mpz_size(k->p)));
}

/* *pz::code_for_get_ui */
static inline
unsigned long mpfq_pz_get_ui(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_src_elt x)
{
        return x[0];
}

/* *pz::code_for_get_mpn */
static inline
void mpfq_pz_get_mpn(mpfq_pz_dst_field k, mp_limb_t * z, mpfq_pz_src_elt x)
{
        mpfq_copy(z, x, mpz_size(k->p));
}

/* *pz::code_for_get_mpz */
static inline
void mpfq_pz_get_mpz(mpfq_pz_dst_field k, mpz_ptr z, mpfq_pz_src_elt x)
{
    int n = mpz_size(k->p);
    mpz_set_ui(z, x[n - 1]);
    for (int i = n - 2; i >= 0; --i) {
        mpz_mul_2exp(z, z, 64);
        mpz_add_ui(z, z, x[i]);
    }
}

/* *pz::code_for_elt_ur_stride */
static inline
ptrdiff_t mpfq_pz_elt_ur_stride(mpfq_pz_dst_field k MAYBE_UNUSED)
{
    return mpz_size(k->bigmul_p) * sizeof(mp_limb_t);
}

/* *pz::code_for_elt_ur_set */
static inline
void mpfq_pz_elt_ur_set(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt_ur x)
{
        mpfq_copy(z, x, mpz_size(k->bigmul_p));
}

/* *pz::code_for_elt_ur_set_elt */
static inline
void mpfq_pz_elt_ur_set_elt(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt x)
{
        mpfq_copy(z, x, mpz_size(k->p));
        mpfq_zero(z + mpz_size(k->p), (mpz_size(k->bigmul_p) - mpz_size(k->p)));
}

/* *pz::code_for_elt_ur_set_zero */
static inline
void mpfq_pz_elt_ur_set_zero(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z)
{
        mpfq_zero(z, (mpz_size(k->bigmul_p)));
}

/* *pz::code_for_elt_ur_set_ui */
static inline
void mpfq_pz_elt_ur_set_ui(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, unsigned long x0)
{
        z[0] = x0;
        mpfq_zero(z + 1, (mpz_size(k->bigmul_p) - 1));
}

/* *simd_pz::code_for_addmul_si_ur */
static inline
void mpfq_pz_addmul_si_ur(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_elt_ur w, mpfq_pz_src_elt u, long v)
{
        mpfq_pz_elt_ur s;
        mpfq_pz_elt vx;
        mpfq_pz_elt_ur_init(K, &s);
        mpfq_pz_init(K, &vx);
        if (v>0) {
            mpfq_pz_set_ui(K, vx, v);
            mpfq_pz_mul_ur(K, s, u, vx);
            mpfq_pz_elt_ur_add(K, w, w, s);
        } else {
            mpfq_pz_set_ui(K, vx, -v);
            mpfq_pz_mul_ur(K, s, u, vx);
            mpfq_pz_elt_ur_sub(K, w, w, s);
        }
        mpfq_pz_clear(K, &vx);
        mpfq_pz_elt_ur_clear(K, &s);
}

/* *pz::code_for_cmp */
static inline
int mpfq_pz_cmp(mpfq_pz_dst_field k, mpfq_pz_src_elt x, mpfq_pz_src_elt y)
{
        return mpn_cmp(x, y, mpz_size(k->p));
}

/* *pz::code_for_cmp_ui */
static inline
int mpfq_pz_cmp_ui(mpfq_pz_dst_field k, mpfq_pz_src_elt x, unsigned long y0)
{
        for (int i = mpz_size(k->p) - 1; i > 0; --i) {
        if (x[i]) return 1;
        }
        return (x[0] > y0) - (x[0] < y0);
}

/* *pz::code_for_is_zero */
static inline
int mpfq_pz_is_zero(mpfq_pz_dst_field k, mpfq_pz_src_elt x)
{
        for (unsigned int i = 0; i < mpz_size(k->p); ++i)
        if (x[i] != 0)
            return 0;
        return 1;
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, pz */
static inline
void mpfq_pz_vec_scal_mul(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec w, mpfq_pz_src_vec u, mpfq_pz_src_elt c, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(K, u, i);
        mpfq_pz_dst_elt y = mpfq_pz_vec_coeff_ptr(K, w, i);
        mpfq_pz_mul(K, y, x, c);
    }
}

/* *pz::code_for_vec_subvec */
static inline
mpfq_pz_dst_vec mpfq_pz_vec_subvec(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec v, long i)
{
    return v + i * mpz_size(K->p);
}

/* *pz::code_for_vec_subvec_const */
static inline
mpfq_pz_src_vec mpfq_pz_vec_subvec_const(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_vec v, long i)
{
    return v + i * mpz_size(K->p);
}

/* *pz::code_for_vec_coeff_ptr */
static inline
mpfq_pz_dst_elt mpfq_pz_vec_coeff_ptr(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec v, long i)
{
    return v + i * mpz_size(K->p);
}

/* *pz::code_for_vec_coeff_ptr_const */
static inline
mpfq_pz_src_elt mpfq_pz_vec_coeff_ptr_const(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_vec v, long i)
{
    return v + i * mpz_size(K->p);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, pz */
static inline
void mpfq_pz_vec_scal_mul_ur(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec u, mpfq_pz_src_elt c, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i++) {
        mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(K, u, i);
        mpfq_pz_dst_elt_ur y = mpfq_pz_vec_ur_coeff_ptr(K, w, i);
        mpfq_pz_mul_ur(K, y, x, c);
    }
}

/* *pz::code_for_vec_ur_subvec */
static inline
mpfq_pz_dst_vec_ur mpfq_pz_vec_ur_subvec(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec_ur v, long i)
{
    return v + i * mpz_size(K->bigmul_p);
}

/* *pz::code_for_vec_ur_subvec_const */
static inline
mpfq_pz_src_vec_ur mpfq_pz_vec_ur_subvec_const(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_vec_ur v, long i)
{
    return v + i * mpz_size(K->bigmul_p);
}

/* *pz::code_for_vec_ur_coeff_ptr */
static inline
mpfq_pz_dst_elt mpfq_pz_vec_ur_coeff_ptr(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec_ur v, long i)
{
    return v + i * mpz_size(K->bigmul_p);
}

/* *pz::code_for_vec_ur_coeff_ptr_const */
static inline
mpfq_pz_src_elt mpfq_pz_vec_ur_coeff_ptr_const(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_vec_ur v, long i)
{
    return v + i * mpz_size(K->bigmul_p);
}

/* *Mpfq::defaults::poly::code_for_poly_init, pz */
static inline
void mpfq_pz_poly_init(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_poly p, unsigned int n)
{
    mpfq_pz_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear, pz */
static inline
void mpfq_pz_poly_clear(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_poly p)
{
    mpfq_pz_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set, pz */
static inline
void mpfq_pz_poly_set(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_pz_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff, pz */
static inline
void mpfq_pz_poly_setcoeff(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_elt x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_pz_vec_set_zero(k, mpfq_pz_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_pz_vec_setcoeff(k, w->c, x, i);
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_setcoeff_ui, pz */
static inline
void mpfq_pz_poly_setcoeff_ui(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, unsigned long x, unsigned int i)
{
    if (w->alloc < (i+1)) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    if (w->size < (i+1)) {
        mpfq_pz_vec_set_zero(k, mpfq_pz_vec_subvec(k, w->c, w->size), (i - w->size));
        w->size = i+1;
    }
    mpfq_pz_vec_setcoeff_ui(k, w->c, x, i);
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_getcoeff, pz */
static inline
void mpfq_pz_poly_getcoeff(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_elt x, mpfq_pz_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       mpfq_pz_set_ui(k,x,0);
    } else {
       mpfq_pz_vec_getcoeff(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg, pz */
static inline
int mpfq_pz_poly_deg(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    mpfq_pz_elt temp;
    mpfq_pz_init(K, &temp);
    mpfq_pz_vec_getcoeff(K, temp, w->c, deg);
    int comp=mpfq_pz_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           mpfq_pz_vec_getcoeff(K, temp, w->c, deg);
           comp=mpfq_pz_cmp_ui(K, temp, 0);
        }
    }
    mpfq_pz_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add, pz */
static inline
void mpfq_pz_poly_add(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, mpfq_pz_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_pz_vec_add(k, w->c, u->c, v->c, su);
        mpfq_pz_vec_set(k, mpfq_pz_vec_subvec(k, w->c, su), mpfq_pz_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_pz_vec_add(k, w->c, u->c, v->c, sv);
        mpfq_pz_vec_set(k, mpfq_pz_vec_subvec(k, w->c, sv), mpfq_pz_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_sub, pz */
static inline
void mpfq_pz_poly_sub(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, mpfq_pz_src_poly v)
{
    unsigned int su = u->size;
    unsigned int sv = v->size;
    unsigned int maxsize = MAX(su, sv);
    if (w->alloc < maxsize) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    w->size = maxsize;
    if (!maxsize) return;
    if (su <= sv) {
        mpfq_pz_vec_sub(k, w->c, u->c, v->c, su);
        mpfq_pz_vec_neg(k, mpfq_pz_vec_subvec(k, w->c, su), mpfq_pz_vec_subvec_const(k, v->c, su), sv-su);
    } else {
        mpfq_pz_vec_sub(k, w->c, u->c, v->c, sv);
        mpfq_pz_vec_set(k, mpfq_pz_vec_subvec(k, w->c, sv), mpfq_pz_vec_subvec_const(k, u->c, sv), su-sv);
    }
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_set_ui, pz */
static inline
void mpfq_pz_poly_set_ui(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, unsigned long x)
{
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_pz_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_pz_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui, pz */
static inline
void mpfq_pz_poly_add_ui(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_pz_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_pz_vec_setcoeff_ui(k, w->c, x, 0);
        w->size = 1;
        w->size = 1 + mpfq_pz_poly_deg(k, w);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_pz_vec_set(k, mpfq_pz_vec_subvec(k, w->c, 1), mpfq_pz_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_pz_add_ui(k, mpfq_pz_vec_coeff_ptr(k, w->c, 0), mpfq_pz_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui, pz */
static inline
void mpfq_pz_poly_sub_ui(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            mpfq_pz_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        mpfq_pz_elt temp;
        mpfq_pz_init(k, &temp);
        mpfq_pz_set_ui(k, temp, x);
        mpfq_pz_neg(k, mpfq_pz_vec_coeff_ptr(k, w->c, 0), temp);
        w->size = mpfq_pz_cmp_ui(k, temp, 0);
        mpfq_pz_clear(k, &temp);
        return;
    }
    if (w->alloc < u->size) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    w->size=u->size;
    mpfq_pz_vec_set(k, mpfq_pz_vec_subvec(k, w->c, 1), mpfq_pz_vec_subvec_const(k, u->c, 1), u->size - 1);
    mpfq_pz_sub_ui(k, mpfq_pz_vec_coeff_ptr(k, w->c, 0), mpfq_pz_vec_coeff_ptr_const(k, u->c, 0), x);
}

/* *Mpfq::defaults::poly::code_for_poly_neg, pz */
static inline
void mpfq_pz_poly_neg(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u)
{
    if (w->alloc < u->size) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    mpfq_pz_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul, pz */
static inline
void mpfq_pz_poly_scal_mul(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, mpfq_pz_src_elt x)
{
    if (mpfq_pz_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_pz_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::poly::code_for_poly_mul, pz */
static inline
void mpfq_pz_poly_mul(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, mpfq_pz_src_poly u, mpfq_pz_src_poly v)
{
    unsigned int usize = mpfq_pz_poly_deg(k, u)+1;
    unsigned int vsize = mpfq_pz_poly_deg(k, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    mpfq_pz_vec_conv(k, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
    w->size = 1 + mpfq_pz_poly_deg(k, w);
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly, pz */
static inline
void mpfq_pz_poly_gcd(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly g, mpfq_pz_src_poly a0, mpfq_pz_src_poly b0)
{
    mpfq_pz_poly a,b,q,r;
    int da0=mpfq_pz_poly_deg(k,a0), db0=mpfq_pz_poly_deg(k,b0);
    if (db0==-1)
     mpfq_pz_poly_set(k,g,a0);
    else {
     mpfq_pz_poly_init(k,a,da0+1);
     mpfq_pz_poly_init(k,b,db0+1);
     mpfq_pz_poly_init(k,q,1);
     mpfq_pz_poly_init(k,r,db0);
     mpfq_pz_poly_set(k,a,a0);
     mpfq_pz_poly_set(k,b,b0);
     while (mpfq_pz_poly_deg(k,b)>=0) {
      mpfq_pz_poly_divmod(k,q,r,a,b);
      mpfq_pz_poly_set(k,a,b);
      mpfq_pz_poly_set(k,b,r); 
     }
     mpfq_pz_poly_setmonic(k,g,a);
    mpfq_pz_poly_clear(k,a);
    mpfq_pz_poly_clear(k,b);
    mpfq_pz_poly_clear(k,q);
    mpfq_pz_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly, pz */
static inline
void mpfq_pz_poly_xgcd(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly g, mpfq_pz_dst_poly u0, mpfq_pz_dst_poly v0, mpfq_pz_src_poly a0, mpfq_pz_src_poly b0)
{
    mpfq_pz_poly a,b,u,v,w,x,q,r;
    mpfq_pz_elt c;
    mpfq_pz_init(k,&c);
    mpfq_pz_set_ui(k,c,0);        /* placate gcc */
    int da0=mpfq_pz_poly_deg(k,a0), db0=mpfq_pz_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      mpfq_pz_poly_set(k,u0,a0);
      mpfq_pz_poly_set(k,v0,b0);
      mpfq_pz_poly_set(k,g,a0);
     } else {
      mpfq_pz_poly_getcoeff(k,c,a0,da0);
      mpfq_pz_inv(k,c,c);
      mpfq_pz_poly_scal_mul(k,g,a0,c);
      mpfq_pz_poly_set(k,v0,b0);
      mpfq_pz_poly_set(k,u0,b0);
      mpfq_pz_poly_setcoeff(k,u0,c,0);
     }
    }
    else {
     mpfq_pz_poly_init(k,a,da0+1);
     mpfq_pz_poly_init(k,b,db0+1);
     mpfq_pz_poly_init(k,q,1);
     mpfq_pz_poly_init(k,r,db0);
     mpfq_pz_poly_set(k,a,a0);
     mpfq_pz_poly_set(k,b,b0);
     mpfq_pz_poly_init(k,u,1);
     mpfq_pz_poly_init(k,v,1);
     mpfq_pz_poly_init(k,w,1);
     mpfq_pz_poly_init(k,x,1);
     mpfq_pz_poly_setcoeff_ui(k,u,1,0);
     mpfq_pz_poly_setcoeff_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (mpfq_pz_poly_deg(k,b)>=0) {
      mpfq_pz_poly_divmod(k,q,r,a,b);
      mpfq_pz_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      mpfq_pz_poly_set(k,b,r);
      mpfq_pz_poly_mul(k,r,q,w);
      mpfq_pz_poly_sub(k,r,u,r);
      mpfq_pz_poly_set(k,u,w);   /* u,w <- w,u-qw */
      mpfq_pz_poly_set(k,w,r);
      mpfq_pz_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      mpfq_pz_poly_sub(k,r,v,r);
      mpfq_pz_poly_set(k,v,x);
      mpfq_pz_poly_set(k,x,r);
     }
     dega=mpfq_pz_poly_deg(k,a);
     ASSERT_FOR_STATIC_ANALYZER(dega >= 0);
     mpfq_pz_poly_getcoeff(k,c,a,dega);
     mpfq_pz_inv(k,c,c);
     mpfq_pz_poly_scal_mul(k,g,a,c);
     mpfq_pz_poly_scal_mul(k,u0,u,c);
     mpfq_pz_poly_scal_mul(k,v0,v,c);
     mpfq_pz_poly_clear(k,a);
     mpfq_pz_poly_clear(k,b);
     mpfq_pz_poly_clear(k,u);
     mpfq_pz_poly_clear(k,v);
     mpfq_pz_poly_clear(k,w);
     mpfq_pz_poly_clear(k,x);
     mpfq_pz_poly_clear(k,q);
     mpfq_pz_poly_clear(k,r);
    }
    mpfq_pz_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_random, pz */
static inline
void mpfq_pz_poly_random(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_pz_vec_random(k, w->c, n,state);
    w->size=n;
    int wdeg = mpfq_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_random2, pz */
static inline
void mpfq_pz_poly_random2(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        mpfq_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    mpfq_pz_vec_random2(k, w->c, n,state);
    w->size=n;
    int wdeg = mpfq_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_cmp, pz */
static inline
int mpfq_pz_poly_cmp(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_src_poly u, mpfq_pz_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return mpfq_pz_vec_cmp(k, u->c, v->c, u->size);
}

/* *Mpfq::defaults::poly::code_for_poly_asprint, pz */
static inline
int mpfq_pz_poly_asprint(mpfq_pz_dst_field k MAYBE_UNUSED, char * * pstr, mpfq_pz_src_poly w)
{
    return mpfq_pz_vec_asprint(k, pstr, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_fprint, pz */
static inline
int mpfq_pz_poly_fprint(mpfq_pz_dst_field k MAYBE_UNUSED, FILE * file, mpfq_pz_src_poly w)
{
    return mpfq_pz_vec_fprint(k, file, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_print, pz */
static inline
int mpfq_pz_poly_print(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_src_poly w)
{
    return mpfq_pz_vec_print(k, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_sscan, pz */
static inline
int mpfq_pz_poly_sscan(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w, const char * str)
{
    int ret;
    ret = mpfq_pz_vec_sscan(k, &(w->c), &(w->alloc), str);
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_fscan, pz */
static inline
int mpfq_pz_poly_fscan(mpfq_pz_dst_field k MAYBE_UNUSED, FILE * file, mpfq_pz_dst_poly w)
{
    int ret;
    ret = mpfq_pz_vec_fscan(k, file, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_scan, pz */
static inline
int mpfq_pz_poly_scan(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_dst_poly w)
{
    int ret;
    ret = mpfq_pz_vec_scan(k, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *simd_pz::code_for_simd_hamming_weight */
static inline
int mpfq_pz_simd_hamming_weight(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_elt p)
{
    return !mpfq_pz_is_zero(K,p);
}

/* *simd_pz::code_for_simd_find_first_set */
static inline
int mpfq_pz_simd_find_first_set(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_elt p)
{
    return mpfq_pz_simd_hamming_weight(K,p) - 1;
}

/* *simd_pz::code_for_simd_get_ui_at */
static inline
unsigned long mpfq_pz_simd_get_ui_at(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_elt p, int k MAYBE_UNUSED)
{
    return mpfq_pz_get_ui(K,p);
}

/* *simd_pz::code_for_simd_set_ui_at */
static inline
void mpfq_pz_simd_set_ui_at(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    mpfq_pz_set_ui(K,p,v);
}

/* *simd_pz::code_for_simd_add_ui_at */
static inline
void mpfq_pz_simd_add_ui_at(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_elt p, mpfq_pz_src_elt p0, int k MAYBE_UNUSED, unsigned long v)
{
    mpfq_pz_add_ui(K,p,p0,v);
}

/* *simd_pz::code_for_simd_set_ui_all */
static inline
void mpfq_pz_simd_set_ui_all(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_elt p, unsigned long v)
{
    mpfq_pz_set_ui(K,p,v);
}

/* Mpfq::engine::oo::oo_field_clear */
/* Triggered by: oo */
static inline
void mpfq_pz_oo_field_clear(mpfq_vbase_ptr f)
{
    mpfq_pz_field_clear((mpfq_pz_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* MPFQ_PZ_H_ */

/* vim:set ft=cpp: */
