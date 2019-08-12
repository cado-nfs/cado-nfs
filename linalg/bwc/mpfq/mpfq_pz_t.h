#ifndef MPFQ_PZ_T_H_
#define MPFQ_PZ_T_H_

/* MPFQ generated file -- do not edit */

#include "mpfq_pz.h"
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


/* Functions operating on the field structure */

/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */

/* Vector functions */

/* Polynomial functions */

/* Functions related to SIMD operation */

/* Member templates related to SIMD operation */

/* Object-oriented interface */
#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_wrapper_add_dotprod(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */

#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_add_dotprod(mpfq_pz_dst_field, mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */

#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_wrapper_addmul_tiny(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */

#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_addmul_tiny(mpfq_pz_dst_field, mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */

#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_wrapper_transpose(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */

#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
void mpfq_pz_pz_transpose(mpfq_pz_dst_field, mpfq_pz_dst_field, mpfq_pz_dst_vec, mpfq_pz_src_vec);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */


#endif  /* MPFQ_PZ_T_H_ */

/* vim:set ft=cpp: */
