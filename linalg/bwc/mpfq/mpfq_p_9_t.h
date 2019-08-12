#ifndef MPFQ_P_9_T_H_
#define MPFQ_P_9_T_H_

/* MPFQ generated file -- do not edit */

#include "mpfq_p_9.h"
/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_9, tag=p_9, }, ],
   fieldtype=prime,
   n=9,
   nn=19,
   opthw=,
   tag=p_9,
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
     [ (?^:mpfq_p_9_elt \*), void *, ],
     [ (?^:mpfq_p_9_src_elt\b), const void *, ],
     [ (?^:mpfq_p_9_elt\b), void *, ],
     [ (?^:mpfq_p_9_dst_elt\b), void *, ],
     [ (?^:mpfq_p_9_elt_ur \*), void *, ],
     [ (?^:mpfq_p_9_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_p_9_elt_ur\b), void *, ],
     [ (?^:mpfq_p_9_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_p_9_vec \*), void *, ],
     [ (?^:mpfq_p_9_src_vec\b), const void *, ],
     [ (?^:mpfq_p_9_vec\b), void *, ],
     [ (?^:mpfq_p_9_dst_vec\b), void *, ],
     [ (?^:mpfq_p_9_vec_ur \*), void *, ],
     [ (?^:mpfq_p_9_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_p_9_vec_ur\b), void *, ],
     [ (?^:mpfq_p_9_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_p_9_poly \*), void *, ],
     [ (?^:mpfq_p_9_src_poly\b), const void *, ],
     [ (?^:mpfq_p_9_poly\b), void *, ],
     [ (?^:mpfq_p_9_dst_poly\b), void *, ],
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
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_wrapper_add_dotprod(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_p_9_dst_vec, mpfq_p_9_src_vec, mpfq_p_9_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */

#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_add_dotprod(mpfq_p_9_dst_field, mpfq_p_9_dst_field, mpfq_p_9_dst_vec, mpfq_p_9_src_vec, mpfq_p_9_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */

#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_wrapper_addmul_tiny(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_p_9_dst_vec, mpfq_p_9_src_vec, mpfq_p_9_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */

#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_addmul_tiny(mpfq_p_9_dst_field, mpfq_p_9_dst_field, mpfq_p_9_dst_vec, mpfq_p_9_src_vec, mpfq_p_9_src_vec, unsigned int);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */

#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_wrapper_transpose(mpfq_vbase_ptr, mpfq_vbase_ptr, mpfq_p_9_dst_vec, mpfq_p_9_src_vec);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */

#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
void mpfq_p_9_p_9_transpose(mpfq_p_9_dst_field, mpfq_p_9_dst_field, mpfq_p_9_dst_vec, mpfq_p_9_src_vec);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */


#endif  /* MPFQ_P_9_T_H_ */

/* vim:set ft=cpp: */
