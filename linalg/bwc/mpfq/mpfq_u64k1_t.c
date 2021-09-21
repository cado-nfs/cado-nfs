#include "cado.h"
/* MPFQ generated file -- do not edit */

#include "mpfq_u64k1_t.h"

#include "binary-dotprods-backends.h"
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
   k=1,
   tag=u64k1,
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
     [ (?^:mpfq_u64k1_elt \*), void *, ],
     [ (?^:mpfq_u64k1_src_elt\b), const void *, ],
     [ (?^:mpfq_u64k1_elt\b), void *, ],
     [ (?^:mpfq_u64k1_dst_elt\b), void *, ],
     [ (?^:mpfq_u64k1_elt_ur \*), void *, ],
     [ (?^:mpfq_u64k1_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_u64k1_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k1_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_u64k1_vec \*), void *, ],
     [ (?^:mpfq_u64k1_src_vec\b), const void *, ],
     [ (?^:mpfq_u64k1_vec\b), void *, ],
     [ (?^:mpfq_u64k1_dst_vec\b), void *, ],
     [ (?^:mpfq_u64k1_vec_ur \*), void *, ],
     [ (?^:mpfq_u64k1_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_u64k1_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k1_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_u64k1_poly \*), void *, ],
     [ (?^:mpfq_u64k1_src_poly\b), const void *, ],
     [ (?^:mpfq_u64k1_poly\b), void *, ],
     [ (?^:mpfq_u64k1_dst_poly\b), void *, ],
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
#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_m128_wrapper_add_dotprod(mpfq_vbase_ptr K0 MAYBE_UNUSED, mpfq_vbase_ptr K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_m128_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    mpfq_u64k1_m128_add_dotprod(K0->obj, K1->obj, xw, xu1, xu0, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* *simd_dotprod::code_for_member_template_add_dotprod */
void mpfq_u64k1_m128_add_dotprod(mpfq_u64k1_dst_field K0 MAYBE_UNUSED, mpfq_m128_dst_field K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_m128_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = (uint64_t *) xw;
    const uint64_t * u0 = (const uint64_t *) xu0;
    const uint64_t * u1 = (const uint64_t *) xu1;
    add_dotprod_64K_128(w,u0,u1,n,1);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k1_wrapper_add_dotprod(mpfq_vbase_ptr K0 MAYBE_UNUSED, mpfq_vbase_ptr K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k1_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    mpfq_u64k1_u64k1_add_dotprod(K0->obj, K1->obj, xw, xu1, xu0, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* *simd_dotprod::code_for_member_template_add_dotprod */
void mpfq_u64k1_u64k1_add_dotprod(mpfq_u64k1_dst_field K0 MAYBE_UNUSED, mpfq_u64k1_dst_field K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k1_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = (uint64_t *) xw;
    const uint64_t * u0 = (const uint64_t *) xu0;
    const uint64_t * u1 = (const uint64_t *) xu1;
    add_dotprod_64K_64(w,u1,u0,n,1);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k2_wrapper_add_dotprod(mpfq_vbase_ptr K0 MAYBE_UNUSED, mpfq_vbase_ptr K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k2_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    mpfq_u64k1_u64k2_add_dotprod(K0->obj, K1->obj, xw, xu1, xu0, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* *simd_dotprod::code_for_member_template_add_dotprod */
void mpfq_u64k1_u64k2_add_dotprod(mpfq_u64k1_dst_field K0 MAYBE_UNUSED, mpfq_u64k2_dst_field K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k2_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = (uint64_t *) xw;
    const uint64_t * u0 = (const uint64_t *) xu0;
    const uint64_t * u1 = (const uint64_t *) xu1;
    add_dotprod_64K_128(w,u0,u1,n,1);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k3_wrapper_add_dotprod(mpfq_vbase_ptr K0 MAYBE_UNUSED, mpfq_vbase_ptr K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k3_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    mpfq_u64k1_u64k3_add_dotprod(K0->obj, K1->obj, xw, xu1, xu0, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* *simd_dotprod::code_for_member_template_add_dotprod */
void mpfq_u64k1_u64k3_add_dotprod(mpfq_u64k1_dst_field K0 MAYBE_UNUSED, mpfq_u64k3_dst_field K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k3_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = (uint64_t *) xw;
    const uint64_t * u0 = (const uint64_t *) xu0;
    const uint64_t * u1 = (const uint64_t *) xu1;
    add_dotprod_64K_64(w,u1,u0,n,3);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k4_wrapper_add_dotprod(mpfq_vbase_ptr K0 MAYBE_UNUSED, mpfq_vbase_ptr K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k4_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    mpfq_u64k1_u64k4_add_dotprod(K0->obj, K1->obj, xw, xu1, xu0, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* *simd_dotprod::code_for_member_template_add_dotprod */
void mpfq_u64k1_u64k4_add_dotprod(mpfq_u64k1_dst_field K0 MAYBE_UNUSED, mpfq_u64k4_dst_field K1 MAYBE_UNUSED, mpfq_u64k1_dst_vec xw, mpfq_u64k4_src_vec xu1, mpfq_u64k1_src_vec xu0, unsigned int n)
{
    uint64_t * w = (uint64_t *) xw;
    const uint64_t * u0 = (const uint64_t *) xu0;
    const uint64_t * u1 = (const uint64_t *) xu1;
    add_dotprod_64K_64(w,u1,u0,n,4);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_m128_wrapper_addmul_tiny(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_u64k1_src_vec u, mpfq_m128_src_vec v, unsigned int n)
{
    mpfq_u64k1_m128_addmul_tiny(K->obj, L->obj, w, u, v, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* *simd_dotprod::code_for_member_template_addmul_tiny */
void mpfq_u64k1_m128_addmul_tiny(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_m128_dst_field L MAYBE_UNUSED, mpfq_m128_dst_vec w, mpfq_u64k1_src_vec u, mpfq_m128_src_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],(const uint64_t*)v[0],n,1,2);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k1_wrapper_addmul_tiny(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k1_src_vec v, unsigned int n)
{
    mpfq_u64k1_u64k1_addmul_tiny(K->obj, L->obj, w, u, v, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* *simd_dotprod::code_for_member_template_addmul_tiny */
void mpfq_u64k1_u64k1_addmul_tiny(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k1_src_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],(const uint64_t*)v[0],n,1,1);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k2_wrapper_addmul_tiny(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k2_src_vec v, unsigned int n)
{
    mpfq_u64k1_u64k2_addmul_tiny(K->obj, L->obj, w, u, v, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* *simd_dotprod::code_for_member_template_addmul_tiny */
void mpfq_u64k1_u64k2_addmul_tiny(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_field L MAYBE_UNUSED, mpfq_u64k2_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k2_src_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],(const uint64_t*)v[0],n,1,2);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k3_wrapper_addmul_tiny(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k3_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k3_src_vec v, unsigned int n)
{
    mpfq_u64k1_u64k3_addmul_tiny(K->obj, L->obj, w, u, v, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* *simd_dotprod::code_for_member_template_addmul_tiny */
void mpfq_u64k1_u64k3_addmul_tiny(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k3_dst_field L MAYBE_UNUSED, mpfq_u64k3_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k3_src_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],(const uint64_t*)v[0],n,1,3);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k4_wrapper_addmul_tiny(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k4_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k4_src_vec v, unsigned int n)
{
    mpfq_u64k1_u64k4_addmul_tiny(K->obj, L->obj, w, u, v, n);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* *simd_dotprod::code_for_member_template_addmul_tiny */
void mpfq_u64k1_u64k4_addmul_tiny(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k4_dst_field L MAYBE_UNUSED, mpfq_u64k4_dst_vec w, mpfq_u64k1_src_vec u, mpfq_u64k4_src_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],(const uint64_t*)v[0],n,1,4);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_m128_wrapper_transpose(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_m128_src_vec u)
{
    mpfq_u64k1_m128_transpose(K->obj, L->obj, w, u);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
/* *simd_dotprod::code_for_member_template_transpose */
void mpfq_u64k1_m128_transpose(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_m128_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_m128_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],1,2);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k1_wrapper_transpose(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k1_src_vec u)
{
    mpfq_u64k1_u64k1_transpose(K->obj, L->obj, w, u);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
/* *simd_dotprod::code_for_member_template_transpose */
void mpfq_u64k1_u64k1_transpose(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k1_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k1_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],1,1);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k2_wrapper_transpose(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k2_src_vec u)
{
    mpfq_u64k1_u64k2_transpose(K->obj, L->obj, w, u);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
/* *simd_dotprod::code_for_member_template_transpose */
void mpfq_u64k1_u64k2_transpose(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k2_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k2_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],1,2);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k3_wrapper_transpose(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k3_src_vec u)
{
    mpfq_u64k1_u64k3_transpose(K->obj, L->obj, w, u);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
/* *simd_dotprod::code_for_member_template_transpose */
void mpfq_u64k1_u64k3_transpose(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k3_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k3_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],1,3);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* Mpfq::engine::handler::create_code */
void mpfq_u64k1_u64k4_wrapper_transpose(mpfq_vbase_ptr K MAYBE_UNUSED, mpfq_vbase_ptr L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k4_src_vec u)
{
    mpfq_u64k1_u64k4_transpose(K->obj, L->obj, w, u);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */

#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
/* *simd_dotprod::code_for_member_template_transpose */
void mpfq_u64k1_u64k4_transpose(mpfq_u64k1_dst_field K MAYBE_UNUSED, mpfq_u64k4_dst_field L MAYBE_UNUSED, mpfq_u64k1_dst_vec w, mpfq_u64k4_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],1,4);
}
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */


/* vim:set ft=cpp: */
