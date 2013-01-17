/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_u64k2_t.h"

#include "binary-dotprods-backends.h"
/* Active handler: simd_u64k */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Options used: w=64 k=2 tag=u64k2 vbase_stuff={
                 'vc:includes' => [
                                    '<stdarg.h>'
                                  ],
                 'member_templates_restrict' => {
                                                  'p_1' => [
                                                             {
                                                               'cpp_ifdef' => 'COMPILE_MPFQ_PRIME_FIELDS',
                                                               'tag' => 'p_1'
                                                             }
                                                           ],
                                                  'p_4' => [
                                                             {
                                                               'cpp_ifdef' => 'COMPILE_MPFQ_PRIME_FIELDS',
                                                               'tag' => 'p_4'
                                                             }
                                                           ],
                                                  'u64k2' => [
                                                               'u64k1',
                                                               'u64k2'
                                                             ],
                                                  'p_3' => [
                                                             {
                                                               'cpp_ifdef' => 'COMPILE_MPFQ_PRIME_FIELDS',
                                                               'tag' => 'p_3'
                                                             }
                                                           ],
                                                  'u64k1' => $vbase_stuff->{'member_templates_restrict'}{'u64k2'}
                                                },
                 'families' => [
                                 $vbase_stuff->{'member_templates_restrict'}{'p_4'},
                                 $vbase_stuff->{'member_templates_restrict'}{'p_1'},
                                 $vbase_stuff->{'member_templates_restrict'}{'p_3'},
                                 $vbase_stuff->{'member_templates_restrict'}{'u64k2'}
                               ],
                 'choose_byfeatures' => sub { "DUMMY" }
               };
 family=[u64k1, u64k2] virtual_base={
                  'filebase' => 'abase_vbase',
                  'substitutions' => [
                                       [
                                         qr/(?^:abase_u64k2_elt \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_src_elt\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_dst_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_elt_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_src_elt_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_dst_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_vec \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_src_vec\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_dst_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_vec_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_src_vec_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_dst_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_poly \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_src_poly\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_poly\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_u64k2_dst_poly\b)/,
                                         'void *'
                                       ]
                                     ],
                  'name' => 'abase_vbase',
                  'global_prefix' => 'abase_'
                };
 */


/* Functions operating on the field structure */

/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */

/* Vector functions */

/* Functions related to SIMD operation */

/* Member templates related to SIMD operation */

/* MPI interface */

/* Object-oriented interface */
/* *simd_dotprod::code_for_member_template_dotprod */
void abase_u64k2_u64k1_dotprod(abase_u64k2_dst_field K0 MAYBE_UNUSED, abase_u64k1_dst_field K1 MAYBE_UNUSED, abase_u64k2_dst_vec xw, abase_u64k1_src_vec xu1, abase_u64k2_src_vec xu0, unsigned int n)
{
    uint64_t * w = xw[0];
    const uint64_t * u0 = xu0[0];
    const uint64_t * u1 = xu1[0];
    dotprod_64K_128(w,u1,u0,n,1);
}

/* *simd_dotprod::code_for_member_template_dotprod */
void abase_u64k2_u64k2_dotprod(abase_u64k2_dst_field K0 MAYBE_UNUSED, abase_u64k2_dst_field K1 MAYBE_UNUSED, abase_u64k2_dst_vec xw, abase_u64k2_src_vec xu1, abase_u64k2_src_vec xu0, unsigned int n)
{
    uint64_t * w = xw[0];
    const uint64_t * u0 = xu0[0];
    const uint64_t * u1 = xu1[0];
    dotprod_64K_128(w,u0,u1,n,2);
}

/* *simd_dotprod::code_for_member_template_addmul_tiny */
void abase_u64k2_u64k1_addmul_tiny(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k1_dst_field L MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k2_src_vec u, abase_u64k1_dst_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const
        uint64_t*)u[0],(const uint64_t*)v[0],n,2,1);
}

/* *simd_dotprod::code_for_member_template_addmul_tiny */
void abase_u64k2_u64k2_addmul_tiny(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_dst_field L MAYBE_UNUSED, abase_u64k2_dst_vec w, abase_u64k2_src_vec u, abase_u64k2_dst_vec v, unsigned int n)
{
    vaddmul_tiny_64K_64L((uint64_t*)w[0],(const
        uint64_t*)u[0],(const uint64_t*)v[0],n,2,2);
}

/* *simd_dotprod::code_for_member_template_transpose */
void abase_u64k2_u64k1_transpose(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k1_dst_field L MAYBE_UNUSED, abase_u64k2_dst_vec w, abase_u64k1_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],2,1);
}

/* *simd_dotprod::code_for_member_template_transpose */
void abase_u64k2_u64k2_transpose(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_dst_field L MAYBE_UNUSED, abase_u64k2_dst_vec w, abase_u64k2_src_vec u)
{
    vtranspose_64K_64L((uint64_t*)w[0],(const uint64_t*)u[0],2,2);
}


/* vim:set ft=cpp: */
