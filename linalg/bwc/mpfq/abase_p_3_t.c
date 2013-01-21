/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_p_3_t.h"

/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used: w=64 fieldtype=prime n=3 nn=7 vtag=p_3 vbase_stuff={
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
 tag=p_3 type=plain virtual_base={
                  'filebase' => 'abase_vbase',
                  'substitutions' => [
                                       [
                                         qr/(?^:abase_p_3_elt \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_src_elt\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_dst_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_elt_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_src_elt_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_dst_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_vec \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_src_vec\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_dst_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_vec_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_src_vec_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_dst_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_poly \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_src_poly\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_poly\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p_3_dst_poly\b)/,
                                         'void *'
                                       ]
                                     ],
                  'name' => 'abase_vbase',
                  'global_prefix' => 'abase_'
                };
 family=[HASH(0x1759fe8)] */


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
/* *simd_gfp::code_for_member_template_dotprod */
void abase_p_3_p_3_dotprod(abase_p_3_dst_field K0 MAYBE_UNUSED, abase_p_3_dst_field K1 MAYBE_UNUSED, abase_p_3_dst_vec xw, abase_p_3_src_vec xu1, abase_p_3_src_vec xu0, unsigned int n)
{
        abase_p_3_elt_ur s,t;
        abase_p_3_elt_ur_init(K0, &s);
        abase_p_3_elt_ur_init(K0, &t);
        abase_p_3_elt_ur_set_zero(K0, s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_p_3_mul_ur(K0, t, xu0[i], xu1[i]);
            abase_p_3_elt_ur_add(K0, s, s, t);
        }
        abase_p_3_reduce(K0, xw[0], s);
        abase_p_3_elt_ur_clear(K0, &s);
        abase_p_3_elt_ur_clear(K0, &t);
}

/* *simd_gfp::code_for_member_template_addmul_tiny */
void abase_p_3_p_3_addmul_tiny(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_field L MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, abase_p_3_dst_vec v, unsigned int n)
{
        abase_p_3_elt s;
        abase_p_3_init(K, &s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_p_3_mul(K, s, u[i], v[0]);
            abase_p_3_add(K, w[i], w[i], s);
        }
        abase_p_3_clear(K, &s);
}

/* *simd_gfp::code_for_member_template_transpose */
void abase_p_3_p_3_transpose(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_field L MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u)
{
    abase_p_3_set(K, w[0], u[0]);
}


/* vim:set ft=cpp: */
