#ifndef ABASE_P_3_H_
#define ABASE_P_3_H_

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
#include "fixmp.h"
#include "mpfq_gfp_common.h"
#include "select_mpi.h"
#include "abase_vbase.h"
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

typedef mpfq_p_field abase_p_3_field;
typedef mpfq_p_dst_field abase_p_3_dst_field;

typedef unsigned long abase_p_3_elt[3];
typedef unsigned long * abase_p_3_dst_elt;
typedef const unsigned long * abase_p_3_src_elt;

typedef unsigned long abase_p_3_elt_ur[7];
typedef unsigned long * abase_p_3_dst_elt_ur;
typedef const unsigned long * abase_p_3_src_elt_ur;

typedef abase_p_3_elt * abase_p_3_vec;
typedef abase_p_3_elt * abase_p_3_dst_vec;
typedef abase_p_3_elt * abase_p_3_src_vec;

typedef abase_p_3_elt_ur * abase_p_3_vec_ur;
typedef abase_p_3_elt_ur * abase_p_3_dst_vec_ur;
typedef abase_p_3_elt_ur * abase_p_3_src_vec_ur;
/* Extra types defined by implementation: */
typedef struct {
  abase_p_3_vec c;
  unsigned int alloc;
  unsigned int size;
} abase_p_3_poly_struct;
typedef abase_p_3_poly_struct abase_p_3_poly [1];
typedef abase_p_3_poly_struct * abase_p_3_dst_poly;
typedef abase_p_3_poly_struct * abase_p_3_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif

/* Functions operating on the field structure */
void abase_p_3_field_characteristic(abase_p_3_dst_field, mpz_t);
/* *Mpfq::gfp::field::code_for_field_degree, Mpfq::gfp */
#define abase_p_3_field_degree(K)	1
static inline
void abase_p_3_field_init(abase_p_3_dst_field);
void abase_p_3_field_clear(abase_p_3_dst_field);
void abase_p_3_field_specify(abase_p_3_dst_field, unsigned long, void *);
/* *Mpfq::gfp::field::code_for_field_setopt, Mpfq::gfp */
#define abase_p_3_field_setopt(f, x, y)	/**/

/* Element allocation functions */
static inline
void abase_p_3_init(abase_p_3_dst_field, abase_p_3_elt *);
static inline
void abase_p_3_clear(abase_p_3_dst_field, abase_p_3_elt *);

/* Elementary assignment functions */
static inline
void abase_p_3_set(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt);
static inline
void abase_p_3_set_ui(abase_p_3_dst_field, abase_p_3_dst_elt, unsigned long);
static inline
void abase_p_3_set_zero(abase_p_3_dst_field, abase_p_3_dst_elt);
static inline
unsigned long abase_p_3_get_ui(abase_p_3_dst_field, abase_p_3_src_elt);
static inline
void abase_p_3_set_mpn(abase_p_3_dst_field, abase_p_3_dst_elt, mp_limb_t *, size_t);
static inline
void abase_p_3_set_mpz(abase_p_3_dst_field, abase_p_3_dst_elt, mpz_t);
static inline
void abase_p_3_get_mpn(abase_p_3_dst_field, mp_limb_t *, abase_p_3_src_elt);
static inline
void abase_p_3_get_mpz(abase_p_3_dst_field, mpz_t, abase_p_3_src_elt);

/* Assignment of random values */
static inline
void abase_p_3_normalize(abase_p_3_dst_field, abase_p_3_dst_elt);
static inline
void abase_p_3_random(abase_p_3_dst_field, abase_p_3_dst_elt);
#define HAVE_abase_p_3_random2
static inline
void abase_p_3_random2(abase_p_3_dst_field, abase_p_3_dst_elt);

/* Arithmetic operations on elements */
static inline
void abase_p_3_add(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, abase_p_3_src_elt);
static inline
void abase_p_3_sub(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, abase_p_3_src_elt);
static inline
void abase_p_3_neg(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt);
static inline
void abase_p_3_mul(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, abase_p_3_src_elt);
static inline
void abase_p_3_sqr(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt);
static inline
int abase_p_3_is_sqr(abase_p_3_dst_field, abase_p_3_src_elt);
int abase_p_3_sqrt(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt);
static inline
void abase_p_3_pow(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, unsigned long *, size_t);
/* *Mpfq::gfp::elt::code_for_frobenius, Mpfq::gfp */
#define abase_p_3_frobenius(k, x, y)	abase_p_3_set(k, x, y)
static inline
void abase_p_3_add_ui(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, unsigned long);
static inline
void abase_p_3_sub_ui(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, unsigned long);
static inline
void abase_p_3_mul_ui(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt, unsigned long);
static inline
int abase_p_3_inv(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_elt);
#define HAVE_abase_p_3_hadamard
static inline
void abase_p_3_hadamard(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_dst_elt, abase_p_3_dst_elt, abase_p_3_dst_elt);

/* Operations involving unreduced elements */
static inline
void abase_p_3_elt_ur_init(abase_p_3_dst_field, abase_p_3_elt_ur *);
static inline
void abase_p_3_elt_ur_clear(abase_p_3_dst_field, abase_p_3_elt_ur *);
static inline
void abase_p_3_elt_ur_set(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt_ur);
static inline
void abase_p_3_elt_ur_set_zero(abase_p_3_dst_field, abase_p_3_dst_elt_ur);
static inline
void abase_p_3_elt_ur_set_ui(abase_p_3_dst_field, abase_p_3_dst_elt_ur, unsigned long);
static inline
void abase_p_3_elt_ur_add(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt_ur, abase_p_3_src_elt_ur);
static inline
void abase_p_3_elt_ur_neg(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt_ur);
static inline
void abase_p_3_elt_ur_sub(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt_ur, abase_p_3_src_elt_ur);
static inline
void abase_p_3_mul_ur(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt, abase_p_3_src_elt);
static inline
void abase_p_3_sqr_ur(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt);
static inline
void abase_p_3_reduce(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_dst_elt_ur);
#define HAVE_abase_p_3_addmul_si_ur
static inline
void abase_p_3_addmul_si_ur(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_elt, long);

/* Comparison functions */
static inline
int abase_p_3_cmp(abase_p_3_dst_field, abase_p_3_src_elt, abase_p_3_src_elt);
static inline
int abase_p_3_cmp_ui(abase_p_3_dst_field, abase_p_3_src_elt, unsigned long);
static inline
int abase_p_3_is_zero(abase_p_3_dst_field, abase_p_3_src_elt);

/* Input/output functions */
void abase_p_3_asprint(abase_p_3_dst_field, char * *, abase_p_3_src_elt);
void abase_p_3_fprint(abase_p_3_dst_field, FILE *, abase_p_3_src_elt);
/* *Mpfq::defaults::code_for_print, Mpfq::gfp */
#define abase_p_3_print(k, x)	abase_p_3_fprint(k,stdout,x)
int abase_p_3_sscan(abase_p_3_dst_field, abase_p_3_dst_elt, const char *);
int abase_p_3_fscan(abase_p_3_dst_field, FILE *, abase_p_3_dst_elt);
/* *Mpfq::gfp::io::code_for_scan, Mpfq::gfp */
#define abase_p_3_scan(k, x)	abase_p_3_fscan(k,stdout,x)

/* Vector functions */
void abase_p_3_vec_init(abase_p_3_dst_field, abase_p_3_vec *, unsigned int);
void abase_p_3_vec_reinit(abase_p_3_dst_field, abase_p_3_vec *, unsigned int, unsigned int);
void abase_p_3_vec_clear(abase_p_3_dst_field, abase_p_3_vec *, unsigned int);
static inline
void abase_p_3_vec_set(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_set_zero(abase_p_3_dst_field, abase_p_3_dst_vec, unsigned int);
static inline
void abase_p_3_vec_setcoef(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_elt, unsigned int);
static inline
void abase_p_3_vec_setcoef_ui(abase_p_3_dst_field, abase_p_3_dst_vec, unsigned long, unsigned int);
static inline
void abase_p_3_vec_getcoef(abase_p_3_dst_field, abase_p_3_dst_elt, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_add(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_neg(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_rev(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_sub(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_scal_mul(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, abase_p_3_src_elt, unsigned int);
static inline
void abase_p_3_vec_conv(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, unsigned int, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_random(abase_p_3_dst_field, abase_p_3_dst_vec, unsigned int);
#define HAVE_abase_p_3_vec_random2
static inline
void abase_p_3_vec_random2(abase_p_3_dst_field, abase_p_3_dst_vec, unsigned int);
static inline
int abase_p_3_vec_cmp(abase_p_3_dst_field, abase_p_3_src_vec, abase_p_3_src_vec, unsigned int);
static inline
int abase_p_3_vec_is_zero(abase_p_3_dst_field, abase_p_3_src_vec, unsigned int);
void abase_p_3_vec_asprint(abase_p_3_dst_field, char * *, abase_p_3_src_vec, unsigned int);
void abase_p_3_vec_fprint(abase_p_3_dst_field, FILE *, abase_p_3_src_vec, unsigned int);
void abase_p_3_vec_print(abase_p_3_dst_field, abase_p_3_src_vec, unsigned int);
int abase_p_3_vec_sscan(abase_p_3_dst_field, abase_p_3_vec *, unsigned int *, const char *);
int abase_p_3_vec_fscan(abase_p_3_dst_field, FILE *, abase_p_3_vec *, unsigned int *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
#define abase_p_3_vec_scan(K, w, n)	abase_p_3_vec_fscan(K,stdout,w,n)
void abase_p_3_vec_ur_init(abase_p_3_dst_field, abase_p_3_vec_ur *, unsigned int);
static inline
void abase_p_3_vec_ur_set_zero(abase_p_3_dst_field, abase_p_3_dst_vec_ur, unsigned int);
void abase_p_3_vec_ur_reinit(abase_p_3_dst_field, abase_p_3_vec_ur *, unsigned int, unsigned int);
void abase_p_3_vec_ur_clear(abase_p_3_dst_field, abase_p_3_vec_ur *, unsigned int);
static inline
void abase_p_3_vec_ur_set(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec_ur, unsigned int);
static inline
void abase_p_3_vec_ur_setcoef(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_elt_ur, unsigned int);
static inline
void abase_p_3_vec_ur_getcoef(abase_p_3_dst_field, abase_p_3_dst_elt_ur, abase_p_3_src_vec_ur, unsigned int);
static inline
void abase_p_3_vec_ur_add(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec_ur, abase_p_3_src_vec_ur, unsigned int);
static inline
void abase_p_3_vec_ur_sub(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec_ur, abase_p_3_src_vec_ur, unsigned int);
static inline
void abase_p_3_vec_scal_mul_ur(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec, abase_p_3_src_elt, unsigned int);
static inline
void abase_p_3_vec_conv_ur_n(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec, abase_p_3_src_vec, unsigned int);
void abase_p_3_vec_conv_ur_ks(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec, unsigned int, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_conv_ur(abase_p_3_dst_field, abase_p_3_dst_vec_ur, abase_p_3_src_vec, unsigned int, abase_p_3_src_vec, unsigned int);
static inline
void abase_p_3_vec_reduce(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_dst_vec_ur, unsigned int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define abase_p_3_vec_elt_stride(K, n)	((n)*sizeof(abase_p_3_elt))

/* Functions related to SIMD operation */
/* *simd_gfp::code_for_groupsize */
#define abase_p_3_groupsize(K)	1
/* *simd_gfp::code_for_offset */
#define abase_p_3_offset(K, n)	n /* TO BE DEPRECATED */
/* *simd_gfp::code_for_stride */
#define abase_p_3_stride(K)	1 /* TO BE DEPRECATED */
static inline
void abase_p_3_set_ui_at(abase_p_3_dst_field, abase_p_3_dst_elt, int, unsigned long);
static inline
void abase_p_3_set_ui_all(abase_p_3_dst_field, abase_p_3_dst_elt, unsigned long);
static inline
void abase_p_3_elt_ur_set_ui_at(abase_p_3_dst_field, abase_p_3_dst_elt, int, unsigned long);
static inline
void abase_p_3_elt_ur_set_ui_all(abase_p_3_dst_field, abase_p_3_dst_elt, unsigned long);
void abase_p_3_dotprod(abase_p_3_dst_field, abase_p_3_dst_vec, abase_p_3_src_vec, abase_p_3_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* MPI interface */
void abase_p_3_mpi_ops_init(abase_p_3_dst_field);
MPI_Datatype abase_p_3_mpi_datatype(abase_p_3_dst_field);
MPI_Datatype abase_p_3_mpi_datatype_ur(abase_p_3_dst_field);
MPI_Op abase_p_3_mpi_addition_op(abase_p_3_dst_field);
MPI_Op abase_p_3_mpi_addition_op_ur(abase_p_3_dst_field);
void abase_p_3_mpi_ops_clear(abase_p_3_dst_field);

/* Object-oriented interface */
#define abase_p_3_oo_impl_name(v)	"p_3"
static inline
void abase_p_3_oo_field_clear(abase_vbase_ptr);
void abase_p_3_oo_field_init(abase_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::gfp::field::code_for_field_init, Mpfq::gfp */
static inline
void abase_p_3_field_init(abase_p_3_dst_field k)
{
    k->p = NULL;
    k->bigmul_p = NULL;
    k->io_base = 10;
    mpz_init(k->factor);
    k->ts_info.e=0;
}

/* *Mpfq::gfp::elt::code_for_init, Mpfq::gfp */
static inline
void abase_p_3_init(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(k->p);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_clear, Mpfq::gfp */
static inline
void abase_p_3_clear(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::defaults::flatdata::code_for_set, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_set(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt r, abase_p_3_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(abase_p_3_elt));
}

/* *Mpfq::gfp::elt::code_for_set_ui, Mpfq::gfp */
static inline
void abase_p_3_set_ui(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt r, unsigned long x)
{
    int i; 
    assert (r);
    r[0] = x;
    for (i = 1; i < 3; ++i)
        r[i] = 0;
}

/* *Mpfq::defaults::flatdata::code_for_set_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_set_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt r)
{
    abase_p_3_vec_set_zero(K,(abase_p_3_dst_vec)r,1);
}

/* *Mpfq::gfp::elt::code_for_get_ui, Mpfq::gfp */
static inline
unsigned long abase_p_3_get_ui(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_src_elt x)
{
    return x[0];
}

/* *Mpfq::gfp::elt::code_for_set_mpn, Mpfq::gfp */
static inline
void abase_p_3_set_mpn(abase_p_3_dst_field k, abase_p_3_dst_elt r, mp_limb_t * x, size_t n)
{
    int i;
    if (n < 3) {
        for (i = 0; i < (int)n; ++i)
            r[i] = x[i];
        for (i = n; i < 3; ++i)
            r[i] = 0;
    } else {
        mp_limb_t tmp[n-3+1];
        mpn_tdiv_qr(tmp, r, 0, x, n, k->p, 3);
    }
}

/* *Mpfq::gfp::elt::code_for_set_mpz, Mpfq::gfp */
static inline
void abase_p_3_set_mpz(abase_p_3_dst_field k, abase_p_3_dst_elt r, mpz_t z)
{
    if (z->_mp_size < 0) {
        abase_p_3_set_mpn(k, r, z->_mp_d, -z->_mp_size);
        abase_p_3_neg(k, r, r);
    } else {
        abase_p_3_set_mpn(k, r, z->_mp_d, z->_mp_size);
    }
}

/* *Mpfq::gfp::elt::code_for_get_mpn, Mpfq::gfp */
static inline
void abase_p_3_get_mpn(abase_p_3_dst_field k MAYBE_UNUSED, mp_limb_t * r, abase_p_3_src_elt x)
{
    int i; 
    assert (r);
    assert (x);
    for (i = 0; i < 3; ++i)
        r[i] = x[i];
}

/* *Mpfq::gfp::elt::code_for_get_mpz, Mpfq::gfp */
static inline
void abase_p_3_get_mpz(abase_p_3_dst_field k MAYBE_UNUSED, mpz_t z, abase_p_3_src_elt y)
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
void abase_p_3_normalize(abase_p_3_dst_field k, abase_p_3_dst_elt x)
{
    if (cmp_3(x,k->p)>=0) {
      mp_limb_t q[3+1];
      abase_p_3_elt r;
      mpn_tdiv_qr(q, r, 0, x, 3, k->p, 3);
      abase_p_3_set(k, x, r);
    }
}

/* *Mpfq::gfp::elt::code_for_random, Mpfq::gfp */
static inline
void abase_p_3_random(abase_p_3_dst_field k, abase_p_3_dst_elt x)
{
    mpn_random(x, 3);
    abase_p_3_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_random2, Mpfq::gfp */
static inline
void abase_p_3_random2(abase_p_3_dst_field k, abase_p_3_dst_elt x)
{
    mpn_random2(x, 3);
    abase_p_3_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_add, Mpfq::gfp */
static inline
void abase_p_3_add(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, abase_p_3_src_elt y)
{
    mp_limb_t cy;
    cy = add_3(z, x, y);
    if (cy || (cmp_3(z, k->p) >= 0))
        sub_3(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_sub, Mpfq::gfp */
static inline
void abase_p_3_sub(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, abase_p_3_src_elt y)
{
    mp_limb_t cy;
    cy = sub_3(z, x, y);
    if (cy) // negative result
        add_3(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_neg, Mpfq::gfp */
static inline
void abase_p_3_neg(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x)
{
    if (cmp_ui_3(x, 0))
        sub_3(z, k->p, x);
    else {
        int i;
        for (i = 0; i < 3; ++i)
            z[i] = 0;
        }
}

/* *Mpfq::gfp::elt::code_for_mul, Mpfq::gfp */
static inline
void abase_p_3_mul(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, abase_p_3_src_elt y)
{
    mp_limb_t tmp[2*3];
    mul_3(tmp, x, y);
    mod_3(z, tmp, k->p);
}

/* *Mpfq::gfp::elt::code_for_sqr, Mpfq::gfp */
static inline
void abase_p_3_sqr(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x)
{
    mp_limb_t tmp[2*3];
    sqr_3(tmp, x);
    mod_3(z, tmp, k->p);
}

/* *Mpfq::gfp::elt::code_for_is_sqr, Mpfq::gfp */
static inline
int abase_p_3_is_sqr(abase_p_3_dst_field k, abase_p_3_src_elt x)
{
    mp_limb_t pp[3];
    abase_p_3_elt y;
    sub_ui_nc_3(pp, k->p, 1);
    rshift_3(pp, 1);
    abase_p_3_init(k, &y);
    abase_p_3_pow(k, y, x, pp, 3);
    int res = cmp_ui_3(y, 1);
    abase_p_3_clear(k, &y);
    if (res == 0)
        return 1;
    else 
        return 0;
}

/* *Mpfq::defaults::pow::code_for_pow, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_pow(abase_p_3_dst_field k, abase_p_3_dst_elt res, abase_p_3_src_elt r, unsigned long * x, size_t n)
{
    abase_p_3_elt u, a;
    long i, j, lead;     /* it is a signed type */
    unsigned long mask;
    
    assert (n>0);
    
    /* get the correct (i,j) position of the most significant bit in x */
    for(i = n-1; i>=0 && x[i]==0; i--)
        ;
    if (i < 0) {
        abase_p_3_set_ui(k, res, 0);
        return;
    }
    j = 64 - 1;
    mask = (1UL<<j);
    for( ; (x[i]&mask)==0 ;j--, mask>>=1)
        ;
    lead = i*64+j;      /* Ensured. */
    
    abase_p_3_init(k, &u);
    abase_p_3_init(k, &a);
    abase_p_3_set(k, a, r);
    for( ; lead > 0; lead--) {
        if (j-- == 0) {
            i--;
            j = 64-1;
            mask = (1UL<<j);
        } else {
            mask >>= 1;
        }
        if (x[i]&mask) {
            abase_p_3_sqr(k, u, a);
            abase_p_3_mul(k, a, u, r);
        } else {
            abase_p_3_sqr(k, a,a);
        }
    }
    abase_p_3_set(k, res, a);
    abase_p_3_clear(k, &u);
    abase_p_3_clear(k, &a);
}

/* *Mpfq::gfp::elt::code_for_add_ui, Mpfq::gfp */
static inline
void abase_p_3_add_ui(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = add_ui_3(z, x, y);
    if (cy || (cmp_3(z, k->p) >= 0))
        sub_3(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_sub_ui, Mpfq::gfp */
static inline
void abase_p_3_sub_ui(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = sub_ui_3(z, x, y);
    if (cy) // negative result
        add_3(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_mul_ui, Mpfq::gfp */
static inline
void abase_p_3_mul_ui(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x, unsigned long y)
{
    abase_p_3_elt yy;
    abase_p_3_init(k, &yy);
    abase_p_3_set_ui(k, yy, y);
    abase_p_3_mul(k, z, x, yy);
    abase_p_3_clear(k, &yy);
}

/* *Mpfq::gfp::elt::code_for_inv, Mpfq::gfp */
static inline
int abase_p_3_inv(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_src_elt x)
{
    int ret=invmod_3(z, x, k->p);
    if (!ret)
        abase_p_3_get_mpz(k, k->factor, z);
    return ret;
}

/* *Mpfq::gfp::elt::code_for_hadamard, Mpfq::gfp */
static inline
void abase_p_3_hadamard(abase_p_3_dst_field k, abase_p_3_dst_elt x, abase_p_3_dst_elt y, abase_p_3_dst_elt z, abase_p_3_dst_elt t)
{
    abase_p_3_elt tmp;
    abase_p_3_init(k, &tmp);
    abase_p_3_add(k, tmp, x, y);
    abase_p_3_sub(k, y, x, y);
    abase_p_3_set(k, x, tmp);
    abase_p_3_add(k, tmp, z, t);
    abase_p_3_sub(k, t, z, t);
    abase_p_3_set(k, z, tmp);
    abase_p_3_sub(k, tmp, x, z);
    abase_p_3_add(k, x, x, z);
    abase_p_3_add(k, z, y, t);
    abase_p_3_sub(k, t, y, t);
    abase_p_3_set(k, y, tmp);
    abase_p_3_clear(k, &tmp); 
}

/* *Mpfq::gfp::elt::code_for_elt_ur_init, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_init(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(k->p);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_clear, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_clear(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_set(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur z, abase_p_3_src_elt_ur x)
{
    int i;
    for (i = 0; i < 7; ++i) 
        z[i] = x[i];
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_set_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt_ur r)
{
    memset(r, 0, sizeof(abase_p_3_elt_ur));
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set_ui, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_set_ui(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur r, unsigned long x)
{
    int i; 
    assert (r); 
    r[0] = x;
    for (i = 1; i < 7; ++i)
        r[i] = 0;
}

/* *Mpfq::gfp::elt::code_for_elt_ur_add, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_add(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur z, abase_p_3_src_elt_ur x, abase_p_3_src_elt_ur y)
{
    mpn_add_n(z, x, y, 7);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_neg, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_neg(abase_p_3_dst_field k, abase_p_3_dst_elt_ur z, abase_p_3_src_elt_ur x)
{
    abase_p_3_elt_ur tmp;
    abase_p_3_elt_ur_init(k, &tmp);
    int i;
    for (i = 0; i < 7; ++i) 
        tmp[i] = 0;
    mpn_sub_n(z, tmp, x, 7);
    abase_p_3_elt_ur_clear(k, &tmp);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_sub, Mpfq::gfp */
static inline
void abase_p_3_elt_ur_sub(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur z, abase_p_3_src_elt_ur x, abase_p_3_src_elt_ur y)
{
    mpn_sub_n(z, x, y, 7);
}

/* *Mpfq::gfp::elt::code_for_mul_ur, Mpfq::gfp */
static inline
void abase_p_3_mul_ur(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur z, abase_p_3_src_elt x, abase_p_3_src_elt y)
{
    mul_3(z, x, y);
    int i;
    for (i = 2*3; i < 7; ++i) {
        z[i] = 0;
    }
}

/* *Mpfq::gfp::elt::code_for_sqr_ur, Mpfq::gfp */
static inline
void abase_p_3_sqr_ur(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_dst_elt_ur z, abase_p_3_src_elt x)
{
    sqr_3(z, x);
    int i;
    for (i = 2*3; i < 7; ++i) {
        z[i] = 0;
    }
}

/* *Mpfq::gfp::elt::code_for_reduce, Mpfq::gfp */
static inline
void abase_p_3_reduce(abase_p_3_dst_field k, abase_p_3_dst_elt z, abase_p_3_dst_elt_ur x)
{
    mp_limb_t q[7+1];
    if (x[7-1]>>(64-1)) {
        // negative number, add bigmul_p to make it positive before reduction
        mpn_add_n(x, x, k->bigmul_p, 7);
    }
    mpn_tdiv_qr(q, z, 0, x, 7, k->p, 3);
}

/* *simd_gfp::code_for_addmul_si_ur */
static inline
void abase_p_3_addmul_si_ur(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt_ur w, abase_p_3_src_elt u, long v)
{
        abase_p_3_elt_ur s;
        abase_p_3_elt vx;
        abase_p_3_elt_ur_init(K, &s);
        abase_p_3_init(K, &vx);
        if (v>0) {
            abase_p_3_set_ui(K, vx, v);
            abase_p_3_mul_ur(K, s, u, vx);
            abase_p_3_elt_ur_add(K, w, w, s);
        } else {
            abase_p_3_set_ui(K, vx, -v);
            abase_p_3_mul_ur(K, s, u, vx);
            abase_p_3_elt_ur_sub(K, w, w, s);
        }
        abase_p_3_clear(K, &vx);
        abase_p_3_elt_ur_clear(K, &s);
}

/* *Mpfq::gfp::elt::code_for_cmp, Mpfq::gfp */
static inline
int abase_p_3_cmp(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_src_elt x, abase_p_3_src_elt y)
{
    return cmp_3(x,y);
}

/* *Mpfq::gfp::elt::code_for_cmp_ui, Mpfq::gfp */
static inline
int abase_p_3_cmp_ui(abase_p_3_dst_field k MAYBE_UNUSED, abase_p_3_src_elt x, unsigned long y)
{
    return cmp_ui_3(x,y);
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
int abase_p_3_is_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(abase_p_3_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_vec_set(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec r, abase_p_3_src_vec s, unsigned int n)
{
    if (r != s) memcpy(r, s, n*sizeof(abase_p_3_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_vec_set_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_p_3_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoef, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_setcoef(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_elt x, unsigned int i)
{
    abase_p_3_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoef_ui, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_setcoef_ui(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, unsigned long x, unsigned int i)
{
    abase_p_3_set_ui(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoef, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_getcoef(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt x, abase_p_3_src_vec w, unsigned int i)
{
    abase_p_3_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_add(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, abase_p_3_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_neg(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_3_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_rev(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_p_3_elt tmp[1];
    abase_p_3_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_p_3_set(K, tmp[0], u[i]);
        abase_p_3_set(K, w[i], u[n-1-i]);
        abase_p_3_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_p_3_set(K, w[nn], u[nn]);
    abase_p_3_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_sub(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, abase_p_3_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_3_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_scal_mul(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, abase_p_3_src_elt x, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_mul(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv, Mpfq::gfp */
static inline
void abase_p_3_vec_conv(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_src_vec u, unsigned int n, abase_p_3_src_vec v, unsigned int m)
{
    abase_p_3_vec_ur tmp;
    abase_p_3_vec_ur_init(K, &tmp, m+n-1);
    abase_p_3_vec_conv_ur(K, tmp, u, n, v, m);
    abase_p_3_vec_reduce(K, w, tmp, m+n-1);
    abase_p_3_vec_ur_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_random(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_3_random(K, w[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random2, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_random2(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_3_random2(K, w[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
int abase_p_3_vec_cmp(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_src_vec u, abase_p_3_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = abase_p_3_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_is_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
int abase_p_3_vec_is_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!abase_p_3_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_set_zero(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_p_3_elt_ur));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_set(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur r, abase_p_3_src_vec_ur s, unsigned int n)
{
    if (r != s) memcpy(r, s, n*sizeof(abase_p_3_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoef, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_setcoef(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_elt_ur x, unsigned int i)
{
    abase_p_3_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoef, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_getcoef(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt_ur x, abase_p_3_src_vec_ur w, unsigned int i)
{
    abase_p_3_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_add(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_vec_ur u, abase_p_3_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_ur_sub(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_vec_ur u, abase_p_3_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_scal_mul_ur(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_vec u, abase_p_3_src_elt x, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_mul_ur(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
static inline
void abase_p_3_vec_conv_ur_n(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_vec u, abase_p_3_src_vec v, unsigned int n)
{
    if (n == 0)
        return;
    if (n == 1) {
        abase_p_3_mul_ur(K, w[0], u[0], v[0]);
        return;
    }
    if (n == 2) {  // Kara 2
        abase_p_3_elt t1, t2;
        abase_p_3_init(K, &t1);
        abase_p_3_init(K, &t2);
        abase_p_3_mul_ur(K, w[0], u[0], v[0]);
        abase_p_3_mul_ur(K, w[2], u[1], v[1]);
        abase_p_3_add(K, t1, u[0], u[1]);
        abase_p_3_add(K, t2, v[0], v[1]);
        abase_p_3_mul_ur(K, w[1], t1, t2);
        abase_p_3_elt_ur_sub(K, w[1], w[1], w[0]);
        abase_p_3_elt_ur_sub(K, w[1], w[1], w[2]);
        abase_p_3_clear(K, &t1);
        abase_p_3_clear(K, &t2);
        return;
    }
    if (n == 3) {  // do it in 6
        abase_p_3_elt t1, t2;
        abase_p_3_elt_ur s;
        abase_p_3_init(K, &t1);
        abase_p_3_init(K, &t2);
        abase_p_3_elt_ur_init(K, &s);
        // a0*b0*(1 - X)
        abase_p_3_mul_ur(K, w[0], u[0], v[0]);
        abase_p_3_elt_ur_neg(K, w[1], w[0]);
        // a1*b1*(-X + 2*X^2 - X^3)
        abase_p_3_mul_ur(K, w[2], u[1], v[1]);
        abase_p_3_elt_ur_neg(K, w[3], w[2]);
        abase_p_3_elt_ur_add(K, w[2], w[2], w[2]);
        abase_p_3_elt_ur_add(K, w[1], w[1], w[3]);
        // a2*b2*(-X^3+X^4)
        abase_p_3_mul_ur(K, w[4], u[2], v[2]);
        abase_p_3_elt_ur_sub(K, w[3], w[3], w[4]);
        // (a0+a1)*(b0+b1)*(X - X^2)
        abase_p_3_add(K, t1, u[0], u[1]);
        abase_p_3_add(K, t2, v[0], v[1]);
        abase_p_3_mul_ur(K, s, t1, t2);
        abase_p_3_elt_ur_add(K, w[1], w[1], s);
        abase_p_3_elt_ur_sub(K, w[2], w[2], s);
        // (a1+a2)*(b1+b2)*(X^3 - X^2)
        abase_p_3_add(K, t1, u[1], u[2]);
        abase_p_3_add(K, t2, v[1], v[2]);
        abase_p_3_mul_ur(K, s, t1, t2);
        abase_p_3_elt_ur_add(K, w[3], w[3], s);
        abase_p_3_elt_ur_sub(K, w[2], w[2], s);
        // (a0+a1+a2)*(b0+b1+b2)* X^2
        abase_p_3_add(K, t1, u[0], t1);
        abase_p_3_add(K, t2, v[0], t2);
        abase_p_3_mul_ur(K, s, t1, t2);
        abase_p_3_elt_ur_add(K, w[2], w[2], s);
        return;
    }
    unsigned int n0, n1;
    n0 = n / 2;
    n1 = n - n0;
    abase_p_3_vec_conv_ur_n(K, w, u, v, n0);
    abase_p_3_vec_conv_ur_n(K, w + 2*n0, u + n0, v + n0, n1);
    abase_p_3_elt_ur_set_ui(K, w[2*n0-1], 0);
    
    abase_p_3_vec tmpu, tmpv;
    abase_p_3_vec_ur tmpw;
    abase_p_3_vec_init(K, &tmpu, n1);
    abase_p_3_vec_init(K, &tmpv, n1);
    abase_p_3_vec_ur_init(K, &tmpw, 2*n1-1);
    
    abase_p_3_vec_set(K, tmpu, u, n0);
    if (n1 != n0) 
        abase_p_3_set_ui(K, tmpu[n0], 0);
    abase_p_3_vec_add(K, tmpu, tmpu, u+n0, n1);
    abase_p_3_vec_set(K, tmpv, v, n0);
    if (n1 != n0) 
        abase_p_3_set_ui(K, tmpv[n0], 0);
    abase_p_3_vec_add(K, tmpv, tmpv, v+n0, n1);
    abase_p_3_vec_conv_ur_n(K, tmpw, tmpu, tmpv, n1);
    abase_p_3_vec_ur_sub(K, tmpw, tmpw, w, 2*n0-1);
    abase_p_3_vec_ur_sub(K, tmpw, tmpw, w + 2*n0, 2*n1-1);
    abase_p_3_vec_ur_add(K, w + n0, w + n0, tmpw, 2*n1-1);
    
    abase_p_3_vec_clear(K, &tmpu, n1);
    abase_p_3_vec_clear(K, &tmpv, n1);
    abase_p_3_vec_ur_clear(K, &tmpw, 2*n1-1);
    return;
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
static inline
void abase_p_3_vec_conv_ur(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec_ur w, abase_p_3_src_vec u, unsigned int n, abase_p_3_src_vec v, unsigned int m)
{
    if ((n > 1) && (m > 1) && (n+m > 15)) {
        abase_p_3_vec_conv_ur_ks(K, w, u, n, v, m);
        return;
    }
    if (n == m) {
        abase_p_3_vec_conv_ur_n(K, w, u, v, n);
        return;
    }
    unsigned int i, j MAYBE_UNUSED, k;
    abase_p_3_elt_ur acc, z;
    abase_p_3_elt_ur_init(K, &acc);
    abase_p_3_elt_ur_init(K, &z);
    // swap pointers to have n <= m
    abase_p_3_src_vec uu, vv;
    if (n <= m) {
        uu = u; vv = v;
    } else {
        uu = v; vv = u;
        unsigned int tmp = n;
        n = m; m = tmp;
    }
    for(k = 0; k < n; ++k) {
        abase_p_3_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i <= k; ++i) {
            abase_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_3_elt_ur_add(K, acc, acc, z);
        }
        abase_p_3_elt_ur_set(K, w[k], acc);
    }
    for(k = n; k < m; ++k) {
        abase_p_3_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i < n; ++i) {
            abase_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_3_elt_ur_add(K, acc, acc, z);
        }
        abase_p_3_elt_ur_set(K, w[k], acc);
    }
    for(k = m; k < n+m-1; ++k) {
        abase_p_3_mul_ur(K, acc, uu[k-m+1], vv[m-1]);
        for(i = k-m+2; i < n; ++i) {
            abase_p_3_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_3_elt_ur_add(K, acc, acc, z);
        }
        abase_p_3_elt_ur_set(K, w[k], acc);
    }
    abase_p_3_elt_ur_clear(K, &acc);
    abase_p_3_elt_ur_clear(K, &z);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec, Mpfq::defaults, Mpfq::gfp */
static inline
void abase_p_3_vec_reduce(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_vec w, abase_p_3_dst_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_3_reduce(K, w[i], u[i]);
}

/* *simd_gfp::code_for_set_ui_at */
static inline
void abase_p_3_set_ui_at(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_p_3_set_ui(K,p,v);
}

/* *simd_gfp::code_for_set_ui_all */
static inline
void abase_p_3_set_ui_all(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt p, unsigned long v)
{
    abase_p_3_set_ui(K,p,v);
}

/* *simd_gfp::code_for_elt_ur_set_ui_at */
static inline
void abase_p_3_elt_ur_set_ui_at(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_p_3_set_ui(K,p,v);
}

/* *simd_gfp::code_for_elt_ur_set_ui_all */
static inline
void abase_p_3_elt_ur_set_ui_all(abase_p_3_dst_field K MAYBE_UNUSED, abase_p_3_dst_elt p, unsigned long v)
{
    abase_p_3_set_ui(K,p,v);
}

static inline
void abase_p_3_oo_field_clear(abase_vbase_ptr f)
{
    abase_p_3_field_clear((abase_p_3_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* ABASE_P_3_H_ */

/* vim:set ft=cpp: */
