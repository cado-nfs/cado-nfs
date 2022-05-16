#include "cado.h"
/* MPFQ generated file -- do not edit */

#include "mpfq_pz.h"

#include <inttypes.h>
#include <limits.h>
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
/* *Mpfq::gfp::field::code_for_field_clear, pz */
void mpfq_pz_field_clear(mpfq_pz_dst_field k)
{
        mpz_clear(k->p);
        mpz_clear(k->bigmul_p);
        if (k->ts_info.e > 0) {
            free(k->ts_info.hh);
            free(k->ts_info.z);
        }
        mpz_clear(k->factor);
}

/* *pz::code_for_field_specify */
void mpfq_pz_field_specify(mpfq_pz_dst_field k, unsigned long dummy MAYBE_UNUSED, const void * vp)
{
        if (dummy == MPFQ_PRIME_MPN) {
            fprintf(stderr, "MPFQ_PRIME_MPN is deprecated\n");
            return;
        } else if (dummy == MPFQ_PRIME_MPZ) {
            mpz_srcptr p = (mpz_srcptr) vp;
            mpz_set(k->p, p);
            if (!(mpz_size(k->p) > 0)) abort();
            {
                /* precompute bigmul_p = largest multiple of p that fits in an
                 * elt_ur: p*Floor( (2^((2*mpz_size(k->p) + 1)*64)-1)/p )
                 */
                mpz_ui_pow_ui(k->bigmul_p, 2, (2*mpz_size(k->p) + 1)*64);
                mpz_sub_ui(k->bigmul_p, k->bigmul_p, 1);
                mpz_fdiv_q(k->bigmul_p, k->bigmul_p, k->p);
                mpz_mul(k->bigmul_p, k->bigmul_p, k->p);
            }
        } else if (dummy == MPFQ_SIMD_GROUPSIZE && *(int*)vp == 1) {
            /* Do nothing, this is an admitted condition */
            return;
        } else {
            return;
        }
}


/* Element allocation functions */
/* *pz::code_for_init */
void mpfq_pz_init(mpfq_pz_dst_field k, mpfq_pz_elt * x)
{
        *x = (mpfq_pz_elt) mpfq_malloc_check(mpz_size(k->p) * sizeof(mp_limb_t));
}

/* *pz::code_for_clear */
void mpfq_pz_clear(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_elt * x)
{
        mpfq_free(*x, mpz_size(k->p) * sizeof(mp_limb_t));
        *x = NULL;
}


/* Elementary assignment functions */
/* *pz::code_for_set_mpn */
void mpfq_pz_set_mpn(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mp_limb_t * x, size_t n)
{
    if (n < mpz_size(k->p)) {
        mpfq_copy(z, x, n);
        mpfq_zero(z + n, (mpz_size(k->p) - n));
    } else {
        mp_limb_t *tmp;
        tmp = (mp_limb_t *) mpfq_malloc_check((n + 1 - mpz_size(k->p)) * sizeof(mp_limb_t));
        mpn_tdiv_qr(tmp, z, 0, x, n, k->p->_mp_d, mpz_size(k->p));
        mpfq_free(tmp, (n + 1 - mpz_size(k->p)) * sizeof(mp_limb_t));
    }
}

/* *pz::code_for_set_mpz */
void mpfq_pz_set_mpz(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpz_srcptr x)
{
    if (mpz_sgn(x) < 0) {
        mpfq_pz_set_mpn(k, z, x->_mp_d, -x->_mp_size);
        mpfq_pz_neg(k, z, z);
    } else {
        mpfq_pz_set_mpn(k, z, x->_mp_d, x->_mp_size);
    }
}


/* Assignment of random values */
/* *pz::code_for_random */
void mpfq_pz_random(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, gmp_randstate_t state)
{
        mpz_t zz;
        mpz_init(zz);
        mpz_urandomb(zz, state, mpz_size(k->p) * GMP_LIMB_BITS);
        mpfq_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
}

/* *pz::code_for_random2 */
void mpfq_pz_random2(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, gmp_randstate_t state)
{
        mpz_t zz;
        mpz_init(zz);
        mpz_rrandomb(zz, state, mpz_size(k->p) * GMP_LIMB_BITS);
        mpfq_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
}


/* Arithmetic operations on elements */
/* *pz::code_for_add */
void mpfq_pz_add(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, mpfq_pz_src_elt y)
{
        mp_limb_t cy;
        cy = mpn_add_n(z, x, y, mpz_size(k->p));
        if (cy || mpn_cmp(z, k->p->_mp_d, mpz_size(k->p)) >= 0)
        mpn_sub_n(z, z, k->p->_mp_d, mpz_size(k->p));
}

/* *pz::code_for_sub */
void mpfq_pz_sub(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, mpfq_pz_src_elt y)
{
        mp_limb_t cy;
        cy = mpn_sub_n(z, x, y, mpz_size(k->p));
        if (cy)
        mpn_add_n(z, z, k->p->_mp_d, mpz_size(k->p));
}

/* *pz::code_for_neg */
void mpfq_pz_neg(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x)
{
        if (mpfq_pz_is_zero(k, x))
        mpfq_pz_set_zero(k, z);
        else
        mpn_sub_n(z, k->p->_mp_d, x, mpz_size(k->p));
}

/* *pz::code_for_mul */
void mpfq_pz_mul(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, mpfq_pz_src_elt y)
{
        mpfq_pz_elt_ur zz;
        mpfq_pz_elt_ur_init(k, &zz);
        mpfq_pz_mul_ur(k, zz, x, y);
        mpfq_pz_reduce(k, z, zz);
        mpfq_pz_elt_ur_clear(k, &zz);
}

/* *pz::code_for_sqr */
void mpfq_pz_sqr(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x)
{
        mpfq_pz_elt_ur zz;
        mpfq_pz_elt_ur_init(k, &zz);
        mpfq_pz_sqr_ur(k, zz, x);
        mpfq_pz_reduce(k, z, zz);
        mpfq_pz_elt_ur_clear(k, &zz);
}

/* *pz::code_for_is_sqr */
int mpfq_pz_is_sqr(mpfq_pz_dst_field k, mpfq_pz_src_elt x)
{
        mpz_t a, p;
        mpz_init(a);
        mpz_init(p);
        mpfq_pz_field_characteristic(k, p);
        mpfq_pz_get_mpz(k, a, x);
        int ret = mpz_jacobi(a, p);
        mpz_clear(a);
        mpz_clear(p);
        return (ret >= 0);
}

/* *pz::code_for_pow */
void mpfq_pz_pow(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, unsigned long * n, size_t nl)
{
        mpz_t p, xx, nn;
        mpz_init(p);
        mpz_init(xx);
        mpfq_pz_field_characteristic(k, p);
        mpfq_pz_get_mpz(k, xx, x);
        nn->_mp_d = n;
        nn->_mp_size = nl;
        nn->_mp_alloc = nl;
        mpz_powm(xx, xx, nn, p);
        mpfq_pz_set_mpz(k, z, xx);
        mpz_clear(p);
        mpz_clear(xx);
}

/* *Mpfq::defaults::pow::code_for_powz, pz */
void mpfq_pz_powz(mpfq_pz_dst_field k, mpfq_pz_dst_elt y, mpfq_pz_src_elt x, mpz_srcptr z)
{
        if (mpz_sgn(z) < 0) {
            mpz_t mz;
            mpz_init(mz);
            mpz_neg(mz, z);
            mpfq_pz_powz(k, y, x, mz);
            mpfq_pz_inv(k, y, y);
            mpz_clear(mz);
        } else if (mpz_sizeinbase(z, 2) > mpfq_pz_field_degree(k) * mpfq_pz_field_characteristic_bits(k)) {
            mpz_t zr;
            mpz_init(zr);
            mpz_t ppz;
            mpz_init(ppz);
            mpfq_pz_field_characteristic(k, ppz);
            mpz_pow_ui(ppz,ppz,mpfq_pz_field_degree(k));
            mpz_sub_ui(ppz,ppz,1);
            mpz_fdiv_r(zr, z, ppz);
            mpfq_pz_powz(k, y, x, zr);
            mpz_clear(ppz);
            mpz_clear(zr);
        } else {
            mpfq_pz_pow(k, y, x, z->_mp_d, mpz_size(z));
        }
}

/* *pz::code_for_add_ui */
void mpfq_pz_add_ui(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, unsigned long y)
{
        mp_limb_t cy;
        cy = mpn_add_1(z, x, mpz_size(k->p), y);
        if (cy || mpn_cmp(z, k->p->_mp_d, mpz_size(k->p)) >= 0)
        mpn_sub_n(z, z, k->p->_mp_d, mpz_size(k->p));
}

/* *pz::code_for_sub_ui */
void mpfq_pz_sub_ui(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, unsigned long y)
{
        mp_limb_t cy;
        cy = mpn_sub_1(z, x, mpz_size(k->p), y);
        if (cy)
        mpn_add_n(z, z, k->p->_mp_d, mpz_size(k->p));
}

/* *pz::code_for_mul_ui */
void mpfq_pz_mul_ui(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x, unsigned long y)
{
        mpz_t xx, p;
        mpz_init(xx);
        mpz_init(p);
        mpfq_pz_field_characteristic(k, p);
        mpfq_pz_get_mpz(k, xx, x);
        mpz_mul_ui(xx, xx, y);
        mpz_tdiv_r(xx, xx, p);
        mpfq_pz_set_mpz(k, z, xx);
        mpz_clear(xx);
        mpz_clear(p);
}

/* *pz::code_for_normalize */
void mpfq_pz_normalize(mpfq_pz_dst_field k, mpfq_pz_dst_elt z)
{
        mp_limb_t tmp;
        mpn_tdiv_qr(&tmp, z, 0, z, mpz_size(k->p), k->p->_mp_d, mpz_size(k->p));
}

/* *pz::code_for_inv */
int mpfq_pz_inv(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_elt x)
{
    mpz_t xx, p;
    mpz_init(xx);
    mpz_init(p);
    mpfq_pz_field_characteristic(k, p);
    mpfq_pz_get_mpz(k, xx, x);
    int ret = mpz_invert(xx, xx, p);
    mpfq_pz_set_mpz(k, z, xx);
    mpz_clear(xx);
    mpz_clear(p);
    if (!ret) {
        if (mpfq_pz_is_zero(k, x)) {
            mpfq_pz_set_ui(k, z, 0);
            return 0;
        } else {
            fprintf(stderr, "Not implemented\n");
            exit(1);
        }
    }
    return 1;
}


/* Operations involving unreduced elements */
/* *pz::code_for_elt_ur_init */
void mpfq_pz_elt_ur_init(mpfq_pz_dst_field k, mpfq_pz_elt_ur * x)
{
        *x = (mpfq_pz_elt_ur) mpfq_malloc_check(mpz_size(k->bigmul_p) * sizeof(mp_limb_t));
}

/* *pz::code_for_elt_ur_clear */
void mpfq_pz_elt_ur_clear(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_elt_ur * x)
{
        mpfq_free(*x, mpz_size(k->bigmul_p) * sizeof(mp_limb_t));
        *x = NULL;
}

/* *pz::code_for_elt_ur_add */
void mpfq_pz_elt_ur_add(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt_ur x, mpfq_pz_src_elt_ur y)
{
        mpn_add_n(z, x, y, mpz_size(k->bigmul_p));
}

/* *pz::code_for_elt_ur_neg */
void mpfq_pz_elt_ur_neg(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt_ur x)
{
#if GMP_VERSION_ATLEAST(5, 0, 0)
        mpn_neg(z, x, mpz_size(k->bigmul_p));
#else
      mp_size_t i;
      mp_size_t n = mpz_size(k->bigmul_p);
    
      for (i = 0; i < n; i++)
        z[i] = ~x[i];
      mpn_add_1 (z, z, n, 1);
#endif
}

/* *pz::code_for_elt_ur_sub */
void mpfq_pz_elt_ur_sub(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt_ur x, mpfq_pz_src_elt_ur y)
{
        mpn_sub_n(z, x, y, mpz_size(k->bigmul_p));
}

/* *pz::code_for_mul_ur */
void mpfq_pz_mul_ur(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt x, mpfq_pz_src_elt y)
{
        mpn_mul_n(z, x, y, mpz_size(k->p));
        z[mpz_size(k->bigmul_p) - 1] = 0;
}

/* *pz::code_for_sqr_ur */
void mpfq_pz_sqr_ur(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_elt x)
{
        mpn_sqr(z, x, mpz_size(k->p));
        z[mpz_size(k->bigmul_p) - 1] = 0;
}

/* *pz::code_for_reduce */
void mpfq_pz_reduce(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_dst_elt_ur x)
{
    //    int neg = 0;
    if (x[mpz_size(k->bigmul_p) - 1] >> (GMP_LIMB_BITS - 1)) {	// negative input
    //        neg = 1;
    //        mpfq_pz_elt_ur_neg(k, x, x);
    mpn_add_n(x, x, k->bigmul_p->_mp_d, mpz_size(k->bigmul_p));
    }
    mp_limb_t *tmp;
    tmp = (mp_limb_t *) mpfq_malloc_check((mpz_size(k->bigmul_p) + 1) * sizeof(mp_limb_t));
    mpn_tdiv_qr(tmp, z, 0, x, mpz_size(k->bigmul_p), k->p->_mp_d, mpz_size(k->p));
    //    if (neg)
    //        mpfq_pz_neg(k, z, z);
    mpfq_free(tmp, (mpz_size(k->bigmul_p) + 1) * sizeof(mp_limb_t));
}


/* Comparison functions */

/* Input/output functions */
/* *pz::code_for_asprint */
int mpfq_pz_asprint(mpfq_pz_dst_field k, char * * pstr, mpfq_pz_src_elt x)
{
        // deal with 0
        if (mpfq_pz_is_zero(k, x)) {
        char * s = *pstr = (char *) malloc(2);
        s[0] = '0';
        s[1] = '\0';
        return 1;
        }
        // allocate enough room for base 2 conversion.
        *pstr = (char *) malloc(((mpz_size(k->p)) * 64 + 1));
        mp_limb_t *tmp;
        tmp = (mp_limb_t *) mpfq_malloc_check((mpz_size(k->p)) * sizeof(mp_limb_t));
        int tl = mpz_size(k->p) - 1;
        while (x[tl] == 0)		// x is non-zero, no need to test tl>0
        tl--;
        for (int i = 0; i <= tl; ++i)
        tmp[i] = x[i];
        int n = mpn_get_str((unsigned char *) (*pstr), k->io_base, tmp, tl + 1);
        mpfq_free(tmp, (mpz_size(k->p)) * sizeof(mp_limb_t));
        for (int i = 0; i < n; ++i)
        (*pstr)[i] += '0';
        (*pstr)[n] = '\0';
        // remove leading 0s
        int shift = 0;
        while (((*pstr)[shift] == '0') && ((*pstr)[shift + 1] != '\0'))
        shift++;
        if (shift) {
            memmove(*pstr, (*pstr) + shift, n + 1 - shift);
            n -= shift;
        }
        return n;
}

/* *pz::code_for_fprint */
int mpfq_pz_fprint(mpfq_pz_dst_field k, FILE * file, mpfq_pz_src_elt x)
{
        char *str;
        int rc;
        mpfq_pz_asprint(k, &str, x);
        rc = fprintf(file, "%s", str);
        free(str);
        return rc;
}

/* *pz::code_for_sscan */
int mpfq_pz_sscan(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, const char * str)
{
        mpz_t zz;
        mpz_init(zz);
        int nread;
        if (gmp_sscanf(str, "%Zd%n", zz, &nread) != 1) {
        mpz_clear(zz);
        return 0;
        }
        mpfq_pz_set_mpz(k, z, zz);
        mpz_clear(zz);
        return nread;
}

/* *pz::code_for_fscan */
int mpfq_pz_fscan(mpfq_pz_dst_field k, FILE * file, mpfq_pz_dst_elt x)
{
        mpz_t zz;
        mpz_init(zz);
        int nread;
        if (gmp_fscanf(file, "%Zd%n", zz, &nread) != 1) {
        mpz_clear(zz);
        return 0;
        }
        mpfq_pz_set_mpz(k, x, zz);
        mpz_clear(zz);
        return nread;
}


/* Vector functions */
/* *pz::code_for_vec_init */
void mpfq_pz_vec_init(mpfq_pz_dst_field k, mpfq_pz_vec * v, unsigned long n)
{
        *v = (mpfq_pz_vec) mpfq_malloc_check((mpz_size(k->p)) * n * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_reinit */
void mpfq_pz_vec_reinit(mpfq_pz_dst_field k, mpfq_pz_vec * v, unsigned long n MAYBE_UNUSED, unsigned long m)
{
        *v = (mpfq_pz_vec) mpfq_realloc_check(*v, (mpz_size(k->p)) * n * sizeof(mp_limb_t), (mpz_size(k->p)) * m * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_clear */
void mpfq_pz_vec_clear(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_vec * v, unsigned long n MAYBE_UNUSED)
{
        mpfq_free(*v, (mpz_size(k->p)) * n * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_set */
void mpfq_pz_vec_set(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec v, unsigned long n)
{
        if (v != w)
        memmove(w, v, n * mpz_size(k->p) * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_set_zero */
void mpfq_pz_vec_set_zero(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, unsigned long n)
{
        if (n) mpfq_zero(w, n * mpz_size(k->p));
}

/* *pz::code_for_vec_setcoeff */
void mpfq_pz_vec_setcoeff(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_elt x, unsigned long i)
{
        mpfq_pz_set(k, w + i * mpz_size(k->p), x);
}

/* *pz::code_for_vec_setcoeff_ui */
void mpfq_pz_vec_setcoeff_ui(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, unsigned long x0, unsigned long i)
{
        mpfq_pz_set_ui(k, w + i * mpz_size(k->p), x0);
}

/* *pz::code_for_vec_getcoeff */
void mpfq_pz_vec_getcoeff(mpfq_pz_dst_field k, mpfq_pz_dst_elt z, mpfq_pz_src_vec v, unsigned long i)
{
        mpfq_pz_set(k, z, v + i * mpz_size(k->p));
}

/* *pz::code_for_vec_add */
void mpfq_pz_vec_add(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec u, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_add(k, w + i * mpz_size(k->p), u + i * mpz_size(k->p), v + i * mpz_size(k->p));
}

/* *pz::code_for_vec_neg */
void mpfq_pz_vec_neg(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_neg(k, w + i * mpz_size(k->p), v + i * mpz_size(k->p));
}

/* *pz::code_for_vec_rev */
void mpfq_pz_vec_rev(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec v, unsigned long n)
{
        unsigned int nn = n >> 1;
        mpfq_pz_elt tmp;
        mpfq_pz_init(k, &tmp);
        for (unsigned int i = 0; i < n - 1 - i; ++i) {
        mpfq_pz_set(k, tmp, v + i * mpz_size(k->p));
        mpfq_pz_set(k, w + i * mpz_size(k->p), v + (n - 1 - i) * mpz_size(k->p));
        mpfq_pz_set(k, w + (n - 1 - i) * mpz_size(k->p), tmp);
        }
        if (n & 1)
        mpfq_pz_set(k, w + nn * mpz_size(k->p), v + nn * mpz_size(k->p));
        mpfq_pz_clear(k, &tmp);
}

/* *pz::code_for_vec_sub */
void mpfq_pz_vec_sub(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec u, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_sub(k, w + i * mpz_size(k->p), u + i * mpz_size(k->p), v + i * mpz_size(k->p));
}

/* *pz::code_for_vec_conv */
void mpfq_pz_vec_conv(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_src_vec u, unsigned long n, mpfq_pz_src_vec v, unsigned long m)
{
        mpfq_pz_vec_ur tmp;
        mpfq_pz_vec_ur_init(k, &tmp, m + n - 1);
        mpfq_pz_vec_conv_ur(k, tmp, u, n, v, m);
        mpfq_pz_vec_reduce(k, w, tmp, m + n - 1);
        mpfq_pz_vec_ur_clear(k, &tmp, m + n - 1);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_random, pz */
void mpfq_pz_vec_random(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, unsigned long n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_dst_elt x = mpfq_pz_vec_coeff_ptr(k, w, i);
        mpfq_pz_random(k, x, state);
        }
}

/* *Mpfq::defaults::vec::generic::code_for_vec_random2, pz */
void mpfq_pz_vec_random2(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, unsigned long n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_dst_elt x = mpfq_pz_vec_coeff_ptr(k, w, i);
        mpfq_pz_random2(k, x, state);
        }
}

/* *Mpfq::defaults::vec::generic::code_for_vec_cmp, pz */
long mpfq_pz_vec_cmp(mpfq_pz_dst_field k, mpfq_pz_src_vec u, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, u, i);
            mpfq_pz_src_elt y = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        int ret = mpfq_pz_cmp(k, x, y);
        if (ret != 0)
            return ret;
        }
        return 0;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_is_zero, pz */
long mpfq_pz_vec_is_zero(mpfq_pz_dst_field k, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_pz_is_zero(k, x))
            return 0;
        }
        return 1;
}

/* *pz::code_for_vec_asprint */
long mpfq_pz_vec_asprint(mpfq_pz_dst_field k, char * * pstr, mpfq_pz_src_vec v, unsigned long n)
{
        if (n == 0) {
        *pstr = (char *) malloc(4);
        sprintf(*pstr, "[ ]");
        return strlen(*pstr);
        }
        int alloc = 100;
        int len = 0;
        *pstr = (char *) malloc(alloc);
        if (!*pstr) abort();
        (*pstr)[len++] = '[';
        (*pstr)[len++] = ' ';
        for (unsigned int i = 0; i < n; ++i) {
        if (i) {
            (*pstr)[len++] = ',';
            (*pstr)[len++] = ' ';
        }
        char *tmp;
        mpfq_pz_asprint(k, &tmp, v + i * mpz_size(k->p));
        int ltmp = strlen(tmp);
        if (len + ltmp + 4 > alloc) {
            alloc = len + ltmp + 100;
            *pstr = (char *) realloc(*pstr, alloc);
        }
            strncpy(*pstr+len, tmp, alloc-len);
        len += ltmp;
        free(tmp);
        }
        (*pstr)[len++] = ' ';
        (*pstr)[len++] = ']';
        (*pstr)[len] = '\0';
        return len;
}

/* *pz::code_for_vec_fprint */
long mpfq_pz_vec_fprint(mpfq_pz_dst_field k, FILE * file, mpfq_pz_src_vec v, unsigned long n)
{
        char *str;
        int rc;
        mpfq_pz_vec_asprint(k, &str, v, n);
        rc = fprintf(file, "%s", str);
        free(str);
        return rc;
}

/* *Mpfq::defaults::vec::io::code_for_vec_print, pz */
long mpfq_pz_vec_print(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_src_vec w, unsigned long n)
{
    return mpfq_pz_vec_fprint(K,stdout,w,n);
}

/* *Mpfq::defaults::vec::io::code_for_vec_sscan, pz */
long mpfq_pz_vec_sscan(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_vec * w, unsigned long * n, const char * str)
{
    // start with a clean vector
    unsigned int nn;
    int len = 0;
    mpfq_pz_vec_reinit(K, w, *n, 0);
    *n = nn = 0;
    while (isspace((int)(unsigned char)str[len]))
        len++;
    if (str[len] != '[')
        return 0;
    len++;
    while (isspace((int)(unsigned char)str[len]))
        len++;
    if (str[len] == ']') {
        len++;
        return len;
    }
    unsigned int i = 0;
    for (;;) {
        if (nn < i+1) {
            mpfq_pz_vec_reinit(K, w, nn, i+1);
            *n = nn = i+1;
        }
        int ret = mpfq_pz_sscan(K, mpfq_pz_vec_coeff_ptr(K, *w, i), str + len);
        if (!ret) {
            *n = 0; /* invalidate data ! */
            return 0;
        }
        i++;
        len += ret;
        while (isspace((int)(unsigned char)str[len]))
            len++;
        if (str[len] == ']') {
            len++;
            break;
        }
        if (str[len] != ',') {
            *n = 0; /* invalidate data ! */
            return 0;
        }
        len++;
        while (isspace((int)(unsigned char)str[len]))
            len++;
    }
    return len;
}

/* *Mpfq::defaults::vec::io::code_for_vec_fscan, pz */
long mpfq_pz_vec_fscan(mpfq_pz_dst_field K MAYBE_UNUSED, FILE * file, mpfq_pz_vec * w, unsigned long * n)
{
    char *tmp;
    int c;
    int allocated, len=0;
    allocated=100;
    tmp = (char *)mpfq_malloc_check(allocated);
    int nest = 0, mnest = 0;
    for(;;) {
        c = fgetc(file);
        if (c==EOF) {
            free(tmp);
            return 0;
        }
        if (c == '[') {
            nest++, mnest++;
        }
        if (len==allocated) {
            allocated = len + 10 + allocated / 4;
            tmp = (char*)realloc(tmp, allocated);
        }
        tmp[len]=c;
        len++;
        if (c == ']') {
            nest--, mnest++;
        }
        if (mnest && nest == 0)
            break;
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated);
    }
    tmp[len]='\0';
    int ret=mpfq_pz_vec_sscan(K,w,n,tmp);
    free(tmp);
    return ret;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_hamming_weight, pz */
long mpfq_pz_vec_hamming_weight(mpfq_pz_dst_field k, mpfq_pz_src_vec v, unsigned long n)
{
        int w = 0;
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        w += !mpfq_pz_is_zero(k, x);
        }
        return w;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_find_first_set, pz */
long mpfq_pz_vec_find_first_set(mpfq_pz_dst_field k, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_pz_is_zero(k, x))
            return i;
        }
        return -1;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_simd_hamming_weight, pz */
long mpfq_pz_vec_simd_hamming_weight(mpfq_pz_dst_field k, mpfq_pz_src_vec v, unsigned long n)
{
        int w = 0;
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        w += mpfq_pz_simd_hamming_weight(k, x);
        }
        return w;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_simd_find_first_set, pz */
long mpfq_pz_vec_simd_find_first_set(mpfq_pz_dst_field k, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_pz_src_elt x = mpfq_pz_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_pz_is_zero(k, x))
            return i * mpfq_pz_simd_groupsize(k) + mpfq_pz_simd_hamming_weight(k, x);
        }
        return -1;
}

/* *pz::code_for_vec_ur_init */
void mpfq_pz_vec_ur_init(mpfq_pz_dst_field k, mpfq_pz_vec_ur * v, unsigned long n)
{
        *v = (mpfq_pz_vec) mpfq_malloc_check((mpz_size(k->bigmul_p)) * n * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_set_zero */
void mpfq_pz_vec_ur_set_zero(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, unsigned long n)
{
        mpfq_zero(w, n * mpz_size(k->bigmul_p));
}

/* *pz::code_for_vec_ur_set_vec */
void mpfq_pz_vec_ur_set_vec(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_elt_ur_set_elt(k, w + i * mpz_size(k->bigmul_p), v + i * mpz_size(k->p));
}

/* *pz::code_for_vec_ur_reinit */
void mpfq_pz_vec_ur_reinit(mpfq_pz_dst_field k, mpfq_pz_vec_ur * v, unsigned long n MAYBE_UNUSED, unsigned long m)
{
        *v = (mpfq_pz_vec) mpfq_realloc_check(*v, (mpz_size(k->bigmul_p)) * n * sizeof(mp_limb_t), (mpz_size(k->bigmul_p)) * m * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_clear */
void mpfq_pz_vec_ur_clear(mpfq_pz_dst_field k MAYBE_UNUSED, mpfq_pz_vec_ur * v, unsigned long n MAYBE_UNUSED)
{
        mpfq_free(*v, (mpz_size(k->bigmul_p)) * n * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_set */
void mpfq_pz_vec_ur_set(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec_ur v, unsigned long n)
{
        if (v != w)
        memmove(w, v, n * mpz_size(k->bigmul_p) * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_ur_setcoeff */
void mpfq_pz_vec_ur_setcoeff(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_elt_ur x, unsigned long i)
{
        mpfq_pz_elt_ur_set(k, w + i * mpz_size(k->bigmul_p), x);
}

/* *pz::code_for_vec_ur_getcoeff */
void mpfq_pz_vec_ur_getcoeff(mpfq_pz_dst_field k, mpfq_pz_dst_elt_ur z, mpfq_pz_src_vec_ur v, unsigned long i)
{
        mpfq_pz_elt_ur_set(k, z, v + i * mpz_size(k->bigmul_p));
}

/* *pz::code_for_vec_ur_add */
void mpfq_pz_vec_ur_add(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec_ur u, mpfq_pz_src_vec_ur v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_elt_ur_add(k, w + i * mpz_size(k->bigmul_p), u + i * mpz_size(k->bigmul_p),
                      v + i * mpz_size(k->bigmul_p));
}

/* *pz::code_for_vec_ur_sub */
void mpfq_pz_vec_ur_sub(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec_ur u, mpfq_pz_src_vec_ur v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_elt_ur_sub(k, w + i * mpz_size(k->bigmul_p), u + i * mpz_size(k->bigmul_p),
                      v + i * mpz_size(k->bigmul_p));
}

/* *pz::code_for_vec_ur_neg */
void mpfq_pz_vec_ur_neg(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec_ur v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_elt_ur_neg(k, w + i * mpz_size(k->bigmul_p), v + i * mpz_size(k->bigmul_p));
}

/* *pz::code_for_vec_ur_rev */
void mpfq_pz_vec_ur_rev(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec_ur v, unsigned long n)
{
        unsigned int nn = n >> 1;
        mpfq_pz_elt_ur tmp;
        mpfq_pz_elt_ur_init(k, &tmp);
        for (unsigned int i = 0; i < nn; ++i) {
        mpfq_pz_elt_ur_set(k, tmp, v + i * mpz_size(k->bigmul_p));
        mpfq_pz_elt_ur_set(k, w + i * mpz_size(k->bigmul_p), v + (n - 1 - i) * mpz_size(k->bigmul_p));
        mpfq_pz_elt_ur_set(k, w + (n - 1 - i) * mpz_size(k->bigmul_p), tmp);
        }
        if (n & 1)
        mpfq_pz_elt_ur_set(k, w + nn * mpz_size(k->bigmul_p), v + nn * mpz_size(k->bigmul_p));
        mpfq_pz_elt_ur_clear(k, &tmp);
}

static void mpfq_pz_vec_conv_ur_ks(mpfq_pz_field, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, unsigned int, mpfq_pz_src_vec, unsigned int);
/* *pz::code_for_vec_conv_ur */
/* Triggered by: vec_conv_ur */
static void mpfq_pz_vec_conv_ur_ks(mpfq_pz_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec u, unsigned int n, mpfq_pz_src_vec v, unsigned int m)
{
        // compute base as a power 2^GMP_NUMB_BITS
        // This is the least number of words that can accomodate
        //     log_2( (p-1)^2 * min(n,m) )
        mpz_t p;
        mpz_init(p);
        mpfq_pz_field_characteristic(k, p);
        mpz_sub_ui(p, p, 1);
        mpz_mul(p, p, p);
        mpz_mul_ui(p, p, MIN(m, n));
    
        long nbits = mpz_sizeinbase(p, 2);
        unsigned long nwords = 1 + ((nbits - 1) / GMP_NUMB_BITS);
        nbits = GMP_NUMB_BITS * nwords;
        mpz_clear(p);
        assert(mpz_size(k->bigmul_p) >= nwords);
    
        // Create big integers
        mp_limb_t *U, *V;
        U = (mp_limb_t *) mpfq_malloc_check(n * nwords * sizeof(mp_limb_t));
        V = (mp_limb_t *) mpfq_malloc_check(m * nwords * sizeof(mp_limb_t));
        mpfq_zero(U, n * nwords);
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_get_mpn(k, U + i * nwords, u + i * mpz_size(k->p));
        mpfq_zero(V, m * nwords);
        for (unsigned int i = 0; i < m; ++i)
        mpfq_pz_get_mpn(k, V + i * nwords, v + i * mpz_size(k->p));
        // Mul !
        mp_limb_t *W;
        W = (mp_limb_t *) mpfq_malloc_check((n + m) * nwords * sizeof(mp_limb_t));
        if (n>=m)
            mpn_mul(W, U, n * nwords, V, m * nwords);
        else
            mpn_mul(W, V, m * nwords, U, n * nwords);
        // Put coeffs in w
        mpfq_zero(w, (n + m - 1) * mpz_size(k->bigmul_p));
        for (unsigned int i = 0; i < n + m - 1; ++i)
        mpfq_copy(w + i * mpz_size(k->bigmul_p), W + i * nwords, nwords);
        mpfq_free(U, n * nwords * sizeof(mp_limb_t));
        mpfq_free(V, m * nwords * sizeof(mp_limb_t));
        mpfq_free(W, (m+n) * nwords * sizeof(mp_limb_t));
}

/* *pz::code_for_vec_conv_ur */
void mpfq_pz_vec_conv_ur(mpfq_pz_dst_field k, mpfq_pz_dst_vec_ur w, mpfq_pz_src_vec u, unsigned long n, mpfq_pz_src_vec v, unsigned long m)
{
        mpfq_pz_vec_conv_ur_ks(k, w, u, n, v, m);
}

/* *pz::code_for_vec_reduce */
void mpfq_pz_vec_reduce(mpfq_pz_dst_field k, mpfq_pz_dst_vec w, mpfq_pz_dst_vec_ur v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i)
        mpfq_pz_reduce(k, w + i * mpz_size(k->p), v + i * mpz_size(k->bigmul_p));
}


/* Polynomial functions */
/* *Mpfq::defaults::poly::code_for_poly_setmonic, pz */
void mpfq_pz_poly_setmonic(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_poly q, mpfq_pz_src_poly p)
{
    long degp = mpfq_pz_poly_deg(K, p);
    if (degp == -1) {
        q->size = 0;
        return;
    }
    if (degp == 0) {
        mpfq_pz_elt aux;
        mpfq_pz_init(K, &aux);
        mpfq_pz_set_ui(K, aux, 1);
        mpfq_pz_poly_setcoeff(K, q, aux, 0);
        mpfq_pz_clear(K, &aux);
        q->size = 1;
        return;
    }
    mpfq_pz_elt lc;	/* spurious uninit warning sometimes */
    mpfq_pz_init(K, &lc);
    mpfq_pz_poly_getcoeff(K, lc, p, degp);
    mpfq_pz_inv(K, lc, lc);
    mpfq_pz_poly_setcoeff_ui(K, q, 1, degp);
    mpfq_pz_vec_scal_mul(K, q->c, p->c, lc, degp);
    q->size = degp+1;
    mpfq_pz_clear(K, &lc);
}

/* *Mpfq::defaults::poly::code_for_poly_divmod, pz */
int mpfq_pz_poly_divmod(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_poly q, mpfq_pz_dst_poly r, mpfq_pz_src_poly a, mpfq_pz_src_poly b)
{
    if (b->size == 0) {
        return 0;
    }
    if (a->size == 0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    int dega = mpfq_pz_poly_deg(K, a);
    if (dega<0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    // Compute deg b and inverse of leading coef
    int degb = mpfq_pz_poly_deg(K, b);
    if (degb<0) {
        return 0;
    }
    if (degb > dega) {
        q->size=0;
        mpfq_pz_poly_set(K, r, a);
        return 1;
    }
    int bmonic;
    mpfq_pz_elt ilb;
    mpfq_pz_init(K, &ilb);
    mpfq_pz_elt temp;
    mpfq_pz_init(K, &temp);
    mpfq_pz_set_zero(K, temp); // silence spurious gcc11 warning :-(
    mpfq_pz_poly_getcoeff(K, temp, b, degb);
    if (mpfq_pz_cmp_ui(K, temp, 1) == 0) {
        mpfq_pz_set_ui(K, ilb, 1);
        bmonic = 1;
    } else {
        mpfq_pz_inv(K, ilb, temp);
        bmonic = 0;
    }
    
    mpfq_pz_poly qq, rr;
    mpfq_pz_poly_init(K, qq, dega-degb+1);
    mpfq_pz_poly_init(K, rr, dega);
    
    mpfq_pz_poly_set(K, rr, a);
    mpfq_pz_elt aux, aux2;
    
    mpfq_pz_init(K, &aux);
    mpfq_pz_init(K, &aux2);
    
    int i;
    int j;
    for (i = dega; i >= (int)degb; --i) {
        mpfq_pz_poly_getcoeff(K, aux, rr, i);
        if (!bmonic) 
            mpfq_pz_mul(K, aux, aux, ilb);
        mpfq_pz_poly_setcoeff(K, qq, aux, i-degb);
        for (j = i-1; j >= (int)(i - degb); --j) {
            mpfq_pz_poly_getcoeff(K, temp, b, j-i+degb);
            mpfq_pz_mul(K, aux2, aux, temp);
            mpfq_pz_poly_getcoeff(K, temp, rr, j);
    
            mpfq_pz_sub(K, temp, temp, aux2);
            mpfq_pz_poly_setcoeff(K, rr, temp, j);
        }
    }    
    
    rr->size = degb;
    int degr = mpfq_pz_poly_deg(K, rr);
    rr->size = degr+1;
    
    if (q != NULL) 
        mpfq_pz_poly_set(K, q, qq);
    if (r != NULL)
        mpfq_pz_poly_set(K, r, rr);
    mpfq_pz_clear(K, &temp);
    mpfq_pz_clear(K, &ilb);
    mpfq_pz_clear(K, &aux);
    mpfq_pz_clear(K, &aux2);
    mpfq_pz_poly_clear(K, rr);
    mpfq_pz_poly_clear(K, qq);
    return 1;
}

static void mpfq_pz_poly_preinv(mpfq_pz_dst_field, mpfq_pz_dst_poly, mpfq_pz_src_poly, unsigned int);
/* *Mpfq::defaults::poly::code_for_poly_precomp_mod, pz */
/* Triggered by: poly_precomp_mod */
static void mpfq_pz_poly_preinv(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_poly q, mpfq_pz_src_poly p, unsigned int n)
{
    // Compute the inverse of p(x) modulo x^n
    // Newton iteration: x_{n+1} = x_n + x_n(1 - a*x_n)
    // Requires p(0) = 1
    // Assume p != q (no alias)
    assert (p != q);
    int m;
    if (n <= 2) {
        mpfq_pz_poly_setcoeff_ui(K, q, 1, 0);
        q->size = 1;
        m = 1;
        if (n == 1)
            return;
    } else {
        // n >= 3: recursive call at prec m = ceil(n/2)
        m = 1 + ((n-1)/2);
        mpfq_pz_poly_preinv(K, q, p, m);
    }
    mpfq_pz_elt temp;	/* spurious uninit warning sometimes */
    mpfq_pz_init(K, &temp);
    mpfq_pz_set_zero(K, temp); // silence spurious gcc11 warning :-(
    mpfq_pz_poly_getcoeff(K, temp, p, 0);//Should be in the assert
    assert( mpfq_pz_cmp_ui(K, temp, 1) == 0);
    // enlarge q if necessary
    if (q->alloc < n) {
        mpfq_pz_vec_reinit(K, &(q->c), q->alloc, n);
        q->alloc = n;
    }
    // refine value
    mpfq_pz_vec tmp;
    mpfq_pz_vec_init(K, &tmp, m+n-1);
    
    mpfq_pz_vec_conv(K, tmp, p->c, MIN(n, p->size), q->c, m);
    int nn = MIN(n, MIN(n, p->size) + m -1);
    mpfq_pz_vec_neg(K, tmp, tmp, nn);
    mpfq_pz_vec_getcoeff(K, temp, tmp, 0);
    mpfq_pz_add_ui(K, temp, temp, 1);
    mpfq_pz_vec_setcoeff(K, tmp, temp, 0);
    mpfq_pz_vec_conv(K, tmp, q->c, m, tmp, nn);
    mpfq_pz_vec_set(K, mpfq_pz_vec_subvec(K, q->c, m), mpfq_pz_vec_subvec(K, tmp, m), n-m);
    q->size = n;
    
    mpfq_pz_clear(K, &temp);
    mpfq_pz_vec_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::poly::code_for_poly_precomp_mod, pz */
void mpfq_pz_poly_precomp_mod(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_poly q, mpfq_pz_src_poly p)
{
    assert(p != q);
    int N = mpfq_pz_poly_deg(K, p);
    if (N < 0) {
        mpfq_pz_poly_set(K, q, p);
        return;
    }
    mpfq_pz_poly rp;
    mpfq_pz_poly_init(K, rp, N+1);
    mpfq_pz_vec_rev(K, rp->c, p->c, N+1);
    rp->size = N+1;
    mpfq_pz_poly_preinv(K, q, rp, N);
    mpfq_pz_poly_clear(K, rp);
}

/* *Mpfq::defaults::poly::code_for_poly_mod_pre, pz */
void mpfq_pz_poly_mod_pre(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_poly r, mpfq_pz_src_poly q, mpfq_pz_src_poly p, mpfq_pz_src_poly irp)
{
    int N = mpfq_pz_poly_deg(K, p);
    int degq = mpfq_pz_poly_deg(K, q);
    if (degq < N) {
        mpfq_pz_poly_set(K, r, q);
        return;
    }
    int m = degq - N;
    assert (degq <= 2*N-2);
    mpfq_pz_poly revq;
    mpfq_pz_poly_init(K, revq, MAX(degq+1, m+1));
    mpfq_pz_vec_rev(K, revq->c, q->c, degq+1);
    revq->size = q->size;
    mpfq_pz_poly_mul(K, revq, revq, irp);
    mpfq_pz_vec_rev(K, revq->c, revq->c, m+1);
    revq->size = m+1;
    
    mpfq_pz_poly_mul(K, revq, revq, p);
    mpfq_pz_poly_sub(K, r, q, revq);
    r->size = mpfq_pz_poly_deg(K, r)+1;
    mpfq_pz_poly_clear(K, revq);
}


/* Functions related to SIMD operation */
/* *simd_pz::code_for_add_dotprod */
void mpfq_pz_add_dotprod(mpfq_pz_dst_field K MAYBE_UNUSED, mpfq_pz_dst_vec xw, mpfq_pz_src_vec xu1, mpfq_pz_src_vec xu0, unsigned int n)
{
        mpfq_pz_elt_ur s,t;
        mpfq_pz_elt_ur_init(K, &s);
        mpfq_pz_elt_ur_init(K, &t);
        mpfq_pz_elt_ur_set_elt(K, s, mpfq_pz_vec_coeff_ptr(K, xw, 0));
        for(unsigned int i = 0 ; i < n ; i++) {
            mpfq_pz_mul_ur(K, t, mpfq_pz_vec_coeff_ptr_const(K, xu0, i), mpfq_pz_vec_coeff_ptr_const(K, xu1, i));
            mpfq_pz_elt_ur_add(K, s, s, t);
        }
        mpfq_pz_reduce(K, mpfq_pz_vec_coeff_ptr(K, xw, 0), s);
        mpfq_pz_elt_ur_clear(K, &s);
        mpfq_pz_elt_ur_clear(K, &t);
}


/* Member templates related to SIMD operation */

/* Object-oriented interface */
static void mpfq_pz_wrapper_oo_field_clear(mpfq_vbase_ptr);
static void mpfq_pz_wrapper_oo_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_pz_oo_field_clear(vbase);
}

static void mpfq_pz_wrapper_oo_field_init(mpfq_vbase_ptr);
static void mpfq_pz_wrapper_oo_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_pz_oo_field_init(vbase);
}

static void mpfq_pz_wrapper_add_dotprod(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned int);
static void mpfq_pz_wrapper_add_dotprod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec xw MAYBE_UNUSED, mpfq_pz_src_vec xu1 MAYBE_UNUSED, mpfq_pz_src_vec xu0 MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_pz_add_dotprod(vbase->obj, xw, xu1, xu0, n);
}

static void mpfq_pz_wrapper_simd_set_ui_all(mpfq_vbase_ptr, mpfq_pz_dst_elt, unsigned long);
static void mpfq_pz_wrapper_simd_set_ui_all(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt p MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_pz_simd_set_ui_all(vbase->obj, p, v);
}

static void mpfq_pz_wrapper_simd_add_ui_at(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, int, unsigned long);
static void mpfq_pz_wrapper_simd_add_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt p MAYBE_UNUSED, mpfq_pz_src_elt p0 MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_pz_simd_add_ui_at(vbase->obj, p, p0, k, v);
}

static void mpfq_pz_wrapper_simd_set_ui_at(mpfq_vbase_ptr, mpfq_pz_dst_elt, int, unsigned long);
static void mpfq_pz_wrapper_simd_set_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_pz_simd_set_ui_at(vbase->obj, p, k, v);
}

static unsigned long mpfq_pz_wrapper_simd_get_ui_at(mpfq_vbase_ptr, mpfq_pz_src_elt, int);
static unsigned long mpfq_pz_wrapper_simd_get_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt p MAYBE_UNUSED, int k MAYBE_UNUSED)
{
    return mpfq_pz_simd_get_ui_at(vbase->obj, p, k);
}

static int mpfq_pz_wrapper_simd_find_first_set(mpfq_vbase_ptr, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_simd_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt p MAYBE_UNUSED)
{
    return mpfq_pz_simd_find_first_set(vbase->obj, p);
}

static int mpfq_pz_wrapper_simd_hamming_weight(mpfq_vbase_ptr, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_simd_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt p MAYBE_UNUSED)
{
    return mpfq_pz_simd_hamming_weight(vbase->obj, p);
}

static int mpfq_pz_wrapper_simd_groupsize(mpfq_vbase_ptr);
static int mpfq_pz_wrapper_simd_groupsize(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_simd_groupsize(vbase->obj);
}

static int mpfq_pz_wrapper_poly_scan(mpfq_vbase_ptr, mpfq_pz_dst_poly);
static int mpfq_pz_wrapper_poly_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_scan(vbase->obj, w);
}

static int mpfq_pz_wrapper_poly_fscan(mpfq_vbase_ptr, FILE *, mpfq_pz_dst_poly);
static int mpfq_pz_wrapper_poly_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_fscan(vbase->obj, file, w);
}

static int mpfq_pz_wrapper_poly_sscan(mpfq_vbase_ptr, mpfq_pz_dst_poly, const char *);
static int mpfq_pz_wrapper_poly_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_pz_poly_sscan(vbase->obj, w, str);
}

static int mpfq_pz_wrapper_poly_print(mpfq_vbase_ptr, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_print(vbase->obj, w);
}

static int mpfq_pz_wrapper_poly_fprint(mpfq_vbase_ptr, FILE *, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_src_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_fprint(vbase->obj, file, w);
}

static int mpfq_pz_wrapper_poly_asprint(mpfq_vbase_ptr, char * *, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_pz_src_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_asprint(vbase->obj, pstr, w);
}

static int mpfq_pz_wrapper_poly_cmp(mpfq_vbase_ptr, mpfq_pz_src_poly, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, mpfq_pz_src_poly v MAYBE_UNUSED)
{
    return mpfq_pz_poly_cmp(vbase->obj, u, v);
}

static void mpfq_pz_wrapper_poly_random2(mpfq_vbase_ptr, mpfq_pz_dst_poly, unsigned int, gmp_randstate_t);
static void mpfq_pz_wrapper_poly_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_poly_random2(vbase->obj, w, n, state);
}

static void mpfq_pz_wrapper_poly_random(mpfq_vbase_ptr, mpfq_pz_dst_poly, unsigned int, gmp_randstate_t);
static void mpfq_pz_wrapper_poly_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_poly_random(vbase->obj, w, n, state);
}

static void mpfq_pz_wrapper_poly_xgcd(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_xgcd(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly g MAYBE_UNUSED, mpfq_pz_dst_poly u0 MAYBE_UNUSED, mpfq_pz_dst_poly v0 MAYBE_UNUSED, mpfq_pz_src_poly a0 MAYBE_UNUSED, mpfq_pz_src_poly b0 MAYBE_UNUSED)
{
    mpfq_pz_poly_xgcd(vbase->obj, g, u0, v0, a0, b0);
}

static void mpfq_pz_wrapper_poly_gcd(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_gcd(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly g MAYBE_UNUSED, mpfq_pz_src_poly a0 MAYBE_UNUSED, mpfq_pz_src_poly b0 MAYBE_UNUSED)
{
    mpfq_pz_poly_gcd(vbase->obj, g, a0, b0);
}

static void mpfq_pz_wrapper_poly_mod_pre(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_mod_pre(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly r MAYBE_UNUSED, mpfq_pz_src_poly q MAYBE_UNUSED, mpfq_pz_src_poly p MAYBE_UNUSED, mpfq_pz_src_poly irp MAYBE_UNUSED)
{
    mpfq_pz_poly_mod_pre(vbase->obj, r, q, p, irp);
}

static void mpfq_pz_wrapper_poly_precomp_mod(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_precomp_mod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly q MAYBE_UNUSED, mpfq_pz_src_poly p MAYBE_UNUSED)
{
    mpfq_pz_poly_precomp_mod(vbase->obj, q, p);
}

static int mpfq_pz_wrapper_poly_divmod(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_divmod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly q MAYBE_UNUSED, mpfq_pz_dst_poly r MAYBE_UNUSED, mpfq_pz_src_poly a MAYBE_UNUSED, mpfq_pz_src_poly b MAYBE_UNUSED)
{
    return mpfq_pz_poly_divmod(vbase->obj, q, r, a, b);
}

static void mpfq_pz_wrapper_poly_mul(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, mpfq_pz_src_poly v MAYBE_UNUSED)
{
    mpfq_pz_poly_mul(vbase->obj, w, u, v);
}

static void mpfq_pz_wrapper_poly_scal_mul(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_poly_scal_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_poly_scal_mul(vbase->obj, w, u, x);
}

static void mpfq_pz_wrapper_poly_neg(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED)
{
    mpfq_pz_poly_neg(vbase->obj, w, u);
}

static void mpfq_pz_wrapper_poly_sub_ui(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, unsigned long);
static void mpfq_pz_wrapper_poly_sub_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_pz_poly_sub_ui(vbase->obj, w, u, x);
}

static void mpfq_pz_wrapper_poly_add_ui(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, unsigned long);
static void mpfq_pz_wrapper_poly_add_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_pz_poly_add_ui(vbase->obj, w, u, x);
}

static void mpfq_pz_wrapper_poly_set_ui(mpfq_vbase_ptr, mpfq_pz_dst_poly, unsigned long);
static void mpfq_pz_wrapper_poly_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_pz_poly_set_ui(vbase->obj, w, x);
}

static void mpfq_pz_wrapper_poly_sub(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, mpfq_pz_src_poly v MAYBE_UNUSED)
{
    mpfq_pz_poly_sub(vbase->obj, w, u, v);
}

static void mpfq_pz_wrapper_poly_add(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED, mpfq_pz_src_poly v MAYBE_UNUSED)
{
    mpfq_pz_poly_add(vbase->obj, w, u, v);
}

static int mpfq_pz_wrapper_poly_deg(mpfq_vbase_ptr, mpfq_pz_src_poly);
static int mpfq_pz_wrapper_poly_deg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_poly w MAYBE_UNUSED)
{
    return mpfq_pz_poly_deg(vbase->obj, w);
}

static void mpfq_pz_wrapper_poly_getcoeff(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_poly, unsigned int);
static void mpfq_pz_wrapper_poly_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt x MAYBE_UNUSED, mpfq_pz_src_poly w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_pz_poly_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_pz_wrapper_poly_setcoeff_ui(mpfq_vbase_ptr, mpfq_pz_dst_poly, unsigned long, unsigned int);
static void mpfq_pz_wrapper_poly_setcoeff_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_pz_poly_setcoeff_ui(vbase->obj, w, x, i);
}

static void mpfq_pz_wrapper_poly_setcoeff(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_elt, unsigned int);
static void mpfq_pz_wrapper_poly_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_pz_poly_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_pz_wrapper_poly_setmonic(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_setmonic(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly q MAYBE_UNUSED, mpfq_pz_src_poly p MAYBE_UNUSED)
{
    mpfq_pz_poly_setmonic(vbase->obj, q, p);
}

static void mpfq_pz_wrapper_poly_set(mpfq_vbase_ptr, mpfq_pz_dst_poly, mpfq_pz_src_poly);
static void mpfq_pz_wrapper_poly_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_poly w MAYBE_UNUSED, mpfq_pz_src_poly u MAYBE_UNUSED)
{
    mpfq_pz_poly_set(vbase->obj, w, u);
}

static void mpfq_pz_wrapper_poly_clear(mpfq_vbase_ptr, mpfq_pz_poly);
static void mpfq_pz_wrapper_poly_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_poly p MAYBE_UNUSED)
{
    mpfq_pz_poly_clear(vbase->obj, p);
}

static void mpfq_pz_wrapper_poly_init(mpfq_vbase_ptr, mpfq_pz_poly, unsigned int);
static void mpfq_pz_wrapper_poly_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_poly p MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_pz_poly_init(vbase->obj, p, n);
}

static ptrdiff_t mpfq_pz_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr, long);
static ptrdiff_t mpfq_pz_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_ur_elt_stride(vbase->obj, n);
}

static ptrdiff_t mpfq_pz_wrapper_vec_elt_stride(mpfq_vbase_ptr, long);
static ptrdiff_t mpfq_pz_wrapper_vec_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_elt_stride(vbase->obj, n);
}

static mpfq_pz_src_elt mpfq_pz_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr, mpfq_pz_src_vec_ur, long);
static mpfq_pz_src_elt mpfq_pz_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_ur_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_pz_dst_elt mpfq_pz_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, long);
static mpfq_pz_dst_elt mpfq_pz_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_ur_coeff_ptr(vbase->obj, v, i);
}

static mpfq_pz_src_vec_ur mpfq_pz_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr, mpfq_pz_src_vec_ur, long);
static mpfq_pz_src_vec_ur mpfq_pz_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_ur_subvec_const(vbase->obj, v, i);
}

static mpfq_pz_dst_vec_ur mpfq_pz_wrapper_vec_ur_subvec(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, long);
static mpfq_pz_dst_vec_ur mpfq_pz_wrapper_vec_ur_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_ur_subvec(vbase->obj, v, i);
}

static void mpfq_pz_wrapper_vec_reduce(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_dst_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_reduce(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_dst_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_reduce(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_conv_ur(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, unsigned long, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_conv_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_pz_vec_conv_ur(vbase->obj, w, u, n, v, m);
}

static void mpfq_pz_wrapper_vec_scal_mul_ur(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_vec_scal_mul_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, mpfq_pz_src_elt c MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_scal_mul_ur(vbase->obj, w, u, c, n);
}

static void mpfq_pz_wrapper_vec_ur_rev(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_rev(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_ur_neg(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_neg(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_ur_sub(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec_ur u MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_sub(vbase->obj, w, u, v, n);
}

static void mpfq_pz_wrapper_vec_ur_add(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec_ur u MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_add(vbase->obj, w, u, v, n);
}

static void mpfq_pz_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_getcoeff(vbase->obj, z, v, i);
}

static void mpfq_pz_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_elt_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_elt_ur x MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_pz_wrapper_vec_ur_set(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_set(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_ur_clear(mpfq_vbase_ptr, mpfq_pz_vec_ur *, unsigned long);
static void mpfq_pz_wrapper_vec_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_clear(vbase->obj, v, n);
}

static void mpfq_pz_wrapper_vec_ur_reinit(mpfq_vbase_ptr, mpfq_pz_vec_ur *, unsigned long, unsigned long);
static void mpfq_pz_wrapper_vec_ur_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_reinit(vbase->obj, v, n, m);
}

static void mpfq_pz_wrapper_vec_ur_set_vec(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_ur_set_vec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_set_vec(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_ur_set_zero(mpfq_vbase_ptr, mpfq_pz_dst_vec_ur, unsigned long);
static void mpfq_pz_wrapper_vec_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec_ur w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_set_zero(vbase->obj, w, n);
}

static void mpfq_pz_wrapper_vec_ur_init(mpfq_vbase_ptr, mpfq_pz_vec_ur *, unsigned long);
static void mpfq_pz_wrapper_vec_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_ur_init(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_simd_find_first_set(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_simd_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_simd_find_first_set(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_simd_hamming_weight(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_simd_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_simd_hamming_weight(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_find_first_set(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_find_first_set(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_hamming_weight(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_hamming_weight(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_scan(mpfq_vbase_ptr, mpfq_pz_vec *, unsigned long *);
static long mpfq_pz_wrapper_vec_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED)
{
    return mpfq_pz_vec_scan(vbase->obj, w, n);
}

static long mpfq_pz_wrapper_vec_fscan(mpfq_vbase_ptr, FILE *, mpfq_pz_vec *, unsigned long *);
static long mpfq_pz_wrapper_vec_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED)
{
    return mpfq_pz_vec_fscan(vbase->obj, file, w, n);
}

static long mpfq_pz_wrapper_vec_sscan(mpfq_vbase_ptr, mpfq_pz_vec *, unsigned long *, const char *);
static long mpfq_pz_wrapper_vec_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_pz_vec_sscan(vbase->obj, w, n, str);
}

static long mpfq_pz_wrapper_vec_print(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_print(vbase->obj, w, n);
}

static long mpfq_pz_wrapper_vec_fprint(mpfq_vbase_ptr, FILE *, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_fprint(vbase->obj, file, v, n);
}

static long mpfq_pz_wrapper_vec_asprint(mpfq_vbase_ptr, char * *, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_asprint(vbase->obj, pstr, v, n);
}

static mpfq_pz_src_elt mpfq_pz_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr, mpfq_pz_src_vec, long);
static mpfq_pz_src_elt mpfq_pz_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_pz_dst_elt mpfq_pz_wrapper_vec_coeff_ptr(mpfq_vbase_ptr, mpfq_pz_dst_vec, long);
static mpfq_pz_dst_elt mpfq_pz_wrapper_vec_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_coeff_ptr(vbase->obj, v, i);
}

static mpfq_pz_src_vec mpfq_pz_wrapper_vec_subvec_const(mpfq_vbase_ptr, mpfq_pz_src_vec, long);
static mpfq_pz_src_vec mpfq_pz_wrapper_vec_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_subvec_const(vbase->obj, v, i);
}

static mpfq_pz_dst_vec mpfq_pz_wrapper_vec_subvec(mpfq_vbase_ptr, mpfq_pz_dst_vec, long);
static mpfq_pz_dst_vec mpfq_pz_wrapper_vec_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_pz_vec_subvec(vbase->obj, v, i);
}

static long mpfq_pz_wrapper_vec_is_zero(mpfq_vbase_ptr, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_is_zero(vbase->obj, v, n);
}

static long mpfq_pz_wrapper_vec_cmp(mpfq_vbase_ptr, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
static long mpfq_pz_wrapper_vec_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_pz_vec_cmp(vbase->obj, u, v, n);
}

static void mpfq_pz_wrapper_vec_random2(mpfq_vbase_ptr, mpfq_pz_dst_vec, unsigned long, gmp_randstate_t);
static void mpfq_pz_wrapper_vec_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_vec_random2(vbase->obj, w, n, state);
}

static void mpfq_pz_wrapper_vec_random(mpfq_vbase_ptr, mpfq_pz_dst_vec, unsigned long, gmp_randstate_t);
static void mpfq_pz_wrapper_vec_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_vec_random(vbase->obj, w, n, state);
}

static void mpfq_pz_wrapper_vec_conv(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_conv(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_pz_vec_conv(vbase->obj, w, u, n, v, m);
}

static void mpfq_pz_wrapper_vec_scal_mul(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_vec_scal_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, mpfq_pz_src_elt c MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_scal_mul(vbase->obj, w, u, c, n);
}

static void mpfq_pz_wrapper_vec_sub(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_sub(vbase->obj, w, u, v, n);
}

static void mpfq_pz_wrapper_vec_rev(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_rev(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_neg(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_neg(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_add(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec u MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_add(vbase->obj, w, u, v, n);
}

static void mpfq_pz_wrapper_vec_getcoeff(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_pz_vec_getcoeff(vbase->obj, z, v, i);
}

static void mpfq_pz_wrapper_vec_setcoeff_ui(mpfq_vbase_ptr, mpfq_pz_dst_vec, unsigned long, unsigned long);
static void mpfq_pz_wrapper_vec_setcoeff_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_pz_vec_setcoeff_ui(vbase->obj, w, x0, i);
}

static void mpfq_pz_wrapper_vec_setcoeff(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_vec_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_pz_vec_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_pz_wrapper_vec_set_zero(mpfq_vbase_ptr, mpfq_pz_dst_vec, unsigned long);
static void mpfq_pz_wrapper_vec_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_set_zero(vbase->obj, w, n);
}

static void mpfq_pz_wrapper_vec_set(mpfq_vbase_ptr, mpfq_pz_dst_vec, mpfq_pz_src_vec, unsigned long);
static void mpfq_pz_wrapper_vec_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_vec w MAYBE_UNUSED, mpfq_pz_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_set(vbase->obj, w, v, n);
}

static void mpfq_pz_wrapper_vec_clear(mpfq_vbase_ptr, mpfq_pz_vec *, unsigned long);
static void mpfq_pz_wrapper_vec_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_clear(vbase->obj, v, n);
}

static void mpfq_pz_wrapper_vec_reinit(mpfq_vbase_ptr, mpfq_pz_vec *, unsigned long, unsigned long);
static void mpfq_pz_wrapper_vec_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_pz_vec_reinit(vbase->obj, v, n, m);
}

static void mpfq_pz_wrapper_vec_init(mpfq_vbase_ptr, mpfq_pz_vec *, unsigned long);
static void mpfq_pz_wrapper_vec_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_pz_vec_init(vbase->obj, v, n);
}

static int mpfq_pz_wrapper_scan(mpfq_vbase_ptr, mpfq_pz_dst_elt);
static int mpfq_pz_wrapper_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt x MAYBE_UNUSED)
{
    return mpfq_pz_scan(vbase->obj, x);
}

static int mpfq_pz_wrapper_fscan(mpfq_vbase_ptr, FILE *, mpfq_pz_dst_elt);
static int mpfq_pz_wrapper_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_dst_elt x MAYBE_UNUSED)
{
    return mpfq_pz_fscan(vbase->obj, file, x);
}

static int mpfq_pz_wrapper_sscan(mpfq_vbase_ptr, mpfq_pz_dst_elt, const char *);
static int mpfq_pz_wrapper_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_pz_sscan(vbase->obj, z, str);
}

static int mpfq_pz_wrapper_print(mpfq_vbase_ptr, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_print(vbase->obj, x);
}

static int mpfq_pz_wrapper_fprint(mpfq_vbase_ptr, FILE *, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_fprint(vbase->obj, file, x);
}

static int mpfq_pz_wrapper_asprint(mpfq_vbase_ptr, char * *, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_asprint(vbase->obj, pstr, x);
}

static int mpfq_pz_wrapper_is_zero(mpfq_vbase_ptr, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_is_zero(vbase->obj, x);
}

static int mpfq_pz_wrapper_cmp_ui(mpfq_vbase_ptr, mpfq_pz_src_elt, unsigned long);
static int mpfq_pz_wrapper_cmp_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long y0 MAYBE_UNUSED)
{
    return mpfq_pz_cmp_ui(vbase->obj, x, y0);
}

static int mpfq_pz_wrapper_cmp(mpfq_vbase_ptr, mpfq_pz_src_elt, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    return mpfq_pz_cmp(vbase->obj, x, y);
}

static void mpfq_pz_wrapper_addmul_si_ur(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt, long);
static void mpfq_pz_wrapper_addmul_si_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur w MAYBE_UNUSED, mpfq_pz_src_elt u MAYBE_UNUSED, long v MAYBE_UNUSED)
{
    mpfq_pz_addmul_si_ur(vbase->obj, w, u, v);
}

static void mpfq_pz_wrapper_reduce(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_dst_elt_ur);
static void mpfq_pz_wrapper_reduce(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_dst_elt_ur x MAYBE_UNUSED)
{
    mpfq_pz_reduce(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_sqr_ur(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_sqr_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_sqr_ur(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_mul_ur(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_mul_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    mpfq_pz_mul_ur(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_elt_ur_sub(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur, mpfq_pz_src_elt_ur);
static void mpfq_pz_wrapper_elt_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt_ur x MAYBE_UNUSED, mpfq_pz_src_elt_ur y MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_sub(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_elt_ur_neg(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur);
static void mpfq_pz_wrapper_elt_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt_ur x MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_neg(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_elt_ur_add(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur, mpfq_pz_src_elt_ur);
static void mpfq_pz_wrapper_elt_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt_ur x MAYBE_UNUSED, mpfq_pz_src_elt_ur y MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_add(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_elt_ur_set_ui(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, unsigned long);
static void mpfq_pz_wrapper_elt_ur_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_set_ui(vbase->obj, z, x0);
}

static void mpfq_pz_wrapper_elt_ur_set_zero(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur);
static void mpfq_pz_wrapper_elt_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_set_zero(vbase->obj, z);
}

static void mpfq_pz_wrapper_elt_ur_set_elt(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_elt_ur_set_elt(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_set_elt(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_elt_ur_set(mpfq_vbase_ptr, mpfq_pz_dst_elt_ur, mpfq_pz_src_elt_ur);
static void mpfq_pz_wrapper_elt_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt_ur z MAYBE_UNUSED, mpfq_pz_src_elt_ur x MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_set(vbase->obj, z, x);
}

static ptrdiff_t mpfq_pz_wrapper_elt_ur_stride(mpfq_vbase_ptr);
static ptrdiff_t mpfq_pz_wrapper_elt_ur_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_elt_ur_stride(vbase->obj);
}

static void mpfq_pz_wrapper_elt_ur_clear(mpfq_vbase_ptr, mpfq_pz_elt_ur *);
static void mpfq_pz_wrapper_elt_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_elt_ur * x MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_clear(vbase->obj, x);
}

static void mpfq_pz_wrapper_elt_ur_init(mpfq_vbase_ptr, mpfq_pz_elt_ur *);
static void mpfq_pz_wrapper_elt_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_elt_ur * x MAYBE_UNUSED)
{
    mpfq_pz_elt_ur_init(vbase->obj, x);
}

static int mpfq_pz_wrapper_inv(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_inv(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_inv(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_normalize(mpfq_vbase_ptr, mpfq_pz_dst_elt);
static void mpfq_pz_wrapper_normalize(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED)
{
    mpfq_pz_normalize(vbase->obj, z);
}

static void mpfq_pz_wrapper_mul_ui(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_mul_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_pz_mul_ui(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_sub_ui(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_sub_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_pz_sub_ui(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_add_ui(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long);
static void mpfq_pz_wrapper_add_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_pz_add_ui(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_frobenius(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_frobenius(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    mpfq_pz_frobenius(vbase->obj, x, y);
}

static void mpfq_pz_wrapper_powz(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpz_srcptr);
static void mpfq_pz_wrapper_powz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt y MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpz_srcptr z MAYBE_UNUSED)
{
    mpfq_pz_powz(vbase->obj, y, x, z);
}

static void mpfq_pz_wrapper_pow(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, unsigned long *, size_t);
static void mpfq_pz_wrapper_pow(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED, size_t nl MAYBE_UNUSED)
{
    mpfq_pz_pow(vbase->obj, z, x, n, nl);
}

static int mpfq_pz_wrapper_is_sqr(mpfq_vbase_ptr, mpfq_pz_src_elt);
static int mpfq_pz_wrapper_is_sqr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_is_sqr(vbase->obj, x);
}

static void mpfq_pz_wrapper_sqr(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_sqr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_sqr(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_mul(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    mpfq_pz_mul(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_neg(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_neg(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_sub(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    mpfq_pz_sub(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_add(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED, mpfq_pz_src_elt y MAYBE_UNUSED)
{
    mpfq_pz_add(vbase->obj, z, x, y);
}

static void mpfq_pz_wrapper_random2(mpfq_vbase_ptr, mpfq_pz_dst_elt, gmp_randstate_t);
static void mpfq_pz_wrapper_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_random2(vbase->obj, z, state);
}

static void mpfq_pz_wrapper_random(mpfq_vbase_ptr, mpfq_pz_dst_elt, gmp_randstate_t);
static void mpfq_pz_wrapper_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_pz_random(vbase->obj, z, state);
}

static void mpfq_pz_wrapper_get_mpz(mpfq_vbase_ptr, mpz_ptr, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_get_mpz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpz_ptr z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_get_mpz(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_get_mpn(mpfq_vbase_ptr, mp_limb_t *, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_get_mpn(mpfq_vbase_ptr vbase MAYBE_UNUSED, mp_limb_t * z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_get_mpn(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_set_mpz(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpz_srcptr);
static void mpfq_pz_wrapper_set_mpz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpz_srcptr x MAYBE_UNUSED)
{
    mpfq_pz_set_mpz(vbase->obj, z, x);
}

static void mpfq_pz_wrapper_set_mpn(mpfq_vbase_ptr, mpfq_pz_dst_elt, mp_limb_t *, size_t);
static void mpfq_pz_wrapper_set_mpn(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mp_limb_t * x MAYBE_UNUSED, size_t n MAYBE_UNUSED)
{
    mpfq_pz_set_mpn(vbase->obj, z, x, n);
}

static unsigned long mpfq_pz_wrapper_get_ui(mpfq_vbase_ptr, mpfq_pz_src_elt);
static unsigned long mpfq_pz_wrapper_get_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    return mpfq_pz_get_ui(vbase->obj, x);
}

static void mpfq_pz_wrapper_set_zero(mpfq_vbase_ptr, mpfq_pz_dst_elt);
static void mpfq_pz_wrapper_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED)
{
    mpfq_pz_set_zero(vbase->obj, z);
}

static void mpfq_pz_wrapper_set_ui(mpfq_vbase_ptr, mpfq_pz_dst_elt, unsigned long);
static void mpfq_pz_wrapper_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, unsigned long x0 MAYBE_UNUSED)
{
    mpfq_pz_set_ui(vbase->obj, z, x0);
}

static void mpfq_pz_wrapper_set(mpfq_vbase_ptr, mpfq_pz_dst_elt, mpfq_pz_src_elt);
static void mpfq_pz_wrapper_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_dst_elt z MAYBE_UNUSED, mpfq_pz_src_elt x MAYBE_UNUSED)
{
    mpfq_pz_set(vbase->obj, z, x);
}

static ptrdiff_t mpfq_pz_wrapper_elt_stride(mpfq_vbase_ptr);
static ptrdiff_t mpfq_pz_wrapper_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_elt_stride(vbase->obj);
}

static void mpfq_pz_wrapper_clear(mpfq_vbase_ptr, mpfq_pz_elt *);
static void mpfq_pz_wrapper_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_elt * x MAYBE_UNUSED)
{
    mpfq_pz_clear(vbase->obj, x);
}

static void mpfq_pz_wrapper_init(mpfq_vbase_ptr, mpfq_pz_elt *);
static void mpfq_pz_wrapper_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_pz_elt * x MAYBE_UNUSED)
{
    mpfq_pz_init(vbase->obj, x);
}

static void mpfq_pz_wrapper_field_setopt(mpfq_vbase_ptr, unsigned long, void *);
static void mpfq_pz_wrapper_field_setopt(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, void * y MAYBE_UNUSED)
{
    mpfq_pz_field_setopt(vbase->obj, x, y);
}

static void mpfq_pz_wrapper_field_specify(mpfq_vbase_ptr, unsigned long, const void *);
static void mpfq_pz_wrapper_field_specify(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long dummy MAYBE_UNUSED, const void * vp MAYBE_UNUSED)
{
    mpfq_pz_field_specify(vbase->obj, dummy, vp);
}

static void mpfq_pz_wrapper_field_clear(mpfq_vbase_ptr);
static void mpfq_pz_wrapper_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_pz_field_clear(vbase->obj);
}

static void mpfq_pz_wrapper_field_init(mpfq_vbase_ptr);
static void mpfq_pz_wrapper_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_pz_field_init(vbase->obj);
}

static int mpfq_pz_wrapper_field_degree(mpfq_vbase_ptr);
static int mpfq_pz_wrapper_field_degree(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_field_degree(vbase->obj);
}

static unsigned long mpfq_pz_wrapper_field_characteristic_bits(mpfq_vbase_ptr);
static unsigned long mpfq_pz_wrapper_field_characteristic_bits(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_field_characteristic_bits(vbase->obj);
}

static mpz_srcptr mpfq_pz_wrapper_field_characteristic_srcptr(mpfq_vbase_ptr);
static mpz_srcptr mpfq_pz_wrapper_field_characteristic_srcptr(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_pz_field_characteristic_srcptr(vbase->obj);
}

static void mpfq_pz_wrapper_field_characteristic(mpfq_vbase_ptr, mpz_ptr);
static void mpfq_pz_wrapper_field_characteristic(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpz_ptr z MAYBE_UNUSED)
{
    mpfq_pz_field_characteristic(vbase->obj, z);
}

static unsigned long mpfq_pz_wrapper_impl_max_degree();
static unsigned long mpfq_pz_wrapper_impl_max_degree()
{
    return mpfq_pz_impl_max_degree();
}

static unsigned long mpfq_pz_wrapper_impl_max_characteristic_bits();
static unsigned long mpfq_pz_wrapper_impl_max_characteristic_bits()
{
    return mpfq_pz_impl_max_characteristic_bits();
}

static const char * mpfq_pz_wrapper_impl_name();
static const char * mpfq_pz_wrapper_impl_name()
{
    return mpfq_pz_impl_name();
}

/* Mpfq::engine::oo::oo_field_init */
/* Triggered by: oo */
void mpfq_pz_oo_field_init(mpfq_vbase_ptr vbase)
{
    memset(vbase, 0, sizeof(struct mpfq_vbase_s));
    vbase->obj = malloc(sizeof(mpfq_pz_field));
    mpfq_pz_field_init((mpfq_pz_dst_field) vbase->obj);
    vbase->impl_name = (const char * (*) ()) mpfq_pz_wrapper_impl_name;
    vbase->impl_max_characteristic_bits = (unsigned long (*) ()) mpfq_pz_wrapper_impl_max_characteristic_bits;
    vbase->impl_max_degree = (unsigned long (*) ()) mpfq_pz_wrapper_impl_max_degree;
    vbase->field_characteristic = (void (*) (mpfq_vbase_ptr, mpz_ptr)) mpfq_pz_wrapper_field_characteristic;
    vbase->field_characteristic_srcptr = (mpz_srcptr (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_field_characteristic_srcptr;
    vbase->field_characteristic_bits = (unsigned long (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_field_characteristic_bits;
    vbase->field_degree = (int (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_field_degree;
    vbase->field_init = (void (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_field_init;
    vbase->field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_field_clear;
    vbase->field_specify = (void (*) (mpfq_vbase_ptr, unsigned long, const void *)) mpfq_pz_wrapper_field_specify;
    vbase->field_setopt = (void (*) (mpfq_vbase_ptr, unsigned long, void *)) mpfq_pz_wrapper_field_setopt;
    vbase->init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_init;
    vbase->clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_clear;
    vbase->elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_elt_stride;
    vbase->set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_set;
    vbase->set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_set_ui;
    vbase->set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_set_zero;
    vbase->get_ui = (unsigned long (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_get_ui;
    vbase->set_mpn = (void (*) (mpfq_vbase_ptr, void *, mp_limb_t *, size_t)) mpfq_pz_wrapper_set_mpn;
    vbase->set_mpz = (void (*) (mpfq_vbase_ptr, void *, mpz_srcptr)) mpfq_pz_wrapper_set_mpz;
    vbase->get_mpn = (void (*) (mpfq_vbase_ptr, mp_limb_t *, const void *)) mpfq_pz_wrapper_get_mpn;
    vbase->get_mpz = (void (*) (mpfq_vbase_ptr, mpz_ptr, const void *)) mpfq_pz_wrapper_get_mpz;
    vbase->random = (void (*) (mpfq_vbase_ptr, void *, gmp_randstate_t)) mpfq_pz_wrapper_random;
    vbase->random2 = (void (*) (mpfq_vbase_ptr, void *, gmp_randstate_t)) mpfq_pz_wrapper_random2;
    vbase->add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_add;
    vbase->sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_sub;
    vbase->neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_neg;
    vbase->mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_mul;
    vbase->sqr = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_sqr;
    vbase->is_sqr = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_is_sqr;
    /* missing sqrt */
    vbase->pow = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long *, size_t)) mpfq_pz_wrapper_pow;
    vbase->powz = (void (*) (mpfq_vbase_ptr, void *, const void *, mpz_srcptr)) mpfq_pz_wrapper_powz;
    vbase->frobenius = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_frobenius;
    vbase->add_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_add_ui;
    vbase->sub_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_sub_ui;
    vbase->mul_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_mul_ui;
    vbase->normalize = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_normalize;
    vbase->inv = (int (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_inv;
    vbase->elt_ur_init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_elt_ur_init;
    vbase->elt_ur_clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_elt_ur_clear;
    vbase->elt_ur_stride = (ptrdiff_t (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_elt_ur_stride;
    vbase->elt_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_elt_ur_set;
    vbase->elt_ur_set_elt = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_elt_ur_set_elt;
    vbase->elt_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_elt_ur_set_zero;
    vbase->elt_ur_set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_elt_ur_set_ui;
    vbase->elt_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_elt_ur_add;
    vbase->elt_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_elt_ur_neg;
    vbase->elt_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_elt_ur_sub;
    vbase->mul_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_mul_ur;
    vbase->sqr_ur = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_sqr_ur;
    vbase->reduce = (void (*) (mpfq_vbase_ptr, void *, void *)) mpfq_pz_wrapper_reduce;
    vbase->addmul_si_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, long)) mpfq_pz_wrapper_addmul_si_ur;
    vbase->cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *)) mpfq_pz_wrapper_cmp;
    vbase->cmp_ui = (int (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_cmp_ui;
    vbase->is_zero = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_is_zero;
    vbase->asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *)) mpfq_pz_wrapper_asprint;
    vbase->fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *)) mpfq_pz_wrapper_fprint;
    vbase->print = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_print;
    vbase->sscan = (int (*) (mpfq_vbase_ptr, void *, const char *)) mpfq_pz_wrapper_sscan;
    vbase->fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *)) mpfq_pz_wrapper_fscan;
    vbase->scan = (int (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_scan;
    /* missing read */
    /* missing importdata */
    /* missing write */
    /* missing exportdata */
    vbase->vec_init = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_init;
    vbase->vec_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_pz_wrapper_vec_reinit;
    vbase->vec_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_clear;
    vbase->vec_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_set;
    vbase->vec_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_set_zero;
    vbase->vec_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_setcoeff;
    vbase->vec_setcoeff_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_pz_wrapper_vec_setcoeff_ui;
    vbase->vec_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_getcoeff;
    vbase->vec_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_add;
    vbase->vec_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_neg;
    vbase->vec_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_rev;
    vbase->vec_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_sub;
    vbase->vec_scal_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_scal_mul;
    vbase->vec_conv = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long)) mpfq_pz_wrapper_vec_conv;
    vbase->vec_random = (void (*) (mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t)) mpfq_pz_wrapper_vec_random;
    vbase->vec_random2 = (void (*) (mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t)) mpfq_pz_wrapper_vec_random2;
    vbase->vec_cmp = (long (*) (mpfq_vbase_ptr, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_cmp;
    vbase->vec_is_zero = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_is_zero;
    vbase->vec_subvec = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_pz_wrapper_vec_subvec;
    vbase->vec_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_pz_wrapper_vec_subvec_const;
    vbase->vec_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_pz_wrapper_vec_coeff_ptr;
    vbase->vec_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_pz_wrapper_vec_coeff_ptr_const;
    vbase->vec_asprint = (long (*) (mpfq_vbase_ptr, char * *, const void *, unsigned long)) mpfq_pz_wrapper_vec_asprint;
    vbase->vec_fprint = (long (*) (mpfq_vbase_ptr, FILE *, const void *, unsigned long)) mpfq_pz_wrapper_vec_fprint;
    vbase->vec_print = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_print;
    vbase->vec_sscan = (long (*) (mpfq_vbase_ptr, void *, unsigned long *, const char *)) mpfq_pz_wrapper_vec_sscan;
    vbase->vec_fscan = (long (*) (mpfq_vbase_ptr, FILE *, void *, unsigned long *)) mpfq_pz_wrapper_vec_fscan;
    vbase->vec_scan = (long (*) (mpfq_vbase_ptr, void *, unsigned long *)) mpfq_pz_wrapper_vec_scan;
    /* missing vec_read */
    /* missing vec_write */
    /* missing vec_import */
    /* missing vec_export */
    vbase->vec_hamming_weight = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_hamming_weight;
    vbase->vec_find_first_set = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_find_first_set;
    vbase->vec_simd_hamming_weight = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_simd_hamming_weight;
    vbase->vec_simd_find_first_set = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_pz_wrapper_vec_simd_find_first_set;
    vbase->vec_ur_init = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_ur_init;
    vbase->vec_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_ur_set_zero;
    vbase->vec_ur_set_vec = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_set_vec;
    vbase->vec_ur_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_pz_wrapper_vec_ur_reinit;
    vbase->vec_ur_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_vec_ur_clear;
    vbase->vec_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_set;
    vbase->vec_ur_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_setcoeff;
    vbase->vec_ur_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_getcoeff;
    vbase->vec_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_add;
    vbase->vec_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_sub;
    vbase->vec_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_neg;
    vbase->vec_ur_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_ur_rev;
    vbase->vec_scal_mul_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_pz_wrapper_vec_scal_mul_ur;
    vbase->vec_conv_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long)) mpfq_pz_wrapper_vec_conv_ur;
    vbase->vec_reduce = (void (*) (mpfq_vbase_ptr, void *, void *, unsigned long)) mpfq_pz_wrapper_vec_reduce;
    vbase->vec_ur_subvec = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_pz_wrapper_vec_ur_subvec;
    vbase->vec_ur_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_pz_wrapper_vec_ur_subvec_const;
    vbase->vec_ur_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_pz_wrapper_vec_ur_coeff_ptr;
    vbase->vec_ur_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_pz_wrapper_vec_ur_coeff_ptr_const;
    vbase->vec_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, long)) mpfq_pz_wrapper_vec_elt_stride;
    vbase->vec_ur_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, long)) mpfq_pz_wrapper_vec_ur_elt_stride;
    vbase->poly_init = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_pz_wrapper_poly_init;
    vbase->poly_clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_poly_clear;
    vbase->poly_set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_poly_set;
    vbase->poly_setmonic = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_poly_setmonic;
    vbase->poly_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_pz_wrapper_poly_setcoeff;
    vbase->poly_setcoeff_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned int)) mpfq_pz_wrapper_poly_setcoeff_ui;
    vbase->poly_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_pz_wrapper_poly_getcoeff;
    vbase->poly_deg = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_poly_deg;
    vbase->poly_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_poly_add;
    vbase->poly_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_poly_sub;
    vbase->poly_set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_poly_set_ui;
    vbase->poly_add_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_poly_add_ui;
    vbase->poly_sub_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_pz_wrapper_poly_sub_ui;
    vbase->poly_neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_poly_neg;
    vbase->poly_scal_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_poly_scal_mul;
    vbase->poly_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_poly_mul;
    vbase->poly_divmod = (int (*) (mpfq_vbase_ptr, void *, void *, const void *, const void *)) mpfq_pz_wrapper_poly_divmod;
    vbase->poly_precomp_mod = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_pz_wrapper_poly_precomp_mod;
    vbase->poly_mod_pre = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, const void *)) mpfq_pz_wrapper_poly_mod_pre;
    vbase->poly_gcd = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_pz_wrapper_poly_gcd;
    vbase->poly_xgcd = (void (*) (mpfq_vbase_ptr, void *, void *, void *, const void *, const void *)) mpfq_pz_wrapper_poly_xgcd;
    vbase->poly_random = (void (*) (mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t)) mpfq_pz_wrapper_poly_random;
    vbase->poly_random2 = (void (*) (mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t)) mpfq_pz_wrapper_poly_random2;
    vbase->poly_cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *)) mpfq_pz_wrapper_poly_cmp;
    vbase->poly_asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *)) mpfq_pz_wrapper_poly_asprint;
    vbase->poly_fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *)) mpfq_pz_wrapper_poly_fprint;
    vbase->poly_print = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_poly_print;
    vbase->poly_sscan = (int (*) (mpfq_vbase_ptr, void *, const char *)) mpfq_pz_wrapper_poly_sscan;
    vbase->poly_fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *)) mpfq_pz_wrapper_poly_fscan;
    vbase->poly_scan = (int (*) (mpfq_vbase_ptr, void *)) mpfq_pz_wrapper_poly_scan;
    vbase->simd_groupsize = (int (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_simd_groupsize;
    vbase->simd_hamming_weight = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_simd_hamming_weight;
    vbase->simd_find_first_set = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_pz_wrapper_simd_find_first_set;
    vbase->simd_get_ui_at = (unsigned long (*) (mpfq_vbase_ptr, const void *, int)) mpfq_pz_wrapper_simd_get_ui_at;
    vbase->simd_set_ui_at = (void (*) (mpfq_vbase_ptr, void *, int, unsigned long)) mpfq_pz_wrapper_simd_set_ui_at;
    vbase->simd_add_ui_at = (void (*) (mpfq_vbase_ptr, void *, const void *, int, unsigned long)) mpfq_pz_wrapper_simd_add_ui_at;
    vbase->simd_set_ui_all = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_pz_wrapper_simd_set_ui_all;
    vbase->add_dotprod = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_pz_wrapper_add_dotprod;
    vbase->oo_field_init = (void (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_oo_field_init;
    vbase->oo_field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_pz_wrapper_oo_field_clear;
}


/* vim:set ft=cpp: */
