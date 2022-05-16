#include "cado.h"
/* MPFQ generated file -- do not edit */

#include "mpfq_p_5.h"

#include <inttypes.h>
/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_5, tag=p_5, }, ],
   fieldtype=prime,
   n=5,
   nn=11,
   opthw=,
   tag=p_5,
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
     [ (?^:mpfq_p_5_elt \*), void *, ],
     [ (?^:mpfq_p_5_src_elt\b), const void *, ],
     [ (?^:mpfq_p_5_elt\b), void *, ],
     [ (?^:mpfq_p_5_dst_elt\b), void *, ],
     [ (?^:mpfq_p_5_elt_ur \*), void *, ],
     [ (?^:mpfq_p_5_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_p_5_elt_ur\b), void *, ],
     [ (?^:mpfq_p_5_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_p_5_vec \*), void *, ],
     [ (?^:mpfq_p_5_src_vec\b), const void *, ],
     [ (?^:mpfq_p_5_vec\b), void *, ],
     [ (?^:mpfq_p_5_dst_vec\b), void *, ],
     [ (?^:mpfq_p_5_vec_ur \*), void *, ],
     [ (?^:mpfq_p_5_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_p_5_vec_ur\b), void *, ],
     [ (?^:mpfq_p_5_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_p_5_poly \*), void *, ],
     [ (?^:mpfq_p_5_src_poly\b), const void *, ],
     [ (?^:mpfq_p_5_poly\b), void *, ],
     [ (?^:mpfq_p_5_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */


/* Functions operating on the field structure */
/* *Mpfq::gfp::field::code_for_field_clear, Mpfq::gfp */
void mpfq_p_5_field_clear(mpfq_p_5_dst_field k)
{
        mpz_clear(k->p);
        mpz_clear(k->bigmul_p);
        if (k->ts_info.e > 0) {
            free(k->ts_info.hh);
            free(k->ts_info.z);
        }
        mpz_clear(k->factor);
}

/* *Mpfq::gfp::field::code_for_field_specify, Mpfq::gfp */
void mpfq_p_5_field_specify(mpfq_p_5_dst_field k, unsigned long dummy MAYBE_UNUSED, const void * vp)
{
        if (dummy == MPFQ_PRIME_MPN) {
            fprintf(stderr, "MPFQ_PRIME_MPN is deprecated\n");
            return;
        } else if (dummy == MPFQ_PRIME_MPZ) {
            mpz_srcptr p = (mpz_srcptr) vp;
            if (!(mpz_size(p) == 5)) {
                fprintf(stderr, "This binary requires the use of a 5-machine words prime number. Here, p spans %zu machine words. Please adapt linalg/bwc/CMakeLists.txt accordingly and re-run\n", mpz_size(p));
                abort();
            }
            mpz_set(k->p, p);
            {
                /* precompute bigmul_p = largest multiple of p that fits in an
                 * elt_ur: p*Floor( (2^(11*64)-1)/p )
                 */
                mpz_ui_pow_ui(k->bigmul_p, 2, 11*64);
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

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */
static void mpfq_p_5_init_ts(mpfq_p_5_dst_field);
/* *Mpfq::gfp::elt::code_for_sqrt, Mpfq::gfp */
/* Triggered by: sqrt */
static void mpfq_p_5_init_ts(mpfq_p_5_dst_field k)
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    
    mp_limb_t pp[5];
    mp_limb_t *ptr = pp;
    mp_limb_t s[5];
    mpfq_fixmp_5_sub_ui_nc(pp, k->p->_mp_d, 1);
    int e = 0;
    for( ; e < 5*64 && *ptr == 0 ; e+=64, ptr++) ;
    if (e >= 5*64) abort();
    int ee;
    ee = mpfq_ctzl(*ptr);
    e += ee;
    if (e < 64) {
        mpfq_fixmp_5_rshift(pp, e);
    } else {
        mpfq_fixmp_5_long_rshift(pp, e/64, e%64);
    }
    s[0] = 1UL;
    int i;
    for (i = 1; i <5; ++i)
        s[i] = 0UL;
    if (e-1 < 64) {
        mpfq_fixmp_5_lshift(s, e-1);
    } else {
        mpfq_fixmp_5_long_lshift(s, (e-1)/64, (e-1)%64);
    }
    k->ts_info.e = e;
    
    k->ts_info.z = mpfq_malloc_check(5*sizeof(mp_limb_t));
    k->ts_info.hh = mpfq_malloc_check(5*sizeof(mp_limb_t));
    
    mpfq_p_5_elt z, r;
    mpfq_p_5_init(k, &z);
    mpfq_p_5_init(k, &r);
    mpfq_p_5_set_ui(k, r, 0);
    do {
        mpfq_p_5_random(k, z, rstate);
        mpfq_p_5_pow(k, z, z, pp, 5);
        mpfq_p_5_pow(k, r, z, s, 5);
        mpfq_p_5_add_ui(k, r, r, 1);
    } while (mpfq_p_5_cmp_ui(k, r, 0)!=0);
    mpfq_p_5_set(k, (mpfq_p_5_dst_elt)k->ts_info.z, z);
    mpfq_p_5_clear(k, &z);
    mpfq_p_5_clear(k, &r);
    
    mpfq_fixmp_5_sub_ui_nc(pp, pp, 1);
    mpfq_fixmp_5_rshift(pp, 1);
    for (i = 0; i < 5; ++i)
        k->ts_info.hh[i] = pp[i];
    gmp_randclear(rstate);
}

/* *Mpfq::gfp::elt::code_for_sqrt, Mpfq::gfp */
int mpfq_p_5_sqrt(mpfq_p_5_dst_field k, mpfq_p_5_dst_elt z, mpfq_p_5_src_elt a)
{
    if (mpfq_p_5_cmp_ui(k, a, 0) == 0) {
        mpfq_p_5_set_ui(k, z, 0);
        return 1;
    }
    if (k->ts_info.e == 0)
        mpfq_p_5_init_ts(k);
    mpfq_p_5_elt b, x, y;
    mpfq_p_5_init(k, &x);
    mpfq_p_5_init(k, &y);
    mpfq_p_5_init(k, &b);
    mp_limb_t r = k->ts_info.e;
    mp_limb_t s; //= (1UL<<(r-1)); not needed...
    mpfq_p_5_set(k, x, a);
    mpfq_p_5_set(k, y, (mpfq_p_5_src_elt)k->ts_info.z);
    
    mpfq_p_5_pow(k, x, a, k->ts_info.hh, 5);
    mpfq_p_5_sqr(k, b, x);
    mpfq_p_5_mul(k, x, x, a);
    mpfq_p_5_mul(k, b, b, a);
    
    mpfq_p_5_elt t;
    mpfq_p_5_init(k, &t);
    mp_limb_t m;
    for(;;) {
        mpfq_p_5_set(k, t, b);
        for(m=0; mpfq_p_5_cmp_ui(k, t, 1)!=0; m++)
            mpfq_p_5_sqr(k, t, t);
        assert(m<=r);
        
        if (m==0 || m==r)
            break;
        
        s = 1UL<<(r-m-1);
        r = m;
        
        mpfq_p_5_pow(k, t, y, &s, 1);
        mpfq_p_5_sqr(k, y, t);
        mpfq_p_5_mul(k, x, x, t);
        mpfq_p_5_mul(k, b, b, y);
    }
    mpfq_p_5_set(k, z, x);
    mpfq_p_5_clear(k, &t);
    mpfq_p_5_clear(k, &x);
    mpfq_p_5_clear(k, &y);
    mpfq_p_5_clear(k, &b);
    return (m==0);
}

/* *Mpfq::defaults::pow::code_for_powz, Mpfq::gfp::elt, Mpfq::gfp */
void mpfq_p_5_powz(mpfq_p_5_dst_field k, mpfq_p_5_dst_elt y, mpfq_p_5_src_elt x, mpz_srcptr z)
{
        if (mpz_sgn(z) < 0) {
            mpz_t mz;
            mpz_init(mz);
            mpz_neg(mz, z);
            mpfq_p_5_powz(k, y, x, mz);
            mpfq_p_5_inv(k, y, y);
            mpz_clear(mz);
        } else if (mpz_sizeinbase(z, 2) > mpfq_p_5_field_degree(k) * mpfq_p_5_field_characteristic_bits(k)) {
            mpz_t zr;
            mpz_init(zr);
            mpz_t ppz;
            mpz_init(ppz);
            mpfq_p_5_field_characteristic(k, ppz);
            mpz_pow_ui(ppz,ppz,mpfq_p_5_field_degree(k));
            mpz_sub_ui(ppz,ppz,1);
            mpz_fdiv_r(zr, z, ppz);
            mpfq_p_5_powz(k, y, x, zr);
            mpz_clear(ppz);
            mpz_clear(zr);
        } else {
            mpfq_p_5_pow(k, y, x, z->_mp_d, mpz_size(z));
        }
}


/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */
/* *Mpfq::gfp::io::code_for_asprint, Mpfq::gfp */
int mpfq_p_5_asprint(mpfq_p_5_dst_field k, char * * pstr, mpfq_p_5_src_elt x)
{
    int i, n;
    mp_size_t size_x;
    i=5-1;
    while ((i>=0)&&(x[i]==0)) {
        i--;
    }
    i++;
    size_x=i;
    mp_limb_t y[size_x];
    for (i = 0; i<size_x; ++i) {
        y[i]=x[i];
    }
    // allocate enough room for base 2 conversion.
    // mpn_get_str may produce one extra byte
    *pstr = (char *)mpfq_malloc_check(size_x * 64 + 2);
    n = mpn_get_str((unsigned char*)(*pstr), k->io_base, (mp_limb_t *) y, size_x);
    for (i = 0; i < n; ++i)
        (*pstr)[i] += '0';
    (*pstr)[n] = '\0';
    // Remove leading 0s
    /* Note that gmp source says: There are no leading zeros on the digits
     * generated at str, but that's not currently a documented feature.
     * This implies that we won't do much here... */
    int shift = 0;
    while (((*pstr)[shift] == '0') && ((*pstr)[shift+1] != '\0')) 
        shift++;
    if (shift>0) {
        memmove(*pstr, (*pstr) + shift, n + 1 - shift);
        n -= shift;
    }
    // Return '0' instead of empty string for zero element
    if ((*pstr)[0] == '\0') {
        (*pstr)[0] = '0';
        (*pstr)[1] = '\0';
        n = 1;
    }
    return n;
}

/* *Mpfq::defaults::code_for_fprint, Mpfq::gfp */
int mpfq_p_5_fprint(mpfq_p_5_dst_field k, FILE * file, mpfq_p_5_src_elt x)
{
    char *str;
    int rc;
    mpfq_p_5_asprint(k,&str,x);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *Mpfq::gfp::io::code_for_sscan, Mpfq::gfp */
int mpfq_p_5_sscan(mpfq_p_5_dst_field k, mpfq_p_5_dst_elt z, const char * str)
{
    mpz_t zz;
    mpz_init(zz);
    int nread;
    if (gmp_sscanf(str, "%Zd%n", zz, &nread) != 1) {
        mpz_clear(zz);
        return 0;
    }
    mpfq_p_5_set_mpz(k, z, zz);
    mpz_clear(zz);
    return nread;
}

/* *Mpfq::gfp::io::code_for_fscan, Mpfq::gfp */
int mpfq_p_5_fscan(mpfq_p_5_dst_field k, FILE * file, mpfq_p_5_dst_elt z)
{
    char *tmp;
    int allocated, len=0;
    int c, start=0;
    allocated=100;
    tmp = (char *)mpfq_malloc_check(allocated);
    for(;;) {
        c = fgetc(file);
        if (c==EOF)
            break;
        if (isspace((int)(unsigned char)c)) {
            if (start==0)
                continue;
            else
                break;
        } else {
            if (len==allocated) {
                allocated+=100 + allocated / 4;
                tmp = (char*)realloc(tmp, allocated);
            }
            tmp[len]=c;
            len++;
            start=1;
        }
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated);
    }
    tmp[len]='\0';
    int ret=mpfq_p_5_sscan(k,z,tmp);
    free(tmp);
    return ret ? len : 0;
}


/* Vector functions */
/* *Mpfq::defaults::vec::alloc::code_for_vec_init, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_init(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec * v, unsigned long n)
{
    unsigned int i;
    *v = (mpfq_p_5_vec) malloc (n*sizeof(mpfq_p_5_elt));
    for(i = 0; i < n; i++)
        mpfq_p_5_init(K, (*v) + i);
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_reinit, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_reinit(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec * v, unsigned long n, unsigned long m)
{
    if (n < m) { // increase size
        unsigned int i;
        *v = (mpfq_p_5_vec) realloc (*v, m * sizeof(mpfq_p_5_elt));
        for(i = n; i < m; i+=1)
            mpfq_p_5_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_p_5_clear(K, (*v) + i);
        *v = (mpfq_p_5_vec) realloc (*v, m * sizeof(mpfq_p_5_elt));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_clear, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_clear(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec * v, unsigned long n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_p_5_clear(K, (*v) + i);
    free(*v);
}

/* *Mpfq::defaults::vec::generic::code_for_vec_random, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_random(mpfq_p_5_dst_field k, mpfq_p_5_dst_vec w, unsigned long n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_dst_elt x = mpfq_p_5_vec_coeff_ptr(k, w, i);
        mpfq_p_5_random(k, x, state);
        }
}

/* *Mpfq::defaults::vec::generic::code_for_vec_random2, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_random2(mpfq_p_5_dst_field k, mpfq_p_5_dst_vec w, unsigned long n, gmp_randstate_t state)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_dst_elt x = mpfq_p_5_vec_coeff_ptr(k, w, i);
        mpfq_p_5_random2(k, x, state);
        }
}

/* *Mpfq::defaults::vec::generic::code_for_vec_cmp, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_cmp(mpfq_p_5_dst_field k, mpfq_p_5_src_vec u, mpfq_p_5_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, u, i);
            mpfq_p_5_src_elt y = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        int ret = mpfq_p_5_cmp(k, x, y);
        if (ret != 0)
            return ret;
        }
        return 0;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_is_zero, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_is_zero(mpfq_p_5_dst_field k, mpfq_p_5_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_p_5_is_zero(k, x))
            return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::io::code_for_vec_asprint, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_asprint(mpfq_p_5_dst_field K MAYBE_UNUSED, char * * pstr, mpfq_p_5_src_vec w, unsigned long n)
{
    if (n == 0) {
        *pstr = (char *)mpfq_malloc_check(4);
        sprintf(*pstr, "[ ]");
        return strlen(*pstr);
    }
    int alloc = 100;
    int len = 0;
    *pstr = (char *)mpfq_malloc_check(alloc);
    (*pstr)[len++] = '[';
    (*pstr)[len++] = ' ';
    unsigned int i;
    for(i = 0; i < n; i+=1) {
        if (i) {
            (*pstr)[len++] = ',';
            (*pstr)[len++] = ' ';
        }
        char *tmp;
        mpfq_p_5_asprint(K, &tmp, w[i]);
        int ltmp = strlen(tmp);
        if (len+ltmp+4 > alloc) {
            alloc = len+ltmp+100 + alloc / 4;
            *pstr = (char *)realloc(*pstr, alloc);
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

/* *Mpfq::defaults::vec::io::code_for_vec_fprint, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_fprint(mpfq_p_5_dst_field K MAYBE_UNUSED, FILE * file, mpfq_p_5_src_vec w, unsigned long n)
{
    char *str;
    int rc;
    mpfq_p_5_vec_asprint(K,&str,w,n);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *Mpfq::defaults::vec::io::code_for_vec_print, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_print(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_src_vec w, unsigned long n)
{
    return mpfq_p_5_vec_fprint(K,stdout,w,n);
}

/* *Mpfq::defaults::vec::io::code_for_vec_sscan, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_sscan(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec * w, unsigned long * n, const char * str)
{
    // start with a clean vector
    unsigned int nn;
    int len = 0;
    mpfq_p_5_vec_reinit(K, w, *n, 0);
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
            mpfq_p_5_vec_reinit(K, w, nn, i+1);
            *n = nn = i+1;
        }
        int ret = mpfq_p_5_sscan(K, mpfq_p_5_vec_coeff_ptr(K, *w, i), str + len);
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

/* *Mpfq::defaults::vec::io::code_for_vec_fscan, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_fscan(mpfq_p_5_dst_field K MAYBE_UNUSED, FILE * file, mpfq_p_5_vec * w, unsigned long * n)
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
    int ret=mpfq_p_5_vec_sscan(K,w,n,tmp);
    free(tmp);
    return ret;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_hamming_weight, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_hamming_weight(mpfq_p_5_dst_field k, mpfq_p_5_src_vec v, unsigned long n)
{
        int w = 0;
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        w += !mpfq_p_5_is_zero(k, x);
        }
        return w;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_find_first_set, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_find_first_set(mpfq_p_5_dst_field k, mpfq_p_5_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_p_5_is_zero(k, x))
            return i;
        }
        return -1;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_simd_hamming_weight, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_simd_hamming_weight(mpfq_p_5_dst_field k, mpfq_p_5_src_vec v, unsigned long n)
{
        int w = 0;
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        w += mpfq_p_5_simd_hamming_weight(k, x);
        }
        return w;
}

/* *Mpfq::defaults::vec::generic::code_for_vec_simd_find_first_set, Mpfq::defaults::vec, Mpfq::gfp */
long mpfq_p_5_vec_simd_find_first_set(mpfq_p_5_dst_field k, mpfq_p_5_src_vec v, unsigned long n)
{
        for (unsigned int i = 0; i < n; ++i) {
            mpfq_p_5_src_elt x = mpfq_p_5_vec_coeff_ptr_const(k, v, i);
        if (!mpfq_p_5_is_zero(k, x))
            return i * mpfq_p_5_simd_groupsize(k) + mpfq_p_5_simd_hamming_weight(k, x);
        }
        return -1;
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_init, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_ur_init(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec_ur * v, unsigned long n)
{
    unsigned int i;
    *v = (mpfq_p_5_vec_ur) malloc (n*sizeof(mpfq_p_5_elt_ur));
    for(i = 0; i < n; i+=1)
        mpfq_p_5_elt_ur_init(K, &( (*v)[i]));
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_reinit, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_ur_reinit(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec_ur * v, unsigned long n, unsigned long m)
{
    if (n < m) { // increase size
        *v = (mpfq_p_5_vec_ur) realloc (*v, m * sizeof(mpfq_p_5_elt_ur));
        unsigned int i;
        for(i = n; i < m; i+=1)
            mpfq_p_5_elt_ur_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_p_5_elt_ur_clear(K, (*v) + i);
        *v = (mpfq_p_5_vec_ur) realloc (*v, m * sizeof(mpfq_p_5_elt_ur));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_clear, Mpfq::defaults::vec, Mpfq::gfp */
void mpfq_p_5_vec_ur_clear(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_vec_ur * v, unsigned long n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_p_5_elt_ur_clear(K, &( (*v)[i]));
    free(*v);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
/* Triggered by: vec_conv_ur */
void mpfq_p_5_vec_conv_ur_ks(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w, mpfq_p_5_src_vec u, unsigned int n, mpfq_p_5_src_vec v, unsigned int m)
{
    // compute base as a power 2^GMP_NUMB_BITS
    // This is the least number of words that can accomodate
    //     log_2( (p-1)^2 * min(n,m) )
    mpz_t p;
    mpz_init(p);
    mpfq_p_5_field_characteristic(K, p);
    mpz_sub_ui(p, p, 1);
    mpz_mul(p, p, p);
    mpz_mul_ui(p, p, MIN(m, n));
    
    long nbits = mpz_sizeinbase(p, 2);
    long nwords = 1 + ((nbits-1) / GMP_NUMB_BITS);
    nbits = GMP_NUMB_BITS*nwords;
    mpz_clear(p);
    
    assert(sizeof(mpfq_p_5_elt_ur) >= nwords*sizeof(unsigned long));
    
    // Create big integers
    mpz_t U, V;
    mpz_init2(U, n*nbits);
    mpz_init2(V, m*nbits);
    memset(U->_mp_d, 0, n*nwords*sizeof(unsigned long));
    memset(V->_mp_d, 0, m*nwords*sizeof(unsigned long));
    unsigned int i;
    assert (U->_mp_alloc == (int) (n*nwords));
    for (i = 0; i < n; ++i)
        mpfq_p_5_get_mpn(K, U->_mp_d + i*nwords, u[i]);
    U->_mp_size = U->_mp_alloc;
    // TODO: in principle one could reduce _mp_size until its true value, 
    // but then one should take care of W->_mp_size as well...
    //while (U->_mp_size > 0 && U->_mp_d[U->_mp_size-1] == 0)
    //    U->_mp_size--;
    assert (V->_mp_alloc == (int) (m*nwords));
    for (i = 0; i < m; ++i)
        mpfq_p_5_get_mpn(K, V->_mp_d + i*nwords, v[i]);
    V->_mp_size = V->_mp_alloc;
    //while (V->_mp_size > 0 && V->_mp_d[V->_mp_size-1] == 0)
    //    V->_mp_size--;
    
    // Multiply
    mpz_t W;
    mpz_init(W);
    mpz_mul(W, U, V);
    mpz_clear(U);
    mpz_clear(V);
    
    // Put coefficients in w
    assert (W->_mp_size >= (int) ((m+n-1)*nwords));
    if (sizeof(mpfq_p_5_elt_ur) == nwords*sizeof(unsigned long)) {
        for (i = 0; i < m+n-1; ++i) 
            mpfq_p_5_elt_ur_set(K, w[i], (mpfq_p_5_src_elt_ur)(W->_mp_d + i*nwords));
    } else {
        for (i = 0; i < m+n-1; ++i) {
            mpfq_p_5_elt_ur_set_ui(K, w[i], 0);
            memcpy(w[i], W->_mp_d + i*nwords, nwords*sizeof(unsigned long));
        }
    }
    
    mpz_clear(W);
}


/* Polynomial functions */
/* *Mpfq::defaults::poly::code_for_poly_setmonic, Mpfq::gfp */
void mpfq_p_5_poly_setmonic(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_poly q, mpfq_p_5_src_poly p)
{
    long degp = mpfq_p_5_poly_deg(K, p);
    if (degp == -1) {
        q->size = 0;
        return;
    }
    if (degp == 0) {
        mpfq_p_5_elt aux;
        mpfq_p_5_init(K, &aux);
        mpfq_p_5_set_ui(K, aux, 1);
        mpfq_p_5_poly_setcoeff(K, q, aux, 0);
        mpfq_p_5_clear(K, &aux);
        q->size = 1;
        return;
    }
    mpfq_p_5_elt lc;	/* spurious uninit warning sometimes */
    mpfq_p_5_init(K, &lc);
    mpfq_p_5_poly_getcoeff(K, lc, p, degp);
    mpfq_p_5_inv(K, lc, lc);
    mpfq_p_5_poly_setcoeff_ui(K, q, 1, degp);
    mpfq_p_5_vec_scal_mul(K, q->c, p->c, lc, degp);
    q->size = degp+1;
    mpfq_p_5_clear(K, &lc);
}

/* *Mpfq::defaults::poly::code_for_poly_divmod, Mpfq::gfp */
int mpfq_p_5_poly_divmod(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_poly q, mpfq_p_5_dst_poly r, mpfq_p_5_src_poly a, mpfq_p_5_src_poly b)
{
    if (b->size == 0) {
        return 0;
    }
    if (a->size == 0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    int dega = mpfq_p_5_poly_deg(K, a);
    if (dega<0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    // Compute deg b and inverse of leading coef
    int degb = mpfq_p_5_poly_deg(K, b);
    if (degb<0) {
        return 0;
    }
    if (degb > dega) {
        q->size=0;
        mpfq_p_5_poly_set(K, r, a);
        return 1;
    }
    int bmonic;
    mpfq_p_5_elt ilb;
    mpfq_p_5_init(K, &ilb);
    mpfq_p_5_elt temp;
    mpfq_p_5_init(K, &temp);
    mpfq_p_5_set_zero(K, temp); // silence spurious gcc11 warning :-(
    mpfq_p_5_poly_getcoeff(K, temp, b, degb);
    if (mpfq_p_5_cmp_ui(K, temp, 1) == 0) {
        mpfq_p_5_set_ui(K, ilb, 1);
        bmonic = 1;
    } else {
        mpfq_p_5_inv(K, ilb, temp);
        bmonic = 0;
    }
    
    mpfq_p_5_poly qq, rr;
    mpfq_p_5_poly_init(K, qq, dega-degb+1);
    mpfq_p_5_poly_init(K, rr, dega);
    
    mpfq_p_5_poly_set(K, rr, a);
    mpfq_p_5_elt aux, aux2;
    
    mpfq_p_5_init(K, &aux);
    mpfq_p_5_init(K, &aux2);
    
    int i;
    int j;
    for (i = dega; i >= (int)degb; --i) {
        mpfq_p_5_poly_getcoeff(K, aux, rr, i);
        if (!bmonic) 
            mpfq_p_5_mul(K, aux, aux, ilb);
        mpfq_p_5_poly_setcoeff(K, qq, aux, i-degb);
        for (j = i-1; j >= (int)(i - degb); --j) {
            mpfq_p_5_poly_getcoeff(K, temp, b, j-i+degb);
            mpfq_p_5_mul(K, aux2, aux, temp);
            mpfq_p_5_poly_getcoeff(K, temp, rr, j);
    
            mpfq_p_5_sub(K, temp, temp, aux2);
            mpfq_p_5_poly_setcoeff(K, rr, temp, j);
        }
    }    
    
    rr->size = degb;
    int degr = mpfq_p_5_poly_deg(K, rr);
    rr->size = degr+1;
    
    if (q != NULL) 
        mpfq_p_5_poly_set(K, q, qq);
    if (r != NULL)
        mpfq_p_5_poly_set(K, r, rr);
    mpfq_p_5_clear(K, &temp);
    mpfq_p_5_clear(K, &ilb);
    mpfq_p_5_clear(K, &aux);
    mpfq_p_5_clear(K, &aux2);
    mpfq_p_5_poly_clear(K, rr);
    mpfq_p_5_poly_clear(K, qq);
    return 1;
}

static void mpfq_p_5_poly_preinv(mpfq_p_5_dst_field, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, unsigned int);
/* *Mpfq::defaults::poly::code_for_poly_precomp_mod, Mpfq::gfp */
/* Triggered by: poly_precomp_mod */
static void mpfq_p_5_poly_preinv(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_poly q, mpfq_p_5_src_poly p, unsigned int n)
{
    // Compute the inverse of p(x) modulo x^n
    // Newton iteration: x_{n+1} = x_n + x_n(1 - a*x_n)
    // Requires p(0) = 1
    // Assume p != q (no alias)
    assert (p != q);
    int m;
    if (n <= 2) {
        mpfq_p_5_poly_setcoeff_ui(K, q, 1, 0);
        q->size = 1;
        m = 1;
        if (n == 1)
            return;
    } else {
        // n >= 3: recursive call at prec m = ceil(n/2)
        m = 1 + ((n-1)/2);
        mpfq_p_5_poly_preinv(K, q, p, m);
    }
    mpfq_p_5_elt temp;	/* spurious uninit warning sometimes */
    mpfq_p_5_init(K, &temp);
    mpfq_p_5_set_zero(K, temp); // silence spurious gcc11 warning :-(
    mpfq_p_5_poly_getcoeff(K, temp, p, 0);//Should be in the assert
    assert( mpfq_p_5_cmp_ui(K, temp, 1) == 0);
    // enlarge q if necessary
    if (q->alloc < n) {
        mpfq_p_5_vec_reinit(K, &(q->c), q->alloc, n);
        q->alloc = n;
    }
    // refine value
    mpfq_p_5_vec tmp;
    mpfq_p_5_vec_init(K, &tmp, m+n-1);
    
    mpfq_p_5_vec_conv(K, tmp, p->c, MIN(n, p->size), q->c, m);
    int nn = MIN(n, MIN(n, p->size) + m -1);
    mpfq_p_5_vec_neg(K, tmp, tmp, nn);
    mpfq_p_5_vec_getcoeff(K, temp, tmp, 0);
    mpfq_p_5_add_ui(K, temp, temp, 1);
    mpfq_p_5_vec_setcoeff(K, tmp, temp, 0);
    mpfq_p_5_vec_conv(K, tmp, q->c, m, tmp, nn);
    mpfq_p_5_vec_set(K, mpfq_p_5_vec_subvec(K, q->c, m), mpfq_p_5_vec_subvec(K, tmp, m), n-m);
    q->size = n;
    
    mpfq_p_5_clear(K, &temp);
    mpfq_p_5_vec_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::poly::code_for_poly_precomp_mod, Mpfq::gfp */
void mpfq_p_5_poly_precomp_mod(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_poly q, mpfq_p_5_src_poly p)
{
    assert(p != q);
    int N = mpfq_p_5_poly_deg(K, p);
    if (N < 0) {
        mpfq_p_5_poly_set(K, q, p);
        return;
    }
    mpfq_p_5_poly rp;
    mpfq_p_5_poly_init(K, rp, N+1);
    mpfq_p_5_vec_rev(K, rp->c, p->c, N+1);
    rp->size = N+1;
    mpfq_p_5_poly_preinv(K, q, rp, N);
    mpfq_p_5_poly_clear(K, rp);
}

/* *Mpfq::defaults::poly::code_for_poly_mod_pre, Mpfq::gfp */
void mpfq_p_5_poly_mod_pre(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_poly r, mpfq_p_5_src_poly q, mpfq_p_5_src_poly p, mpfq_p_5_src_poly irp)
{
    int N = mpfq_p_5_poly_deg(K, p);
    int degq = mpfq_p_5_poly_deg(K, q);
    if (degq < N) {
        mpfq_p_5_poly_set(K, r, q);
        return;
    }
    int m = degq - N;
    assert (degq <= 2*N-2);
    mpfq_p_5_poly revq;
    mpfq_p_5_poly_init(K, revq, MAX(degq+1, m+1));
    mpfq_p_5_vec_rev(K, revq->c, q->c, degq+1);
    revq->size = q->size;
    mpfq_p_5_poly_mul(K, revq, revq, irp);
    mpfq_p_5_vec_rev(K, revq->c, revq->c, m+1);
    revq->size = m+1;
    
    mpfq_p_5_poly_mul(K, revq, revq, p);
    mpfq_p_5_poly_sub(K, r, q, revq);
    r->size = mpfq_p_5_poly_deg(K, r)+1;
    mpfq_p_5_poly_clear(K, revq);
}


/* Functions related to SIMD operation */
/* *simd_gfp::code_for_add_dotprod */
void mpfq_p_5_add_dotprod(mpfq_p_5_dst_field K MAYBE_UNUSED, mpfq_p_5_dst_vec xw, mpfq_p_5_src_vec xu1, mpfq_p_5_src_vec xu0, unsigned int n)
{
        mpfq_p_5_elt_ur s,t;
        mpfq_p_5_elt_ur_init(K, &s);
        mpfq_p_5_elt_ur_init(K, &t);
        mpfq_p_5_elt_ur_set_elt(K, s, xw[0]);
        for(unsigned int i = 0 ; i < n ; i++) {
            mpfq_p_5_mul_ur(K, t, xu0[i], xu1[i]);
            mpfq_p_5_elt_ur_add(K, s, s, t);
        }
        mpfq_p_5_reduce(K, xw[0], s);
        mpfq_p_5_elt_ur_clear(K, &s);
        mpfq_p_5_elt_ur_clear(K, &t);
}


/* Member templates related to SIMD operation */

/* Object-oriented interface */
static void mpfq_p_5_wrapper_oo_field_clear(mpfq_vbase_ptr);
static void mpfq_p_5_wrapper_oo_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_p_5_oo_field_clear(vbase);
}

static void mpfq_p_5_wrapper_oo_field_init(mpfq_vbase_ptr);
static void mpfq_p_5_wrapper_oo_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_p_5_oo_field_init(vbase);
}

static void mpfq_p_5_wrapper_add_dotprod(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, mpfq_p_5_src_vec, unsigned int);
static void mpfq_p_5_wrapper_add_dotprod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec xw MAYBE_UNUSED, mpfq_p_5_src_vec xu1 MAYBE_UNUSED, mpfq_p_5_src_vec xu0 MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_p_5_add_dotprod(vbase->obj, xw, xu1, xu0, n);
}

static void mpfq_p_5_wrapper_simd_set_ui_all(mpfq_vbase_ptr, mpfq_p_5_dst_elt, unsigned long);
static void mpfq_p_5_wrapper_simd_set_ui_all(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt p MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_p_5_simd_set_ui_all(vbase->obj, p, v);
}

static void mpfq_p_5_wrapper_simd_add_ui_at(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, int, unsigned long);
static void mpfq_p_5_wrapper_simd_add_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt p MAYBE_UNUSED, mpfq_p_5_src_elt p0 MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_p_5_simd_add_ui_at(vbase->obj, p, p0, k, v);
}

static void mpfq_p_5_wrapper_simd_set_ui_at(mpfq_vbase_ptr, mpfq_p_5_dst_elt, int, unsigned long);
static void mpfq_p_5_wrapper_simd_set_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    mpfq_p_5_simd_set_ui_at(vbase->obj, p, k, v);
}

static unsigned long mpfq_p_5_wrapper_simd_get_ui_at(mpfq_vbase_ptr, mpfq_p_5_src_elt, int);
static unsigned long mpfq_p_5_wrapper_simd_get_ui_at(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt p MAYBE_UNUSED, int k MAYBE_UNUSED)
{
    return mpfq_p_5_simd_get_ui_at(vbase->obj, p, k);
}

static int mpfq_p_5_wrapper_simd_find_first_set(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_simd_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt p MAYBE_UNUSED)
{
    return mpfq_p_5_simd_find_first_set(vbase->obj, p);
}

static int mpfq_p_5_wrapper_simd_hamming_weight(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_simd_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt p MAYBE_UNUSED)
{
    return mpfq_p_5_simd_hamming_weight(vbase->obj, p);
}

static int mpfq_p_5_wrapper_simd_groupsize(mpfq_vbase_ptr);
static int mpfq_p_5_wrapper_simd_groupsize(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_simd_groupsize(vbase->obj);
}

static int mpfq_p_5_wrapper_poly_scan(mpfq_vbase_ptr, mpfq_p_5_dst_poly);
static int mpfq_p_5_wrapper_poly_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_scan(vbase->obj, w);
}

static int mpfq_p_5_wrapper_poly_fscan(mpfq_vbase_ptr, FILE *, mpfq_p_5_dst_poly);
static int mpfq_p_5_wrapper_poly_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_fscan(vbase->obj, file, w);
}

static int mpfq_p_5_wrapper_poly_sscan(mpfq_vbase_ptr, mpfq_p_5_dst_poly, const char *);
static int mpfq_p_5_wrapper_poly_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_p_5_poly_sscan(vbase->obj, w, str);
}

static int mpfq_p_5_wrapper_poly_print(mpfq_vbase_ptr, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_print(vbase->obj, w);
}

static int mpfq_p_5_wrapper_poly_fprint(mpfq_vbase_ptr, FILE *, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_src_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_fprint(vbase->obj, file, w);
}

static int mpfq_p_5_wrapper_poly_asprint(mpfq_vbase_ptr, char * *, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_p_5_src_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_asprint(vbase->obj, pstr, w);
}

static int mpfq_p_5_wrapper_poly_cmp(mpfq_vbase_ptr, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, mpfq_p_5_src_poly v MAYBE_UNUSED)
{
    return mpfq_p_5_poly_cmp(vbase->obj, u, v);
}

static void mpfq_p_5_wrapper_poly_random2(mpfq_vbase_ptr, mpfq_p_5_dst_poly, unsigned int, gmp_randstate_t);
static void mpfq_p_5_wrapper_poly_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_poly_random2(vbase->obj, w, n, state);
}

static void mpfq_p_5_wrapper_poly_random(mpfq_vbase_ptr, mpfq_p_5_dst_poly, unsigned int, gmp_randstate_t);
static void mpfq_p_5_wrapper_poly_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_poly_random(vbase->obj, w, n, state);
}

static void mpfq_p_5_wrapper_poly_xgcd(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_dst_poly, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_xgcd(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly g MAYBE_UNUSED, mpfq_p_5_dst_poly u0 MAYBE_UNUSED, mpfq_p_5_dst_poly v0 MAYBE_UNUSED, mpfq_p_5_src_poly a0 MAYBE_UNUSED, mpfq_p_5_src_poly b0 MAYBE_UNUSED)
{
    mpfq_p_5_poly_xgcd(vbase->obj, g, u0, v0, a0, b0);
}

static void mpfq_p_5_wrapper_poly_gcd(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_gcd(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly g MAYBE_UNUSED, mpfq_p_5_src_poly a0 MAYBE_UNUSED, mpfq_p_5_src_poly b0 MAYBE_UNUSED)
{
    mpfq_p_5_poly_gcd(vbase->obj, g, a0, b0);
}

static void mpfq_p_5_wrapper_poly_mod_pre(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_mod_pre(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly r MAYBE_UNUSED, mpfq_p_5_src_poly q MAYBE_UNUSED, mpfq_p_5_src_poly p MAYBE_UNUSED, mpfq_p_5_src_poly irp MAYBE_UNUSED)
{
    mpfq_p_5_poly_mod_pre(vbase->obj, r, q, p, irp);
}

static void mpfq_p_5_wrapper_poly_precomp_mod(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_precomp_mod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly q MAYBE_UNUSED, mpfq_p_5_src_poly p MAYBE_UNUSED)
{
    mpfq_p_5_poly_precomp_mod(vbase->obj, q, p);
}

static int mpfq_p_5_wrapper_poly_divmod(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_divmod(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly q MAYBE_UNUSED, mpfq_p_5_dst_poly r MAYBE_UNUSED, mpfq_p_5_src_poly a MAYBE_UNUSED, mpfq_p_5_src_poly b MAYBE_UNUSED)
{
    return mpfq_p_5_poly_divmod(vbase->obj, q, r, a, b);
}

static void mpfq_p_5_wrapper_poly_mul(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, mpfq_p_5_src_poly v MAYBE_UNUSED)
{
    mpfq_p_5_poly_mul(vbase->obj, w, u, v);
}

static void mpfq_p_5_wrapper_poly_scal_mul(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_poly_scal_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    mpfq_p_5_poly_scal_mul(vbase->obj, w, u, x);
}

static void mpfq_p_5_wrapper_poly_neg(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED)
{
    mpfq_p_5_poly_neg(vbase->obj, w, u);
}

static void mpfq_p_5_wrapper_poly_sub_ui(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, unsigned long);
static void mpfq_p_5_wrapper_poly_sub_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_p_5_poly_sub_ui(vbase->obj, w, u, x);
}

static void mpfq_p_5_wrapper_poly_add_ui(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, unsigned long);
static void mpfq_p_5_wrapper_poly_add_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_p_5_poly_add_ui(vbase->obj, w, u, x);
}

static void mpfq_p_5_wrapper_poly_set_ui(mpfq_vbase_ptr, mpfq_p_5_dst_poly, unsigned long);
static void mpfq_p_5_wrapper_poly_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_p_5_poly_set_ui(vbase->obj, w, x);
}

static void mpfq_p_5_wrapper_poly_sub(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, mpfq_p_5_src_poly v MAYBE_UNUSED)
{
    mpfq_p_5_poly_sub(vbase->obj, w, u, v);
}

static void mpfq_p_5_wrapper_poly_add(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED, mpfq_p_5_src_poly v MAYBE_UNUSED)
{
    mpfq_p_5_poly_add(vbase->obj, w, u, v);
}

static int mpfq_p_5_wrapper_poly_deg(mpfq_vbase_ptr, mpfq_p_5_src_poly);
static int mpfq_p_5_wrapper_poly_deg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_poly w MAYBE_UNUSED)
{
    return mpfq_p_5_poly_deg(vbase->obj, w);
}

static void mpfq_p_5_wrapper_poly_getcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_poly, unsigned int);
static void mpfq_p_5_wrapper_poly_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, mpfq_p_5_src_poly w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_p_5_poly_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_p_5_wrapper_poly_setcoeff_ui(mpfq_vbase_ptr, mpfq_p_5_dst_poly, unsigned long, unsigned int);
static void mpfq_p_5_wrapper_poly_setcoeff_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_p_5_poly_setcoeff_ui(vbase->obj, w, x, i);
}

static void mpfq_p_5_wrapper_poly_setcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_elt, unsigned int);
static void mpfq_p_5_wrapper_poly_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    mpfq_p_5_poly_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_p_5_wrapper_poly_setmonic(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_setmonic(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly q MAYBE_UNUSED, mpfq_p_5_src_poly p MAYBE_UNUSED)
{
    mpfq_p_5_poly_setmonic(vbase->obj, q, p);
}

static void mpfq_p_5_wrapper_poly_set(mpfq_vbase_ptr, mpfq_p_5_dst_poly, mpfq_p_5_src_poly);
static void mpfq_p_5_wrapper_poly_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_poly w MAYBE_UNUSED, mpfq_p_5_src_poly u MAYBE_UNUSED)
{
    mpfq_p_5_poly_set(vbase->obj, w, u);
}

static void mpfq_p_5_wrapper_poly_clear(mpfq_vbase_ptr, mpfq_p_5_poly);
static void mpfq_p_5_wrapper_poly_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_poly p MAYBE_UNUSED)
{
    mpfq_p_5_poly_clear(vbase->obj, p);
}

static void mpfq_p_5_wrapper_poly_init(mpfq_vbase_ptr, mpfq_p_5_poly, unsigned int);
static void mpfq_p_5_wrapper_poly_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_poly p MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    mpfq_p_5_poly_init(vbase->obj, p, n);
}

static ptrdiff_t mpfq_p_5_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr, long);
static ptrdiff_t mpfq_p_5_wrapper_vec_ur_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_ur_elt_stride(vbase->obj, n);
}

static ptrdiff_t mpfq_p_5_wrapper_vec_elt_stride(mpfq_vbase_ptr, long);
static ptrdiff_t mpfq_p_5_wrapper_vec_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED, long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_elt_stride(vbase->obj, n);
}

static mpfq_p_5_src_elt mpfq_p_5_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr, mpfq_p_5_src_vec_ur, long);
static mpfq_p_5_src_elt mpfq_p_5_wrapper_vec_ur_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_ur_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_p_5_dst_elt mpfq_p_5_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, long);
static mpfq_p_5_dst_elt mpfq_p_5_wrapper_vec_ur_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_ur_coeff_ptr(vbase->obj, v, i);
}

static mpfq_p_5_src_vec_ur mpfq_p_5_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr, mpfq_p_5_src_vec_ur, long);
static mpfq_p_5_src_vec_ur mpfq_p_5_wrapper_vec_ur_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_ur_subvec_const(vbase->obj, v, i);
}

static mpfq_p_5_dst_vec_ur mpfq_p_5_wrapper_vec_ur_subvec(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, long);
static mpfq_p_5_dst_vec_ur mpfq_p_5_wrapper_vec_ur_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_ur_subvec(vbase->obj, v, i);
}

static void mpfq_p_5_wrapper_vec_reduce(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_dst_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_reduce(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_dst_vec_ur u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_reduce(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_conv_ur(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec, unsigned long, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_conv_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_p_5_vec_conv_ur(vbase->obj, w, u, n, v, m);
}

static void mpfq_p_5_wrapper_vec_scal_mul_ur(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_vec_scal_mul_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, mpfq_p_5_src_elt c MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_scal_mul_ur(vbase->obj, w, u, c, n);
}

static void mpfq_p_5_wrapper_vec_ur_rev(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec_ur u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_rev(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_ur_neg(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec_ur u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_neg(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_ur_sub(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec_ur u MAYBE_UNUSED, mpfq_p_5_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_sub(vbase->obj, w, u, v, n);
}

static void mpfq_p_5_wrapper_vec_ur_add(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec_ur u MAYBE_UNUSED, mpfq_p_5_src_vec_ur v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_add(vbase->obj, w, u, v, n);
}

static void mpfq_p_5_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur x MAYBE_UNUSED, mpfq_p_5_src_vec_ur w MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_p_5_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_elt_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_elt_ur x MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_p_5_wrapper_vec_ur_set(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur r MAYBE_UNUSED, mpfq_p_5_src_vec_ur s MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_set(vbase->obj, r, s, n);
}

static void mpfq_p_5_wrapper_vec_ur_clear(mpfq_vbase_ptr, mpfq_p_5_vec_ur *, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_clear(vbase->obj, v, n);
}

static void mpfq_p_5_wrapper_vec_ur_reinit(mpfq_vbase_ptr, mpfq_p_5_vec_ur *, unsigned long, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_reinit(vbase->obj, v, n, m);
}

static void mpfq_p_5_wrapper_vec_ur_set_vec(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_set_vec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_set_vec(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_ur_set_zero(mpfq_vbase_ptr, mpfq_p_5_dst_vec_ur, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec_ur r MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_set_zero(vbase->obj, r, n);
}

static void mpfq_p_5_wrapper_vec_ur_init(mpfq_vbase_ptr, mpfq_p_5_vec_ur *, unsigned long);
static void mpfq_p_5_wrapper_vec_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec_ur * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_ur_init(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_simd_find_first_set(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_simd_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_simd_find_first_set(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_simd_hamming_weight(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_simd_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_simd_hamming_weight(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_find_first_set(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_find_first_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_find_first_set(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_hamming_weight(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_hamming_weight(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_hamming_weight(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_scan(mpfq_vbase_ptr, mpfq_p_5_vec *, unsigned long *);
static long mpfq_p_5_wrapper_vec_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_scan(vbase->obj, w, n);
}

static long mpfq_p_5_wrapper_vec_fscan(mpfq_vbase_ptr, FILE *, mpfq_p_5_vec *, unsigned long *);
static long mpfq_p_5_wrapper_vec_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_fscan(vbase->obj, file, w, n);
}

static long mpfq_p_5_wrapper_vec_sscan(mpfq_vbase_ptr, mpfq_p_5_vec *, unsigned long *, const char *);
static long mpfq_p_5_wrapper_vec_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec * w MAYBE_UNUSED, unsigned long * n MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_p_5_vec_sscan(vbase->obj, w, n, str);
}

static long mpfq_p_5_wrapper_vec_print(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_print(vbase->obj, w, n);
}

static long mpfq_p_5_wrapper_vec_fprint(mpfq_vbase_ptr, FILE *, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_src_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_fprint(vbase->obj, file, w, n);
}

static long mpfq_p_5_wrapper_vec_asprint(mpfq_vbase_ptr, char * *, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_p_5_src_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_asprint(vbase->obj, pstr, w, n);
}

static mpfq_p_5_src_elt mpfq_p_5_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr, mpfq_p_5_src_vec, long);
static mpfq_p_5_src_elt mpfq_p_5_wrapper_vec_coeff_ptr_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_coeff_ptr_const(vbase->obj, v, i);
}

static mpfq_p_5_dst_elt mpfq_p_5_wrapper_vec_coeff_ptr(mpfq_vbase_ptr, mpfq_p_5_dst_vec, long);
static mpfq_p_5_dst_elt mpfq_p_5_wrapper_vec_coeff_ptr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_coeff_ptr(vbase->obj, v, i);
}

static mpfq_p_5_src_vec mpfq_p_5_wrapper_vec_subvec_const(mpfq_vbase_ptr, mpfq_p_5_src_vec, long);
static mpfq_p_5_src_vec mpfq_p_5_wrapper_vec_subvec_const(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_subvec_const(vbase->obj, v, i);
}

static mpfq_p_5_dst_vec mpfq_p_5_wrapper_vec_subvec(mpfq_vbase_ptr, mpfq_p_5_dst_vec, long);
static mpfq_p_5_dst_vec mpfq_p_5_wrapper_vec_subvec(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec v MAYBE_UNUSED, long i MAYBE_UNUSED)
{
    return mpfq_p_5_vec_subvec(vbase->obj, v, i);
}

static long mpfq_p_5_wrapper_vec_is_zero(mpfq_vbase_ptr, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_is_zero(vbase->obj, v, n);
}

static long mpfq_p_5_wrapper_vec_cmp(mpfq_vbase_ptr, mpfq_p_5_src_vec, mpfq_p_5_src_vec, unsigned long);
static long mpfq_p_5_wrapper_vec_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    return mpfq_p_5_vec_cmp(vbase->obj, u, v, n);
}

static void mpfq_p_5_wrapper_vec_random2(mpfq_vbase_ptr, mpfq_p_5_dst_vec, unsigned long, gmp_randstate_t);
static void mpfq_p_5_wrapper_vec_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_vec_random2(vbase->obj, w, n, state);
}

static void mpfq_p_5_wrapper_vec_random(mpfq_vbase_ptr, mpfq_p_5_dst_vec, unsigned long, gmp_randstate_t);
static void mpfq_p_5_wrapper_vec_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_vec_random(vbase->obj, w, n, state);
}

static void mpfq_p_5_wrapper_vec_conv(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, unsigned long, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_conv(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_p_5_vec_conv(vbase->obj, w, u, n, v, m);
}

static void mpfq_p_5_wrapper_vec_scal_mul(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_vec_scal_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, mpfq_p_5_src_elt c MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_scal_mul(vbase->obj, w, u, c, n);
}

static void mpfq_p_5_wrapper_vec_sub(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_sub(vbase->obj, w, u, v, n);
}

static void mpfq_p_5_wrapper_vec_rev(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_rev(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_rev(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_neg(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_neg(vbase->obj, w, u, n);
}

static void mpfq_p_5_wrapper_vec_add(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_vec u MAYBE_UNUSED, mpfq_p_5_src_vec v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_add(vbase->obj, w, u, v, n);
}

static void mpfq_p_5_wrapper_vec_getcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_getcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, mpfq_p_5_src_vec w MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_p_5_vec_getcoeff(vbase->obj, x, w, i);
}

static void mpfq_p_5_wrapper_vec_setcoeff_ui(mpfq_vbase_ptr, mpfq_p_5_dst_vec, unsigned long, unsigned long);
static void mpfq_p_5_wrapper_vec_setcoeff_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_p_5_vec_setcoeff_ui(vbase->obj, w, x, i);
}

static void mpfq_p_5_wrapper_vec_setcoeff(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_vec_setcoeff(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec w MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned long i MAYBE_UNUSED)
{
    mpfq_p_5_vec_setcoeff(vbase->obj, w, x, i);
}

static void mpfq_p_5_wrapper_vec_set_zero(mpfq_vbase_ptr, mpfq_p_5_dst_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec r MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_set_zero(vbase->obj, r, n);
}

static void mpfq_p_5_wrapper_vec_set(mpfq_vbase_ptr, mpfq_p_5_dst_vec, mpfq_p_5_src_vec, unsigned long);
static void mpfq_p_5_wrapper_vec_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_vec r MAYBE_UNUSED, mpfq_p_5_src_vec s MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_set(vbase->obj, r, s, n);
}

static void mpfq_p_5_wrapper_vec_clear(mpfq_vbase_ptr, mpfq_p_5_vec *, unsigned long);
static void mpfq_p_5_wrapper_vec_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_clear(vbase->obj, v, n);
}

static void mpfq_p_5_wrapper_vec_reinit(mpfq_vbase_ptr, mpfq_p_5_vec *, unsigned long, unsigned long);
static void mpfq_p_5_wrapper_vec_reinit(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED, unsigned long m MAYBE_UNUSED)
{
    mpfq_p_5_vec_reinit(vbase->obj, v, n, m);
}

static void mpfq_p_5_wrapper_vec_init(mpfq_vbase_ptr, mpfq_p_5_vec *, unsigned long);
static void mpfq_p_5_wrapper_vec_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_vec * v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
    mpfq_p_5_vec_init(vbase->obj, v, n);
}

static int mpfq_p_5_wrapper_scan(mpfq_vbase_ptr, mpfq_p_5_dst_elt);
static int mpfq_p_5_wrapper_scan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_scan(vbase->obj, x);
}

static int mpfq_p_5_wrapper_fscan(mpfq_vbase_ptr, FILE *, mpfq_p_5_dst_elt);
static int mpfq_p_5_wrapper_fscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED)
{
    return mpfq_p_5_fscan(vbase->obj, file, z);
}

static int mpfq_p_5_wrapper_sscan(mpfq_vbase_ptr, mpfq_p_5_dst_elt, const char *);
static int mpfq_p_5_wrapper_sscan(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return mpfq_p_5_sscan(vbase->obj, z, str);
}

static int mpfq_p_5_wrapper_print(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_print(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_print(vbase->obj, x);
}

static int mpfq_p_5_wrapper_fprint(mpfq_vbase_ptr, FILE *, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_fprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_fprint(vbase->obj, file, x);
}

static int mpfq_p_5_wrapper_asprint(mpfq_vbase_ptr, char * *, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_asprint(mpfq_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_asprint(vbase->obj, pstr, x);
}

static int mpfq_p_5_wrapper_is_zero(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_is_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt r MAYBE_UNUSED)
{
    return mpfq_p_5_is_zero(vbase->obj, r);
}

static int mpfq_p_5_wrapper_cmp_ui(mpfq_vbase_ptr, mpfq_p_5_src_elt, unsigned long);
static int mpfq_p_5_wrapper_cmp_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    return mpfq_p_5_cmp_ui(vbase->obj, x, y);
}

static int mpfq_p_5_wrapper_cmp(mpfq_vbase_ptr, mpfq_p_5_src_elt, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_cmp(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    return mpfq_p_5_cmp(vbase->obj, x, y);
}

static void mpfq_p_5_wrapper_addmul_si_ur(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt, long);
static void mpfq_p_5_wrapper_addmul_si_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur w MAYBE_UNUSED, mpfq_p_5_src_elt u MAYBE_UNUSED, long v MAYBE_UNUSED)
{
    mpfq_p_5_addmul_si_ur(vbase->obj, w, u, v);
}

static void mpfq_p_5_wrapper_reduce(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_dst_elt_ur);
static void mpfq_p_5_wrapper_reduce(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_dst_elt_ur x MAYBE_UNUSED)
{
    mpfq_p_5_reduce(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_sqr_ur(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_sqr_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    mpfq_p_5_sqr_ur(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_mul_ur(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_mul_ur(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_mul_ur(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_elt_ur_sub(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt_ur, mpfq_p_5_src_elt_ur);
static void mpfq_p_5_wrapper_elt_ur_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt_ur x MAYBE_UNUSED, mpfq_p_5_src_elt_ur y MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_sub(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_elt_ur_neg(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt_ur);
static void mpfq_p_5_wrapper_elt_ur_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt_ur x MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_neg(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_elt_ur_add(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt_ur, mpfq_p_5_src_elt_ur);
static void mpfq_p_5_wrapper_elt_ur_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt_ur x MAYBE_UNUSED, mpfq_p_5_src_elt_ur y MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_add(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_elt_ur_set_ui(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, unsigned long);
static void mpfq_p_5_wrapper_elt_ur_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur r MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_set_ui(vbase->obj, r, x);
}

static void mpfq_p_5_wrapper_elt_ur_set_zero(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur);
static void mpfq_p_5_wrapper_elt_ur_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur r MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_set_zero(vbase->obj, r);
}

static void mpfq_p_5_wrapper_elt_ur_set_elt(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_elt_ur_set_elt(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur r MAYBE_UNUSED, mpfq_p_5_src_elt s MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_set_elt(vbase->obj, r, s);
}

static void mpfq_p_5_wrapper_elt_ur_set(mpfq_vbase_ptr, mpfq_p_5_dst_elt_ur, mpfq_p_5_src_elt_ur);
static void mpfq_p_5_wrapper_elt_ur_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt_ur z MAYBE_UNUSED, mpfq_p_5_src_elt_ur x MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_set(vbase->obj, z, x);
}

static ptrdiff_t mpfq_p_5_wrapper_elt_ur_stride(mpfq_vbase_ptr);
static ptrdiff_t mpfq_p_5_wrapper_elt_ur_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_elt_ur_stride(vbase->obj);
}

static void mpfq_p_5_wrapper_elt_ur_clear(mpfq_vbase_ptr, mpfq_p_5_elt_ur *);
static void mpfq_p_5_wrapper_elt_ur_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_elt_ur * x MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_clear(vbase->obj, x);
}

static void mpfq_p_5_wrapper_elt_ur_init(mpfq_vbase_ptr, mpfq_p_5_elt_ur *);
static void mpfq_p_5_wrapper_elt_ur_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_elt_ur * x MAYBE_UNUSED)
{
    mpfq_p_5_elt_ur_init(vbase->obj, x);
}

static void mpfq_p_5_wrapper_hadamard(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_dst_elt, mpfq_p_5_dst_elt, mpfq_p_5_dst_elt);
static void mpfq_p_5_wrapper_hadamard(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, mpfq_p_5_dst_elt y MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_dst_elt t MAYBE_UNUSED)
{
    mpfq_p_5_hadamard(vbase->obj, x, y, z, t);
}

static int mpfq_p_5_wrapper_inv(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_inv(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_inv(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_normalize(mpfq_vbase_ptr, mpfq_p_5_dst_elt);
static void mpfq_p_5_wrapper_normalize(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED)
{
    mpfq_p_5_normalize(vbase->obj, x);
}

static void mpfq_p_5_wrapper_mul_ui(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_mul_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_p_5_mul_ui(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_sub_ui(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_sub_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_p_5_sub_ui(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_add_ui(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, unsigned long);
static void mpfq_p_5_wrapper_add_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, unsigned long y MAYBE_UNUSED)
{
    mpfq_p_5_add_ui(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_frobenius(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_frobenius(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_frobenius(vbase->obj, x, y);
}

static void mpfq_p_5_wrapper_powz(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, mpz_srcptr);
static void mpfq_p_5_wrapper_powz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt y MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpz_srcptr z MAYBE_UNUSED)
{
    mpfq_p_5_powz(vbase->obj, y, x, z);
}

static void mpfq_p_5_wrapper_pow(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, unsigned long *, size_t);
static void mpfq_p_5_wrapper_pow(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt res MAYBE_UNUSED, mpfq_p_5_src_elt r MAYBE_UNUSED, unsigned long * x MAYBE_UNUSED, size_t n MAYBE_UNUSED)
{
    mpfq_p_5_pow(vbase->obj, res, r, x, n);
}

static int mpfq_p_5_wrapper_sqrt(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_sqrt(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt a MAYBE_UNUSED)
{
    return mpfq_p_5_sqrt(vbase->obj, z, a);
}

static int mpfq_p_5_wrapper_is_sqr(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static int mpfq_p_5_wrapper_is_sqr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_is_sqr(vbase->obj, x);
}

static void mpfq_p_5_wrapper_sqr(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_sqr(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    mpfq_p_5_sqr(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_mul(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_mul(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_mul(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_neg(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_neg(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    mpfq_p_5_neg(vbase->obj, z, x);
}

static void mpfq_p_5_wrapper_sub(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_sub(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_sub(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_add(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_add(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt z MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_add(vbase->obj, z, x, y);
}

static void mpfq_p_5_wrapper_random2(mpfq_vbase_ptr, mpfq_p_5_dst_elt, gmp_randstate_t);
static void mpfq_p_5_wrapper_random2(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_random2(vbase->obj, x, state);
}

static void mpfq_p_5_wrapper_random(mpfq_vbase_ptr, mpfq_p_5_dst_elt, gmp_randstate_t);
static void mpfq_p_5_wrapper_random(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt x MAYBE_UNUSED, gmp_randstate_t state MAYBE_UNUSED)
{
    mpfq_p_5_random(vbase->obj, x, state);
}

static void mpfq_p_5_wrapper_get_mpz(mpfq_vbase_ptr, mpz_ptr, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_get_mpz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpz_ptr z MAYBE_UNUSED, mpfq_p_5_src_elt y MAYBE_UNUSED)
{
    mpfq_p_5_get_mpz(vbase->obj, z, y);
}

static void mpfq_p_5_wrapper_get_mpn(mpfq_vbase_ptr, mp_limb_t *, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_get_mpn(mpfq_vbase_ptr vbase MAYBE_UNUSED, mp_limb_t * r MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    mpfq_p_5_get_mpn(vbase->obj, r, x);
}

static void mpfq_p_5_wrapper_set_mpz(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpz_srcptr);
static void mpfq_p_5_wrapper_set_mpz(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt r MAYBE_UNUSED, mpz_srcptr z MAYBE_UNUSED)
{
    mpfq_p_5_set_mpz(vbase->obj, r, z);
}

static void mpfq_p_5_wrapper_set_mpn(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mp_limb_t *, size_t);
static void mpfq_p_5_wrapper_set_mpn(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt r MAYBE_UNUSED, mp_limb_t * x MAYBE_UNUSED, size_t n MAYBE_UNUSED)
{
    mpfq_p_5_set_mpn(vbase->obj, r, x, n);
}

static unsigned long mpfq_p_5_wrapper_get_ui(mpfq_vbase_ptr, mpfq_p_5_src_elt);
static unsigned long mpfq_p_5_wrapper_get_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_src_elt x MAYBE_UNUSED)
{
    return mpfq_p_5_get_ui(vbase->obj, x);
}

static void mpfq_p_5_wrapper_set_zero(mpfq_vbase_ptr, mpfq_p_5_dst_elt);
static void mpfq_p_5_wrapper_set_zero(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt r MAYBE_UNUSED)
{
    mpfq_p_5_set_zero(vbase->obj, r);
}

static void mpfq_p_5_wrapper_set_ui(mpfq_vbase_ptr, mpfq_p_5_dst_elt, unsigned long);
static void mpfq_p_5_wrapper_set_ui(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt r MAYBE_UNUSED, unsigned long x MAYBE_UNUSED)
{
    mpfq_p_5_set_ui(vbase->obj, r, x);
}

static void mpfq_p_5_wrapper_set(mpfq_vbase_ptr, mpfq_p_5_dst_elt, mpfq_p_5_src_elt);
static void mpfq_p_5_wrapper_set(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_dst_elt r MAYBE_UNUSED, mpfq_p_5_src_elt s MAYBE_UNUSED)
{
    mpfq_p_5_set(vbase->obj, r, s);
}

static ptrdiff_t mpfq_p_5_wrapper_elt_stride(mpfq_vbase_ptr);
static ptrdiff_t mpfq_p_5_wrapper_elt_stride(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_elt_stride(vbase->obj);
}

static void mpfq_p_5_wrapper_clear(mpfq_vbase_ptr, mpfq_p_5_elt *);
static void mpfq_p_5_wrapper_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_elt * x MAYBE_UNUSED)
{
    mpfq_p_5_clear(vbase->obj, x);
}

static void mpfq_p_5_wrapper_init(mpfq_vbase_ptr, mpfq_p_5_elt *);
static void mpfq_p_5_wrapper_init(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpfq_p_5_elt * x MAYBE_UNUSED)
{
    mpfq_p_5_init(vbase->obj, x);
}

static void mpfq_p_5_wrapper_field_setopt(mpfq_vbase_ptr, unsigned long, void *);
static void mpfq_p_5_wrapper_field_setopt(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, void * y MAYBE_UNUSED)
{
    mpfq_p_5_field_setopt(vbase->obj, x, y);
}

static void mpfq_p_5_wrapper_field_specify(mpfq_vbase_ptr, unsigned long, const void *);
static void mpfq_p_5_wrapper_field_specify(mpfq_vbase_ptr vbase MAYBE_UNUSED, unsigned long dummy MAYBE_UNUSED, const void * vp MAYBE_UNUSED)
{
    mpfq_p_5_field_specify(vbase->obj, dummy, vp);
}

static void mpfq_p_5_wrapper_field_clear(mpfq_vbase_ptr);
static void mpfq_p_5_wrapper_field_clear(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_p_5_field_clear(vbase->obj);
}

static void mpfq_p_5_wrapper_field_init(mpfq_vbase_ptr);
static void mpfq_p_5_wrapper_field_init(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    mpfq_p_5_field_init(vbase->obj);
}

static int mpfq_p_5_wrapper_field_degree(mpfq_vbase_ptr);
static int mpfq_p_5_wrapper_field_degree(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_field_degree(vbase->obj);
}

static unsigned long mpfq_p_5_wrapper_field_characteristic_bits(mpfq_vbase_ptr);
static unsigned long mpfq_p_5_wrapper_field_characteristic_bits(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_field_characteristic_bits(vbase->obj);
}

static mpz_srcptr mpfq_p_5_wrapper_field_characteristic_srcptr(mpfq_vbase_ptr);
static mpz_srcptr mpfq_p_5_wrapper_field_characteristic_srcptr(mpfq_vbase_ptr vbase MAYBE_UNUSED)
{
    return mpfq_p_5_field_characteristic_srcptr(vbase->obj);
}

static void mpfq_p_5_wrapper_field_characteristic(mpfq_vbase_ptr, mpz_ptr);
static void mpfq_p_5_wrapper_field_characteristic(mpfq_vbase_ptr vbase MAYBE_UNUSED, mpz_ptr z MAYBE_UNUSED)
{
    mpfq_p_5_field_characteristic(vbase->obj, z);
}

static unsigned long mpfq_p_5_wrapper_impl_max_degree();
static unsigned long mpfq_p_5_wrapper_impl_max_degree()
{
    return mpfq_p_5_impl_max_degree();
}

static unsigned long mpfq_p_5_wrapper_impl_max_characteristic_bits();
static unsigned long mpfq_p_5_wrapper_impl_max_characteristic_bits()
{
    return mpfq_p_5_impl_max_characteristic_bits();
}

static const char * mpfq_p_5_wrapper_impl_name();
static const char * mpfq_p_5_wrapper_impl_name()
{
    return mpfq_p_5_impl_name();
}

/* Mpfq::engine::oo::oo_field_init */
/* Triggered by: oo */
void mpfq_p_5_oo_field_init(mpfq_vbase_ptr vbase)
{
    memset(vbase, 0, sizeof(struct mpfq_vbase_s));
    vbase->obj = malloc(sizeof(mpfq_p_5_field));
    mpfq_p_5_field_init((mpfq_p_5_dst_field) vbase->obj);
    vbase->impl_name = (const char * (*) ()) mpfq_p_5_wrapper_impl_name;
    vbase->impl_max_characteristic_bits = (unsigned long (*) ()) mpfq_p_5_wrapper_impl_max_characteristic_bits;
    vbase->impl_max_degree = (unsigned long (*) ()) mpfq_p_5_wrapper_impl_max_degree;
    vbase->field_characteristic = (void (*) (mpfq_vbase_ptr, mpz_ptr)) mpfq_p_5_wrapper_field_characteristic;
    vbase->field_characteristic_srcptr = (mpz_srcptr (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_field_characteristic_srcptr;
    vbase->field_characteristic_bits = (unsigned long (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_field_characteristic_bits;
    vbase->field_degree = (int (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_field_degree;
    vbase->field_init = (void (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_field_init;
    vbase->field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_field_clear;
    vbase->field_specify = (void (*) (mpfq_vbase_ptr, unsigned long, const void *)) mpfq_p_5_wrapper_field_specify;
    vbase->field_setopt = (void (*) (mpfq_vbase_ptr, unsigned long, void *)) mpfq_p_5_wrapper_field_setopt;
    vbase->init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_init;
    vbase->clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_clear;
    vbase->elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_elt_stride;
    vbase->set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_set;
    vbase->set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_set_ui;
    vbase->set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_set_zero;
    vbase->get_ui = (unsigned long (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_get_ui;
    vbase->set_mpn = (void (*) (mpfq_vbase_ptr, void *, mp_limb_t *, size_t)) mpfq_p_5_wrapper_set_mpn;
    vbase->set_mpz = (void (*) (mpfq_vbase_ptr, void *, mpz_srcptr)) mpfq_p_5_wrapper_set_mpz;
    vbase->get_mpn = (void (*) (mpfq_vbase_ptr, mp_limb_t *, const void *)) mpfq_p_5_wrapper_get_mpn;
    vbase->get_mpz = (void (*) (mpfq_vbase_ptr, mpz_ptr, const void *)) mpfq_p_5_wrapper_get_mpz;
    vbase->random = (void (*) (mpfq_vbase_ptr, void *, gmp_randstate_t)) mpfq_p_5_wrapper_random;
    vbase->random2 = (void (*) (mpfq_vbase_ptr, void *, gmp_randstate_t)) mpfq_p_5_wrapper_random2;
    vbase->add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_add;
    vbase->sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_sub;
    vbase->neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_neg;
    vbase->mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_mul;
    vbase->sqr = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_sqr;
    vbase->is_sqr = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_is_sqr;
    vbase->sqrt = (int (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_sqrt;
    vbase->pow = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long *, size_t)) mpfq_p_5_wrapper_pow;
    vbase->powz = (void (*) (mpfq_vbase_ptr, void *, const void *, mpz_srcptr)) mpfq_p_5_wrapper_powz;
    vbase->frobenius = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_frobenius;
    vbase->add_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_add_ui;
    vbase->sub_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_sub_ui;
    vbase->mul_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_mul_ui;
    vbase->normalize = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_normalize;
    vbase->inv = (int (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_inv;
    vbase->hadamard = (void (*) (mpfq_vbase_ptr, void *, void *, void *, void *)) mpfq_p_5_wrapper_hadamard;
    vbase->elt_ur_init = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_elt_ur_init;
    vbase->elt_ur_clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_elt_ur_clear;
    vbase->elt_ur_stride = (ptrdiff_t (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_elt_ur_stride;
    vbase->elt_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_elt_ur_set;
    vbase->elt_ur_set_elt = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_elt_ur_set_elt;
    vbase->elt_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_elt_ur_set_zero;
    vbase->elt_ur_set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_elt_ur_set_ui;
    vbase->elt_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_elt_ur_add;
    vbase->elt_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_elt_ur_neg;
    vbase->elt_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_elt_ur_sub;
    vbase->mul_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_mul_ur;
    vbase->sqr_ur = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_sqr_ur;
    vbase->reduce = (void (*) (mpfq_vbase_ptr, void *, void *)) mpfq_p_5_wrapper_reduce;
    vbase->addmul_si_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, long)) mpfq_p_5_wrapper_addmul_si_ur;
    vbase->cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *)) mpfq_p_5_wrapper_cmp;
    vbase->cmp_ui = (int (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_cmp_ui;
    vbase->is_zero = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_is_zero;
    vbase->asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *)) mpfq_p_5_wrapper_asprint;
    vbase->fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *)) mpfq_p_5_wrapper_fprint;
    vbase->print = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_print;
    vbase->sscan = (int (*) (mpfq_vbase_ptr, void *, const char *)) mpfq_p_5_wrapper_sscan;
    vbase->fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *)) mpfq_p_5_wrapper_fscan;
    vbase->scan = (int (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_scan;
    /* missing read */
    /* missing importdata */
    /* missing write */
    /* missing exportdata */
    vbase->vec_init = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_init;
    vbase->vec_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_p_5_wrapper_vec_reinit;
    vbase->vec_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_clear;
    vbase->vec_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_set;
    vbase->vec_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_set_zero;
    vbase->vec_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_setcoeff;
    vbase->vec_setcoeff_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_p_5_wrapper_vec_setcoeff_ui;
    vbase->vec_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_getcoeff;
    vbase->vec_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_add;
    vbase->vec_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_neg;
    vbase->vec_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_rev;
    vbase->vec_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_sub;
    vbase->vec_scal_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_scal_mul;
    vbase->vec_conv = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long)) mpfq_p_5_wrapper_vec_conv;
    vbase->vec_random = (void (*) (mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t)) mpfq_p_5_wrapper_vec_random;
    vbase->vec_random2 = (void (*) (mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t)) mpfq_p_5_wrapper_vec_random2;
    vbase->vec_cmp = (long (*) (mpfq_vbase_ptr, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_cmp;
    vbase->vec_is_zero = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_is_zero;
    vbase->vec_subvec = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_p_5_wrapper_vec_subvec;
    vbase->vec_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_p_5_wrapper_vec_subvec_const;
    vbase->vec_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_p_5_wrapper_vec_coeff_ptr;
    vbase->vec_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_p_5_wrapper_vec_coeff_ptr_const;
    vbase->vec_asprint = (long (*) (mpfq_vbase_ptr, char * *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_asprint;
    vbase->vec_fprint = (long (*) (mpfq_vbase_ptr, FILE *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_fprint;
    vbase->vec_print = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_print;
    vbase->vec_sscan = (long (*) (mpfq_vbase_ptr, void *, unsigned long *, const char *)) mpfq_p_5_wrapper_vec_sscan;
    vbase->vec_fscan = (long (*) (mpfq_vbase_ptr, FILE *, void *, unsigned long *)) mpfq_p_5_wrapper_vec_fscan;
    vbase->vec_scan = (long (*) (mpfq_vbase_ptr, void *, unsigned long *)) mpfq_p_5_wrapper_vec_scan;
    /* missing vec_read */
    /* missing vec_write */
    /* missing vec_import */
    /* missing vec_export */
    vbase->vec_hamming_weight = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_hamming_weight;
    vbase->vec_find_first_set = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_find_first_set;
    vbase->vec_simd_hamming_weight = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_simd_hamming_weight;
    vbase->vec_simd_find_first_set = (long (*) (mpfq_vbase_ptr, const void *, unsigned long)) mpfq_p_5_wrapper_vec_simd_find_first_set;
    vbase->vec_ur_init = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_init;
    vbase->vec_ur_set_zero = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_set_zero;
    vbase->vec_ur_set_vec = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_set_vec;
    vbase->vec_ur_reinit = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned long)) mpfq_p_5_wrapper_vec_ur_reinit;
    vbase->vec_ur_clear = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_clear;
    vbase->vec_ur_set = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_set;
    vbase->vec_ur_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_setcoeff;
    vbase->vec_ur_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_getcoeff;
    vbase->vec_ur_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_add;
    vbase->vec_ur_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_sub;
    vbase->vec_ur_neg = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_neg;
    vbase->vec_ur_rev = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_ur_rev;
    vbase->vec_scal_mul_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned long)) mpfq_p_5_wrapper_vec_scal_mul_ur;
    vbase->vec_conv_ur = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long)) mpfq_p_5_wrapper_vec_conv_ur;
    vbase->vec_reduce = (void (*) (mpfq_vbase_ptr, void *, void *, unsigned long)) mpfq_p_5_wrapper_vec_reduce;
    vbase->vec_ur_subvec = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_p_5_wrapper_vec_ur_subvec;
    vbase->vec_ur_subvec_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_p_5_wrapper_vec_ur_subvec_const;
    vbase->vec_ur_coeff_ptr = (void * (*) (mpfq_vbase_ptr, void *, long)) mpfq_p_5_wrapper_vec_ur_coeff_ptr;
    vbase->vec_ur_coeff_ptr_const = (const void * (*) (mpfq_vbase_ptr, const void *, long)) mpfq_p_5_wrapper_vec_ur_coeff_ptr_const;
    vbase->vec_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, long)) mpfq_p_5_wrapper_vec_elt_stride;
    vbase->vec_ur_elt_stride = (ptrdiff_t (*) (mpfq_vbase_ptr, long)) mpfq_p_5_wrapper_vec_ur_elt_stride;
    vbase->poly_init = (void (*) (mpfq_vbase_ptr, void *, unsigned int)) mpfq_p_5_wrapper_poly_init;
    vbase->poly_clear = (void (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_poly_clear;
    vbase->poly_set = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_poly_set;
    vbase->poly_setmonic = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_poly_setmonic;
    vbase->poly_setcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_p_5_wrapper_poly_setcoeff;
    vbase->poly_setcoeff_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long, unsigned int)) mpfq_p_5_wrapper_poly_setcoeff_ui;
    vbase->poly_getcoeff = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned int)) mpfq_p_5_wrapper_poly_getcoeff;
    vbase->poly_deg = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_poly_deg;
    vbase->poly_add = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_add;
    vbase->poly_sub = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_sub;
    vbase->poly_set_ui = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_poly_set_ui;
    vbase->poly_add_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_poly_add_ui;
    vbase->poly_sub_ui = (void (*) (mpfq_vbase_ptr, void *, const void *, unsigned long)) mpfq_p_5_wrapper_poly_sub_ui;
    vbase->poly_neg = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_poly_neg;
    vbase->poly_scal_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_scal_mul;
    vbase->poly_mul = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_mul;
    vbase->poly_divmod = (int (*) (mpfq_vbase_ptr, void *, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_divmod;
    vbase->poly_precomp_mod = (void (*) (mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_wrapper_poly_precomp_mod;
    vbase->poly_mod_pre = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, const void *)) mpfq_p_5_wrapper_poly_mod_pre;
    vbase->poly_gcd = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_gcd;
    vbase->poly_xgcd = (void (*) (mpfq_vbase_ptr, void *, void *, void *, const void *, const void *)) mpfq_p_5_wrapper_poly_xgcd;
    vbase->poly_random = (void (*) (mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t)) mpfq_p_5_wrapper_poly_random;
    vbase->poly_random2 = (void (*) (mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t)) mpfq_p_5_wrapper_poly_random2;
    vbase->poly_cmp = (int (*) (mpfq_vbase_ptr, const void *, const void *)) mpfq_p_5_wrapper_poly_cmp;
    vbase->poly_asprint = (int (*) (mpfq_vbase_ptr, char * *, const void *)) mpfq_p_5_wrapper_poly_asprint;
    vbase->poly_fprint = (int (*) (mpfq_vbase_ptr, FILE *, const void *)) mpfq_p_5_wrapper_poly_fprint;
    vbase->poly_print = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_poly_print;
    vbase->poly_sscan = (int (*) (mpfq_vbase_ptr, void *, const char *)) mpfq_p_5_wrapper_poly_sscan;
    vbase->poly_fscan = (int (*) (mpfq_vbase_ptr, FILE *, void *)) mpfq_p_5_wrapper_poly_fscan;
    vbase->poly_scan = (int (*) (mpfq_vbase_ptr, void *)) mpfq_p_5_wrapper_poly_scan;
    vbase->simd_groupsize = (int (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_simd_groupsize;
    vbase->simd_hamming_weight = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_simd_hamming_weight;
    vbase->simd_find_first_set = (int (*) (mpfq_vbase_ptr, const void *)) mpfq_p_5_wrapper_simd_find_first_set;
    vbase->simd_get_ui_at = (unsigned long (*) (mpfq_vbase_ptr, const void *, int)) mpfq_p_5_wrapper_simd_get_ui_at;
    vbase->simd_set_ui_at = (void (*) (mpfq_vbase_ptr, void *, int, unsigned long)) mpfq_p_5_wrapper_simd_set_ui_at;
    vbase->simd_add_ui_at = (void (*) (mpfq_vbase_ptr, void *, const void *, int, unsigned long)) mpfq_p_5_wrapper_simd_add_ui_at;
    vbase->simd_set_ui_all = (void (*) (mpfq_vbase_ptr, void *, unsigned long)) mpfq_p_5_wrapper_simd_set_ui_all;
    vbase->add_dotprod = (void (*) (mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_5_wrapper_add_dotprod;
    vbase->oo_field_init = (void (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_oo_field_init;
    vbase->oo_field_clear = (void (*) (mpfq_vbase_ptr)) mpfq_p_5_wrapper_oo_field_clear;
}


/* vim:set ft=cpp: */
