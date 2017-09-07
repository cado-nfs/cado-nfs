/* MPFQ generated file -- do not edit */

#include "mpfq_2_64.h"

/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::gf2n::field */
/* Automatically generated code for GF(2^64) */
/* Definition polynomial P = X^64 + X^4 + X^3 + X + 1 */
/* Active handler: Mpfq::gf2n::trivialities */
/* Active handler: Mpfq::gf2n::io */
/* Active handler: Mpfq::gf2n::linearops */
/* Active handler: Mpfq::gf2n::inversion */
/* Active handler: Mpfq::gf2n::reduction */
/* Active handler: Mpfq::gf2n::mul */
/* Active handler: Mpfq::defaults::poly */
/* Options used:{
   coeffs=[ 64, 4, 3, 1, 0, ],
   helper=/tmp/mpfq-cado/gf2n/helper/helper,
   n=64,
   no_gmp=1,
   output_path=i386,
   slice=4,
   table=/tmp/mpfq-cado/gf2x/wizard.table,
   tag=2_64,
   w=32,
   } */


/* Functions operating on the field structure */

/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */
/* missing powz */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */
/* missing fprint */
/* missing print */
/* missing scan */

/* Vector functions */
/* *Mpfq::defaults::vec::alloc::code_for_vec_init, Mpfq::defaults::vec */
void mpfq_2_64_vec_init(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec * v, unsigned int n)
{
    unsigned int i;
    *v = (mpfq_2_64_vec) malloc (n*sizeof(mpfq_2_64_elt));
    for(i = 0; i < n; i++)
        mpfq_2_64_init(K, (*v) + i);
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_reinit, Mpfq::defaults::vec */
void mpfq_2_64_vec_reinit(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        unsigned int i;
        *v = (mpfq_2_64_vec) realloc (*v, m * sizeof(mpfq_2_64_elt));
        for(i = n; i < m; i+=1)
            mpfq_2_64_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_2_64_clear(K, (*v) + i);
        *v = (mpfq_2_64_vec) realloc (*v, m * sizeof(mpfq_2_64_elt));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_clear, Mpfq::defaults::vec */
void mpfq_2_64_vec_clear(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec * v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_clear(K, (*v) + i);
    free(*v);
}

/* missing vec_random */
/* missing vec_random2 */
/* missing vec_asprint */
/* missing vec_fprint */
/* missing vec_print */
/* missing vec_sscan */
/* missing vec_fscan */
/* missing vec_scan */
/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_init, Mpfq::defaults::vec */
void mpfq_2_64_vec_ur_init(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec_ur * v, unsigned int n)
{
    unsigned int i;
    *v = (mpfq_2_64_vec_ur) malloc (n*sizeof(mpfq_2_64_elt_ur));
    for(i = 0; i < n; i+=1)
        mpfq_2_64_elt_ur_init(K, &( (*v)[i]));
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_reinit, Mpfq::defaults::vec */
void mpfq_2_64_vec_ur_reinit(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec_ur * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        *v = (mpfq_2_64_vec_ur) realloc (*v, m * sizeof(mpfq_2_64_elt_ur));
        unsigned int i;
        for(i = n; i < m; i+=1)
            mpfq_2_64_elt_ur_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            mpfq_2_64_elt_ur_clear(K, (*v) + i);
        *v = (mpfq_2_64_vec_ur) realloc (*v, m * sizeof(mpfq_2_64_elt_ur));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_clear, Mpfq::defaults::vec */
void mpfq_2_64_vec_ur_clear(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec_ur * v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        mpfq_2_64_elt_ur_clear(K, &( (*v)[i]));
    free(*v);
}


/* Polynomial functions */
/* *Mpfq::defaults::poly::code_for_poly_setmonic */
void mpfq_2_64_poly_setmonic(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_poly q, mpfq_2_64_src_poly p)
{
    long degp = mpfq_2_64_poly_deg(K, p);
    if (degp == -1) {
        q->size = 0;
        return;
    }
    if (degp == 0) {
        mpfq_2_64_elt aux;
        mpfq_2_64_init(K, &aux);
        mpfq_2_64_set_ui(K, aux, 1);
        mpfq_2_64_poly_setcoeff(K, q, aux, 0);
        mpfq_2_64_clear(K, &aux);
        q->size = 1;
        return;
    }
    mpfq_2_64_elt lc;	/* spurious uninit warning sometimes */
    mpfq_2_64_init(K, &lc);
    mpfq_2_64_poly_getcoeff(K, lc, p, degp);
    mpfq_2_64_inv(K, lc, lc);
    mpfq_2_64_poly_setcoeff_ui(K, q, 1, degp);
    mpfq_2_64_vec_scal_mul(K, q->c, p->c, lc, degp);
    q->size = degp+1;
    mpfq_2_64_clear(K, &lc);
}

/* *Mpfq::defaults::poly::code_for_poly_divmod */
int mpfq_2_64_poly_divmod(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_poly q, mpfq_2_64_dst_poly r, mpfq_2_64_src_poly a, mpfq_2_64_src_poly b)
{
    if (b->size == 0) {
        return 0;
    }
    if (a->size == 0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    int dega = mpfq_2_64_poly_deg(K, a);
    if (dega<0) {
        q->size = 0; r->size = 0;
        return 1;
    }
    // Compute deg b and inverse of leading coef
    int degb = mpfq_2_64_poly_deg(K, b);
    if (degb<0) {
        return 0;
    }
    if (degb > dega) {
        q->size=0;
        mpfq_2_64_poly_set(K, r, a);
        return 1;
    }
    int bmonic;
    mpfq_2_64_elt ilb;
    mpfq_2_64_init(K, &ilb);
    mpfq_2_64_elt temp;
    mpfq_2_64_init(K, &temp);
    mpfq_2_64_poly_getcoeff(K, temp, b, degb);
    if (mpfq_2_64_cmp_ui(K, temp, 1) == 0) {
        mpfq_2_64_set_ui(K, ilb, 1);
        bmonic = 1;
    } else {
        mpfq_2_64_inv(K, ilb, temp);
        bmonic = 0;
    }
    
    mpfq_2_64_poly qq, rr;
    mpfq_2_64_poly_init(K, qq, dega-degb+1);
    mpfq_2_64_poly_init(K, rr, dega);
    
    mpfq_2_64_poly_set(K, rr, a);
    mpfq_2_64_elt aux, aux2;
    
    mpfq_2_64_init(K, &aux);
    mpfq_2_64_init(K, &aux2);
    
    int i;
    int j;
    for (i = dega; i >= (int)degb; --i) {
        mpfq_2_64_poly_getcoeff(K, aux, rr, i);
        if (!bmonic) 
            mpfq_2_64_mul(K, aux, aux, ilb);
        mpfq_2_64_poly_setcoeff(K, qq, aux, i-degb);
        for (j = i-1; j >= (int)(i - degb); --j) {
            mpfq_2_64_poly_getcoeff(K, temp, b, j-i+degb);
            mpfq_2_64_mul(K, aux2, aux, temp);
            mpfq_2_64_poly_getcoeff(K, temp, rr, j);
    
            mpfq_2_64_sub(K, temp, temp, aux2);
            mpfq_2_64_poly_setcoeff(K, rr, temp, j);
        }
    }    
    
    rr->size = degb;
    int degr = mpfq_2_64_poly_deg(K, rr);
    rr->size = degr+1;
    
    if (q != NULL) 
        mpfq_2_64_poly_set(K, q, qq);
    if (r != NULL)
        mpfq_2_64_poly_set(K, r, rr);
    mpfq_2_64_clear(K, &temp);
    mpfq_2_64_clear(K, &ilb);
    mpfq_2_64_clear(K, &aux);
    mpfq_2_64_clear(K, &aux2);
    mpfq_2_64_poly_clear(K, rr);
    mpfq_2_64_poly_clear(K, qq);
    return 1;
}

static void mpfq_2_64_poly_preinv(mpfq_2_64_dst_field, mpfq_2_64_dst_poly, mpfq_2_64_src_poly, unsigned int);
/* *Mpfq::defaults::poly::code_for_poly_precomp_mod */
/* Triggered by: poly_precomp_mod */
static void mpfq_2_64_poly_preinv(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_poly q, mpfq_2_64_src_poly p, unsigned int n)
{
    // Compute the inverse of p(x) modulo x^n
    // Newton iteration: x_{n+1} = x_n + x_n(1 - a*x_n)
    // Requires p(0) = 1
    // Assume p != q (no alias)
    mpfq_2_64_elt temp;	/* spurious uninit warning sometimes */
    mpfq_2_64_init(K, &temp);
    mpfq_2_64_poly_getcoeff(K, temp, p, 0);//Should be in the assert
    assert( mpfq_2_64_cmp_ui(K, temp, 1) == 0);
    assert (p != q);
    int m;
    if (n <= 2) {
        mpfq_2_64_poly_setcoeff_ui(K, q, 1, 0);
        q->size = 1;
        m = 1;
        if (n == 1)
            return;
    } else {
        // n >= 3: recursive call at prec m = ceil(n/2)
        m = 1 + ((n-1)/2);
        mpfq_2_64_poly_preinv(K, q, p, m);
    }
    // enlarge q if necessary
    if (q->alloc < n) {
        mpfq_2_64_vec_reinit(K, &(q->c), q->alloc, n);
        q->alloc = n;
    }
    // refine value
    mpfq_2_64_vec tmp;
    mpfq_2_64_vec_init(K, &tmp, m+n-1);
    
    mpfq_2_64_vec_conv(K, tmp, p->c, MIN(n, p->size), q->c, m);
    int nn = MIN(n, MIN(n, p->size) + m -1);
    mpfq_2_64_vec_neg(K, tmp, tmp, nn);
    mpfq_2_64_vec_getcoeff(K, temp, tmp, 0);
    mpfq_2_64_add_ui(K, temp, temp, 1);
    mpfq_2_64_vec_setcoeff(K, tmp, temp, 0);
    mpfq_2_64_vec_conv(K, tmp, q->c, m, tmp, nn);
    mpfq_2_64_vec_set(K, mpfq_2_64_vec_subvec(K, q->c, m), mpfq_2_64_vec_subvec(K, tmp, m), n-m);
    q->size = n;
    
    mpfq_2_64_clear(K, &temp);
    mpfq_2_64_vec_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::poly::code_for_poly_precomp_mod */
void mpfq_2_64_poly_precomp_mod(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_poly q, mpfq_2_64_src_poly p)
{
    assert(p != q);
    int N = mpfq_2_64_poly_deg(K, p);
    mpfq_2_64_poly rp;
    mpfq_2_64_poly_init(K, rp, N+1);
    mpfq_2_64_vec_rev(K, rp->c, p->c, N+1);
    rp->size = N+1;
    mpfq_2_64_poly_preinv(K, q, rp, N);
    mpfq_2_64_poly_clear(K, rp);
}

/* *Mpfq::defaults::poly::code_for_poly_mod_pre */
void mpfq_2_64_poly_mod_pre(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_dst_poly r, mpfq_2_64_src_poly q, mpfq_2_64_src_poly p, mpfq_2_64_src_poly irp)
{
    int N = mpfq_2_64_poly_deg(K, p);
    int degq = mpfq_2_64_poly_deg(K, q);
    if (degq < N) {
        mpfq_2_64_poly_set(K, r, q);
        return;
    }
    int m = degq - N;
    assert (degq <= 2*N-2);
    mpfq_2_64_poly revq;
    mpfq_2_64_poly_init(K, revq, MAX(degq+1, m+1));
    mpfq_2_64_vec_rev(K, revq->c, q->c, degq+1);
    revq->size = q->size;
    mpfq_2_64_poly_mul(K, revq, revq, irp);
    mpfq_2_64_vec_rev(K, revq->c, revq->c, m+1);
    revq->size = m+1;
    
    mpfq_2_64_poly_mul(K, revq, revq, p);
    mpfq_2_64_poly_sub(K, r, q, revq);
    r->size = mpfq_2_64_poly_deg(K, r)+1;
    mpfq_2_64_poly_clear(K, revq);
}

/* missing poly_random */
/* missing poly_random2 */
/* missing poly_asprint */
/* missing poly_fprint */
/* missing poly_print */
/* missing poly_sscan */
/* missing poly_fscan */
/* missing poly_scan */

/* vim:set ft=cpp: */
