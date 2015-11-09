/* MPFQ generated file -- do not edit */

#include "mpfq_2_64.h"

#ifndef GMP_LIMB_BITS
#error "Please arrange so that GMP_LIMB_BITS is defined before including this file"
#endif

#if !(GMP_LIMB_BITS == 32)
#error "Constraints not met for this file: GMP_LIMB_BITS == 32"
#endif
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
/* *Mpfq::defaults::pow::code_for_powz, Mpfq::gf2n::trivialities */
void mpfq_2_64_powz(mpfq_2_64_dst_field k, mpfq_2_64_dst_elt y, mpfq_2_64_src_elt x, mpz_srcptr z)
{
        if (mpz_sgn(z) < 0) {
            mpz_t mz;
            mpz_init(mz);
            mpz_neg(mz, z);
            mpfq_2_64_powz(k, y, x, mz);
            mpfq_2_64_inv(k, y, y);
            mpz_clear(mz);
        } else if (mpz_sizeinbase(z, 2) > mpfq_2_64_field_degree(k) * mpfq_2_64_field_characteristic_bits(k)) {
            mpz_t zr;
            mpz_init(zr);
            mpz_t ppz;
            mpz_init(ppz);
            mpfq_2_64_field_characteristic(k, ppz);
            mpz_pow_ui(ppz,ppz,mpfq_2_64_field_degree(k));
            mpz_sub_ui(ppz,ppz,1);
            mpz_fdiv_r(zr, z, ppz);
            mpfq_2_64_powz(k, y, x, zr);
            mpz_clear(ppz);
            mpz_clear(zr);
        } else {
            mpfq_2_64_pow(k, y, x, z->_mp_d, mpz_size(z));
        }
}


/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */
/* *Mpfq::gf2n::io::code_for_asprint */
int mpfq_2_64_asprint(mpfq_2_64_dst_field k, char * * pstr, mpfq_2_64_src_elt x)
{
    int type = k->io_type;
    int i, n; 
    
    // Numerical io.
    if (type <= 16) {
        // allocate enough room for base 2 conversion.
        *pstr = (char *)mpfq_malloc_check((64+1));
    
        mp_limb_t tmp[2 + 1];
        for (i = 0; i < 2; ++i)
            tmp[i] = x[i];
    
        // mpn_get_str() needs a non-zero most significant limb
        int msl = 2 - 1;
        while ((msl > 0) && (tmp[msl] == 0))
            msl--;
        msl++;
        if ((msl == 1) && (tmp[0] == 0)) {
            (*pstr)[0] = '0';
            (*pstr)[1] = '\0';
            return 1;
        }
        n = mpn_get_str((unsigned char*)(*pstr), type, tmp, msl);
        for (i = 0; i < n; ++i) {
            char c = (*pstr)[i] + '0';
            if (c > '9')
                c = c-'0'+'a'-10;
            (*pstr)[i] = c;
        }
        (*pstr)[n] = '\0';
    
        // Remove leading 0s
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
    // Polynomial io.
    else {
        char c = (char)type;
        // allocate (more than) enough room for polynomial conversion.
        // Warning: this is for exponent that fit in 3 digits
        *pstr = (char *)mpfq_malloc_check((8*64+1));
        {
            unsigned int j;
            int sth = 0;
            char *ptr = *pstr;
            for(j = 0 ; j < 64 ; j++) {
                if (x[j/32] >> (j % 32) & 1UL) {
                	if (sth) {
                        *ptr++ = ' ';
                        *ptr++ = '+';
                        *ptr++ = ' ';
                    }
                	sth = 1;
                	if (j == 0) {
                        *ptr++ = '1';      
                	} else if (j == 1) {
                        *ptr++ = c;      
                	} else {
                        int ret = sprintf(ptr,"%c^%d",c,j);
                        ptr += ret;
                	}
                }
            }
            if (!sth) {
                *ptr++ = '0';
            }
            *ptr = '\0';
            return ptr - *pstr;
        }
    }
}

/* *Mpfq::defaults::code_for_fprint */
int mpfq_2_64_fprint(mpfq_2_64_dst_field k, FILE * file, mpfq_2_64_src_elt x)
{
    char *str;
    int rc;
    mpfq_2_64_asprint(k,&str,x);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *Mpfq::gf2n::io::code_for_sscan */
int mpfq_2_64_sscan(mpfq_2_64_dst_field k, mpfq_2_64_dst_elt z, const char * str)
{
    if (k->io_type <= 16) {
        unsigned char *tmp;
        int len = strlen(str);
        tmp = (unsigned char *)mpfq_malloc_check(len+1);
        int i;
        for (i = 0; i < len; ++i) {
            if (str[i] >= '0' && str[i] <= '9' && str[i] < '0' + k->io_type) {
                tmp[i] = str[i] - '0';
            } else if (str[i] >= 'a' && str[i] < 'a' + k->io_type - 10) {
                tmp[i] = str[i] + 10 - 'a';
            } else if (str[i] >= 'A' && str[i] < 'A' + k->io_type - 10) {
                tmp[i] = str[i] + 10 - 'A';
            } else {
                break;
            }
        }
        if (i == 0)
            return 0;
        len = i;
        tmp[len]='\0';
        mp_limb_t *zz;
        // Allocate one limb per byte... very conservative.
        zz = (mp_limb_t *)mpfq_malloc_check(len*sizeof(mp_limb_t));
        int ret = mpn_set_str(zz, tmp, len, k->io_type);
        free(tmp);
        if (ret > 2) {
            free(zz);
            return 0;
        }
        mpfq_copy(z, zz, ret);
        mpfq_zero(z + ret, 2 - ret);
        free(zz);
        return len;
    } else {
        fprintf(stderr, "Polynomial io not implemented for reading\n");
        return 0;
    }
}

/* *Mpfq::gf2n::io::code_for_fscan */
int mpfq_2_64_fscan(mpfq_2_64_dst_field k, FILE * file, mpfq_2_64_dst_elt z)
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
                allocated = len + allocated / 4;
                tmp = (char*)realloc(tmp, allocated);
            }
            tmp[len]=c;
            len++;
        }
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated);
    }
    tmp[len]='\0';
    int ret=mpfq_2_64_sscan(k,z,tmp);
    free(tmp);
    return ret ? len : 0;
}


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

/* *Mpfq::defaults::vec::io::code_for_vec_asprint, Mpfq::defaults::vec */
int mpfq_2_64_vec_asprint(mpfq_2_64_dst_field K MAYBE_UNUSED, char * * pstr, mpfq_2_64_src_vec w, unsigned int n)
{
    if (n == 0) {
        *pstr = (char *)mpfq_malloc_check(4);
        sprintf(*pstr, "[ ]");
        return strlen(*pstr);
    }
    int alloc = 100;
    int len = 0;
    *pstr = (char *)mpfq_malloc_check(alloc);
    char *str = *pstr;
    *str++ = '[';
    *str++ = ' ';
    len = 2;
    unsigned int i;
    for(i = 0; i < n; i+=1) {
        if (i) {
            (*pstr)[len++] = ',';
            (*pstr)[len++] = ' ';
        }
        char *tmp;
        mpfq_2_64_asprint(K, &tmp, w[i]);
        int ltmp = strlen(tmp);
        if (len+ltmp+4 > alloc) {
            alloc = len+ltmp+100 + alloc / 4;
            *pstr = (char *)realloc(*pstr, alloc);
        }
        strncpy(*pstr+len, tmp, ltmp+4);
        len += ltmp;
        free(tmp);
    }
    (*pstr)[len++] = ' ';
    (*pstr)[len++] = ']';
    (*pstr)[len] = '\0';
    return len;
}

/* *Mpfq::defaults::vec::io::code_for_vec_fprint, Mpfq::defaults::vec */
int mpfq_2_64_vec_fprint(mpfq_2_64_dst_field K MAYBE_UNUSED, FILE * file, mpfq_2_64_src_vec w, unsigned int n)
{
    char *str;
    int rc;
    mpfq_2_64_vec_asprint(K,&str,w,n);
    rc = fprintf(file,"%s",str);
    free(str);
    return rc;
}

/* *Mpfq::defaults::vec::io::code_for_vec_print, Mpfq::defaults::vec */
int mpfq_2_64_vec_print(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_src_vec w, unsigned int n)
{
    return mpfq_2_64_vec_fprint(K,stdout,w,n);
}

/* *Mpfq::defaults::vec::io::code_for_vec_sscan, Mpfq::defaults::vec */
int mpfq_2_64_vec_sscan(mpfq_2_64_dst_field K MAYBE_UNUSED, mpfq_2_64_vec * w, unsigned int * n, const char * str)
{
    // start with a clean vector
    unsigned int nn;
    int len = 0;
    mpfq_2_64_vec_reinit(K, w, *n, 0);
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
            mpfq_2_64_vec_reinit(K, w, nn, i+1);
            *n = nn = i+1;
        }
        int ret = mpfq_2_64_sscan(K, mpfq_2_64_vec_coeff_ptr(K, *w, i), str + len);
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

/* *Mpfq::defaults::vec::io::code_for_vec_fscan, Mpfq::defaults::vec */
int mpfq_2_64_vec_fscan(mpfq_2_64_dst_field K MAYBE_UNUSED, FILE * file, mpfq_2_64_vec * w, unsigned int * n)
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
    int ret=mpfq_2_64_vec_sscan(K,w,n,tmp);
    free(tmp);
    return ret;
}

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


/* vim:set ft=cpp: */
