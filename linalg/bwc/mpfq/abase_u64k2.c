/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_u64k2.h"

#include "binary-dotprods-backends.h"
#include <inttypes.h>
static int abase_u64k2_impl_mpi_attr;     /* for MPI functions */
static MPI_Datatype abase_u64k2_impl_mpi_datatype;
static MPI_Datatype abase_u64k2_impl_mpi_datatype_ur;
static MPI_Op abase_u64k2_impl_mpi_addition_op;
static MPI_Op abase_u64k2_impl_mpi_addition_op_ur;
static int abase_u64k2_impl_mpi_use_count;   /* several stacked init()/clear() pairs are supported */
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
/* *simd_u64k::code_for_field_specify */
void abase_u64k2_field_specify(abase_u64k2_dst_field K MAYBE_UNUSED, unsigned long tag, void * x MAYBE_UNUSED)
{
        if (tag == MPFQ_GROUPSIZE) {
            assert(*(int*)x == 128);
        } else if (tag == MPFQ_PRIME_MPZ) {
            assert(mpz_cmp_ui((mpz_srcptr)x, 2) == 0);
        } else {
            fprintf(stderr, "Unsupported field_specify tag %ld\n", tag);
        }
}


/* Element allocation functions */

/* Elementary assignment functions */
/* missing set_ui */
/* missing get_ui */
/* missing set_mpn */
/* missing set_mpz */
/* missing get_mpn */
/* missing get_mpz */

/* Assignment of random values */

/* Arithmetic operations on elements */
/* missing mul */
/* missing sqr */
/* missing is_sqr */
/* missing sqrt */
/* missing pow */
/* missing frobenius */
/* missing add_ui */
/* missing sub_ui */
/* missing mul_ui */
/* missing inv */

/* Operations involving unreduced elements */
/* missing elt_ur_set_ui */
/* missing mul_ur */
/* missing sqr_ur */
/* missing reduce */

/* Comparison functions */
/* missing cmp_ui */

/* Input/output functions */
/* *io::code_for_asprint */
void abase_u64k2_asprint(abase_u64k2_dst_field K MAYBE_UNUSED, char * * ps, abase_u64k2_src_elt x)
{
        const uint64_t * y = x;
        const unsigned int stride = abase_u64k2_vec_elt_stride(K,1)/sizeof(uint64_t);
        *ps = malloc(stride * 16 + 1);
        if (!*ps) MALLOC_FAILED();
        for(unsigned int i = 0 ; i < stride ; i++) {
            snprintf((*ps) + i * 16, 17, "%" PRIx64, y[i]);
        }
}

/* *io::code_for_fprint */
void abase_u64k2_fprint(abase_u64k2_dst_field k, FILE * file, abase_u64k2_src_elt x)
{
    char *str;
    abase_u64k2_asprint(k,&str,x);
    fprintf(file,"%s",str);
    free(str);
}

/* *io::code_for_sscan */
int abase_u64k2_sscan(abase_u64k2_dst_field k MAYBE_UNUSED, abase_u64k2_dst_elt z, const char * str)
{
        char tmp[17];
        uint64_t * y = z;
        const unsigned int stride = abase_u64k2_vec_elt_stride(K,1)/sizeof(uint64_t);
        assert(strlen(str) >= 1 * 16);
        int r = 0;
        for(unsigned int i = 0 ; i < stride ; i++) {
            memcpy(tmp, str + i * 16, 16);
            tmp[16]=0;
            if (sscanf(tmp, "%" SCNx64, &(y[i])) == 1) {
                r+=16;
            } else {
                return r;
            }
        }
        return r;
}

/* *io::code_for_fscan */
int abase_u64k2_fscan(abase_u64k2_dst_field k, FILE * file, abase_u64k2_dst_elt z)
{
    char *tmp;
    int allocated, len=0;
    int c, start=0;
    allocated=100;
    tmp = (char *)malloc(allocated*sizeof(char));
    if (!tmp)
        MALLOC_FAILED();
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
                allocated+=100;
                tmp = (char*)realloc(tmp, allocated*sizeof(char));
            }
            tmp[len]=c;
            len++;
            start=1;
        }
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated*sizeof(char));
    }
    tmp[len]='\0';
    int ret=abase_u64k2_sscan(k,z,tmp);
    free(tmp);
    return ret;
}


/* Vector functions */
/* *Mpfq::defaults::vec::alloc::code_for_vec_init, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_init(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec * v, unsigned int n)
{
    unsigned int i;
    *v = (abase_u64k2_vec) malloc (n*sizeof(abase_u64k2_elt));
    for(i = 0; i < n; i+=1)
        abase_u64k2_init(K, (*v) + i);
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_reinit, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_reinit(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        unsigned int i;
        *v = (abase_u64k2_vec) realloc (*v, m * sizeof(abase_u64k2_elt));
        for(i = n; i < m; i+=1)
            abase_u64k2_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            abase_u64k2_clear(K, (*v) + i);
        *v = (abase_u64k2_vec) realloc (*v, m * sizeof(abase_u64k2_elt));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_clear, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_clear(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec * v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k2_clear(K, (*v) + i);
    free(*v);
}

/* missing vec_setcoef_ui */
/* missing vec_scal_mul */
/* missing vec_conv */
/* *Mpfq::defaults::vec::io::code_for_vec_asprint, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_asprint(abase_u64k2_dst_field K MAYBE_UNUSED, char * * pstr, abase_u64k2_src_vec w, unsigned int n)
{
    if (n == 0) {
        *pstr = (char *)malloc(4*sizeof(char));
        sprintf(*pstr, "[ ]");
        return;
    }
    int alloc = 100;
    int len = 0;
    *pstr = (char *)malloc(alloc*sizeof(char));
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
        abase_u64k2_asprint(K, &tmp, w[i]);
        int ltmp = strlen(tmp);
        if (len+ltmp+4 > alloc) {
            alloc = len+ltmp+100;
            *pstr = (char *)realloc(*pstr, alloc*sizeof(char));
        }
        strncpy(*pstr+len, tmp, ltmp+4);
        len += ltmp;
        free(tmp);
    }
    (*pstr)[len++] = ' ';
    (*pstr)[len++] = ']';
    (*pstr)[len] = '\0';
}

/* *Mpfq::defaults::vec::io::code_for_vec_fprint, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_fprint(abase_u64k2_dst_field K MAYBE_UNUSED, FILE * file, abase_u64k2_src_vec w, unsigned int n)
{
    char *str;
    abase_u64k2_vec_asprint(K,&str,w,n);
    fprintf(file,"%s",str);
    free(str);
}

/* *Mpfq::defaults::vec::io::code_for_vec_print, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_print(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_src_vec w, unsigned int n)
{
    abase_u64k2_vec_fprint(K,stdout,w,n);
}

/* *Mpfq::defaults::vec::io::code_for_vec_sscan, Mpfq::defaults::vec, Mpfq::defaults */
int abase_u64k2_vec_sscan(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec * w, unsigned int * n, const char * str)
{
    // start with a clean vector
    abase_u64k2_vec_reinit(K, w, *n, 0);
    *n = 0;
    while (isspace((int)(unsigned char)str[0]))
        str++;
    if (str[0] != '[')
        return 0;
    str++;
    if (str[0] != ' ')
        return 0;
    str++;
    if (str[0] == ']') {
        return 1;
    }
    unsigned int i = 0;
    for (;;) {
        if (*n < i+1) {
            abase_u64k2_vec_reinit(K, w, *n, i+1);
            *n = i+1;
        }
        int ret;
        ret = abase_u64k2_sscan(K, (*w)[i], str);
        if (!ret) {
            return 0;
        }
        i++;
        while (isdigit((int)(unsigned char)str[0]))
            str++;
        while (isspace((int)(unsigned char)str[0]))
            str++;
        if (str[0] == ']')
            break;
        if (str[0] != ',')
            return 0;
        str++;
        while (isspace((int)(unsigned char)str[0]))
            str++;
    }
    return 1;
}

/* *Mpfq::defaults::vec::io::code_for_vec_fscan, Mpfq::defaults::vec, Mpfq::defaults */
int abase_u64k2_vec_fscan(abase_u64k2_dst_field K MAYBE_UNUSED, FILE * file, abase_u64k2_vec * w, unsigned int * n)
{
    char *tmp;
    int c;
    int allocated, len=0;
    allocated=100;
    tmp = (char *)malloc(allocated*sizeof(char));
    if (!tmp)
        MALLOC_FAILED();
    for(;;) {
        c = fgetc(file);
        if (c==EOF)
            return 0;
        if (len==allocated) {
            allocated+=100;
            tmp = (char*)realloc(tmp, allocated*sizeof(char));
        }
        tmp[len]=c;
        len++;
        if (c==']')
            break;
    }
    if (len==allocated) {
        allocated+=1;
        tmp = (char*)realloc(tmp, allocated*sizeof(char));
    }
    tmp[len]='\0';
    int ret=abase_u64k2_vec_sscan(K,w,n,tmp);
    free(tmp);
    return ret;
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_init, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_ur_init(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec_ur * v, unsigned int n)
{
    unsigned int i;
    *v = (abase_u64k2_vec_ur) malloc (n*sizeof(abase_u64k2_elt_ur));
    for(i = 0; i < n; i+=1)
        abase_u64k2_elt_ur_init(K, &( (*v)[i]));
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_reinit, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_ur_reinit(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec_ur * v, unsigned int n, unsigned int m)
{
    if (n < m) { // increase size
        *v = (abase_u64k2_vec_ur) realloc (*v, m * sizeof(abase_u64k2_elt_ur));
        unsigned int i;
        for(i = n; i < m; i+=1)
            abase_u64k2_elt_ur_init(K, (*v) + i);
    } else if (m < n) { // decrease size
        unsigned int i;
        for(i = m; i < n; i+=1)
            abase_u64k2_elt_ur_clear(K, (*v) + i);
        *v = (abase_u64k2_vec_ur) realloc (*v, m * sizeof(abase_u64k2_elt_ur));
    }
}

/* *Mpfq::defaults::vec::alloc::code_for_vec_ur_clear, Mpfq::defaults::vec, Mpfq::defaults */
void abase_u64k2_vec_ur_clear(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_vec_ur * v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k2_elt_ur_clear(K, &( (*v)[i]));
    free(*v);
}

/* missing vec_scal_mul_ur */
/* missing vec_conv_ur */
/* missing vec_reduce */

/* Functions related to SIMD operation */
/* *simd_dotprod::code_for_dotprod */
void abase_u64k2_dotprod(abase_u64k2_dst_field K MAYBE_UNUSED, abase_u64k2_dst_vec xw, abase_u64k2_src_vec xu1, abase_u64k2_src_vec xu0, unsigned int n)
{
    uint64_t * w = xw[0];
    const uint64_t * u0 = xu0[0];
    const uint64_t * u1 = xu1[0];
    dotprod_64K_128(w,u0,u1,n,2);
}


/* Member templates related to SIMD operation */

/* MPI interface */
static void abase_u64k2_mpi_op_inner(void *, void *, int *, MPI_Datatype *);
static /* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_u64k2_mpi_op_inner(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    abase_u64k2_dst_field K;
    MPI_Type_get_attr(*datatype, abase_u64k2_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    abase_u64k2_vec_add(K, inoutvec, inoutvec, invec, *len);
}

static void abase_u64k2_mpi_op_inner_ur(void *, void *, int *, MPI_Datatype *);
static /* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_u64k2_mpi_op_inner_ur(void * invec, void * inoutvec, int * len, MPI_Datatype * datatype)
{
    int got_it;
    abase_u64k2_dst_field K;
    MPI_Type_get_attr(*datatype, abase_u64k2_impl_mpi_attr, (void*) &K, &got_it);
    assert(got_it);
    abase_u64k2_vec_ur_add(K, inoutvec, inoutvec, invec, *len);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_init */
void abase_u64k2_mpi_ops_init(abase_u64k2_dst_field K MAYBE_UNUSED)
{
        if (abase_u64k2_impl_mpi_use_count++) return;
    MPI_Type_create_keyval(MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN, &abase_u64k2_impl_mpi_attr, NULL);
    MPI_Type_contiguous(sizeof(abase_u64k2_elt), MPI_BYTE, &abase_u64k2_impl_mpi_datatype);
    MPI_Type_commit(&abase_u64k2_impl_mpi_datatype);
    MPI_Type_contiguous(sizeof(abase_u64k2_elt_ur), MPI_BYTE, &abase_u64k2_impl_mpi_datatype_ur);
    MPI_Type_commit(&abase_u64k2_impl_mpi_datatype_ur);
    MPI_Type_set_attr(abase_u64k2_impl_mpi_datatype, abase_u64k2_impl_mpi_attr, K);
    MPI_Type_set_attr(abase_u64k2_impl_mpi_datatype_ur, abase_u64k2_impl_mpi_attr, K);
    /* 1 here indicates that our operation is always taken to be
     * commutative */
    MPI_Op_create(&abase_u64k2_mpi_op_inner, 1, &abase_u64k2_impl_mpi_addition_op);
    MPI_Op_create(&abase_u64k2_mpi_op_inner_ur, 1, &abase_u64k2_impl_mpi_addition_op_ur);
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype */
MPI_Datatype abase_u64k2_mpi_datatype(abase_u64k2_dst_field K MAYBE_UNUSED)
{
    return abase_u64k2_impl_mpi_datatype;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_datatype_ur */
MPI_Datatype abase_u64k2_mpi_datatype_ur(abase_u64k2_dst_field K MAYBE_UNUSED)
{
    return abase_u64k2_impl_mpi_datatype_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op */
MPI_Op abase_u64k2_mpi_addition_op(abase_u64k2_dst_field K MAYBE_UNUSED)
{
    return abase_u64k2_impl_mpi_addition_op;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_addition_op_ur */
MPI_Op abase_u64k2_mpi_addition_op_ur(abase_u64k2_dst_field K MAYBE_UNUSED)
{
    return abase_u64k2_impl_mpi_addition_op_ur;
}

/* *Mpfq::defaults::mpi_flat::code_for_mpi_ops_clear */
void abase_u64k2_mpi_ops_clear(abase_u64k2_dst_field K MAYBE_UNUSED)
{
        if (--abase_u64k2_impl_mpi_use_count) return;
    MPI_Op_free(&abase_u64k2_impl_mpi_addition_op);
    MPI_Op_free(&abase_u64k2_impl_mpi_addition_op_ur);
    MPI_Type_delete_attr(abase_u64k2_impl_mpi_datatype, abase_u64k2_impl_mpi_attr);
    MPI_Type_delete_attr(abase_u64k2_impl_mpi_datatype_ur, abase_u64k2_impl_mpi_attr);
    MPI_Type_free(&abase_u64k2_impl_mpi_datatype);
    MPI_Type_free(&abase_u64k2_impl_mpi_datatype_ur);
    MPI_Type_free_keyval(&abase_u64k2_impl_mpi_attr);
}


/* Object-oriented interface */
static void abase_u64k2_wrapper_field_characteristic(abase_vbase_ptr, mpz_t);
static void abase_u64k2_wrapper_field_characteristic(abase_vbase_ptr vbase MAYBE_UNUSED, mpz_t z MAYBE_UNUSED)
{
    abase_u64k2_field_characteristic(vbase->obj, z);
}

static int abase_u64k2_wrapper_field_degree(abase_vbase_ptr);
static int abase_u64k2_wrapper_field_degree(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_field_degree(vbase->obj);
}

static void abase_u64k2_wrapper_field_init(abase_vbase_ptr);
static void abase_u64k2_wrapper_field_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_field_init(vbase->obj);
}

static void abase_u64k2_wrapper_field_clear(abase_vbase_ptr);
static void abase_u64k2_wrapper_field_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_field_clear(vbase->obj);
}

static void abase_u64k2_wrapper_field_specify(abase_vbase_ptr, unsigned long, void *);
static void abase_u64k2_wrapper_field_specify(abase_vbase_ptr vbase MAYBE_UNUSED, unsigned long tag MAYBE_UNUSED, void * x MAYBE_UNUSED)
{
    abase_u64k2_field_specify(vbase->obj, tag, x);
}

static void abase_u64k2_wrapper_field_setopt(abase_vbase_ptr, unsigned long, void *);
static void abase_u64k2_wrapper_field_setopt(abase_vbase_ptr vbase MAYBE_UNUSED, unsigned long x MAYBE_UNUSED, void * y MAYBE_UNUSED)
{
    abase_u64k2_field_setopt(vbase->obj, x, y);
}

static void abase_u64k2_wrapper_init(abase_vbase_ptr, abase_u64k2_elt *);
static void abase_u64k2_wrapper_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_elt * px MAYBE_UNUSED)
{
    abase_u64k2_init(vbase->obj, px);
}

static void abase_u64k2_wrapper_clear(abase_vbase_ptr, abase_u64k2_elt *);
static void abase_u64k2_wrapper_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_elt * px MAYBE_UNUSED)
{
    abase_u64k2_clear(vbase->obj, px);
}

static void abase_u64k2_wrapper_set(abase_vbase_ptr, abase_u64k2_dst_elt, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, abase_u64k2_src_elt s MAYBE_UNUSED)
{
    abase_u64k2_set(vbase->obj, r, s);
}

static void abase_u64k2_wrapper_set_zero(abase_vbase_ptr, abase_u64k2_dst_elt);
static void abase_u64k2_wrapper_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED)
{
    abase_u64k2_set_zero(vbase->obj, r);
}

static void abase_u64k2_wrapper_random(abase_vbase_ptr, abase_u64k2_dst_elt);
static void abase_u64k2_wrapper_random(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED)
{
    abase_u64k2_random(vbase->obj, r);
}

static void abase_u64k2_wrapper_add(abase_vbase_ptr, abase_u64k2_dst_elt, abase_u64k2_src_elt, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, abase_u64k2_src_elt s1 MAYBE_UNUSED, abase_u64k2_src_elt s2 MAYBE_UNUSED)
{
    abase_u64k2_add(vbase->obj, r, s1, s2);
}

static void abase_u64k2_wrapper_sub(abase_vbase_ptr, abase_u64k2_dst_elt, abase_u64k2_src_elt, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, abase_u64k2_src_elt s1 MAYBE_UNUSED, abase_u64k2_src_elt s2 MAYBE_UNUSED)
{
    abase_u64k2_sub(vbase->obj, r, s1, s2);
}

static void abase_u64k2_wrapper_neg(abase_vbase_ptr, abase_u64k2_dst_elt, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, abase_u64k2_src_elt s MAYBE_UNUSED)
{
    abase_u64k2_neg(vbase->obj, r, s);
}

static void abase_u64k2_wrapper_elt_ur_init(abase_vbase_ptr, abase_u64k2_elt_ur *);
static void abase_u64k2_wrapper_elt_ur_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_elt_ur * px MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_init(vbase->obj, px);
}

static void abase_u64k2_wrapper_elt_ur_clear(abase_vbase_ptr, abase_u64k2_elt_ur *);
static void abase_u64k2_wrapper_elt_ur_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_elt_ur * px MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_clear(vbase->obj, px);
}

static void abase_u64k2_wrapper_elt_ur_set(abase_vbase_ptr, abase_u64k2_dst_elt_ur, abase_u64k2_src_elt_ur);
static void abase_u64k2_wrapper_elt_ur_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur r MAYBE_UNUSED, abase_u64k2_src_elt_ur s MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_set(vbase->obj, r, s);
}

static void abase_u64k2_wrapper_elt_ur_set_zero(abase_vbase_ptr, abase_u64k2_dst_elt_ur);
static void abase_u64k2_wrapper_elt_ur_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur r MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_set_zero(vbase->obj, r);
}

static void abase_u64k2_wrapper_elt_ur_add(abase_vbase_ptr, abase_u64k2_dst_elt_ur, abase_u64k2_src_elt_ur, abase_u64k2_src_elt_ur);
static void abase_u64k2_wrapper_elt_ur_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur r MAYBE_UNUSED, abase_u64k2_src_elt_ur s1 MAYBE_UNUSED, abase_u64k2_src_elt_ur s2 MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_add(vbase->obj, r, s1, s2);
}

static void abase_u64k2_wrapper_elt_ur_neg(abase_vbase_ptr, abase_u64k2_dst_elt_ur, abase_u64k2_src_elt_ur);
static void abase_u64k2_wrapper_elt_ur_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur r MAYBE_UNUSED, abase_u64k2_src_elt_ur s MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_neg(vbase->obj, r, s);
}

static void abase_u64k2_wrapper_elt_ur_sub(abase_vbase_ptr, abase_u64k2_dst_elt_ur, abase_u64k2_src_elt_ur, abase_u64k2_src_elt_ur);
static void abase_u64k2_wrapper_elt_ur_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur r MAYBE_UNUSED, abase_u64k2_src_elt_ur s1 MAYBE_UNUSED, abase_u64k2_src_elt_ur s2 MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_sub(vbase->obj, r, s1, s2);
}

static int abase_u64k2_wrapper_cmp(abase_vbase_ptr, abase_u64k2_src_elt, abase_u64k2_src_elt);
static int abase_u64k2_wrapper_cmp(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_elt r MAYBE_UNUSED, abase_u64k2_src_elt s MAYBE_UNUSED)
{
    return abase_u64k2_cmp(vbase->obj, r, s);
}

static int abase_u64k2_wrapper_is_zero(abase_vbase_ptr, abase_u64k2_src_elt);
static int abase_u64k2_wrapper_is_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_elt r MAYBE_UNUSED)
{
    return abase_u64k2_is_zero(vbase->obj, r);
}

static void abase_u64k2_wrapper_asprint(abase_vbase_ptr, char * *, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_asprint(abase_vbase_ptr vbase MAYBE_UNUSED, char * * ps MAYBE_UNUSED, abase_u64k2_src_elt x MAYBE_UNUSED)
{
    abase_u64k2_asprint(vbase->obj, ps, x);
}

static void abase_u64k2_wrapper_fprint(abase_vbase_ptr, FILE *, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_fprint(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_u64k2_src_elt x MAYBE_UNUSED)
{
    abase_u64k2_fprint(vbase->obj, file, x);
}

static void abase_u64k2_wrapper_print(abase_vbase_ptr, abase_u64k2_src_elt);
static void abase_u64k2_wrapper_print(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_elt x MAYBE_UNUSED)
{
    abase_u64k2_print(vbase->obj, x);
}

static int abase_u64k2_wrapper_sscan(abase_vbase_ptr, abase_u64k2_dst_elt, const char *);
static int abase_u64k2_wrapper_sscan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt z MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return abase_u64k2_sscan(vbase->obj, z, str);
}

static int abase_u64k2_wrapper_fscan(abase_vbase_ptr, FILE *, abase_u64k2_dst_elt);
static int abase_u64k2_wrapper_fscan(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_u64k2_dst_elt z MAYBE_UNUSED)
{
    return abase_u64k2_fscan(vbase->obj, file, z);
}

static int abase_u64k2_wrapper_scan(abase_vbase_ptr, abase_u64k2_dst_elt);
static int abase_u64k2_wrapper_scan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt x MAYBE_UNUSED)
{
    return abase_u64k2_scan(vbase->obj, x);
}

static void abase_u64k2_wrapper_vec_init(abase_vbase_ptr, abase_u64k2_vec *, unsigned int);
static void abase_u64k2_wrapper_vec_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_init(vbase->obj, v, n);
}

static void abase_u64k2_wrapper_vec_reinit(abase_vbase_ptr, abase_u64k2_vec *, unsigned int, unsigned int);
static void abase_u64k2_wrapper_vec_reinit(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_u64k2_vec_reinit(vbase->obj, v, n, m);
}

static void abase_u64k2_wrapper_vec_clear(abase_vbase_ptr, abase_u64k2_vec *, unsigned int);
static void abase_u64k2_wrapper_vec_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_clear(vbase->obj, v, n);
}

static void abase_u64k2_wrapper_vec_set(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec r MAYBE_UNUSED, abase_u64k2_src_vec s MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_set(vbase->obj, r, s, n);
}

static void abase_u64k2_wrapper_vec_set_zero(abase_vbase_ptr, abase_u64k2_dst_vec, unsigned int);
static void abase_u64k2_wrapper_vec_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_set_zero(vbase->obj, r, n);
}

static void abase_u64k2_wrapper_vec_setcoef(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_elt, unsigned int);
static void abase_u64k2_wrapper_vec_setcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, abase_u64k2_src_elt x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_u64k2_vec_setcoef(vbase->obj, w, x, i);
}

static void abase_u64k2_wrapper_vec_getcoef(abase_vbase_ptr, abase_u64k2_dst_elt, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_getcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt x MAYBE_UNUSED, abase_u64k2_src_vec w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_u64k2_vec_getcoef(vbase->obj, x, w, i);
}

static void abase_u64k2_wrapper_vec_add(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, abase_u64k2_src_vec u MAYBE_UNUSED, abase_u64k2_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_add(vbase->obj, w, u, v, n);
}

static void abase_u64k2_wrapper_vec_neg(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_neg(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, abase_u64k2_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_neg(vbase->obj, w, u, n);
}

static void abase_u64k2_wrapper_vec_rev(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_rev(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, abase_u64k2_src_vec u MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_rev(vbase->obj, w, u, n);
}

static void abase_u64k2_wrapper_vec_sub(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, abase_u64k2_src_vec u MAYBE_UNUSED, abase_u64k2_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_sub(vbase->obj, w, u, v, n);
}

static void abase_u64k2_wrapper_vec_random(abase_vbase_ptr, abase_u64k2_dst_vec, unsigned int);
static void abase_u64k2_wrapper_vec_random(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_random(vbase->obj, w, n);
}

static int abase_u64k2_wrapper_vec_cmp(abase_vbase_ptr, abase_u64k2_src_vec, abase_u64k2_src_vec, unsigned int);
static int abase_u64k2_wrapper_vec_cmp(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_vec u MAYBE_UNUSED, abase_u64k2_src_vec v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return abase_u64k2_vec_cmp(vbase->obj, u, v, n);
}

static int abase_u64k2_wrapper_vec_is_zero(abase_vbase_ptr, abase_u64k2_src_vec, unsigned int);
static int abase_u64k2_wrapper_vec_is_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_vec r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    return abase_u64k2_vec_is_zero(vbase->obj, r, n);
}

static void abase_u64k2_wrapper_vec_asprint(abase_vbase_ptr, char * *, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_asprint(abase_vbase_ptr vbase MAYBE_UNUSED, char * * pstr MAYBE_UNUSED, abase_u64k2_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_asprint(vbase->obj, pstr, w, n);
}

static void abase_u64k2_wrapper_vec_fprint(abase_vbase_ptr, FILE *, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_fprint(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_u64k2_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_fprint(vbase->obj, file, w, n);
}

static void abase_u64k2_wrapper_vec_print(abase_vbase_ptr, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_vec_print(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_src_vec w MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_print(vbase->obj, w, n);
}

static int abase_u64k2_wrapper_vec_sscan(abase_vbase_ptr, abase_u64k2_vec *, unsigned int *, const char *);
static int abase_u64k2_wrapper_vec_sscan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED, const char * str MAYBE_UNUSED)
{
    return abase_u64k2_vec_sscan(vbase->obj, w, n, str);
}

static int abase_u64k2_wrapper_vec_fscan(abase_vbase_ptr, FILE *, abase_u64k2_vec *, unsigned int *);
static int abase_u64k2_wrapper_vec_fscan(abase_vbase_ptr vbase MAYBE_UNUSED, FILE * file MAYBE_UNUSED, abase_u64k2_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED)
{
    return abase_u64k2_vec_fscan(vbase->obj, file, w, n);
}

static int abase_u64k2_wrapper_vec_scan(abase_vbase_ptr, abase_u64k2_vec *, unsigned int *);
static int abase_u64k2_wrapper_vec_scan(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec * w MAYBE_UNUSED, unsigned int * n MAYBE_UNUSED)
{
    return abase_u64k2_vec_scan(vbase->obj, w, n);
}

static void abase_u64k2_wrapper_vec_ur_init(abase_vbase_ptr, abase_u64k2_vec_ur *, unsigned int);
static void abase_u64k2_wrapper_vec_ur_init(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_init(vbase->obj, v, n);
}

static void abase_u64k2_wrapper_vec_ur_set_zero(abase_vbase_ptr, abase_u64k2_dst_vec_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_set_zero(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec_ur r MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_set_zero(vbase->obj, r, n);
}

static void abase_u64k2_wrapper_vec_ur_reinit(abase_vbase_ptr, abase_u64k2_vec_ur *, unsigned int, unsigned int);
static void abase_u64k2_wrapper_vec_ur_reinit(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, unsigned int m MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_reinit(vbase->obj, v, n, m);
}

static void abase_u64k2_wrapper_vec_ur_clear(abase_vbase_ptr, abase_u64k2_vec_ur *, unsigned int);
static void abase_u64k2_wrapper_vec_ur_clear(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_vec_ur * v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_clear(vbase->obj, v, n);
}

static void abase_u64k2_wrapper_vec_ur_set(abase_vbase_ptr, abase_u64k2_dst_vec_ur, abase_u64k2_src_vec_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_set(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec_ur r MAYBE_UNUSED, abase_u64k2_src_vec_ur s MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_set(vbase->obj, r, s, n);
}

static void abase_u64k2_wrapper_vec_ur_setcoef(abase_vbase_ptr, abase_u64k2_dst_vec_ur, abase_u64k2_src_elt_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_setcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec_ur w MAYBE_UNUSED, abase_u64k2_src_elt_ur x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_setcoef(vbase->obj, w, x, i);
}

static void abase_u64k2_wrapper_vec_ur_getcoef(abase_vbase_ptr, abase_u64k2_dst_elt_ur, abase_u64k2_src_vec_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_getcoef(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt_ur x MAYBE_UNUSED, abase_u64k2_src_vec_ur w MAYBE_UNUSED, unsigned int i MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_getcoef(vbase->obj, x, w, i);
}

static void abase_u64k2_wrapper_vec_ur_add(abase_vbase_ptr, abase_u64k2_dst_vec_ur, abase_u64k2_src_vec_ur, abase_u64k2_src_vec_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_add(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec_ur w MAYBE_UNUSED, abase_u64k2_src_vec_ur u MAYBE_UNUSED, abase_u64k2_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_add(vbase->obj, w, u, v, n);
}

static void abase_u64k2_wrapper_vec_ur_sub(abase_vbase_ptr, abase_u64k2_dst_vec_ur, abase_u64k2_src_vec_ur, abase_u64k2_src_vec_ur, unsigned int);
static void abase_u64k2_wrapper_vec_ur_sub(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec_ur w MAYBE_UNUSED, abase_u64k2_src_vec_ur u MAYBE_UNUSED, abase_u64k2_src_vec_ur v MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_vec_ur_sub(vbase->obj, w, u, v, n);
}

static ptrdiff_t abase_u64k2_wrapper_vec_elt_stride(abase_vbase_ptr, int);
static ptrdiff_t abase_u64k2_wrapper_vec_elt_stride(abase_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return abase_u64k2_vec_elt_stride(vbase->obj, n);
}

static int abase_u64k2_wrapper_groupsize(abase_vbase_ptr);
static int abase_u64k2_wrapper_groupsize(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_groupsize(vbase->obj);
}

static int abase_u64k2_wrapper_offset(abase_vbase_ptr, int);
static int abase_u64k2_wrapper_offset(abase_vbase_ptr vbase MAYBE_UNUSED, int n MAYBE_UNUSED)
{
    return abase_u64k2_offset(vbase->obj, n);
}

static int abase_u64k2_wrapper_stride(abase_vbase_ptr);
static int abase_u64k2_wrapper_stride(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_stride(vbase->obj);
}

static void abase_u64k2_wrapper_set_ui_at(abase_vbase_ptr, abase_u64k2_dst_elt, int, unsigned long);
static void abase_u64k2_wrapper_set_ui_at(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_u64k2_set_ui_at(vbase->obj, p, k, v);
}

static void abase_u64k2_wrapper_set_ui_all(abase_vbase_ptr, abase_u64k2_dst_elt, unsigned long);
static void abase_u64k2_wrapper_set_ui_all(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_u64k2_set_ui_all(vbase->obj, r, v);
}

static void abase_u64k2_wrapper_elt_ur_set_ui_at(abase_vbase_ptr, abase_u64k2_dst_elt, int, unsigned long);
static void abase_u64k2_wrapper_elt_ur_set_ui_at(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt p MAYBE_UNUSED, int k MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_set_ui_at(vbase->obj, p, k, v);
}

static void abase_u64k2_wrapper_elt_ur_set_ui_all(abase_vbase_ptr, abase_u64k2_dst_elt, unsigned long);
static void abase_u64k2_wrapper_elt_ur_set_ui_all(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_elt r MAYBE_UNUSED, unsigned long v MAYBE_UNUSED)
{
    abase_u64k2_elt_ur_set_ui_all(vbase->obj, r, v);
}

static void abase_u64k2_wrapper_dotprod(abase_vbase_ptr, abase_u64k2_dst_vec, abase_u64k2_src_vec, abase_u64k2_src_vec, unsigned int);
static void abase_u64k2_wrapper_dotprod(abase_vbase_ptr vbase MAYBE_UNUSED, abase_u64k2_dst_vec xw MAYBE_UNUSED, abase_u64k2_src_vec xu1 MAYBE_UNUSED, abase_u64k2_src_vec xu0 MAYBE_UNUSED, unsigned int n MAYBE_UNUSED)
{
    abase_u64k2_dotprod(vbase->obj, xw, xu1, xu0, n);
}

static void abase_u64k2_wrapper_mpi_ops_init(abase_vbase_ptr);
static void abase_u64k2_wrapper_mpi_ops_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_mpi_ops_init(vbase->obj);
}

static MPI_Datatype abase_u64k2_wrapper_mpi_datatype(abase_vbase_ptr);
static MPI_Datatype abase_u64k2_wrapper_mpi_datatype(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_mpi_datatype(vbase->obj);
}

static MPI_Datatype abase_u64k2_wrapper_mpi_datatype_ur(abase_vbase_ptr);
static MPI_Datatype abase_u64k2_wrapper_mpi_datatype_ur(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_mpi_datatype_ur(vbase->obj);
}

static MPI_Op abase_u64k2_wrapper_mpi_addition_op(abase_vbase_ptr);
static MPI_Op abase_u64k2_wrapper_mpi_addition_op(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_mpi_addition_op(vbase->obj);
}

static MPI_Op abase_u64k2_wrapper_mpi_addition_op_ur(abase_vbase_ptr);
static MPI_Op abase_u64k2_wrapper_mpi_addition_op_ur(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_mpi_addition_op_ur(vbase->obj);
}

static void abase_u64k2_wrapper_mpi_ops_clear(abase_vbase_ptr);
static void abase_u64k2_wrapper_mpi_ops_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_mpi_ops_clear(vbase->obj);
}

static const char * abase_u64k2_wrapper_oo_impl_name(abase_vbase_ptr);
static const char * abase_u64k2_wrapper_oo_impl_name(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    return abase_u64k2_oo_impl_name(vbase);
}

static void abase_u64k2_wrapper_oo_field_init(abase_vbase_ptr);
static void abase_u64k2_wrapper_oo_field_init(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_oo_field_init(vbase);
}

static void abase_u64k2_wrapper_oo_field_clear(abase_vbase_ptr);
static void abase_u64k2_wrapper_oo_field_clear(abase_vbase_ptr vbase MAYBE_UNUSED)
{
    abase_u64k2_oo_field_clear(vbase);
}

void abase_u64k2_oo_field_init(abase_vbase_ptr vbase)
{
    memset(vbase, 0, sizeof(struct abase_vbase_s));
    vbase->obj = malloc(sizeof(abase_u64k2_field));
    abase_u64k2_field_init((abase_u64k2_dst_field) vbase->obj);
    vbase->field_characteristic = (void (*) (abase_vbase_ptr, mpz_t)) abase_u64k2_wrapper_field_characteristic;
    vbase->field_degree = (int (*) (abase_vbase_ptr)) abase_u64k2_wrapper_field_degree;
    vbase->field_init = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_field_init;
    vbase->field_clear = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_field_clear;
    vbase->field_specify = (void (*) (abase_vbase_ptr, unsigned long, void *)) abase_u64k2_wrapper_field_specify;
    vbase->field_setopt = (void (*) (abase_vbase_ptr, unsigned long, void *)) abase_u64k2_wrapper_field_setopt;
    vbase->init = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_init;
    vbase->clear = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_clear;
    vbase->set = (void (*) (abase_vbase_ptr, void *, const void *)) abase_u64k2_wrapper_set;
    /* missing set_ui */
    vbase->set_zero = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_set_zero;
    /* missing get_ui */
    /* missing set_mpn */
    /* missing set_mpz */
    /* missing get_mpn */
    /* missing get_mpz */
    vbase->random = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_random;
    vbase->add = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_u64k2_wrapper_add;
    vbase->sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_u64k2_wrapper_sub;
    vbase->neg = (void (*) (abase_vbase_ptr, void *, const void *)) abase_u64k2_wrapper_neg;
    /* missing mul */
    /* missing sqr */
    /* missing is_sqr */
    /* missing sqrt */
    /* missing pow */
    /* missing frobenius */
    /* missing add_ui */
    /* missing sub_ui */
    /* missing mul_ui */
    /* missing inv */
    vbase->elt_ur_init = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_elt_ur_init;
    vbase->elt_ur_clear = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_elt_ur_clear;
    vbase->elt_ur_set = (void (*) (abase_vbase_ptr, void *, const void *)) abase_u64k2_wrapper_elt_ur_set;
    vbase->elt_ur_set_zero = (void (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_elt_ur_set_zero;
    /* missing elt_ur_set_ui */
    vbase->elt_ur_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_u64k2_wrapper_elt_ur_add;
    vbase->elt_ur_neg = (void (*) (abase_vbase_ptr, void *, const void *)) abase_u64k2_wrapper_elt_ur_neg;
    vbase->elt_ur_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *)) abase_u64k2_wrapper_elt_ur_sub;
    /* missing mul_ur */
    /* missing sqr_ur */
    /* missing reduce */
    vbase->cmp = (int (*) (abase_vbase_ptr, const void *, const void *)) abase_u64k2_wrapper_cmp;
    /* missing cmp_ui */
    vbase->is_zero = (int (*) (abase_vbase_ptr, const void *)) abase_u64k2_wrapper_is_zero;
    vbase->asprint = (void (*) (abase_vbase_ptr, char * *, const void *)) abase_u64k2_wrapper_asprint;
    vbase->fprint = (void (*) (abase_vbase_ptr, FILE *, const void *)) abase_u64k2_wrapper_fprint;
    vbase->print = (void (*) (abase_vbase_ptr, const void *)) abase_u64k2_wrapper_print;
    vbase->sscan = (int (*) (abase_vbase_ptr, void *, const char *)) abase_u64k2_wrapper_sscan;
    vbase->fscan = (int (*) (abase_vbase_ptr, FILE *, void *)) abase_u64k2_wrapper_fscan;
    vbase->scan = (int (*) (abase_vbase_ptr, void *)) abase_u64k2_wrapper_scan;
    vbase->vec_init = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_init;
    vbase->vec_reinit = (void (*) (abase_vbase_ptr, void *, unsigned int, unsigned int)) abase_u64k2_wrapper_vec_reinit;
    vbase->vec_clear = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_clear;
    vbase->vec_set = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_set;
    vbase->vec_set_zero = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_set_zero;
    vbase->vec_setcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_setcoef;
    /* missing vec_setcoef_ui */
    vbase->vec_getcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_getcoef;
    vbase->vec_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_add;
    vbase->vec_neg = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_neg;
    vbase->vec_rev = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_rev;
    vbase->vec_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_sub;
    /* missing vec_scal_mul */
    /* missing vec_conv */
    vbase->vec_random = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_random;
    vbase->vec_cmp = (int (*) (abase_vbase_ptr, const void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_cmp;
    vbase->vec_is_zero = (int (*) (abase_vbase_ptr, const void *, unsigned int)) abase_u64k2_wrapper_vec_is_zero;
    vbase->vec_asprint = (void (*) (abase_vbase_ptr, char * *, const void *, unsigned int)) abase_u64k2_wrapper_vec_asprint;
    vbase->vec_fprint = (void (*) (abase_vbase_ptr, FILE *, const void *, unsigned int)) abase_u64k2_wrapper_vec_fprint;
    vbase->vec_print = (void (*) (abase_vbase_ptr, const void *, unsigned int)) abase_u64k2_wrapper_vec_print;
    vbase->vec_sscan = (int (*) (abase_vbase_ptr, void *, unsigned int *, const char *)) abase_u64k2_wrapper_vec_sscan;
    vbase->vec_fscan = (int (*) (abase_vbase_ptr, FILE *, void *, unsigned int *)) abase_u64k2_wrapper_vec_fscan;
    vbase->vec_scan = (int (*) (abase_vbase_ptr, void *, unsigned int *)) abase_u64k2_wrapper_vec_scan;
    vbase->vec_ur_init = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_ur_init;
    vbase->vec_ur_set_zero = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_ur_set_zero;
    vbase->vec_ur_reinit = (void (*) (abase_vbase_ptr, void *, unsigned int, unsigned int)) abase_u64k2_wrapper_vec_ur_reinit;
    vbase->vec_ur_clear = (void (*) (abase_vbase_ptr, void *, unsigned int)) abase_u64k2_wrapper_vec_ur_clear;
    vbase->vec_ur_set = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_ur_set;
    vbase->vec_ur_setcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_ur_setcoef;
    vbase->vec_ur_getcoef = (void (*) (abase_vbase_ptr, void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_ur_getcoef;
    vbase->vec_ur_add = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_ur_add;
    vbase->vec_ur_sub = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_wrapper_vec_ur_sub;
    /* missing vec_scal_mul_ur */
    /* missing vec_conv_ur */
    /* missing vec_reduce */
    vbase->vec_elt_stride = (ptrdiff_t (*) (abase_vbase_ptr, int)) abase_u64k2_wrapper_vec_elt_stride;
    vbase->groupsize = (int (*) (abase_vbase_ptr)) abase_u64k2_wrapper_groupsize;
    vbase->offset = (int (*) (abase_vbase_ptr, int)) abase_u64k2_wrapper_offset;
    vbase->stride = (int (*) (abase_vbase_ptr)) abase_u64k2_wrapper_stride;
    vbase->set_ui_at = (void (*) (abase_vbase_ptr, void *, int, unsigned long)) abase_u64k2_wrapper_set_ui_at;
    vbase->set_ui_all = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_u64k2_wrapper_set_ui_all;
    vbase->elt_ur_set_ui_at = (void (*) (abase_vbase_ptr, void *, int, unsigned long)) abase_u64k2_wrapper_elt_ur_set_ui_at;
    vbase->elt_ur_set_ui_all = (void (*) (abase_vbase_ptr, void *, unsigned long)) abase_u64k2_wrapper_elt_ur_set_ui_all;
    vbase->dotprod = (void (*) (abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_wrapper_dotprod;
    vbase->mpi_ops_init = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_ops_init;
    vbase->mpi_datatype = (MPI_Datatype (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_datatype;
    vbase->mpi_datatype_ur = (MPI_Datatype (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_datatype_ur;
    vbase->mpi_addition_op = (MPI_Op (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_addition_op;
    vbase->mpi_addition_op_ur = (MPI_Op (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_addition_op_ur;
    vbase->mpi_ops_clear = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_mpi_ops_clear;
    vbase->oo_impl_name = (const char * (*) (abase_vbase_ptr)) abase_u64k2_wrapper_oo_impl_name;
    vbase->oo_field_init = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_oo_field_init;
    vbase->oo_field_clear = (void (*) (abase_vbase_ptr)) abase_u64k2_wrapper_oo_field_clear;
}


/* vim:set ft=cpp: */
