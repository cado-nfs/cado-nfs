#include "cado.h"       /* HAVE_* macros ! */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

#include "portability.h"
#include "macros.h"
#include "matops.h"

/*{{{ random data...*/
#if 0
static int urandom_fd = -1;

static inline uint64_t rand64()
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    uint64_t r;
    size_t s = read(urandom_fd, &r, sizeof(uint64_t));
    if (s != sizeof(uint64_t)) abort();
    return r;
}

static inline void memfill_random(uint64_t * ptr, size_t z)
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    size_t s = read(urandom_fd, ptr, z);
    if (s != z) abort();
    // close(urandom_fd);
}
#else
static inline uint64_t rand64()
{
    uint64_t r0 = rand();
    uint64_t r1 = rand(); r1 <<= 22;
    uint64_t r2 = rand(); r2 <<= 44;
    return r0 ^r1 ^ r2;
}

static inline void memfill_random(void *p, size_t s)
{
    size_t w = sizeof(uint64_t);
    uint64_t * pw = (uint64_t *)p;
    for( ; s >= w ; pw++, s -= w) *pw = rand64();
    char * pc = (char*) pw;
    for( ; s ; pc++, s--) *pc = rand();
}

#endif/*}}}*/

int test_accel = 1;

#define t_and_unit_from_clock__bare(t, unit, t1, j)                     \
    double t = t1;							\
    t /= CLOCKS_PER_SEC;						\
    if (j) t /= j; else t = 0;						\
    const char * unit = "s";						\
    if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			        \
    } else if (t < 1.0e-4) { unit = "micros"; t *= 1.0e6;		\
    } else if (t < 1.0e-1) { unit = "ms"; t *= 1.0e3; }                 \
    do { } while (0)

#define TIME1__bare(maxtime, what, args) 		        	\
    clock_t measuring_time = maxtime * CLOCKS_PER_SEC / test_accel;	\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    for (j = 0; ; j++) {						\
        what args;							\
        t1 = clock() - t0;						\
        if (j && t1 > measuring_time)					\
            break;							\
    }									\
    t_and_unit_from_clock__bare(t, unit, t1, j);

#define TIME1(maxtime, what, args) do {			        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what " \t%d times in %.4f %s each\n",       		\
            j, t, unit);		                        	\
} while (0)

#define TIME1N(maxtime, what, args) do {		        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what "(n=%d) \t%d times in %.4f %s each\n", n,     	        \
            j, t, unit);		                        	\
} while (0)

#define TIME1N_spins(rexpr, maxtime, what, args, spinexpr, spinmax) do {        \
    clock_t ts[spinmax];						\
    int ns[spinmax];                                                    \
    for(int s = 0 ; s < spinmax ; s++) ts[s] = ns[s] = 0;		\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    clock_t t = t0;                                                     \
    clock_t fence = t0 + maxtime * CLOCKS_PER_SEC / test_accel;		\
    for (j = 0; ; j++) {						\
        rexpr;                                                          \
        int ret = what args;						\
        int s = spinexpr;                                               \
        ts[s] += (t1 = clock()) - t;                                    \
        ns[s] ++;                                                       \
        t = t1;                                                         \
        if (j && t1 > fence)					        \
            break;							\
    }									\
    int nch=0;                                                          \
    for(int s = 0 ; s < spinmax ; s++) {				\
        if (s == 0) nch=printf(#what"(n=%d)", n);			\
        else for(int k = nch ; k-- ; putchar(' '));		        \
        t_and_unit_from_clock__bare(t, unit, ts[s], ns[s]);		\
        printf(" \t[%d] %d times in %.4f %s each\n",s,ns[s],t,unit);	\
    }									\
} while (0)



/* {{{ test routines */
struct l1_data_s {
    uint64_t * xr;
    uint64_t * r;
    uint64_t * a;
    uint64_t * w;
    uint64_t * wt;
#ifdef  HAVE_M4RI
    mzd_t *R;
    mzd_t *A;
    mzd_t *W;
    mzd_t *WT;
#endif  /* HAVE_M4RI */
    unsigned int n;
};

typedef struct l1_data_s l1_data[1];
typedef struct l1_data_s * l1_data_ptr;
typedef const struct l1_data_s * l1_data_srcptr;

void l1_data_init_set(l1_data_ptr D, unsigned int n)
{
    D->n = n;
    D->xr = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->r = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->a = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->w = (uint64_t *) malloc(64 * sizeof(uint64_t));
    D->wt = (uint64_t *) malloc(64 * sizeof(uint64_t));
#ifdef  HAVE_M4RI
    D->R = mzd_init(n, 64);
    D->A = mzd_init(n, 64);
    D->W = mzd_init(64, 64);
    D->WT = mzd_init(64, 64);
#endif  /* HAVE_M4RI */

    memfill_random(D->a, (n) * sizeof(uint64_t));
    memfill_random(D->w, (64) * sizeof(uint64_t));
    transp_6464(D->wt, D->w);
#ifdef  HAVE_M4RI
    mzd_set_mem(D->A, D->a, n);
    mzd_set_mem(D->W, D->w, 64);
    mzd_set_mem(D->WT, D->wt, 64);
#endif  /* HAVE_M4RI */
}

void l1_data_clear(l1_data_ptr D)
{
    free(D->xr);
    free(D->r);
    free(D->a);
    free(D->w);
    free(D->wt);
#ifdef  HAVE_M4RI
    mzd_free(D->R);
    mzd_free(D->A);
    mzd_free(D->W);
    mzd_free(D->WT);
#endif  /* HAVE_M4RI */
}

#define SCOPE_L1_DATA_MEMBERS_STANDARD(D)				\
    unsigned int n __attribute__((unused)) = D->n;			\
    uint64_t * xr __attribute__((unused)) = D->xr;      		\
    uint64_t * r __attribute__((unused)) = D->r;			\
    const uint64_t * a __attribute__((unused)) = D->a;			\
    const uint64_t * w __attribute__((unused)) = D->w;			\
    const uint64_t * wt __attribute__((unused)) = D->wt

#define SCOPE_L1_DATA_MEMBERS_M4RI(D)   				\
    mzd_t * R __attribute__((unused)) = D->R;				\
    mzd_t * A __attribute__((unused)) = D->A;				\
    mzd_t * W __attribute__((unused)) = D->W;				\
    mzd_t * WT __attribute__((unused)) = D->WT

#ifdef  HAVE_M4RI
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D);      \
        SCOPE_L1_DATA_MEMBERS_M4RI(D)
#else  /* HAVE_M4RI */
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D)
#endif  /* HAVE_M4RI */

void level1_basic_tests()
{
    l1_data D;
    l1_data_init_set(D, 64);
    SCOPE_L1_DATA_MEMBERS(D);

    /* transposition */
    transp_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    transp_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, transp_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_transpose(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_transpose, (R, A));
#endif  /* HAVE_M4RI */

    /* copy */
    copy_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    copy_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, copy_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_copy(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_copy, (R, A));
#endif  /* HAVE_M4RI */

    /* add */
    add_6464_6464(r, a, w); memcpy(xr, r, 64 * sizeof(uint64_t));

    add_6464_6464_C(r, a, w);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, add_6464_6464_C, (r,a,w));

#ifdef  HAVE_M4RI
    mzd_add(R, A, W);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_add, (R, A, W));
#endif  /* HAVE_M4RI */

    /*******/

    l1_data_clear(D);
}

void level1_mul_tests_N_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);

    mul_N64_6464_vec(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_vec, (r,a,w,n));

    mul_N64_6464_transB(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_transB, (r,a,w,n));

#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_N64_6464_sse(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_sse, (r,a,w,n));
#endif

    mul_N64_6464_lookup4(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_lookup4, (r,a,w,n));

    mul_N64_6464_lookup8(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_lookup8, (r,a,w,n));

    mul_N64_T6464_vec(r, a, wt, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_T6464_vec, (r,a,w,n));

    mul_N64_T6464_transB(r, a, wt, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_T6464_transB, (r,a,w,n));

#ifdef  HAVE_M4RI
    mzd_mul_naive(R, A, W);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_naive, (R, A, W));

    _mzd_mul_naive(R, A, WT, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_naive, (R, A, WT, 1));

    mzd_mul_m4rm(R, A, W, 0);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_m4rm, (R, A, W, 0));

    _mzd_mul_m4rm(R, A, W, 0, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_m4rm, (R, A, W, 0, 1));
#endif  /* HAVE_M4RI */
}


void level1_mul_tests_N(unsigned int n)
{
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_N64_6464_vec(D->r, D->a, D->w, D->n);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}

void level1_mul_tests_1_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 1);

    mul_o64_6464_C_lsb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_lsb, (r, *a, w));

    mul_o64_6464_C_msb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_msb, (r, *a, w));

    mul_o64_T6464_C_parity(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity, (r, *a, wt));

    mul_o64_T6464_C_parity3(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity3, (r, *a, wt));
}

void level1_mul_tests_1()
{
    unsigned int n = 1;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_o64_6464(D->r, *D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_1_list(D);
    /* Functions which can do any n can also do n=1 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}
void level1_mul_tests_64_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 64);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_6464_6464_sse(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_sse, (r, a, w));
#endif

    mul_6464_6464_v2(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_v2, (r, a, w));
}

void level1_mul_tests_64()
{
    unsigned int n = 64;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplicate of two 64x64 matrices */
    mul_6464_6464(D->r, D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));
    
    level1_mul_tests_64_list(D);
    /* Functions which can do any n can also do n=64 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}
/* }}} */

/* PLUQ -- well we're not computing exactly PLUQ 
 * PLUQ says: Any m*n matrix A with rank r , can be written A = P*L*U*Q
 * where P and Q are two permutation matrices, of dimension respectively
 * m*m and n*n, L is m*r unit lower triangular and U is r*n upper
 * triangular.
 *
 * Here we compute p,l,u,q such that p*l*a*transpose(q) = an upper
 * triangular matrix (and u is l*a).
 *
 * p*l is lower triangular.
 */
void check_pluq(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m, int n)
{
    mat64 pm[(n/64)*(n/64)];
    pmat_get_matrix(pm, p);

    pmat qt;
    pmat_init(qt, n);
    pmat_transpose(qt, q);

    mat64 qmt[(n/64)*(n/64)];
    pmat_get_matrix(qmt, qt);

    /* compute p*u*transpose(q) */
    mat64 pu[(n/64)*(n/64)];
    memset(pu, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pm[i*(n/64)+k], u[k*(n/64)+j]);
        add_6464_6464(pu[i*(n/64)+j], pu[i*(n/64)+j], tmp);
    }

    mat64 puq[(n/64)*(n/64)];
    memset(puq, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pu[i*(n/64)+k], qmt[k*(n/64)+j]);
        add_6464_6464(puq[i*(n/64)+j], puq[i*(n/64)+j], tmp);
    }
    
    /* at this point puq = p*u*transpose(q) should be a upper triangular,
     * with normalized diagonal. */
    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_uppertriangular(puq[i*(n/64)+i]));
    }

    mat64 lm[(n/64)*(n/64)];
    memset(lm, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k <= i ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, l[i*(n/64)+k], m[k*(n/64)+j]);
        add_6464_6464(lm[i*(n/64)+j], lm[i*(n/64)+j], tmp);
    }

    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_lowertriangular(l[i*(n/64)+i]));
        ASSERT_ALWAYS(mat64_triangular_is_unit(l[i*(n/64)+i]));
        for(int j = 0 ; j < (n/64) ; j++ ) {
            ASSERT_ALWAYS(mat64_eq(lm[i*(n/64)+j], u[i*(n/64)+j]));
        }
    }
    pmat_clear(qt);
}
/*{{{ more test routines */
void level3_gauss_tests_N(int n __attribute__((unused)))
{
#ifdef  HAVE_M4RI
    mzd_t * M;
    mzd_t * LU;
    mzp_t * P, *Q;
    M = mzd_init(n, n);
#if 0
    mzd_set_mem(M, m, n);
    uint64_t * m = (uint64_t *) malloc((n*n/64)*sizeof(uint64_t));
    memfill_random(m, (64) * sizeof(uint64_t));
    free(m);
#else
    my_mzd_randomize(M);
#endif
    LU = mzd_init(n,n);
    P = mzp_init(n);
    Q = mzp_init(n);
    TIME1N(2, mzd_mypluq, (LU, M, P, Q, 0));
    TIME1N(2, mzd_myechelonize_m4ri, (LU, M, 0, 0));
    TIME1N(2, mzd_myechelonize_pluq, (LU, M, 0));
    mzd_free(M);
    mzd_free(LU);
    mzp_free(P);
    mzp_free(Q);
#endif
}
/*}}}*/

int main(int argc, char * argv[])
{
    unsigned int n = 2 * 1000 * 1000;
    int seed = 0;
    argc--,argv++;
    for( ; argc ; argc--, argv++) {
        if (strcmp(argv[0], "--test-fast") == 0 && argc >= 2) {
            test_accel = atoi(argv[1]);
            argc--,argv++;
        } else if (strcmp(argv[0], "--n") == 0 && argc >= 2) {
            n = atoi(argv[1]);
            argc--,argv++;
        } else if (strcmp(argv[0], "-seed") == 0 && argc >= 2) {
            seed = atoi(argv[1]);
            argc--,argv++;
        } else {
            fprintf(stderr, "arguments: [--test-fast <x>] [--n <n>]\n");
            exit(EXIT_FAILURE);
        }
    }
    if (!seed) seed = rand();
    printf("# seeding with seed %d\n", seed);
#if !defined(HAVE_USUAL_SRAND_DETERMINISTIC_BEHAVIOR) && defined(HAVE_SRAND_DETERMINISTIC)
    if (seed) srand_deterministic(seed);
#else
    /* won't be deterministic if
     * !defined(HAVE_USUAL_SRAND_DETERMINISTIC_BEHAVIOR), but so be it.
     */
    if (seed) srand(seed);
#endif
    printf("## compile-time features\n");
#ifdef HAVE_M4RI
    printf("## HAVE_M4RI\n");
#endif /* HAVE_M4RI */
#ifdef HAVE_M4RIE
    printf("## HAVE_M4RIE\n");
#endif /* HAVE_M4RIE */
#ifdef HAVE_PCLMUL
    printf("## HAVE_PCLMUL\n");
#endif /* HAVE_PCLMUL */
#ifdef HAVE_SSE2
    printf("## HAVE_SSE2\n");
#endif /* HAVE_SSE2 */
#ifdef HAVE_SSE41
    printf("## HAVE_SSE41\n");
#endif /* HAVE_SSE41 */
#ifdef VALGRIND
    printf("## VALGRIND\n");
#endif /* VALGRIND */
    printf("## ULONG_BITS=%d\n", ULONG_BITS);

    if (1) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t a = rand64();
        uint64_t w = rand64();

        printf("-- level-1 benches --\n");
        level1_basic_tests();
        level1_mul_tests_1();
        printf("-- level-1 benches, n=64 --\n");
        level1_mul_tests_64();
        printf("-- level-1 benches, misc --\n");
        TIME1(1, addmul_To64_o64_lsb, (r, a, w));
        TIME1(1, addmul_To64_o64_msb, (r, a, w));
        TIME1(1, addmul_To64_o64_lsb_packof2, (r, a, w));
#if defined(HAVE_SSE2) && ULONG_BITS == 64
        TIME1(1, addmul_To64_o64_lsb_sse_v1, (r, a, w));
#endif
        TIME1(1, addmul_To64_o64, (r, a, w));
        free(r);
    }

    if (1) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t * a = (uint64_t *) malloc(n * sizeof(uint64_t));
        uint64_t * w = (uint64_t *) malloc(n * sizeof(uint64_t));
        memfill_random(a, (n) * sizeof(uint64_t));
        memfill_random(w, (n) * sizeof(uint64_t));

        printf("-- level-2 benches (N=%u) --\n", n);
        level1_mul_tests_N(n);
        TIME1N(1, mul_64N_N64_addmul, (r,a,w,n));
        TIME1N(5, mul_TN32_N64_C, (r,(uint32_t*)a,w,n));
        TIME1N(5, mul_TN64_N64_C, (r,a,w,n));

        free(r); free(a); free(w);
    }

    if (1) {
        mat64 m[4],l[4],u[4];
        {
            pmat p,q;
            pmat_init(p, 64);
            pmat_init(q, 64);
            memfill_random(m, sizeof(mat64));
            PLUQ64(p,l,u,q,m);
            check_pluq(p,l,u,q,m,64);
            pmat_clear(p);
            pmat_clear(q);
        }

        {
            pmat p,q;
            pmat_init(p, 128);
            pmat_init(q, 128);
            memfill_random(m, 4 * sizeof(mat64));
            PLUQ128(p,l,u,q,m);
            check_pluq(p,l,u,q,m,128);
            pmat_clear(p);
            pmat_clear(q);
        }

        printf("PLUQ128 stub executed ok\n");


            // pmat_6464(m[0]); printf("\n");
            // pmat_6464(u); printf("\n");
            // pmat_6464(l); printf("\n");
            // pmat_6464(p); printf("\n");
            // pmat_6464(t); printf("\n");
            

            // mul_N64_T6464(t,u,p,64);
            // pmat_6464(t); printf("\n");

            
#if 0
            memset(e,0,sizeof(e));
            memfill_random(m[0], (64) * sizeof(uint32_t));
            memfill_random(m[1], (64) * sizeof(uint32_t));
            memfill_random(m[2], (64) * sizeof(uint32_t));
            memfill_random(m[3], (64) * sizeof(uint32_t));
            // pmat_6464(m[0]); printf("\n");
            full_echelon_6464_imm(mm,e[0],m[0]);

            /* 
            mul_6464

            pmat_6464(mm); printf("\n");
            pmat_6464(e); printf("\n"); printf("\n");
            */
#endif
    }

    if (1) {
        mat64 m;
        mat64 e;
        mat64 mm;
        mat64 l,u,p;
        memfill_random(m, sizeof(mat64));
        memfill_random(e, sizeof(mat64));
        memfill_random(mm, sizeof(mat64));
        mat64 m4[4];
        mat64 u4[4];
        memfill_random(m4, 4 * sizeof(mat64));
        // printf("-- for reference: best matrix mult, 64x64 --\n");
        // TIME1(2, mul_6464_6464, (mm,e,m));
        // TIME1(2, mul_N64_T6464, (mm,e,m,64));
        printf("-- level-3 (reduction) benches, n=64 --\n");
        // TIME1(2, gauss_6464_C, (mm,e,m));
        // TIME1(2, gauss_6464_imm, (mm,e,m));
        // TIME1(2, PLUQ64_inner, (NULL,l,u,m,0));
        int phi[128];
        {
            pmat p,q;
            mat64 m[4],l[4],u[4];
            memfill_random(m, 4 * sizeof(mat64));
            pmat_init(p, 128);
            pmat_init(q, 128);
            TIME1(2, PLUQ128, (p,l,u,q,m));
            pmat_clear(p);
            pmat_clear(q);
        }
        int n=2;
        TIME1N(2, memfill_random(m4, n*sizeof(mat64)), );
        TIME1N_spins(, 2, PLUQ64_n, (phi,l,u4,m4,64*n), ret/64,n+1);
        TIME1N_spins(memfill_random(m4, n*sizeof(mat64)), 2, PLUQ64_n, (phi,l,u4,m4,64*n), ret/64,n+1);
        TIME1(2, LUP64_imm, (l,u,p,m));
        TIME1(2, full_echelon_6464_imm, (mm,e,m));
        TIME1(2, gauss_128128_C, (m));
        level3_gauss_tests_N(64);
        level3_gauss_tests_N(128);
        level3_gauss_tests_N(256);
        level3_gauss_tests_N(512);
        level3_gauss_tests_N(1024);
    }

    if (1) {
        size_t n = 64;
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 *n * sizeof(mat64));
        memfill_random(A, n * sizeof(mat64));
        memfill_random(B, n * sizeof(mat64));
        printf("-- polynomials (N=%zu) --\n", n);
        TIME1(5, m64pol_mul, (C,A,B,n,n));
        TIME1(5, m64pol_mul_kara, (C,A,B,n,n));
        free(A);
        free(B);
        free(C);
    }

    if (1) {
        size_t n = 64;
        unsigned int K = 2;
        mat64 * A = (mat64 *) malloc(K * K * n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(K * K * n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(K * K * 2 *n * sizeof(mat64));
        memfill_random(A, K * K * n * sizeof(mat64));
        memfill_random(B, K * K * n * sizeof(mat64));
        printf("-- polynomials, larger matrices (K=%u, N=%zu) --\n", K, n);
        TIME1(5, m64polblock_mul, (C,A,B,n,n,2));
        TIME1(5, m64polblock_mul_kara, (C,A,B,n,n,K));
        free(A);
        free(B);
        free(C);
    }

    if (1) {
        size_t n = 128;
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(n * sizeof(mat64));
        uint64_t * Al = (uint64_t *) A;
        uint64_t * Bl = (uint64_t *) B;
        uint64_t * Cl = (uint64_t *) C;
        memfill_random(A, n * sizeof(mat64));
        memfill_random(B, n * sizeof(mat64));
        printf("-- 64x64 matrices over GF(2^64) --\n");
        TIME1(5, m64pol_mul_gf2_64_bitslice, (C,A,B));
        TIME1(5, m64pol_mul_gf2_64_nobitslice, (Cl,Al,Bl));
        printf("-- 64x64 matrices over GF(2^128) --\n");
        TIME1(5, m64pol_mul_gf2_128_bitslice, (C,A,B));
        TIME1(5, m64pol_mul_gf2_128_nobitslice, (Cl,Al,Bl));
        free(A);
        free(B);
        free(C);
        /* On Core i5 (magret), it's almost a tie between the two
         * options... */
#if 0
-- 64x64 matrices over GF(2^64) --
m64pol_mul_gf2_64_bitslice       6351 times in 0.7889 ms each
m64pol_mul_gf2_64_nobitslice    4773 times in 1.0497 ms each
-- 64x64 matrices over GF(2^128) --
m64pol_mul_gf2_128_bitslice      2067 times in 2.4238 ms each
m64pol_mul_gf2_128_nobitslice   1521 times in 3.2939 ms each
#endif
        /* Without pclmul, of course the situation is more clear (truffe,
         * Core2 Duo U9400 */
#if 0
-- 64x64 matrices over GF(2^64) --
m64pol_mul_gf2_64_bitslice       2695 times in 1.8590 ms each
m64pol_mul_gf2_64_nobitslice    435 times in 11.5172 ms each
-- 64x64 matrices over GF(2^128) --
m64pol_mul_gf2_128_bitslice      871 times in 5.7520 ms each
m64pol_mul_gf2_128_nobitslice   158 times in 31.8354 ms each

#endif
    }


    if (1) {
        /* Now multiplication by a scalar. We'll do both GF(2^64) and
         * GF(2^128), so let's allocate room for both */
        size_t n = 128;
        /* random values with average hamming weight. */
        uint64_t scalar[2] = { UINT64_C(0x8d5511cbd7f0d885), UINT64_C(0x2073a477a8b5dd8a) };
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        uint64_t * Al = (uint64_t *) A;
        uint64_t * Bl = (uint64_t *) B;
        memfill_random(A, n * sizeof(mat64));
        memfill_random(B, n * sizeof(mat64));
        printf("-- 64x64 matrix over GF(2^64), multiplication by scalar --\n");
        TIME1(5, m64pol_scalmul_gf2_64_bitslice, (B,A,scalar));
        TIME1(5, m64pol_scalmul_gf2_64_bitslice2, (B,A,scalar));
        TIME1(5, m64pol_scalmul_gf2_64_nobitslice, (Bl,Al,scalar));
        printf("-- 64x64 matrix over GF(2^128), multiplication by scalar --\n");
        TIME1(5, m64pol_scalmul_gf2_128_bitslice, (B,A,scalar));
        TIME1(5, m64pol_scalmul_gf2_128_nobitslice, (Bl,Al,scalar));
        free(A);
        free(B);
        /* The bitsliced version sucks. Really.
         * TODO: See if we can do something. Abandon L1 cache focus, and
         * be content with L2 ? */
    }
#ifdef HAVE_M4RIE
    if (1) {
        printf("-- 64x64 matrices over GF(2^64) using M4RIE --\n");
        /* Now try to see if m4rie can improve these timings */
        /* Unfortunately as of version 20111203, m4rie supporst only
         * GF(2^n) up until n==10. Which cleary won't do, for our
         * objectives. So the following code aborts with segmentation
         * fault. */
        GFqDom<int> GF = GFqDom<int>(2,64);
        FiniteField *F = (FiniteField*)&GF;
        gf2e * ff = gf2e_init_givgfq(F);
        mzed_t *Az = mzed_init(ff, 64, 64);
        mzed_t *Bz = mzed_init(ff, 64, 64);
        mzed_t *Cz = mzed_init(ff, 64, 64);
        mzed_randomize(Az);
        mzed_randomize(Bz);
        TIME1(5, mzed_mul, (Cz, Az, Bz));
        mzed_free(Az);
        mzed_free(Bz);
        mzed_free(Cz);
        gf2e_free(ff);
    }
#endif /* HAVE_M4RIE */

    return 0;
}
