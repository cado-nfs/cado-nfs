#include <limits.h>
#ifndef BIGMATPOLY_FT_H_
#define BIGMATPOLY_FT_H_

#include "mpfq_layer.h"
#include "lingen-matpoly-ft.h"
#include "lingen-bigmatpoly.h"
#include "select_mpi.h"
struct fft_transform_info; // IWYU pragma: keep

/* This defines an MPI-shared polynomial matrix type */

struct bigmatpoly_ft_s {
    /* XXX the first four fields must be compatible with bigmatpoly_s */
    unsigned int m1;      /* number of block rows, index i */
    unsigned int n1;      /* number of block cols, index j */
    MPI_Comm com[3];      /* [0]: MPI_COMM_WORLD, reordered.
                             [1]: row-wise       ; size == n1
                             [2]: column-wise    ; size == m1
                           */

    unsigned int m;     /* total number of rows */
    unsigned int n;     /* total number of cols */
    /* The following three are also in cells */
    unsigned int m0;      /* number of rows per block */
    unsigned int n0;      /* number of cols per block */
    matpoly_ft * cells;
};
typedef struct bigmatpoly_ft_s bigmatpoly_ft[1];
typedef struct bigmatpoly_ft_s * bigmatpoly_ft_ptr;
typedef const struct bigmatpoly_ft_s * bigmatpoly_ft_srcptr;


#ifdef __cplusplus
extern "C" {
#endif

    /* This is only a query-replace from the bigmatpoly interface. Not
     * all methods are relevant, I suppose */

void bigmatpoly_ft_init(abdst_field ab, bigmatpoly_ft_ptr p, bigmatpoly_ft_srcptr model, unsigned int m, unsigned int n, struct fft_transform_info * fti);
void bigmatpoly_ft_finish_init(abdst_field ab, bigmatpoly_ft_ptr p, unsigned int m, unsigned int n, struct fft_transform_info * fti);
void bigmatpoly_ft_init_model(bigmatpoly_ft_ptr model, MPI_Comm * comm, unsigned int m, unsigned int n);
int bigmatpoly_ft_check_pre_init(bigmatpoly_ft_srcptr p);
matpoly_ft_ptr bigmatpoly_ft_my_cell(bigmatpoly_ft_ptr p);
// void bigmatpoly_ft_realloc(bigmatpoly_ft_ptr p, int newalloc, struct
// fft_transform_info * fti);
void bigmatpoly_ft_zero(abdst_field ab, bigmatpoly_ft_ptr p, struct fft_transform_info * fti);
void bigmatpoly_ft_clear(abdst_field ab, bigmatpoly_ft_ptr p, struct fft_transform_info * fti);
void bigmatpoly_ft_clear_model(bigmatpoly_ft_ptr p);

void bigmatpoly_ft_swap(bigmatpoly_ft_ptr a, bigmatpoly_ft_ptr b);

double bigmatpoly_ft_mul(abdst_field ab, bigmatpoly_ft_ptr c, bigmatpoly_ft_ptr a, bigmatpoly_ft_ptr b, struct fft_transform_info * fti, int draft);
double bigmatpoly_ft_mul2(abdst_field ab, bigmatpoly_ft_ptr c, bigmatpoly_ft_ptr a, bigmatpoly_ft_ptr b, struct fft_transform_info * fti, int draft);

double bigmatpoly_ft_dft(abdst_field ab, bigmatpoly_ft_ptr ta, bigmatpoly_ptr a, struct fft_transform_info * fti, int draft);
double bigmatpoly_ft_ift(abdst_field ab, bigmatpoly_ptr a, bigmatpoly_ft_ptr ta, struct fft_transform_info * fti, int draft);
double bigmatpoly_ft_ift_mp(abdst_field ab, bigmatpoly_ptr a, bigmatpoly_ft_ptr ta, unsigned int shift, struct fft_transform_info * fti, int draft);


/* {{{ access interface for bigmatpoly_ft */
static inline matpoly_ft * bigmatpoly_ft_part(bigmatpoly_ft_ptr p, unsigned int i, unsigned int j) {
    return p->cells+i*p->n1+j;
}
static inline matpoly_ft_ptr bigmatpoly_ft_cell(bigmatpoly_ft_ptr p, unsigned int i, unsigned int j) {
    return *bigmatpoly_ft_part(p,i,j);
}
/* }}} */

double bigmatpoly_mul_caching_adj(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b, unsigned int adj, int draft);

double bigmatpoly_mp_caching_adj(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b, unsigned int adj, int draft);

static inline double bigmatpoly_mul_caching(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b, int draft) { return bigmatpoly_mul_caching_adj(ab, c, a, b, UINT_MAX, draft); }

static inline double bigmatpoly_mp_caching(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b, int draft) { return bigmatpoly_mp_caching_adj(ab, c, a, b, UINT_MAX, draft); }

#ifdef __cplusplus
}
#endif

#endif	/* BIGMATPOLY_FT_H_ */
