#include <limits.h>
#ifndef LINGEN_MATPOLY_FT_H_
#define LINGEN_MATPOLY_FT_H_

#include "lingen-matpoly.h"
#include "flint-fft/fft.h"

struct matpoly_ft_s {
    unsigned int m;
    unsigned int n;
    void * data;
};

typedef struct matpoly_ft_s matpoly_ft[1];
typedef struct matpoly_ft_s * matpoly_ft_ptr;
typedef const struct matpoly_ft_s * matpoly_ft_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void matpoly_ft_init(abdst_field ab, matpoly_ft_ptr t, unsigned int m, unsigned int n, struct fft_transform_info * fti);
void matpoly_ft_zero(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_export(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_import(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_clear(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
double matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr p, struct fft_transform_info * fti, int draft);
void matpoly_ft_add(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
*/
double matpoly_ft_mul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti, int draft);
double matpoly_ft_addmul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti, int draft);
double matpoly_ft_ift(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, struct fft_transform_info * fti, int draft);
double matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti, int draft);


/* In a way, this is the only real API exported by this module */
double matpoly_mul_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, int draft);
static inline double matpoly_mul_caching(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft) { return matpoly_mul_caching_adj(ab, c, a, b, UINT_MAX, draft); }

double matpoly_mp_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, int draft);
static inline double matpoly_mp_caching(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft) { return matpoly_mp_caching_adj(ab, c, a, b, UINT_MAX, draft); }


#ifdef __cplusplus
}
#endif


#endif	/* LINGEN_MATPOLY_FT_H_ */
