#ifndef LINGEN_MATPOLY_FT_H_
#define LINGEN_MATPOLY_FT_H_

#include "lingen-matpoly.hpp"
#include "flint-fft/fft.h"
#include "lingen-substep-schedule.h"

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

void matpoly_ft_init(abdst_field ab, matpoly_ft_ptr t, unsigned int m, unsigned int n, const struct fft_transform_info * fti);
void matpoly_ft_zero(abdst_field ab, matpoly_ft_ptr t, const struct fft_transform_info * fti);
void matpoly_ft_export(abdst_field ab, matpoly_ft_ptr t, const struct fft_transform_info * fti);
void matpoly_ft_import(abdst_field ab, matpoly_ft_ptr t, const struct fft_transform_info * fti);
void matpoly_ft_clear(abdst_field ab, matpoly_ft_ptr t, const struct fft_transform_info * fti);
void matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr p, const struct fft_transform_info * fti);
void matpoly_ft_add(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
*/
void matpoly_ft_mul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
void matpoly_ft_addmul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
void matpoly_ft_ift(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, const struct fft_transform_info * fti);
void matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, unsigned int shift, const struct fft_transform_info * fti);


/* In a way, this is the only real API exported by this module */
void matpoly_mul_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void matpoly_mul_caching(abdst_field ab, matpoly c, matpoly a, matpoly b, const struct lingen_substep_schedule * S) { return matpoly_mul_caching_adj(ab, c, a, b, UINT_MAX, S); }

void matpoly_mp_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void matpoly_mp_caching(abdst_field ab, matpoly c, matpoly a, matpoly b, const struct lingen_substep_schedule * S) { return matpoly_mp_caching_adj(ab, c, a, b, UINT_MAX, S); }


#ifdef __cplusplus
}
#endif


#endif	/* LINGEN_MATPOLY_FT_H_ */
