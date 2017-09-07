#ifndef MATOPS_H_
#define MATOPS_H_

#include <stdint.h>
#include "macros.h"

typedef uint64_t mat64[64] ATTRIBUTE((aligned(64)));
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

struct pmat_s {
    int * v;
    int n;
};
typedef struct pmat_s pmat[1];
typedef struct pmat_s * pmat_ptr;
typedef const struct pmat_s * pmat_srcptr;

typedef uint64_t (*m64pol_ptr)[64];
typedef uint64_t (*const m64pol_srcptr)[64];

#ifdef __cplusplus
extern "C" {
#endif


extern void pmat_init(pmat_ptr x, int n);
extern void pmat_clear(pmat_ptr x);
extern void pmat_transpose(pmat_ptr x, pmat_srcptr y);
extern void pmat_get_matrix(mat64 * qm, pmat_ptr qp);
extern void pmattab_complete(int * phi, uint64_t * bits, int nbits);
extern void pqperms_from_phi(pmat_ptr p, pmat_ptr q, int * phi, int m, int n);
static inline int pmat_get(pmat_srcptr x, int k) { return x->v[k]; }
static inline void pmat_set(pmat_srcptr x, int k, int w) { x->v[k]=w; }

extern void mat64_set_identity(mat64_ptr m);
extern int mat64_is_uppertriangular(mat64_srcptr u);
extern int mat64_is_lowertriangular(mat64_srcptr u);
extern int mat64_triangular_is_unit(mat64_srcptr u);
extern int mat64_eq(mat64_srcptr a, mat64_srcptr b);
extern void mat64_copy(mat64_ptr a, mat64_srcptr b);

/* These are the only important entry points. These functions
 * are meant to be accessors for tuned implementations (ok, for now
 * tuning is static).
 */
extern void transp_6464(mat64_ptr dst, mat64_srcptr src);
extern void copy_6464(mat64_ptr dst, mat64_srcptr src);
extern void mul_6464_6464(mat64 C, mat64 A, mat64 B);
extern void add_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
extern void mul_o64_6464(uint64_t *r, uint64_t a, mat64_srcptr w);
extern void mul_N64_6464(uint64_t *C, const uint64_t *A,
 	 const uint64_t *B, size_t m);
extern void addmul_N64_6464(uint64_t *C, const uint64_t *A,
 	 const uint64_t *B, size_t m);
extern void mul_N64_T6464(uint64_t *C, const uint64_t *A,
            const uint64_t *B, size_t m);
extern void addmul_To64_o64(uint64_t * r, uint64_t a, uint64_t w);
extern void mul_o64_6464(uint64_t * r, uint64_t a, mat64_srcptr w);
extern void mul_o64_T6464(uint64_t * w, uint64_t a, mat64_srcptr b);
extern void addmul_TN64_N64(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol);
extern void mul_TN64_N64(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol);

/* Here come the various implementations of each */
extern void add_6464_6464_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);

extern void mul_N64_6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
extern void mul_N64_6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
#if defined(HAVE_SSE2) && ULONG_BITS == 64
extern void mul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, size_t m);
#endif
extern void mul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
extern void mul_N64_6464_lookup8(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
extern void mul_N64_T6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
extern void mul_N64_T6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m);
extern void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64_srcptr w);
extern void mul_o64_6464_C_msb(uint64_t *r, uint64_t a, mat64_srcptr w);
extern void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64_srcptr b);
extern void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64_srcptr b);
#if defined(HAVE_SSE2) && ULONG_BITS == 64
extern void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
#endif
extern void mul_6464_6464_v2(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);


extern void addmul_To64_o64_lsb(uint64_t * r, uint64_t a, uint64_t w);
extern void addmul_To64_o64_msb(uint64_t * r, uint64_t a, uint64_t w);
extern void addmul_To64_o64_lsb_packof2(uint64_t * r, uint64_t a, uint64_t w);
extern void addmul_To64_o64_lsb_sse_v1(uint64_t * r, uint64_t a, uint64_t w);
extern void mul_64N_N64_addmul(uint64_t *r, uint64_t *a, uint64_t *w, size_t n);
extern void mul_TN32_N64_C(uint64_t * b, uint32_t * A, uint64_t * x, unsigned int ncol);
extern void addmul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol);
extern void mul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol);

extern int PLUQ64(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m);
extern int PLUQ128(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m);
extern int PLUQ64_n(int * phi, mat64 l, mat64 * u, mat64 * a, int n);
extern int LUP64_imm(mat64 l, mat64 u, mat64 p, mat64 a);
extern int full_echelon_6464_imm(mat64 mm, mat64 e, mat64 m);
extern int gauss_128128_C(uint64_t * m);

extern void m64pol_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);
extern void m64pol_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2);
extern void m64polblock_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K);
extern void m64polblock_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K);

extern void m64pol_mul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2);
extern void m64pol_mul_gf2_64_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2);
extern void m64pol_mul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2);
extern void m64pol_mul_gf2_128_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2);
extern void m64pol_scalmul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);
extern void m64pol_scalmul_gf2_64_bitslice2(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);
extern void m64pol_scalmul_gf2_64_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar);
extern void m64pol_scalmul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s);
extern void m64pol_scalmul_gf2_128_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar);

#ifdef __cplusplus
}
#endif

#endif	/* MATOPS_H_ */
