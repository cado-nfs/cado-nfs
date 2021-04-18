#ifndef PERM_MATRIX_HPP_
#define PERM_MATRIX_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include <stdint.h>         // for uint64_t
#include "bblas_mat64.hpp"

/* We probably want this interface to disappear. A vector of int's should
 * do just equally well.
 */
struct perm_matrix_s {
    int * v;
    int n;
};
typedef struct perm_matrix_s perm_matrix[1];
typedef struct perm_matrix_s * perm_matrix_ptr;
typedef const struct perm_matrix_s * perm_matrix_srcptr;

 
extern void perm_matrix_init(perm_matrix_ptr x, int n);
extern void perm_matrix_clear(perm_matrix_ptr x);
extern void perm_matrix_transpose(perm_matrix_ptr x, perm_matrix_srcptr y);
extern void perm_matrix_get_matrix(mat64 * qm, perm_matrix_ptr qp);
extern void perm_matrixtab_complete(int * phi, uint64_t * bits, int nbits);
extern void pqperms_from_phi(perm_matrix_ptr p, perm_matrix_ptr q, int * phi, int m, int n);
static inline int perm_matrix_get(perm_matrix_srcptr x, int k) { return x->v[k]; }
static inline void perm_matrix_set(perm_matrix_srcptr x, int k, int w) { x->v[k]=w; }

#endif	/* PERM_MATRIX_HPP_ */
