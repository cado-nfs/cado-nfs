#ifndef CADO_LINALG_GAUSS_H_
#define CADO_LINALG_GAUSS_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute the right nullspace of the nrows times ncols matrix given by
 * mat. The matrix is given as a flat array of mp_limb_t, with the
 * specified number of limbs per row. The kernel is written to the array
 * of arrays given by the ker argument, where each member is expected to
 * hold space for at least limbs_per_col mp_limb_t values. Caution leads
 * to allocate as many as ncols pointers in the ker array.
 *
 * The dimension of the kernel is given by the return value. If ker ==
 * NULL, this is the only thing computed (and limbs_per_col is unused).
 *
 * limbs_per_row (and accordingly limbs_per_col) must of course being
 * larger than or equal to ceiling(nrows/GMP_LIMB_BITS). We allow this
 * value to be exceeded so as to allow some padding.
 *
 * In case you wonder, this function is not reentrant at all. Sorry.
 */
extern int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
		  int limbs_per_row, int limbs_per_col);
/* This is the dual function. It returns into lmat an extraction matrix
 * such that the first rows of lmat*mat are linearly independent, while
 * the rest is zero. lmat is full rank. Of course the number of
 * independent rows matches the rank of mat, and is actually returned by
 * the function.
 */
extern int spanned_basis(mp_limb_t * lmat, mp_limb_t * mat, int nrows, int ncols,
        int limbs_per_row, int limbs_per_col, mp_limb_t * elim_table);
#ifdef __cplusplus
}
#endif

#endif	/* CADO_LINALG_GAUSS_H_ */
