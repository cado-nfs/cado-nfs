#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "purgedfile.h"
#include "typedefs.h"

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

#ifndef FOR_DL
#define typerow_t index_t
#define cmp_typerow_t cmp_index
#else
#define typerow_t ideal_merge_t
#define cmp_typerow_t cmp_ideal_merge
#endif

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int verbose;         /* verbose level */
  uint64_t nrows;
  uint64_t ncols;
  uint64_t rem_nrows;  /* number of remaining rows */
  uint64_t rem_ncols;  /* number of remaining columns, including the buried
                          columns */
  typerow_t **rows;    /* rows[i][k] contains indices of an ideal of row[i] 
                          with 1 <= k <= rows[i][0] */
                       /* FOR_DL: struct containing also the exponent */
  unsigned char *wt;   /* weight w of column j, if w <= cwmax,
                          else 0 for a deleted column */
                       /* 8 bits is sufficient as we only want precise weight
                          for column of low weight. If the weight exceeds
                          255, we saturate at 255 */
  uint64_t skip;       /* number of buried/skipped columns of smaller index */
  uint64_t weight;     /* number of non-zero coefficients in the active part */
  uint64_t tot_weight; /* Initial total number of non-zero coefficients */
  int cwmax;           /* bound on weight of j to enter the SWAR structure */
  index_t **R;         /* R[j][k] contains the indices of the rows containing
                          the ideal of index j, 0 <= j < ncols,
                          1 <= k <= R[j][0], for weight(j) <= cwmax.
                          R[j][k] = UMAX(index_t) if the corresponding row has
                                                                  been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
  index_t *MKZQ;         /* priority queue for Markowitz stuff */
  index_t *MKZA;         /* MKZA[j] gives u s.t. MKZQ[2*u] = j and
                            MKZQ[2*u+1] is the Markowitz cost of column j,
                            otherwise it is MKZ_INF if the column is inactive
                            (either too heavy initially or deleted). It is safe
                            to have a 32-bit table as long as the Markowitz
                            queue has less than 2^32 entries. */
} filter_matrix_t;

#ifdef __cplusplus
extern "C" {
#endif

#define compute_WN(mat) ((double) (mat)->rem_nrows * (double) (mat)->weight)
#define compute_WoverN(mat) (((double)(mat)->weight)/((double)(mat)->rem_nrows))

void initMat(filter_matrix_t *, uint32_t);
void clearMat (filter_matrix_t *mat);

void print_row(filter_matrix_t *mat, index_t i);

#define matCell(mat, i, k) rowCell(mat->rows[i], k)
#define rowLength(row, i) rowCell(row[i], 0)
#define matLengthRow(mat, i) matCell(mat, i, 0)
#define isRowNull(mat, i) ((mat)->rows[(i)] == NULL)

#define SPARSE_ITERATE(mat, i, k) for((k)=1; (k)<=lengthRow((mat),(i)); (k)++)

int weightSum(filter_matrix_t *mat, index_t i1, index_t i2, index_t j);

#ifdef __cplusplus
}
#endif

#endif	/* SPARSE_MAT_H_ */
