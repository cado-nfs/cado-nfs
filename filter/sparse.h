#ifndef CADO_SPARSE_H_
#define CADO_SPARSE_H_

#include <stdio.h>
#include <stdint.h> /* for int32_t */
#include "typedefs.h"   // index_t
#include "merge_replay_matrix.h"        // typerow_t

#ifdef FOR_DL
#define rowCell(row, k) (row)[k].id
#define rowFullCell(row, k) (row)[k]
#define setRawCell(row, k, v, c) (row)[k] = (ideal_merge_t) {.id = (v), .e = (c)}
#define setCell(row, k, v, c) (row)[k] = (ideal_merge_t) {.id = (v), .e = (c)}
#define compressRow(row, buf, n) memcpy(row,buf,(n+1)*sizeof (typerow_t))
#else
#define setRawCell(row, k, v, c) (row)[k] = (v)
#define rowFullCell(row, k) rowCell((row), (k))
#if SIZEOF_INDEX == 4 || SIZEOF_INDEX == 8
#define rowCell(row, k) (row)[k]
#define setCell(row, k, v, c) (row)[k] = (v)
#define compressRow(row, buf, n) memcpy((row),buf,((n)+1)*sizeof (typerow_t))
#else /* experimental code for 5 <= SIZEOF_INDEX <= 7, for factorization */
#ifdef __cplusplus
extern "C" {
#endif
index_t rowCell (index_t *row, int k);
void setCell (index_t *row, int k, index_t j, exponent_t e);
void compressRow (index_t *row, index_t *buf, int n);
#ifdef __cplusplus
}
#endif
#endif
#endif

// Structures for the data to create the index file.
// This is an array of relation-sets.
// One relation-set is an array of pairs (relation,multiplicity).
//
// Note: the reason why it's here is because the addRows function is
// shared by merge and replay, and needs to know about it.

typedef struct {
    /* we assume the number of rows at the end of purge is < 2^32 */
    uint32_t ind_row;
#ifdef FOR_DL
    int32_t e;
#endif
} multirel_t;

typedef struct {
    unsigned int n;
    multirel_t * rels;
} relset_t;

typedef relset_t * index_data_t;


// Protos.

#ifdef __cplusplus
extern "C" {
#endif
extern void fprintRow(FILE *file, typerow_t *row);

extern int parse_hisfile_line (index_signed_t *ind, const char *t, index_t *j);

extern void addRowsUpdateIndex(typerow_t **rows, index_data_t index_data_t,
        index_t i1, index_t i2, index_t j);

static inline void addRows(typerow_t **rows, index_t i1, index_t i2, index_t j);

typerow_t* mallocRow (uint32_t nb);
typerow_t* reallocRow (typerow_t* row, uint32_t nb);
#ifdef __cplusplus
}
#endif

// The following version without index updating, is for merge.
static inline void addRows(typerow_t **rows, index_t i1, index_t i2, index_t j) {
    addRowsUpdateIndex(rows, NULL, i1, i2, j);
}

#endif  /* CADO_SPARSE_H_ */
