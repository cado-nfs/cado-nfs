#ifndef SQRT_ABFILES_HPP_
#define SQRT_ABFILES_HPP_

#include <cstddef>
#include <cstdio>
#include <cstdint>

#include "select_mpi.h"

struct ab_source_s {
    const char * fname0;
    char * sname;
    size_t sname_len;
    size_t nab;
    char * prefix;
    int depnum;
    // int cado; // use !numfiles instead.

    /* this relates to rough estimations based on the size(s) of the
     * file(s). Note however that apparently the guess logic isn't too
     * good at taking the leading coefficient into account, so it's not
     * meant to be used before it gets fixed. Doing the accurate
     * evaluation via complex embeddings is ultra-fast anyway */
    size_t nab_estim;
    // size_t digitbytes_estim;

    /* This relates to the different files (if several), their sizes, and
     * the current position. */
    size_t * file_bases;
    int nfiles; // 0 for cado format.
    size_t totalsize;

    FILE * f;
    int c;
    size_t cpos;        // position within current file.
    size_t tpos;        // position within totality.
};

typedef struct ab_source_s ab_source[1];
typedef struct ab_source_s * ab_source_ptr;

extern void ab_source_init(ab_source_ptr ab, const char * fname, int rank, int root, MPI_Comm comm);
extern void ab_source_rewind(ab_source_ptr ab);
extern void ab_source_init_set(ab_source_ptr ab, ab_source_ptr ab0);
extern void ab_source_clear(ab_source_ptr ab);
extern int ab_source_next(ab_source_ptr ab, int64_t * a, uint64_t * b);
extern void ab_source_move_afterpos(ab_source_ptr ab, size_t offset);





#endif	/* SQRT_ABFILES_HPP_ */
