#ifndef BLOCKMATRIX_H_
#define BLOCKMATRIX_H_

#include "bblas.hpp"
#include "macros.h"

struct blockmatrix {
    mat64 * mb;
    unsigned int nrblocks;   // row blocks
    unsigned int ncblocks;   // col blocks
    unsigned int nrows;      // only for convenience
    unsigned int ncols;      // only for convenience
    unsigned int stride;     // may differ from nrblocks if owner==false
    bool owner;
    blockmatrix(unsigned int nrows, unsigned int ncols);
    ~blockmatrix();
    blockmatrix(blockmatrix const &) = delete;
    blockmatrix(blockmatrix & k, int i0, int j0, unsigned int nrows, unsigned int ncols);
    blockmatrix(blockmatrix && x) {
        mb = x.mb; x.mb = NULL;
        nrblocks = x.nrblocks;
        ncblocks = x.ncblocks;
        nrows = x.nrows;
        ncols = x.ncols;
        stride = x.stride;
        owner = x.owner; x.owner = false;
    }
    blockmatrix& operator=(blockmatrix && x) {
        std::swap(mb, x.mb);
        std::swap(nrblocks, x.nrblocks);
        std::swap(ncblocks, x.ncblocks);
        std::swap(nrows, x.nrows);
        std::swap(ncols, x.ncols);
        std::swap(stride, x.stride);
        std::swap(owner, x.owner);
        return *this;
    }

    void set_zero();
    void set_identity();
    uint64_t * subrow_ptr(int i, int j);
    void copy_colrange(blockmatrix const & A, int j0, int j1);

    void mul_Ta_b(blockmatrix const & a, blockmatrix const & b);
    void mul_smallb(
            blockmatrix const & a,
            blockmatrix const & b);
    void copy_to_flat(uint64_t * tiny, unsigned int stride,
            int i0, int j0) const;
    void copy_transpose_to_flat(uint64_t * tiny, unsigned int stride,
            int i0, int j0) const;
    void copy_transpose_from_flat(uint64_t * tiny, unsigned int stride, int i0, int j0);
    void copy_from_flat(uint64_t * tiny, unsigned int stride, int i0, int j0);
    void read_from_flat_file(int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols);
    void write_to_flat_file(const char * name, int i0, int j0, unsigned int fnrows, unsigned int fncols) const;
    void transpose(blockmatrix const & a);
    void reverse_columns();
    void reverse_rows();
    void read_transpose_from_flat_file(int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols);
    void swap(blockmatrix & A);
    static void swap_words_if_needed (uint64_t *v, unsigned long n);
    void print(const char *vname) const;
};


/* Use this macro to allocate flat matrix areas with proper readahead
 * padding. In some situations, it is also necessary to use zero out the
 * padding data as well, so this macro must also be used in the memset()
 * calls.
 */
#define FLAT_BYTES_WITH_READAHEAD(nr, nc) \
    iceildiv((nr), 64) * iceildiv((nc), 64) * sizeof(mat64)

#endif	/* BLOCKMATRIX_H_ */
