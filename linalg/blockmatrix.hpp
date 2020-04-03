#ifndef BLOCKMATRIX_HPP_
#define BLOCKMATRIX_HPP_

#include "macros.h"
#include "bblas.hpp"
#include "submatrix_range.hpp"

// template<typename matrix>
struct blockmatrix {
    typedef mat64 matrix;
    private:
    mat64 * mb;
    unsigned int nrblocks;   // row blocks
    unsigned int ncblocks;   // col blocks
    unsigned int _nrows;      // only for convenience
    unsigned int _ncols;      // only for convenience
    unsigned int stride;     // may differ from nrblocks if owner==false
    bool owner;
    public:
    inline unsigned int nrows() const { return _nrows; }
    inline unsigned int ncols() const { return _ncols; }
    inline unsigned int nrows_padded() const { return nrblocks * matrix::width; }
    inline unsigned int ncols_padded() const { return ncblocks * matrix::width; }
    blockmatrix(unsigned int nrows, unsigned int ncols);
    ~blockmatrix();
    blockmatrix(blockmatrix const &) = delete;
    blockmatrix(blockmatrix & k, int i0, int j0, unsigned int nrows, unsigned int ncols);
    blockmatrix(blockmatrix && x) {
        mb = x.mb; x.mb = NULL;
        nrblocks = x.nrblocks;
        ncblocks = x.ncblocks;
        _nrows = x._nrows;
        _ncols = x._ncols;
        stride = x.stride;
        owner = x.owner; x.owner = false;
    }
    blockmatrix& operator=(blockmatrix && x) {
        std::swap(mb, x.mb);
        std::swap(nrblocks, x.nrblocks);
        std::swap(ncblocks, x.ncblocks);
        std::swap(_nrows, x._nrows);
        std::swap(_ncols, x._ncols);
        std::swap(stride, x.stride);
        std::swap(owner, x.owner);
        return *this;
    }

    struct view_t : public submatrix_range {/*{{{*/
        blockmatrix & M;
        view_t(blockmatrix & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(blockmatrix & M) : submatrix_range(M), M(M) {}
        inline matrix::datatype & part(unsigned int ii, unsigned int jj) {
            ii += i0;
            jj += j0;
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
        inline matrix::datatype const & part(unsigned int ii, unsigned int jj) const {
            ii += i0;
            jj += j0;
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
        void zero();
    };
/*}}}*/
    struct const_view_t : public submatrix_range {/*{{{*/
        blockmatrix const & M;
        const_view_t(blockmatrix const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(blockmatrix const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline matrix::datatype const & part(unsigned int ii, unsigned int jj) const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
    };
/*}}}*/
    struct row_view_t {/*{{{*/
        blockmatrix & M;
        unsigned int ii;
        row_view_t(blockmatrix & M, unsigned int ii) : M(M), ii(ii) {}
        inline matrix::datatype & operator[](unsigned int jj) {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
        inline matrix::datatype const & operator[](unsigned int jj)  const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
    };
/*}}}*/
    struct const_row_view_t {/*{{{*/
        blockmatrix const & M;
        unsigned int ii;
        const_row_view_t(blockmatrix const & M, unsigned int ii) : M(M), ii(ii) {}
        inline matrix::datatype const & operator[](unsigned int jj)  const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.mb[(ii / B) + (jj / B) * M.stride][ii % B];
        }
    };
/*}}}*/
    row_view_t operator[](unsigned int ii) { return row_view_t(*this, ii); }
    const_row_view_t operator[](unsigned int ii) const { return const_row_view_t(*this, ii); }

    void set_zero();
    void set_identity();

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

#endif	/* BLOCKMATRIX_HPP_ */
