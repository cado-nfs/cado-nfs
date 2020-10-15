#ifndef BLOCKMATRIX_HPP_
#define BLOCKMATRIX_HPP_

#include <stddef.h>             // for NULL
#include <stdint.h>             // for uint64_t
#include <utility>              // for swap
#include "bblas.hpp"      // for mat64
#include "macros.h"
#include "submatrix_range.hpp"

/* This interface is old. We most probably want to use bpack<mat64>
 * instead, as provided in bblas/bpack.hpp ; however we have more things
 * here. We might find it interesting to augment the (primitive) bpack
 * interface to include the features from here.
 */
// template<typename matrix>
struct blockmatrix {
    typedef mat64 matrix;
    private:
    mat64 * mb;
    unsigned int nrblocks;   // row blocks
    unsigned int ncblocks;   // col blocks
    unsigned int _nrows;      // only for convenience
    unsigned int _ncols;      // only for convenience
    public:
    inline unsigned int rowstride() const { return 1; }
    inline unsigned int colstride() const { return nrblocks; }
    inline mat64 & getblock(unsigned int bi, unsigned int bj) {
        return mb[bi * rowstride() + bj * colstride()];
    }
    inline mat64 const & getblock(unsigned int bi, unsigned int bj) const {
        return mb[bi * rowstride() + bj * colstride()];
    }

    inline unsigned int nrows() const { return _nrows; }
    inline unsigned int ncols() const { return _ncols; }
    inline unsigned int nrows_padded() const { return nrblocks * matrix::width; }
    inline unsigned int ncols_padded() const { return ncblocks * matrix::width; }
    blockmatrix(unsigned int nrows, unsigned int ncols);
    ~blockmatrix();
    blockmatrix(blockmatrix const &) = delete;
    blockmatrix(blockmatrix && x) {
        mb = x.mb; x.mb = NULL;
        nrblocks = x.nrblocks;
        ncblocks = x.ncblocks;
        _nrows = x._nrows;
        _ncols = x._ncols;
    }
    blockmatrix& operator=(blockmatrix && x) {
        std::swap(mb, x.mb);
        std::swap(nrblocks, x.nrblocks);
        std::swap(ncblocks, x.ncblocks);
        std::swap(_nrows, x._nrows);
        std::swap(_ncols, x._ncols);
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
            return M.getblock(ii/B, jj/B)[ii % B];
        }
        inline matrix::datatype part(unsigned int ii, unsigned int jj) const {
            ii += i0;
            jj += j0;
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.getblock(ii/B, jj/B)[ii % B];
        }
        inline unsigned int nrblocks() const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(i0 % B == 0);
            return iceildiv(nrows(), B);
        }
        inline unsigned int ncblocks() const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(j0 % B == 0);
            return iceildiv(ncols(), B);
        }
        inline mat64 & getblock(unsigned int bi, unsigned int bj) {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(i0 % B == 0);
            ASSERT_ALWAYS(j0 % B == 0);
            return M.getblock(i0 / B + bi, j0 / B + bj);
        }
        inline mat64 const & getblock(unsigned int bi, unsigned int bj) const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(i0 % B == 0);
            ASSERT_ALWAYS(j0 % B == 0);
            return M.getblock(i0 / B + bi, j0 / B + bj);
        }
        void set_zero();
    };
/*}}}*/
    struct const_view_t : public submatrix_range {/*{{{*/
        blockmatrix const & M;
        const_view_t(blockmatrix const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(blockmatrix const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline matrix::datatype part(unsigned int ii, unsigned int jj) const {
            ii += i0;
            jj += j0;
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.getblock(ii/B, jj/B)[ii % B];
        }
        inline unsigned int nrblocks() const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(i0 % B == 0);
            return iceildiv(nrows(), B);
        }
        inline unsigned int ncblocks() const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(j0 % B == 0);
            return iceildiv(ncols(), B);
        }
        inline mat64 const & getblock(unsigned int bi, unsigned int bj) const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(i0 % B == 0);
            ASSERT_ALWAYS(j0 % B == 0);
            return M.getblock(i0 / B + bi, j0 / B + bj);
        }
    };
/*}}}*/
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }
    // operator view_t() { return view(); }
    // operator const_view_t() const { return view(); }

    struct row_view_t {/*{{{*/
        blockmatrix & M;
        unsigned int ii;
        row_view_t(blockmatrix & M, unsigned int ii) : M(M), ii(ii) {}
        inline matrix::datatype & operator[](unsigned int jj) {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.getblock(ii/B, jj/B)[ii % B];
        }
        inline matrix::datatype operator[](unsigned int jj)  const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.getblock(ii/B, jj/B)[ii % B];
        }
    };
/*}}}*/
    struct const_row_view_t {/*{{{*/
        blockmatrix const & M;
        unsigned int ii;
        const_row_view_t(blockmatrix const & M, unsigned int ii) : M(M), ii(ii) {}
        inline matrix::datatype operator[](unsigned int jj)  const {
            constexpr const unsigned int B = matrix::width;
            ASSERT_ALWAYS(jj % B == 0);
            return M.getblock(ii/B, jj/B)[ii % B];
        }
    };
/*}}}*/
    row_view_t operator[](unsigned int ii) { return row_view_t(*this, ii); }
    const_row_view_t operator[](unsigned int ii) const { return const_row_view_t(*this, ii); }

    void set_zero();
    void set_identity();

    void copy_colrange(blockmatrix const & A, unsigned int j0, unsigned int j1);

    static void mul_Ta_b(view_t c, const_view_t a, const_view_t b);

    // replaces a=*this by a*b, in place.
    static void mul_smallb(view_t a, const_view_t b);
    static void mul_smallb(blockmatrix & c, const_view_t a, const_view_t b);

    struct flat_area {/*{{{*/
        /* a flat area has the matrix stored row-major starting from
         * pointer p, but with a possibly non-trivial stride from one row
         * to the next.
         */
        matrix::datatype * p;
        unsigned int row_stride;
        flat_area(matrix::datatype * p, unsigned int row_stride) : p(p), row_stride(row_stride) {}
    };/*}}}*/
    struct const_flat_area {/*{{{*/
        const matrix::datatype * p;
        unsigned int row_stride;
        const_flat_area(const matrix::datatype * p, unsigned int row_stride) : p(p), row_stride(row_stride) {}
        const_flat_area(flat_area F) : p(F.p), row_stride(F.row_stride) {}
    };/*}}}*/

    static void copy_to_flat(flat_area F, const_view_t V);
    static void copy_transpose_to_flat(flat_area F, const_view_t V);
    static void copy_transpose_from_flat(view_t V, const_flat_area F);
    static void copy_from_flat(view_t V, const_flat_area F);

    static void copy_to_flat(matrix::datatype * flat, unsigned int flat_row_stride, const_view_t V) {/*{{{*/
        copy_to_flat(flat_area(flat, flat_row_stride), V);
    }/*}}}*/
    static void copy_transpose_to_flat(matrix::datatype * flat, unsigned int flat_row_stride, const_view_t V) {/*{{{*/
        copy_transpose_to_flat(flat_area(flat, flat_row_stride), V);
    }/*}}}*/
    static void copy_transpose_from_flat(view_t V, const matrix::datatype * flat, unsigned int flat_row_stride) {/*{{{*/
        copy_transpose_from_flat(V, const_flat_area(flat, flat_row_stride));
    }/*}}}*/
    static void copy_from_flat(view_t V, const matrix::datatype * flat, unsigned int flat_row_stride) {/*{{{*/
        copy_from_flat(V, const_flat_area(flat, flat_row_stride));
    }/*}}}*/

    /* This reads/writes the blockmatrix from a data file that is stored
     * in flat format, with a row stride of exactly ceil(ncols/64).
     */
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
