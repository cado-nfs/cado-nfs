#ifndef MATRIX_FILE_HPP_
#define MATRIX_FILE_HPP_

#include <stdint.h>
#include <stddef.h>
#include <string>
#include <utility>
#include <array>
#include "parallelizing_info.hpp"
#include "params.h"

/* 202207. The cado-nfs matrix format changes.
 *
 * A matrix is actually a collection of files, all sharing the same base
 * name (say $X below).
 * - $X.rows contains the column positions of the non-zero coefficients
 *   in each row, thus in row major order. Each coefficient position is
 *   stored as a 32-bit unsigned little-endian integer. Therefore if the
 *   matrix has a total of N non-zero coefficients, $X.rows cwishas exactly
 *   N*4 bytes long.
 * - $X.cols is the same with row indices in column major order
 * - $X.row_offsets is a list of 32-bit unsigned little endian integers.
 *   The $i$-th of these integers is the number of non-zero coefficients
 *   in rows $0$ to $i$, inclusive. The number of integers in the
 *   row_offsets file is exactly the number of rows of the matrix.
 * - $X.col_offsets is the same for columns.
 *
 * - $X.row_coeffs exists if and only if the matrix is defined over a
 *   field different from GF(2). It has exactly the same size as $X.rows,
 *   and contains 32-bit *signed* little-endian integers. The k-th
 *   integer in the file $X.rows encodes that some non-zero coefficient
 *   is at a given (r,c) position within the matrix, and the value of
 *   this coefficient is of course the k-th integer in $X.row_coeffs
 *
 * - $X.col_coeffs is the same for columns.
 *
 * It is never sufficient to have only one of these files, and it is
 * hardly ever useful to have all of them (and anyway not all of them get
 * created by default).
 *
 * The following combinations are valid. Note that the row_offsets and
 * col_offsets must always be present.
 *
 * - $X.rows + $X.row_offsets + $X.col_offsets :
 *   This is the most common data set for binary matrices.
 *
 * - $X.rows + $X.row_offsets + $X.col_offsets + $X.row_coeffs
 *   This extends the above scheme for matrices over GF(p).
 *
 * - $X.cols + $X.row_offsets + $X.col_offsets :
 *   This is an alternative representation. It makes sense to have this
 *   if the matrix is to be read in transposed form.
 *
 * - $X.cols + $X.row_offsets + $X.col_offsets + $X.col_coeffs
 *   This is an alternative representation.
 *
 * Also, we may have:
 *
 * - $X.rows + $X.row_offsets + $X.col_offsets + $X.cols :
 *   This allows speedy reading of the matrix data, either transposed or
 *   not.
 *   (likewise with coefficients added)
 *
 * Conceptually, $X.rows + $X.row_offsets (or the same with columns) is
 * enough to reconstruct a complete matrix file set. Of course, if the
 * matrix is not over GF(2), the coefficients files must be present as
 * well.
 *
 * The tools to convert to/from other matrix formats and to complete an
 * existing matrix file set into another one are not written yet.
 */
struct matrix_file : public std::vector<uint32_t> {
    // static inline const char * rowcol(int d) { const char * x[2] = { "row", "col", }; return x[d != 0]; }
    static constexpr const char * rowcol[2] { "row", "col" };
    private:
    static inline std::string dotrowcols(int d) {
        return std::string(".") + rowcol[d] + "s";
    }
    static inline std::string dotrowcol_coeffs(int d) {
        return std::string(".") + rowcol[d] + "_coeffs";
    }
    static inline std::string dotrowcol_offsets(int d) {
        return std::string(".") + rowcol[d] + "_offsets";
    }
    
    /* These private fields are filled by lookup() */
    bool has_rowcols[2];
    bool has_rowcol_offsets[2];
    bool has_rowcol_coeffs[2];
    std::string lookup_diagnostic;

    std::string mfile;  // This is the base of the matrix name. The
                        // following file names are possibly looked up:
                        //      mfile + dotrowcols(0)
                        //      mfile + dotrowcols(1)
                        //      mfile + dotrowcol_offsets(0)
                        //      mfile + dotrowcol_offsets(1)
                        //      mfile + dotrowcol_coeffs(0)
                        //      mfile + dotrowcol_coeffs(1)

    public:

    const std::string cmdline_filename;
                          // matrix file name. This is the name that ends
                          // in the command line. It's only used as a
                          // public field for external access, and
                          // remains consistent with the file name that
                          // was used to construct the object. A priori,
                          // it can be either of the form .rows or .cols,
                          // it does not matter.  Note that a random
                          // matrix has no file name!
                      
    std::string bfile;    // balancing file name ; empty means auto-detect.
                          // The balancing file name is ignored in a
                          // single read call.

    matrix_file(std::string const & mfile, std::string const & bfile = {})
        : cmdline_filename(mfile)
        , bfile(bfile)
    {
        for(int i = 0 ; i < 2 ; i++) {
            std::string dotXs = dotrowcols(i);
            std::string::size_type d = mfile.rfind(dotXs);
            if (d + dotXs.size() == mfile.size()) {
                this->direction = i;
                this->mfile = mfile.substr(0, d);
                return;
            }
        }
        throw std::runtime_error("Matrix filename must end in .rows or .cols");
    }

    /* The 'direction' flag just tells that the in-memory representation
     * of the sparse matrix below is given by columns and not by rows.
     * This is sometimes the format that is preferred by some MM layers,
     * and cado-nfs decides to present the local matrices to the MM
     * layers precisely in their preferred format.
     *
     * direction=0: data in *p is row-major
     * direction=1: data in *p is col-major
     *
     * This takes some default value when the matrix_file object is
     * created, as well as when the lookup is performed. However this
     * default value is only an indication of what was found on disk.
     *
     * The moment where this gets important is when the read() function
     * gets called.
     */
    bool direction;

    /* Does this matrix have coefficients? This is determined when
     * lookup() is performed.
     */
    int withcoeffs = -1;        /* -1 means auto-detect */

    /* lookup() can be either a collective call or a single call.
     *
     * The collective call requires the balancing file to be present. By
     * default, it is looked up as $X.${N0}x${N1}.bin, where N0 and N1
     * are the dimensions of the thread mesh.
     *
     */
    int lookup();
    int lookup(parallelizing_info_ptr);
    std::string get_lookup_diagnostic() const { return lookup_diagnostic; }

    /* filled by lookup() */
    std::array<uint32_t, 2> nrowcols = {};
    inline uint32_t nrows() const { return nrowcols[0]; }
    inline uint32_t ncols() const { return nrowcols[1]; }
    uint64_t ncoeffs = 0;

    /* The in-memory representation, which is obtained as a result of
     * the read() calls, is not exactly the same as the representation on
     * disk. Maybe it's a terrible idea. For the moment, what we have in
     * memory is what we used to have as an on-disk format, namely: a
     * long string of 32-bit little endian integers, each row being of
     * the form
     *          (row length) ((col index)(coeff value))*.
     *
     * with the following important modifications:
     * - of course, if read() was called with transform=1, then the above
     *   is to be read with columns instead of rows.
     * - in a collective read, each thread gets only its own fragment of
     *   the matrix (and also, this is done according to the balancing
     *   file).
     *
     * The steal() and release() functions are meant for the ultimate
     * consumers of the matrix_file type. They may either decide to steal
     * the memory area, and play tricks with it, or read it and then
     * dospose of it when they're done.
     */
    void steal(std::vector<uint32_t> & x) { std::swap(x, *this); }
    void release() { clear(); }

    /*
     * The two read() calls are actually vastly different procedures, and
     * implemented in different places.
     */

    /* defined in matrix_file.cpp */
    /* This reads the matrix data from mfile into p, arranging it
     * according to what is documented above regarding p.
     * If sanity_check_vector is not the empty string, this is understood
     * as a filename where the result of the product of the matrix times
     * a specific vector, which can be checked by the bwc/dispatch
     * program later on.  Note that this sanity check thing is specific
     * to GF(2) for the moment.
     *
     * This read() call is _NOT_ a collective call *BUT* the read may be
     * done in parallel by several openmp threads.
     */
    void read(int direction, std::string const & = {});

    /* This is defined in balancing_workhorse.cpp
     *
     * This is a collective call, and the resulting matrix that ends up
     * in the struct is the submatrix that this job/thread uses. The
     * splitting is done according to the balancing file bfile.
     */
    void read(parallelizing_info_ptr pi, int direction, std::string const & sanity_check_vector = {});

    /* The three dump() functions dump the data in different formats.
     * dump_mixed is the legacy format, and currently this also matches
     * the in-memory format.
     * If direction==0:
     *  dump_data outputs the contents of the .rows file
     *  dump_offsets outputs the contents of the .row_offsets file
     *  dump_coeffs outputs the contents of the .row_coeffs file
     */
    void dump_data(std::ostream&) const;
    void dump_offsets(std::ostream&) const;
    void dump_coeffs(std::ostream&) const;
    void dump_mixed(std::ostream&) const;

    /* This is defined in random_matrix.cpp */
    static matrix_file build_random(parallelizing_info_ptr pi, cxx_param_list & pl, unsigned int nrows, unsigned int ncols, int coeff_bound);
    static matrix_file build_random(cxx_param_list & pl, unsigned int nrows, unsigned int ncols, int coeff_bound);

    private:

    void require_lookup() {
        if (withcoeffs == -1) {
            if (!lookup())
                throw std::runtime_error(lookup_diagnostic);
        }
        /* This is set by lookup() */
        ASSERT_ALWAYS(withcoeffs != -1);
    }
};

#endif	/* MATRIX_FILE_HPP_ */
