#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>

#include <string>
#include <memory>

#include <sys/stat.h>

#include "fmt/base.h"

#include "utils_cxx.hpp"

#include "matrix_u32.hpp"
#include "misc.h"
#include "macros.h"

static size_t number_of_words(std::string const & filename)
{
    struct stat sbuf[1];

    int const rc = stat(filename.c_str(), sbuf);
    if (rc < 0) {
        fmt::print(stderr, "stat({}): {}\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return sbuf->st_size / sizeof(uint32_t);
}

/* This used to be matrix_read_passed (now gone), called from
 * mf_prepare_matrix_u32 (now gone as well)
 */
matrix_u32::with_dimensions
matrix_u32::from_file(std::string const & mfile,
        transpose_option const & transpose,
        withcoeffs_option const & withcoeffs)
{
    matrix_u32 m(/* transpose, */ withcoeffs);
    m.mfile = mfile;
    std::string const rwfile = std::unique_ptr<char, free_delete<char>>(
            derived_filename(mfile.c_str(), "rw", ".bin")).get();
    std::string const cwfile = std::unique_ptr<char, free_delete<char>>(
            derived_filename(mfile.c_str(), "cw", ".bin")).get();

    size_t const msize = number_of_words(mfile);

    /* if the matrix is stored column major, then the rw/cw file names
     * are swapped, becaused I prefer to _always_ have the simple rules
     * that msize = rwsize + (1 + withcoeffs) * withcoeffs; but in
     * effect, in that case (column major), rwsize is the number of
     * columns.
     */
    unsigned int const rwsize = number_of_words(rwfile);
    unsigned int const cwsize = number_of_words(cwfile);

    ASSERT_ALWAYS(!withcoeffs || (msize - rwsize) % 2);
    size_t const ncoeffs = (msize - rwsize) / (1 + m.withcoeffs);
    m.p.assign(msize, 0);

    fmt::print("Reading {}\n", mfile);
    std::unique_ptr<FILE, delete_FILE> const f(fopen(mfile.c_str(), "rb"));
    ASSERT_ALWAYS(f);
    size_t const nread = fread(m.p.data(), sizeof(uint32_t), msize, f.get());
    if (nread < msize) {
        fmt::print(stderr, "{}: short read ({} < {})\n", mfile, nread, msize);
        exit(EXIT_FAILURE);
    }


    /* Beware. We really have dim[0] = nrows and dim[1] = ncols as far as
     * this matrix bit is concerned. However, when m->transpose == 1
     * the submatrix file which gets written on disk is transposed. This
     * is the reason why in this case we need to swap indices, and work
     * de facto with a transposed matrix.
     *
     * Let's state it a second time. A bit matrix of size 100M * 99M, if
     * split over a processor grid of size 100x100, and with the
     * implementation expecting transposed matrix data, will (with
     * save_submatrices=1) save submatrices of size 100k * 99k, but in
     * column major order. In turn,
     * if we use the bench program on such submatrices, we will (as per
     * matrix_read_pass above) read them in row major order to count a
     * ``number of rows'' (in fact, columns), and a ``number of columns''
     * (in fact, rows).
     *
     * XXX These fields used to be set from mf_prepare_matrix_u32, in the
     * matmul_ptr field. Now this step no longer takes place here, so
     * this all needs a new home.
     */

    if (!transpose) {
        return { m, { rwsize, cwsize}, ncoeffs };
    } else {
        return { m, { cwsize, rwsize}, ncoeffs };
    }
}
