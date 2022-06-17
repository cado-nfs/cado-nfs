#include "cado.h" // IWYU pragma: keep
#include <stdint.h>    // for uint32_t
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/stat.h>
#include "matmul-mf.hpp"
#include "mf.h"
#include "portability.h" // strdup // IWYU pragma: keep
#include "macros.h"

/* This is some scotch tape between the balancing layer and the
 * lower-level matmul layer. It also uses some code from mf. This, in
 * effect, bypasses matmul_common_read_stupid_data, which is no longer
 * used.
 */
void mf_prepare_matrix_u32(matmul_ptr mm, matrix_u32_ptr m, const char * file, int withcoeffs)
{
    struct mf_io_file mf[1];
    struct mf_io_file rw[1];
    struct mf_io_file cw[1];
    memset(mf, 0, sizeof(mf));
    memset(rw, 0, sizeof(rw));
    memset(cw, 0, sizeof(cw));
    struct stat sbuf[1];
    int rc = stat(file, sbuf);
    if (rc < 0) {
        fprintf(stderr, "stat(%s): %s\n", file, strerror(errno));
        exit(1);
    }
    mf->alloc = mf->size = sbuf->st_size / sizeof(uint32_t);
    mf->p = (uint32_t *) malloc(sbuf->st_size);
    ASSERT_ALWAYS(mf->p);
    FILE * f = fopen(file, "rb");
    ASSERT_ALWAYS(f);
    int nread = fread(mf->p, sizeof(uint32_t), mf->size, f);
    if (nread < (int) mf->size) {
        fprintf(stderr, "%s: short read (%d < %" PRIu64 ")\n", file, nread, mf->size);
        exit(1);
    }
    fclose(f);
    matrix_read_pass(mf, NULL, rw, cw, 0, 0, 1, withcoeffs);

    memset(m, 0, sizeof(matrix_u32));
    m->mfile = strdup(file);
    m->bfile = NULL;
    m->transpose = mm->store_transposed;
    m->size = sbuf->st_size / sizeof(uint32_t);
    m->p = mf->p;

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
     */
    mm->dim[m->transpose] = rw->size;
    mm->dim[!m->transpose] = cw->size;
}

