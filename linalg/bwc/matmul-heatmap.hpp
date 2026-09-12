#ifndef CADO_MATMUL_HEATMAP_HPP
#define CADO_MATMUL_HEATMAP_HPP

#include <cstdint>

#include <array>
#include <string>
#include <vector>

/* Where the time of the matrix times vector product goes, block by
 * block, as the matmul implementation understands blocks. Only
 * mm_impl=bucket fills this in for the moment, when the
 * mm_bucket_heatmap parameter is set. The result is written as a json
 * file, to be examined with scripts/bwc-heatmap/heatmap.html
 *
 * matmul_top also collects the blocks of all the submatrices of the
 * mpi/thread grid into a single file, which is the picture of the whole
 * matrix. Hence a header that lives outside of the backends.
 *
 * Nothing here is used unless the parameter is set.
 */

struct heatmap_block {
    /* Always (rows, columns) of the matrix, whatever the storage
     * ordering of the implementation that produced them. */
    uint32_t i0, i1;    /* row range */
    uint32_t j0, j1;    /* column range */
    uint64_t ncoeffs;
    uint16_t type;      /* index in heatmap_info::type_names */
    /* Both times below are wct seconds, summed over all the iterations
     * that were timed. The writer divides by the iteration count.
     *
     * "dispatch" is measured for this block. "combine" is only non-zero
     * for the vertical staircase code, where it accounts for the second
     * pass over the coefficients ; it is not measured per block, but
     * attributed (see matmul-bucket.cpp).
     */
    double dispatch;
    double combine;
};

struct heatmap_info {
    unsigned int nrows = 0;
    unsigned int ncols = 0;
    /* Sum over the blocks below, and not matmul_public::ncoeffs, so that
     * the file is self-consistent. (The two do not always agree:
     * bench_matcache and build_matcache pre-set matmul_public::ncoeffs
     * before build_cache, which adds to it.) */
    uint64_t ncoeffs = 0;
    std::array<int, 2> iterations { 0, 0 };
    /* wct of the whole cpu-bound loop, summed over all the iterations
     * that were timed, and over the submatrices when there are several */
    double total = 0;
    /* wct of the busiest submatrix alone. Equal to total when there is
     * only one, and a measure of the load imbalance otherwise. */
    double total_max = 0;
    std::string matrix; /* name of the matrix, for reference */
    std::vector<std::string> type_names;

    /* Only set when the blocks of a whole mpi/thread grid have been
     * collected: the shape of the grid, and the dimensions of one
     * submatrix. nrows and ncols are then the padded dimensions of the
     * whole matrix, which is what the grid really works on. */
    unsigned int nh = 0, nv = 0;
    unsigned int submatrix_nrows = 0, submatrix_ncols = 0;
};

/* Write the json file. filename is the user-provided name. A non-empty
 * tag identifies one particular piece of the matrix among the ones that
 * the mpi/thread grid deals with, and is inserted before the extension
 * so that concurrent instances do not fight over the same file; the
 * collected picture of the whole matrix goes to the plain name.
 */
void heatmap_dump(std::string const & filename,
        std::string const & tag,
        heatmap_info const & info,
        std::vector<heatmap_block> const & blocks);

#endif	/* CADO_MATMUL_HEATMAP_HPP */
