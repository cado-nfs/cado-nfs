#ifndef CADO_MATMUL_BUCKET_HEATMAP_HPP
#define CADO_MATMUL_BUCKET_HEATMAP_HPP

#include <cstdint>

#include <array>
#include <string>
#include <vector>

/* Data collected by mm_impl=bucket when the mm_bucket_heatmap parameter
 * is set. One record per block of the matrix, as the bucket
 * implementation understands blocks. The result is written as a json
 * file, to be examined with scripts/bwc-heatmap/heatmap.html
 *
 * Nothing here is used unless the parameter is set.
 */

struct heatmap_block {
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
    double total = 0;   /* wct of the whole cpu-bound loop, all iterations */
    std::string matrix; /* name of the local matrix, for reference */
    std::vector<std::string> type_names;
};

/* Write the json file. filename is the user-provided name, and tag
 * identifies this particular piece of the matrix among the ones that the
 * mpi/thread grid deals with: it is inserted before the extension, so
 * that concurrent instances do not fight over the same file.
 */
void heatmap_dump(std::string const & filename,
        std::string const & tag,
        heatmap_info const & info,
        std::vector<heatmap_block> const & blocks);

#endif	/* CADO_MATMUL_BUCKET_HEATMAP_HPP */
