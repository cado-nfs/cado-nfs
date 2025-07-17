#ifndef CADO_MATRIX_U32_HPP
#define CADO_MATRIX_U32_HPP

#include "params.h"
#include "parallelizing_info.hpp"

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <tuple>
#include <array>

struct random_matrix_ddata;

struct matrix_u32 {
    // input arguments.
    std::string mfile;    // matrix file name
    std::string bfile;    // balancing file name ; NULL will mean auto-detect

    /* The 'transpose' flag just tells that the in-memory representation
     * of the sparse matrix below is given by columns and not by rows.
     * (perhaps it should be called "column_major" instead).
     *
     * This is sometimes the format that is preferred by some MM layers,
     * and cado-nfs decides to present the local matrices to the MM
     * layers precisely in their preferred format.
     */

    struct transpose_option {
        bool value = false;
        explicit operator bool() const { return value; }
    };

    struct withcoeffs_option {
        bool value = false;
        explicit operator bool() const { return value; }
    };

    // it's not absolutely obvious that we need to keep the [transpose]
    // field (even though we do use it in some of the named ctors).
    // Oftentimes (always?) the mm layer is around, and tells us the
    // expected layout for the matrix.

    // bool transpose;
    bool withcoeffs;

    // output arguments.
    std::vector<uint32_t> p;

    typedef std::tuple<matrix_u32, std::array<unsigned int, 2>, uint64_t>
        with_dimensions;

    static with_dimensions from_file(
            std::string const & mfile,
            transpose_option const & transpose = { false },
            withcoeffs_option const & withcoeffs = { false });

    /* We mark it as deprecated because this function has only two users,
     * currently, and it's quite telling that those are precisely missing
     * the info that we don't know where to put, about the dimensions of
     * the matrix that we just read. We must find a way.
     */
    /*
    explicit matrix_u32(std::string const & mfile,
            // transpose_option const & transpose = { false },
            withcoeffs_option const & withcoeffs = { false }) __attribute__((deprecated));
            */

    /* A few interfaces are allowed to call the default ctor */
    friend struct random_matrix_ddata;
    friend matrix_u32 balancing_get_matrix_u32(
            parallelizing_info_ptr pi,
            cxx_param_list & pl,
            std::string const & mfile,
            std::string const & bfile,
            bool withcoeffs,
            bool transpose_while_dispatching);

    private:
    explicit matrix_u32(
            // transpose_option const & transpose = { false },
            withcoeffs_option const & withcoeffs = { false })
        : withcoeffs(withcoeffs)
        {}
};

#endif	/* MATRIX_U32_HPP_ */
