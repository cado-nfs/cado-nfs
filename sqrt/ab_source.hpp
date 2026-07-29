#ifndef SQRT_AB_SOURCE_HPP_
#define SQRT_AB_SOURCE_HPP_

#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

#include "select_mpi.h"

struct ab_source {
    std::string fname0;
    std::string sname;
    size_t sname_len = 0;
    size_t nab = 0;
    std::string prefix;
    int depnum = 0;

    size_t nab_estim = 0;

    std::vector<size_t> file_bases;
    int nfiles = 0; // 0 for cado format.
    size_t totalsize = 0;

    FILE * f = nullptr;
    int c = 0;
    size_t cpos = 0;        // position within current file.
    size_t tpos = 0;        // position within totality.

    ab_source() = default;
    ab_source(std::string const & fname, int rank, int root, MPI_Comm comm);

    ab_source(ab_source const & o);
    ab_source & operator=(ab_source const & o);

    ab_source(ab_source && o) noexcept;
    ab_source & operator=(ab_source && o) noexcept;

    ~ab_source();

    void init(std::string const & fname, int rank, int root, MPI_Comm comm);
    void rewind();
    int next(int64_t & a, uint64_t & b);
    void move_afterpos(size_t offset);

private:
    int openfile_internal();
};

#endif	/* SQRT_AB_SOURCE_HPP_ */

