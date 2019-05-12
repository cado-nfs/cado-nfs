#ifndef BIGMATPOLY_HPP_
#define BIGMATPOLY_HPP_

#include <vector>

#include "mpfq_layer.h"
#include "lingen-matpoly.hpp"

#include "select_mpi.h"

/* This defines an MPI-shared polynomial matrix type */

struct bigmatpoly_model {
    unsigned int m1 = 0;      /* number of block rows, index i */
    unsigned int n1 = 0;      /* number of block cols, index j */
    MPI_Comm com[3] = {};     /* [0]: MPI_COMM_WORLD, reordered.
                                 [1]: row-wise       ; size == n1
                                 [2]: column-wise    ; size == m1
                                 */
    bigmatpoly_model() = default;
    bigmatpoly_model(MPI_Comm * comm, unsigned int m, unsigned int n);

    int rank() const;
    int irank() const;
    int jrank() const;
};

struct bigmatpoly : public bigmatpoly_model {
    abdst_field ab = NULL;
    unsigned int m = 0;     /* total number of rows */
    unsigned int n = 0;     /* total number of cols */
    /* The following three are also in cells */
    unsigned int m0 = 0;      /* number of rows per block */
    unsigned int n0 = 0;      /* number of cols per block */
    size_t size = 0;
    std::vector<matpoly> cells;

    bigmatpoly(bigmatpoly_model const &);
    bigmatpoly(abdst_field, bigmatpoly_model const &, unsigned int m, unsigned int n, int len);
    bigmatpoly(bigmatpoly const&) = delete;
    bigmatpoly& operator=(bigmatpoly const&) = delete;
    bigmatpoly(bigmatpoly &&);
    bigmatpoly& operator=(bigmatpoly &&);
    matpoly & my_cell() { return cell(irank(), jrank()); }
    matpoly const & my_cell() const { return cell(irank(), jrank()); }
    /* {{{ access interface for bigmatpoly */
    inline matpoly & cell(unsigned int i, unsigned int j) {
        return cells[i*n1+j];
    }
    inline matpoly const & cell(unsigned int i, unsigned int j) const {
        return cells[i*n1+j];
    }
    /* }}} */

    void finish_init(abdst_field ab, unsigned int m, unsigned int n, int len);
    int check_pre_init() const;
    // void realloc(int newalloc);
    void zero();

    void set_size(size_t size);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero_loc(unsigned int k);
    int bigmatpoly_coeff_is_zero(abdst_field ab, bigmatpoly const & pi, unsigned int k);
    
    // void swap(bigmatpoly & b);

    void truncate_loc(bigmatpoly & src, unsigned int size);
    void rshift(bigmatpoly & src, unsigned int k);

    void mul(bigmatpoly & a, bigmatpoly & b);
    void mp(bigmatpoly & a, bigmatpoly & b);

    void scatter_mat(matpoly const & src);
    void gather_mat(matpoly & dst) const;

    void scatter_mat_partial(matpoly const & src, size_t offset, size_t length);
    void gather_mat_partial(matpoly & dst, size_t offset, size_t length) const;

    private:
    int provisioned() const;
    void provision_row();
    void provision_col();
    void allgather_row();
    void allgather_col();
};

#endif	/* BIGMATPOLY_HPP_ */
