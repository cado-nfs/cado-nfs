#ifndef BIGMATPOLY_HPP_
#define BIGMATPOLY_HPP_

#include <vector>

#ifdef SELECT_MPFQ_LAYER_u64k1
#include "lingen_matpoly_binary.hpp"
#else
#include "lingen_matpoly.hpp"
#endif

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
    bigmatpoly_model & get_model() { return *this; }
    bigmatpoly_model const & get_model() const { return *this; }

    bool is_square() const { return m1 == n1; }
    bool operator==(bigmatpoly_model const & o) const {
        return memcmp(this, &o, sizeof(*this)) == 0;
    }
    int rank() const;
    int irank() const;
    int jrank() const;
};

/* TODO: it now seems that the entire cells vector can be dropped in
 * favor of just storing the local cell. This would simplify the
 * implementation significantly.
 */
struct bigmatpoly : public bigmatpoly_model {
    abdst_field ab = NULL;
    unsigned int m = 0;     /* total number of rows */
    unsigned int n = 0;     /* total number of cols */
    /* The following three are also in cells */
    unsigned int m0 = 0;      /* upper bound on the number of rows per block */
    unsigned int n0 = 0;      /* upper bound on the number of cols per block */
    /* This gives the *relevant* entries. On some nodes, du to rounding,
     * not all entries in the local cell are relevant */
    unsigned int m0r() const { return subdivision(m, m1).nth_block_size(irank()); }
    unsigned int n0r() const { return subdivision(n, n1).nth_block_size(jrank()); }
    private:
    size_t size = 0;
    std::vector<matpoly> cells;
    public:

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
    bool check_pre_init() const { return size == 0; }

    // void realloc(int newalloc);
    inline void shrink_to_fit() { my_cell().shrink_to_fit(); }
    void zero();

    inline size_t get_size() const { return size; }
    void set_size(size_t size);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero_loc(unsigned int k);
    int bigmatpoly_coeff_is_zero(abdst_field ab, bigmatpoly const & pi, unsigned int k);
    
    // void swap(bigmatpoly & b);

    void truncate_loc(bigmatpoly & src, unsigned int size);
    void rshift(bigmatpoly & src, unsigned int k);

#if 0
    void mul(bigmatpoly & a, bigmatpoly & b);
    void mp(bigmatpoly & a, bigmatpoly & b);
#endif

    void scatter_mat(matpoly const & src);
    void gather_mat(matpoly & dst) const;

    void scatter_mat_partial(matpoly const & src, size_t offset, size_t length);
    void gather_mat_partial(matpoly & dst, size_t offset, size_t length) const;

    bigmatpoly truncate_and_rshift(unsigned int truncated_size, unsigned int shiftcount);

#if 0
    private:
    int provisioned() const;
    void provision_row();
    void provision_col();
    void allgather_row();
    void allgather_col();
#endif
};

#endif	/* BIGMATPOLY_HPP_ */
