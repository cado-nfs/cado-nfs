#ifndef CADO_LINGEN_BIGMATPOLY_HPP
#define CADO_LINGEN_BIGMATPOLY_HPP

#include <cstring>

#include <vector>

#include "lingen_matpoly_select.hpp"
#include "lingen_call_companion.hpp"
#include "subdivision.hpp"
#include "select_mpi.h"
#include "macros.h"

class tree_stats;

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
template<bool is_binary_arg>
class bigmatpoly : public bigmatpoly_model {
    public:
    static constexpr bool is_binary = is_binary_arg;
    typedef matpoly<is_binary> matpoly_type;
    typedef typename matpoly_type::arith_hard arith_hard_t;

    arith_hard_t * ab = nullptr;
    unsigned int m = 0;     /* total number of rows */
    unsigned int n = 0;     /* total number of cols */
    /* The following three are also in cells */
    unsigned int m0 = 0;      /* upper bound on the number of rows per block */
    unsigned int n0 = 0;      /* upper bound on the number of cols per block */
    /* This gives the *relevant* entries. On some nodes, du to rounding,
     * not all entries in the local cell are relevant */
    unsigned int m0r() const { return subdivision(m, m1).nth_block_size(irank()); }
    unsigned int n0r() const { return subdivision(n, n1).nth_block_size(jrank()); }
    unsigned int nrows() const { return m; }
    unsigned int ncols() const { return n; }
    private:
    size_t size = 0;
    std::vector<matpoly_type> cells;
    public:

    explicit bigmatpoly(bigmatpoly_model const &);
    bigmatpoly(arith_hard_t *, bigmatpoly_model const &, unsigned int m, unsigned int n, int len);
    bigmatpoly similar_shell() const { return { ab, get_model(), m, n, 0 }; }
    bigmatpoly(bigmatpoly const&) = delete;
    bigmatpoly& operator=(bigmatpoly const&) = delete;
    bigmatpoly(bigmatpoly &&) noexcept;
    bigmatpoly& operator=(bigmatpoly &&) noexcept;
    ~bigmatpoly() = default;
    matpoly_type & my_cell() { return cell(irank(), jrank()); }
    matpoly_type const & my_cell() const { return cell(irank(), jrank()); }
    /* {{{ access interface for bigmatpoly */
    private:
    matpoly_type & cell(unsigned int i, unsigned int j) {
        return cells[i*n1+j];
    }
    matpoly_type const & cell(unsigned int i, unsigned int j) const {
        return cells[i*n1+j];
    }
    public:
    /* }}} */

    void finish_init(arith_hard_t * ab, unsigned int m, unsigned int n, int len);
    bool check_pre_init() const { return size == 0; }

    // void realloc(int newalloc);
    void shrink_to_fit() { my_cell().shrink_to_fit(); }
    void zero();
    void clear() { *this = bigmatpoly(get_model()); }

    size_t get_size() const { return size; }
    void set_size(size_t size);
    void zero_pad(size_t size);
    void zero_with_size(size_t size) { set_size(0); zero_pad(size); }
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero_loc(unsigned int k);
    int bigmatpoly_coeff_is_zero(arith_hard_t * ab, bigmatpoly const & pi, unsigned int k);
    /* not to be confused with the former. the following two are in fact
     * relevant only to the binary interface. They're noops in the prime
     * field case. Here we're just agnostic, so we'll pass on the action
     * to the underlying layer.
     */
    bool high_word_is_clear() const;
    void clear_high_word();
    
    // void swap(bigmatpoly & b);

    void truncate(bigmatpoly const & src, unsigned int size);
    void rshift(bigmatpoly & src, unsigned int k);

    static bigmatpoly mul(bigmatpoly & a, bigmatpoly & b) ATTRIBUTE_DEPRECATED;
    static bigmatpoly mp(bigmatpoly & a, bigmatpoly & b) ATTRIBUTE_DEPRECATED;

    static bigmatpoly mp(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);
    static bigmatpoly mul(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M);

    void scatter_mat(matpoly_type const & src);
    void gather_mat(matpoly_type & dst) const;

    void scatter_mat_partial(matpoly_type const & src, size_t src_k, size_t offset, size_t length);
    void gather_mat_partial(matpoly_type & dst, size_t dst_k, size_t offset, size_t length) const;

    bigmatpoly truncate_and_rshift(unsigned int truncated_size, unsigned int shiftcount);

    private:
#if 0
    int provisioned() const;
#endif
    void provision_row();
    void provision_col();
    void allgather_row();
    void allgather_col();
};

#endif	/* LINGEN_BIGMATPOLY_HPP_ */
