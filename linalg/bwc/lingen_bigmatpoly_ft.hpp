#ifndef BIGMATPOLY_FT_HPP_
#define BIGMATPOLY_FT_HPP_

#include "mpfq_layer.h"
#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly.hpp"
#include "flint-fft/fft.h"

#include "select_mpi.h"
#include "lingen_substep_schedule.h"
#include "tree_stats.hpp"

/* This defines an MPI-shared polynomial matrix type */

#if 0
struct bigmatpoly_ft : public bigmatpoly_model {
    abdst_field ab;
    unsigned int m;     /* total number of rows */
    unsigned int n;     /* total number of cols */
    const struct fft_transform_info * fti;
    /* The following three are also in cells */
    unsigned int m0;      /* number of rows per block */
    unsigned int n0;      /* number of cols per block */
    std::vector<matpoly_ft> cells;

    bigmatpoly_ft(bigmatpoly_model const & model);
    bigmatpoly_ft(abdst_field ab, bigmatpoly_model const & model, unsigned int m, unsigned int n, const struct fft_transform_info * fti);
    matpoly_ft & my_cell() { return cell(irank(), jrank()); }
    matpoly_ft const & my_cell() const { return cell(irank(), jrank()); }
    /* {{{ access interface for bigmatpoly_ft */
    inline matpoly_ft & cell(unsigned int i, unsigned int j) {
        return cells[i*n1+j];
    }
    inline matpoly_ft const & cell(unsigned int i, unsigned int j) const {
        return cells[i*n1+j];
    }
    /* }}} */

    void finish_init(abdst_field ab, unsigned int m, unsigned int n, const struct fft_transform_info * fti);
    int check_pre_init() const;
    void zero();

    // void swap(bigmatpoly_ft &);

    void mul2(bigmatpoly_ft & a, bigmatpoly_ft & b);
    void bigmatpoly_ft_mul2(bigmatpoly_ft & c, bigmatpoly_ft & a, bigmatpoly_ft & b);
    void dft(bigmatpoly const & a);
    void ift(bigmatpoly & a);
    void ift_mp(bigmatpoly & a, unsigned int shift);
private:
    /* all deprecated */
    void mul(bigmatpoly_ft & a, bigmatpoly_ft & b);
    void allgather_row();
    void allgather_col();
    void provision_row();
    void provision_col();
};
#endif

void bigmatpoly_mul_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

void bigmatpoly_mp_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S);

static inline void bigmatpoly_mul_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, const struct lingen_substep_schedule * S) { return bigmatpoly_mul_caching_adj(t, c, a, b, UINT_MAX, S); }

static inline void bigmatpoly_mp_caching(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, const struct lingen_substep_schedule * S) { return bigmatpoly_mp_caching_adj(t, c, a, b, UINT_MAX, S); }

#endif	/* BIGMATPOLY_FT_HPP_ */
