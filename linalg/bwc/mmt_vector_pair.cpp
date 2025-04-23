#include "cado.h" // IWYU pragma: keep

#include "mmt_vector_pair.hpp"
#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"


mmt_vector_pair::mmt_vector_pair(matmul_top_data & mmt, int dir)
{
    /* we allocate as many vectors as we have matrices, plus one if the
     * number of matrices is odd (so we always have an even number of
     * vectors). If the number of matrices is odd, then
     * the first vector may be shared.  Otherwise, I believe it cannot
     * (but I'm not really sure)
     *
     * Storage for vectors need actually not be present at all times.
     * This could be improved.
     *
     * Interleaving could defined twice as many interleaved levels as we
     * have matrices. It is probably not relevant.
     */

    auto const nmats_odd = mmt.matrices.size() & 1;

    matmul_top_matrix * mptr = & mmt.matrices[dir ? (mmt.matrices.size() - 1) : 0];
    for(size_t i = 0 ; i < mmt.matrices.size() ; i++) {
        int const shared = (i == 0) && nmats_odd;
        emplace_back(mmt, nullptr, nullptr, dir ^ (i&1), shared, mptr->n[dir]);
        mmt_full_vec_set_zero(back());

        mptr += dir ? -1 : 1;
    }
    if (nmats_odd) {
        emplace_back(mmt, nullptr, nullptr, !dir, 0, mmt.matrices[0].n[dir]);
        mmt_full_vec_set_zero(back());
    }
}
