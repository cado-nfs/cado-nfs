#include "cado.h"

#include "mmt_vector_pair.hpp"


mmt_vector_pair::mmt_vector_pair(matmul_top_data_ptr mmt, int dir)
    : std::vector<mmt_vec_s>(mmt->nmatrices + (mmt->nmatrices & 1))
    , mmt(mmt)
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

    int nmats_odd = mmt->nmatrices & 1;

    matmul_top_matrix_ptr mptr;
    mptr = (matmul_top_matrix_ptr) mmt->matrices + (dir ? (mmt->nmatrices - 1) : 0);
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        int shared = (i == 0) && nmats_odd;
        mmt_vec_init(mmt,0,0, (*this)[i], dir ^ (i&1), shared, mptr->n[dir]);
        mmt_full_vec_set_zero((*this)[i]);

        mptr += dir ? -1 : 1;
    }
    if (nmats_odd) {
        mmt_vec_init(mmt,0,0, (*this)[mmt->nmatrices], !dir, 0, mmt->matrices[0]->n[dir]);
        mmt_full_vec_set_zero((*this)[mmt->nmatrices]);
    }
}

mmt_vector_pair::~mmt_vector_pair()
{
    for(auto & v : *this) {
        mmt_vec_clear(mmt, &v);
    }
}
