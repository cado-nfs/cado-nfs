#ifndef MMT_VECTOR_PAIR_HPP_
#define MMT_VECTOR_PAIR_HPP_

#include "matmul_top.hpp"

class mmt_vector_pair : public std::vector<mmt_vec_s> {
    matmul_top_data_ptr mmt;
public:
    inline mmt_vec_ptr operator[](int i) { return &((std::vector<mmt_vec_s> &)*this)[i]; }
    // inline mmt_vec_srcptr operator[](int i) const { return &((std::vector<mmt_vec_s> const &)*this)[i]; }

    mmt_vector_pair(matmul_top_data_ptr, int);
    ~mmt_vector_pair();
    mmt_vector_pair(mmt_vector_pair const &) = delete;
    mmt_vector_pair(mmt_vector_pair &&) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair const &) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair &&) = delete;
    // mmt_vec * vectors() { return data(); }
    mmt_vec * vectors() { return reinterpret_cast<mmt_vec*>(data()); }


};

#endif	/* MMT_VECTOR_PAIR_HPP_ */
