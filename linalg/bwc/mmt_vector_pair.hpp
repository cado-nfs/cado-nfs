#ifndef MMT_VECTOR_PAIR_HPP_
#define MMT_VECTOR_PAIR_HPP_

#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"

class mmt_vector_pair : public std::vector<mmt_vec> {
public:
    inline mmt_vec & operator[](int i) { return ((std::vector<mmt_vec> &)*this)[i]; }
    inline mmt_vec const & operator[](int i) const { return ((std::vector<mmt_vec> const &)*this)[i]; }

    mmt_vector_pair(matmul_top_data &, int);
    mmt_vector_pair(mmt_vector_pair const &) = delete;
    mmt_vector_pair(mmt_vector_pair &&) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair const &) = delete;
    mmt_vector_pair& operator=(mmt_vector_pair &&) = delete;
    // mmt_vec * vectors() { return data(); }
    mmt_vec * vectors() { return reinterpret_cast<mmt_vec*>(data()); }


};

#endif	/* MMT_VECTOR_PAIR_HPP_ */
