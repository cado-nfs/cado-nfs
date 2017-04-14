#ifndef MINIMUM_SPANNING_TREE_HPP_
#define MINIMUM_SPANNING_TREE_HPP_

#include <vector>
#include <cstddef>

template<typename T>
struct matrix : private std::vector<T>
{
private:
    typedef std::vector<T> super;
public:
    using typename super::size_type;
    using typename super::reference;
    using typename super::const_reference;
    using typename super::pointer;
    using typename super::const_pointer;
    using typename super::value_type;
private:
    size_type dims[2];
public:
    reference operator()(size_type i, size_type j) {
        return ((super&)*this)[i*dims[0]+j];
    }
    const_reference operator()(size_type i, size_type j) const {
        return ((super const&)*this)[i*dims[0]+j];
    }
    matrix(size_type m, size_type n) : super(m*n) {}
};

std::pair<int, std::vector<std::pair<int, int>>>
minimalSpanningTreePrimNaive (std::vector<int> weights, size_t m);

#endif	/* MINIMUM_SPANNING_TREE_HPP_ */
