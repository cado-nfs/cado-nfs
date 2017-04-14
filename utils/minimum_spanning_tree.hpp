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
    inline size_type dimension(int i) const { return dims[i]; }
    reference operator()(size_type i, size_type j) {
        return ((super&)*this)[i*dims[0]+j];
    }
    const_reference operator()(size_type i, size_type j) const {
        return ((super const&)*this)[i*dims[0]+j];
    }
    matrix(size_type m, size_type n, value_type v = value_type()) : super(m*n, v) { dims[0]=m; dims[1]=n; }
};

/* Of course we could templatize the MST code. No need for that presently
 */
std::pair<int, std::vector<std::pair<int, int>>> minimum_spanning_tree(matrix<int> weights);

#endif	/* MINIMUM_SPANNING_TREE_HPP_ */
