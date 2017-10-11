#ifndef GET_SUCCESSIVE_MINIMA_HPP_
#define GET_SUCCESSIVE_MINIMA_HPP_

#include <iterator>
#include <functional>
#include <algorithm>
#include <vector>
#include <queue>
#include <type_traits>

/* This function returns the indices of the n smallest (in the sense
 * given by the comparison function) elements in the range [a,b).
 *
 * a and b must be input iterators whose value type is pair<T, U>.
 * Comparison is done on the U member, and the T indices are returned. As
 * such, it is well suited to an iterator to an std::map. 
 *
 * The cost is typically O(b-a) + O(n*log(n)), and might go up to
 * O((b-a)*log(n)) for very degenerate cases.
 *
 */
template<
    typename InputIterator,
    typename ValueCompare = std::less<typename InputIterator::value_type::second_type>,
    typename KeyCompare = std::less<typename InputIterator::value_type::first_type>
>
std::vector<typename std::remove_const<typename InputIterator::value_type::first_type>::type>
    get_successive_minima(
            InputIterator const & a,
            InputIterator const & b,
            size_t n,
            ValueCompare const& comp2 = ValueCompare(),
            KeyCompare const& comp1 = KeyCompare())
{
    typedef std::pair<typename std::remove_const<typename InputIterator::value_type::first_type>::type,typename std::remove_const<typename InputIterator::value_type::second_type>::type> data_t;
    typedef std::vector<typename std::remove_const<typename InputIterator::value_type::first_type>::type> return_type;
    if (!n) return return_type();
    /* We'll maintain a priority queue of the largest values among the
     * best score tables. Whenever we find an entry which beats it, then
     * it deserves to enter the table (and then we kick the previous
     * largest out).
     * So the default priority_queue behaviour which is a max heap, with
     * respect to our comparison operator, is exactly what we need.
     */
    struct data_comp_t {
        KeyCompare comp1;
        ValueCompare comp2;
        data_comp_t(ValueCompare const& comp2 = ValueCompare(),
                            KeyCompare const& comp1 = KeyCompare())
            : comp1(comp1), comp2(comp2) {}
        bool operator()(data_t const& a, data_t const& b) const {
            if (comp2(a.second, b.second)) return true;
            if (comp2(b.second, a.second)) return false;
            if (comp1(a.first, b.first)) return true;
            if (comp1(b.first, a.first)) return false;
            return false;
        }
    } data_comp(comp2, comp1);
    std::priority_queue<data_t, std::vector<data_t>, data_comp_t> Q(data_comp);
    for(InputIterator x = a ; x != b ; ++x) {
        data_t P(*x);
        if (Q.size() < n || data_comp(P, Q.top())) {
            Q.push(P);
            if (Q.size() > n)
                Q.pop();
        }
    }
    return_type res(Q.size(), 0);
    for(size_t d = Q.size(); d-- ; Q.pop())
        res[d] = Q.top().first;
    return res;
}
template<
    typename T,
    typename Compare = std::less<T>
>
std::vector<size_t>
    get_successive_minima(
            std::vector<T> const & v,
            size_t n,
            Compare const& inner_comp = Compare())
{
    typedef typename std::vector<T>::const_iterator unindexed;
    struct indexed_iterator : public unindexed {
      typedef std::pair<size_t, T> value_type;
      std::vector<T> const& V;
      indexed_iterator(std::vector<T> const& V, unindexed const & x) : unindexed(x), V(V) {}
      value_type operator*() {
          unindexed& s(*this);
          return std::make_pair(s - V.begin(), *s);
      }
    };
    indexed_iterator v0(v, v.begin());
    indexed_iterator v1(v, v.end());
    return get_successive_minima(v0, v1, n, inner_comp);
}

#endif	/* GET_SUCCESSIVE_MINIMA_HPP_ */
