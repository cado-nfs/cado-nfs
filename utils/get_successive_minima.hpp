#ifndef GET_SUCCESSIVE_MINIMA_HPP_
#define GET_SUCCESSIVE_MINIMA_HPP_

#include <iterator>
#include <functional>
#include <algorithm>
#include <vector>
#include <queue>

/* This function returns the indices of the n smallest (in the sense
 * given by the comparison function) elements in the range [a,b).
 *
 * The cost is typically O(b-a) + O(n*log(n)), and might go up to
 * O((b-a)*log(n)) for very degenerate cases.
 *
 */
template<
    typename RandomAccessIterator,
    typename Compare = std::less<typename RandomAccessIterator::value_type>
>
std::vector<size_t>
    get_successive_minima(
            RandomAccessIterator const & a,
            RandomAccessIterator const & b,
            size_t n,
            Compare const& inner_comp = Compare())
{
    typedef typename RandomAccessIterator::value_type data_t;
    typedef std::pair<size_t, data_t> pair_t;
    /* We'll maintain a priority queue of the largest values among the
     * best score tables. Whenever we find an entry which beats it, then
     * it deserves to enter the table (and then we kick the previous
     * largest out).
     * So the default priority_queue behaviour which is a max heap, with
     * respect to our comparison operator, is exactly what we need.
     */
    struct pair_comp_t {
        Compare inner_comp;
        pair_comp_t(Compare const& c) : inner_comp(c) {}
        bool operator()(pair_t const& a, pair_t const& b) const {
            if (inner_comp(a.second, b.second)) return true;
            if (inner_comp(b.second, a.second)) return false;
            return a.first < b.first;
        }
    } pair_comp(inner_comp);
    std::priority_queue<pair_t, std::vector<pair_t>, pair_comp_t> Q(pair_comp);
    for(RandomAccessIterator x = a ; x != b ; ++x) {
        pair_t P(x-a, *x);
        if (Q.size() < n || pair_comp(P, Q.top())) {
            Q.push(P);
            if (Q.size() > n)
                Q.pop();
        }
    }
    std::vector<size_t> res(Q.size(), 0);
    for(size_t d = Q.size(); d-- ;) {
        res[d] = Q.top().first;
        Q.pop();
    }
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
    return get_successive_minima(v.begin(), v.end(), n, inner_comp);
}

#endif	/* GET_SUCCESSIVE_MINIMA_HPP_ */
