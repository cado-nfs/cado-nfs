#ifndef CADO_UTILS_INTEGER_PARTITIONS_HPP_
#define CADO_UTILS_INTEGER_PARTITIONS_HPP_

#include <vector>
#include "macros.h"

/* This provides a range-iterable abstractions for the partitions of a
 * positive integer n. Dereferencing the iterator gives an STL vector.
 */
struct integer_partitions {
    unsigned int n;
    explicit integer_partitions(unsigned int n)
        : n(n)
    {
        ASSERT_ALWAYS(n);
    }
    integer_partitions(integer_partitions const &) = default;
    integer_partitions(integer_partitions &&) = default;
    integer_partitions& operator=(integer_partitions const &) noexcept = default;
    integer_partitions& operator=(integer_partitions &&) noexcept = default;
    ~integer_partitions() = default;

    struct const_iterator :
        /* an array, sorted in decreasing order */
        std::vector<unsigned int>
    {
        typedef std::vector<unsigned int> super;
        // pre-increment
        const_iterator& operator++()
        {
            unsigned int s = size();
            for( ; s && (*this)[s-1] == 1 ; s--) ;
            if (s == 0) {
                clear();
            } else {
                unsigned int r = size() - s;
                erase(super::begin() + s, super::end());
                back()--;
                r++;
                for( ; back() < r ; ) {
                    push_back(back());
                    r -= back();
                }
                if (r)
                    push_back(r);
            }
            return *this;
        }

        // post-increment
        const_iterator operator++(int)
        {
            const_iterator v = *this;
            ++*this;
            return v;
        }

        super const & operator*() const { return (super const &) *this; }
    };

    const_iterator begin() const { const_iterator ret; ret.assign(1, n); return ret; }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    const_iterator end() const { const_iterator ret; return ret; }
};

#endif	/* CADO_UTILS_INTEGER_PARTITIONS_HPP_ */
