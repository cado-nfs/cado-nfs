#ifndef CADO_UTILS_INTEGER_PARTITIONS_HPP_
#define CADO_UTILS_INTEGER_PARTITIONS_HPP_

#include <vector>
#include "macros.h"

/* This provides a range-iterable abstractions for the partitions of a
 * positive integer n. Dereferencing the iterator gives an STL vector.
 */
struct integer_partitions {
    unsigned int n;
    unsigned int min_part;
    explicit integer_partitions(unsigned int n, unsigned int min_part=1)
        : n(n)
        , min_part(min_part)
    {
        ASSERT_ALWAYS(n > 0);
        ASSERT_ALWAYS(min_part > 0);
    }
    integer_partitions(integer_partitions const &) = default;
    integer_partitions(integer_partitions &&) = default;
    integer_partitions& operator=(integer_partitions const &) noexcept = default;
    integer_partitions& operator=(integer_partitions &&) noexcept = default;
    ~integer_partitions() = default;

    /* return the "first" (in lexicographically decreasing order)
     * partition of n that has all summands between min and max (both
     * inclusive).
     * An empty vector is returned to indicate that no such partition
     * exists.
     */
    static std::vector<unsigned int>
    first_partition(unsigned int n, unsigned int min, unsigned int max)
    {
        ASSERT_ALWAYS(n);
        ASSERT_ALWAYS(min >= 1);
        ASSERT_ALWAYS(max >= min);
        if (min == 1) {
            std::vector<unsigned int> res(n / max, max);
            if (n % max)
                res.push_back(n % max);
            return res;
        }
        if (n < min)
            return {};
        if (min <= n && n <= max)
            return { n };
        for( ; max >= min ; max--) {
            auto t = first_partition(n - max, min, max);
            if (!t.empty()) {
                t.insert(t.begin(), max);
                return t;
            }
        }
        return {};
    }

    struct const_iterator :
        /* an array, sorted in decreasing order */
        std::vector<unsigned int>
    {
        unsigned int min_part = 1;
        explicit const_iterator(unsigned int min_part)
            : min_part(min_part)
        {}
        typedef std::vector<unsigned int> super;
        // pre-increment
        const_iterator& operator++()
        {
            unsigned int r = 0;
            for(;;) {
                unsigned int s = size();
                for( ; s && (*this)[s-1] == min_part ; s--) ;
                if (!s) {
                    clear();
                    return *this;
                }
                r += (size() - s) * min_part;
                erase(super::begin() + s, super::end());
                back()--;
                r++;

                /* How can we spill r units, while retaining the property
                 * that all summands are more than min_part ? Try to
                 * see if it is at all possible, otherwise reduce the
                 * stem of the partition even more in a further turn.
                 */
                auto tail = first_partition(r, min_part, back());
                if (!tail.empty()) {
                    insert(end(), tail.begin(), tail.end());
                    return *this;
                }
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

    const_iterator begin() const {
        const_iterator ret(min_part);
        if (n >= min_part)
            ret.assign(1, n);
        return ret;
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    const_iterator end() const {
        const_iterator ret(min_part);
        return ret;
    }
};

#endif	/* CADO_UTILS_INTEGER_PARTITIONS_HPP_ */
