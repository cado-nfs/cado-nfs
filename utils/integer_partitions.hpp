#ifndef CADO_UTILS_INTEGER_PARTITIONS_HPP
#define CADO_UTILS_INTEGER_PARTITIONS_HPP

#include <vector>
#include <algorithm>

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
        /* n == 0 is supported */
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
        typedef std::vector<unsigned int> super;
        bool operator==(const_iterator & o) const {
            if (is_end || o.is_end) return is_end && o.is_end;
            return size() == o.size() && std::ranges::equal(*this, o);
        }
        bool operator!=(const_iterator & o) const {
            return !operator==(o);
        }
        unsigned int min_part = 1;
        bool is_end = false;
        explicit const_iterator(unsigned int min_part)
            : min_part(min_part)
        {}
        // pre-increment
        const_iterator& operator++()
        {
            unsigned int r = 0;
            for(;;) {
                unsigned int s = size();
                for( ; s && (*this)[s-1] == min_part ; s--) ;
                if (!s) {
                    clear();
                    is_end = true;
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
        else if (n)
            ret.is_end = true;
        return ret;
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    const_iterator end() const {
        const_iterator ret(min_part);
        ret.is_end = true;
        return ret;
    }
};

struct integer_partitions_in_k_parts {
    unsigned int n;
    unsigned int k;
    unsigned int min_part;
    integer_partitions_in_k_parts(unsigned int n, unsigned int k, unsigned int min_part=1)
        : n(n)
        , k(k)
        , min_part(min_part)
    {
        /* n == 0 is supported */
        ASSERT_ALWAYS(min_part > 0);
    }
    integer_partitions_in_k_parts(integer_partitions_in_k_parts const &) = default;
    integer_partitions_in_k_parts(integer_partitions_in_k_parts &&) = default;
    integer_partitions_in_k_parts& operator=(integer_partitions_in_k_parts const &) noexcept = default;
    integer_partitions_in_k_parts& operator=(integer_partitions_in_k_parts &&) noexcept = default;
    ~integer_partitions_in_k_parts() = default;

    struct const_iterator :
        /* an array, sorted in decreasing order */
        std::vector<unsigned int>
    {
        typedef std::vector<unsigned int> super;
        bool operator==(const_iterator & o) const {
            if (is_end || o.is_end) return is_end && o.is_end;
            return size() == o.size() && std::ranges::equal(*this, o);
        }
        bool operator!=(const_iterator & o) const {
            return !operator==(o);
        }
        unsigned int min_part = 1;
        bool is_end = false;
        explicit const_iterator(unsigned int min_part)
            : min_part(min_part)
        {}
        // pre-increment
        const_iterator& operator++()
        {
            const unsigned int k = size();
            unsigned int t;
            super & l = *this;
            unsigned int spill = 0;
            for (t = k - 1; t > 0; t--) {
                /* try to increase l[t-1] and decrease l[t] */
                if (l[t-1] + 1 <= l[t])
                {
                    l[t-1] ++;
                    unsigned int s = 1;
                    unsigned int dec = 0;
                    for (unsigned int u = t; u < k; u++)
                    {
                        dec += l[u];
                        l[u] = l[u-1];
                        s += l[u];
                    }
                    /* recover spillover from previous iteration */
                    s += spill;
                    /* we also lost some */
                    if (dec <= s) {
                        s -= dec;
                        dec = 0;
                    } else {
                        dec -= s;
                        s = 0;
                    }
                    /* the total has increased by s-dec */
                    if (l[k-1] + dec >= s) {
                        l[k-1] += dec;
                        l[k-1] -= s;
                        spill = 0;
                        if (l[k-2] <= l[k-1])
                            break;
                    } else {
                        spill = s - (l[k-1] + dec);
                        l[k-1] = 0;
                    }
                }
            }
            is_end = t == 0;
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
        if (n < k * min_part) {
            ret.is_end = true;
        } else {
            ret.assign(k, min_part);
            ret[k-1] += n - k * min_part;
        }
        return ret;
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    const_iterator end() const {
        const_iterator ret(min_part);
        ret.is_end = true;
        return ret;
    }
};

#endif	/* CADO_UTILS_INTEGER_PARTITIONS_HPP_ */
