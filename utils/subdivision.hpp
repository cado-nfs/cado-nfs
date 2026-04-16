#ifndef CADO_SUBDIVISION_HPP
#define CADO_SUBDIVISION_HPP

#include <algorithm>
#include <tuple>
#include <concepts>

#include "macros.h"

template<std::unsigned_integral T = unsigned int>
class subdivision {
    T n = 0;
    T k = 0;
    T q = 0;
    T r = 0;
    T scale = 1;
    public:
    T total_size() const { return n; }
    T nblocks() const { return k; }
    subdivision() = default;
    subdivision(T n, T k, T scale = 1)
        : n(n/scale), k(k), q(n/scale/k), r((n/scale)%k), scale(scale)
    {}
    T nth_block_size(T i) const
    {
        return (q + (i < r)) * scale;
    }
    T nth_block_start(T i) const {
        return (i * q + std::min(i, r)) * scale;
    }
    T nth_block_end(T i) const {
        return nth_block_start(i) + nth_block_size(i);
    }
    std::tuple<T, T> nth_block(T i) const
    {
        const T i0 = nth_block_start(i);
        const T i1 = nth_block_end(i);
        return std::make_tuple(i0, i1);
    }
    T block_size_upper_bound() const {
        return (q + (r != 0)) * scale;
    }
    T flatten(T idx, T pos) const {
        return (idx * q + std::min(idx, r)) * scale + pos;
    }
    static subdivision by_block_size(T n, T b) {
        ASSERT_ALWAYS(n > 0);
        return { n, iceildiv(n, b) };
    }
    subdivision operator*(T x) const {
        subdivision res = *this;
        res.scale = scale * x;
        return res;
    }
};


#endif	/* CADO_SUBDIVISION_HPP */
