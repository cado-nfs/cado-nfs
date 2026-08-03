#ifndef CADO_SUBDIVISION_HPP
#define CADO_SUBDIVISION_HPP

#include <algorithm>
#include <concepts>

#include "macros.h"

template<std::unsigned_integral T = unsigned int>
class subdivision {
    T total_elements = 0;
    T k = 0;
    T q = 0;
    T r = 0;
    T scale = 1;
    public:
    struct Interval { T start; T end; };

    constexpr T total_size() const noexcept { return total_elements; }
    constexpr T nblocks() const noexcept { return k; }
    constexpr subdivision() noexcept = default;
    constexpr subdivision(T n, T k, T scale = 1) noexcept
        : total_elements(n)
        , k(k)
        , scale(scale)
    {
        if (k > 0) {
            T scaled_n = n / scale;
            q = scaled_n / k;
            r = scaled_n % k;
        }
    }
    constexpr T nth_block_size(T i) const noexcept
    {
        return (q + (i < r)) * scale;
    }
    constexpr T nth_block_start(T i) const noexcept
    {
        return (i * q + std::min(i, r)) * scale;
    }
    constexpr T nth_block_end(T i) const noexcept
    {
        return nth_block_start(i) + nth_block_size(i);
    }
    constexpr Interval nth_block(T i) const noexcept {
        return { nth_block_start(i), nth_block_end(i) };
    }
    constexpr T block_size_upper_bound() const noexcept
    {
        return (q + (r != 0)) * scale;
    }
    constexpr T flatten(T idx, T pos) const noexcept
    {
        return (idx * q + std::min(idx, r)) * scale + pos;
    }
    static subdivision by_block_size(T n, T b) {
        ASSERT_ALWAYS(n > 0);
        return { n, iceildiv(n, b) };
    }
    constexpr subdivision operator*(T x) const noexcept {
        subdivision res = *this;
        res.scale = scale * x;
        return res;
    }
    struct iterator {
        const subdivision* self = nullptr;
        T idx = 0;

        // Custom struct or std::pair enables structured bindings [i0, i1]
        constexpr auto operator*() const noexcept { 
            return self->nth_block(idx); 
        }

        constexpr iterator& operator++() noexcept { 
            ++idx; 
            return *this; 
        }

        constexpr bool operator!=(const iterator& other) const noexcept { 
            return idx != other.idx; 
        }
    };

    constexpr iterator begin() const noexcept { return { this, 0 }; }
    constexpr iterator end()   const noexcept { return { this, k }; }
};


#endif	/* CADO_SUBDIVISION_HPP */
