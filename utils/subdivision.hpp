#ifndef SUBDIVISION_HPP_
#define SUBDIVISION_HPP_

#include <tuple>
#include "macros.h"

class subdivision {
    unsigned int n = 0;
    unsigned int k = 0;
    unsigned int q = 0;
    unsigned int r = 0;
    unsigned int scale = 1;
    public:
    unsigned int total_size() const { return n; }
    unsigned int nblocks() const { return k; }
    subdivision() {}
    subdivision(unsigned int n, unsigned int k, unsigned int scale = 1)
        : n(n/scale), k(k), q(n/scale/k), r((n/scale)%k), scale(scale)
    {}
    inline unsigned int nth_block_size(unsigned int i) const
    {
        return (q + (i < r)) * scale;
    }
    inline unsigned int nth_block_start(unsigned int i) const {
        return (i * q + std::min(i, r)) * scale;
    }
    inline unsigned int nth_block_end(unsigned int i) const {
        return nth_block_start(i) + nth_block_size(i);
    }
    inline std::tuple<unsigned int, unsigned int> nth_block(unsigned int i) const
    {
        unsigned int i0 = nth_block_start(i);
        unsigned int i1 = nth_block_end(i);
        return std::make_tuple(i0, i1);
    }
    inline unsigned int block_size_upper_bound() const {
        return (q + (r != 0)) * scale;
    }
    inline unsigned int flatten(unsigned int idx, unsigned int pos) const {
        return (idx * q + std::min(idx, r)) * scale + pos;
    }
    static subdivision by_block_size(unsigned int n, unsigned int b) {
        ASSERT_ALWAYS(n > 0);
        return subdivision(n, iceildiv(n, b));
    }
    subdivision operator*(unsigned int x) const {
        subdivision res = *this;
        res.scale = scale * x;
        return res;
    }
};


#endif	/* SUBDIVISION_HPP_ */
