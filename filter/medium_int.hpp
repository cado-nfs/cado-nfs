#ifndef MEDIUM_INT_HPP_
#define MEDIUM_INT_HPP_

#include <cstdint>
#include <algorithm>

/* Something which is presumably larger than uint32_t, but less than
 * uint64_t */

template<int width>
struct uint64_fits {
    bool operator()(uint64_t a) const { return !(a >> (width * 8)); }
};
template<>
struct uint64_fits<8> {
    bool operator()(uint64_t) const { return true; }
};

template<int width>
struct medium_int {
    uint8_t x[width];
    static_assert(width <= 8, "medium_int must be smaller than uint64_t");

    medium_int() {}

    medium_int(uint64_t a) {
	ASSERT_ALWAYS(uint64_fits<width>()(a));
	memcpy(x, &a, width);
    }
    operator uint64_t () const {
	uint64_t r = 0;
	memcpy(&r, x, width);
	return r;
    }
    template<int otherwidth>
    medium_int operator=(medium_int<otherwidth> const & b) {
        if (otherwidth < width) {
            memset(x,0,width);
            memcpy(x,b.x,otherwidth);
        } else {
            memcpy(x,b.x,width);
        }
        return *this;
    }
    medium_int operator+(medium_int const &b) const {
	return medium_int((uint64_t) b + (uint64_t) *this);
    }
    medium_int & operator+=(medium_int const &b) {
	uint64_t a = (uint64_t) b + (uint64_t) *this;
	ASSERT_ALWAYS(uint64_fits<width>()(a));
	memcpy(x, &a, width);
	return *this;
    }
    medium_int operator-(medium_int const &b) const {
	return medium_int((uint64_t) b - (uint64_t) *this);
    }
    medium_int & operator-=(medium_int const &b) {
	uint64_t a = (uint64_t) b - (uint64_t) *this;
	ASSERT_ALWAYS(uint64_fits<width>()(a));
	memcpy(x, &a, width);
	return *this;
    }
};

#endif	/* MEDIUM_INT_HPP_ */
