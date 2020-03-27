#ifndef BBLAS_MAT64_HPP_
#define BBLAS_MAT64_HPP_

#include <cstdint>

typedef uint64_t mat64[64] ATTRIBUTE((aligned(64)));
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

#endif	/* BBLAS_MAT64_HPP_ */
