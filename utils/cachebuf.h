#ifndef CACHEBUF_H_
#define CACHEBUF_H_

// pragma no prototypes

#ifdef HAVE_SSE2

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <xmmintrin.h>
#include "macros.h"

#ifndef CACHELINESIZE
#define CACHELINESIZE 64
#endif

#define DECLARE_CACHE_BUFFER(__type, __indices)	\
struct __type##_##__indices##_cachebuffer_s {	\
  __type buf[__indices][CACHELINESIZE / sizeof(__type)] ATTRIBUTE((aligned (CACHELINESIZE)));	\
  unsigned char buf_full[__indices];	\
  __type **mem; \
  size_t idx_used; \
};	\
typedef struct __type##_##__indices##_cachebuffer_s __type##_##__indices##_cachebuffer[1];	\
\
static inline void \
init_##__type##_##__indices##_cachebuffer(__type##_##__indices##_cachebuffer buf, __type **mem, const size_t idx_used) \
{ \
  buf->mem = mem; \
  buf->idx_used = idx_used; \
  for (size_t i = 0; i < idx_used; i++) { \
    const uintptr_t uptr = (uintptr_t)(mem[i]); \
    ASSERT (uptr % sizeof(__type) == 0); \
    const size_t offset = (uptr % CACHELINESIZE) / sizeof(__type);  \
    buf->buf_full[i] = offset; \
    mem[i] -= offset; \
    memcpy(buf->buf[i], mem[i], offset * sizeof(__type)); \
  } \
} \
\
static inline void	\
add_##__type##_##__indices##_to_cachebuffer(__type##_##__indices##_cachebuffer buf, 	\
                    const size_t idx, const __type data)	\
{	\
    const size_t elem_per_cl = CACHELINESIZE / sizeof(__type);	\
    buf->buf[idx][buf->buf_full[idx]++] = data;	\
    if (buf->buf_full[idx] == elem_per_cl) {	\
        __m128i *ptr = (__m128i *)buf->mem[idx];	\
        __builtin_ia32_movntdq(ptr, ((__m128i *)buf->buf[idx])[0]);	\
        __builtin_ia32_movntdq(ptr + 1, ((__m128i *)buf->buf[idx])[1]);	\
        __builtin_ia32_movntdq(ptr + 2, ((__m128i *)buf->buf[idx])[2]);	\
        __builtin_ia32_movntdq(ptr + 3, ((__m128i *)buf->buf[idx])[3]);	\
        /* __builtin_ia32_clflush(ptr); */	\
        buf->buf_full[idx] = 0;	\
        buf->mem[idx] = (__type *)ptr + elem_per_cl;	\
    }	\
} \
\
static inline void \
flush_##__type##_##__indices##_cachebuffer(__type##_##__indices##_cachebuffer buf) \
{ \
  for (size_t i = 0; i < buf->idx_used; i++) { \
    ASSERT(buf->buf_full[i] < buf->idx_used); \
    memcpy(buf->mem[i], buf->buf[i], buf->buf_full[i] * sizeof(__type)); \
    buf->mem[i] += buf->buf_full[i]; \
    buf->buf_full[i] = 0; \
  } \
}

#endif /* HAVE_SSE2 */ 
#endif /* CACHEBUF_H_ */
