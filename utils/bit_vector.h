#ifndef BIT_VECTOR_H_
#define BIT_VECTOR_H_

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

#define BV_BITS 64      // since we're using uint64_t's
#define LN2_BV_BITS 6   // 2^^LN2_BV_BITS = BV_BITS
/* Changing bv_t to something else is possibly dangerous and should not
 * be taken lightly. Some code down the line may make indirect
 * assumptions on bv_t being uint64_t */
typedef uint64_t bv_t;

struct bit_vector_s {
    bv_t *p;
    size_t n;
};
typedef struct bit_vector_s bit_vector[1];
typedef struct bit_vector_s * bit_vector_ptr;
typedef const struct bit_vector_s * bit_vector_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void bit_vector_init(bit_vector_ptr b, size_t n);
extern void bit_vector_init_set(bit_vector_ptr b, size_t n, int s);
extern void bit_vector_set(bit_vector_ptr b, int s);
extern void bit_vector_clear(bit_vector_ptr b);
extern void bit_vector_neg(bit_vector_ptr b, bit_vector_srcptr c);

extern int bit_vector_getbit(bit_vector_srcptr b, size_t pos);
/* The value returned by the two following functions reflect the _old_
 * value of the flag, which is overwritten by the function */
extern int bit_vector_setbit(bit_vector_ptr b, size_t pos);
extern int bit_vector_clearbit(bit_vector_ptr b, size_t pos);

/* In contrast, this returns the new value */
extern int bit_vector_flipbit(bit_vector_ptr b, size_t pos);

extern size_t bit_vector_popcount(bit_vector_ptr b);
extern size_t bit_vector_memory_footprint(bit_vector_srcptr b);

extern void bit_vector_read_from_stream(bit_vector_ptr b, FILE * f);
extern void bit_vector_write_to_stream(bit_vector_srcptr b, FILE * f);

#ifdef __cplusplus
}
#endif

#endif	/* BIT_VECTOR_H_ */
