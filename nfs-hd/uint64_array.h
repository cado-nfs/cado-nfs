#ifndef UINT64_ARRAY_H
#define UINT64_ARRAY_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  uint64_t length;
  uint64_t * array;
} s_uint64_array_t;

typedef s_uint64_array_t uint64_array_t[1];
typedef s_uint64_array_t * uint64_array_ptr;
typedef const s_uint64_array_t * uint64_array_srcptr;

/*
 * Init an uint64_array.
 */
void uint64_array_init(uint64_array_ptr array, uint64_t length);

/*
 * Set coeff at index in array.
 */
void uint64_array_set_coeff(uint64_array_ptr array, uint64_t index,
    uint64_t coeff);

/*
 * Decrease the size of an array.
 */
void uint64_array_realloc(uint64_array_ptr array, uint64_t number);

/*
 * Clear an array.
 */
void uint64_array_clear(uint64_array_ptr array);

/*
 * Print an array.
 */
void uint64_array_fprintf(FILE * file, uint64_array_srcptr array);


#ifdef __cplusplus
}
#endif
#endif /* UINT64_ARRAY_H */
