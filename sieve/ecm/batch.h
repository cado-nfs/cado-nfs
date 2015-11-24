#ifndef COFAC_LIST_H
#define COFAC_LIST_H

#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include "facul.h"
#include "facul_doit.h"
#include "utils.h"

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1

/* TODO: should replace by mpz_array_t from utils/mpz_array.h */
typedef struct {
  mpz_t *l;
  size_t alloc;
  size_t size;
} mpz_list_t;
typedef mpz_list_t mpz_list[1];

typedef struct {
  int64_t *a;
  uint64_t *b;
  mpz_t *R;    /* cofactor on side 0 */
  mpz_t *A;    /* cofactor on side 1 */
  mpz_t *R0;   /* initial cofactor on side 0 */
  mpz_t *A0;   /* initial cofactor on side 1 */
  mpz_t *sq;   /* special-q */
  size_t alloc;
  size_t size;
  uint32_t *perm; /* permutation to access elements */
} cofac_list_t;
typedef cofac_list_t cofac_list[1];

#ifdef __cplusplus
extern "C" {
#endif

void cofac_list_init (cofac_list);
void cofac_list_realloc (cofac_list, size_t);
void cofac_list_clear (cofac_list);
void cofac_list_add (cofac_list, long, unsigned long, mpz_t, mpz_t, mpz_t);
unsigned long prime_product (mpz_t, prime_info, unsigned long, unsigned long);
unsigned long find_smooth (cofac_list, int[2], unsigned long[2], mpz_t[2],
                           FILE*);
void factor (cofac_list, unsigned long, cado_poly, int, int, FILE*);
void create_batch_file (const char*, mpz_t, unsigned long, unsigned long,
                        mpz_poly_t, FILE*);

#ifdef __cplusplus
}
#endif

#endif /* COFAC_LIST_H */
