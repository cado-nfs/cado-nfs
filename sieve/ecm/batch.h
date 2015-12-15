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

typedef struct {
  unsigned long *l;
  size_t alloc;
  size_t size;
} ulong_list_t;
typedef ulong_list_t ulong_list[1];

/* structure to compute on-line a product tree, avoiding to first compute a
   list of mpz_t (which might take too much memory) */
typedef struct {
  mpz_t *l;     /* the value stored is l[0] * l[1] * ... * l[size-1],
                   where l[0] is the product of n[0] elements, l[1] is
                   the product of n[1] elements, ..., with n[0]=0 or 1,
                   n[1]=0 or 2, ..., n[k]=0 or 2^k */
  unsigned long *n;
  size_t size;
} mpz_product_tree_t;
typedef mpz_product_tree_t mpz_product_tree[1];

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
unsigned long find_smooth (cofac_list l, mpz_t batchP[2],
    FILE *out, int nthreads MAYBE_UNUSED);
void factor (cofac_list, unsigned long, cado_poly, int[], int, FILE*, int);
void create_batch_file (const char *f, mpz_t P, unsigned long L,
    mpz_poly_t pol, FILE *out, int nthreads MAYBE_UNUSED);

#ifdef __cplusplus
}
#endif

#endif /* COFAC_LIST_H */
