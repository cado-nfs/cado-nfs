#ifndef MPZ_VECTOR_H_
#define MPZ_VECTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

    /* TODO: This should be an std::vector<cxx_mpz>.
     */
typedef struct 
{
  unsigned int dim; /* dimension of the vector */
  mpz_t *c;         /* its coordinates */
} mpz_vector_struct_t;

typedef mpz_vector_struct_t mpz_vector_t[1];
typedef mpz_vector_struct_t * mpz_vector_ptr;
typedef const mpz_vector_struct_t * mpz_vector_srcptr;

/* Management of the structure: init, clear, set and swap. */
void mpz_vector_init (mpz_vector_t, unsigned int);
void mpz_vector_clear(mpz_vector_t);
void mpz_vector_swap (mpz_vector_t, mpz_vector_t);
void mpz_vector_set (mpz_vector_t, mpz_vector_t);
void mpz_vector_setcoordinate (mpz_vector_t, unsigned int, mpz_t);
void mpz_vector_setcoordinate_ui (mpz_vector_t, unsigned int, unsigned int);
void mpz_vector_setcoordinate_si (mpz_vector_t, unsigned int, int);
void mpz_vector_setcoordinate_uint64 (mpz_vector_t, unsigned int, uint64_t);
void mpz_vector_setcoordinate_int64 (mpz_vector_t, unsigned int, int64_t);
int mpz_vector_is_coordinate_zero (mpz_vector_t, unsigned int);
/*
  Return 0 if a and b are equal,
  -1 if a is smaller and 1 if a is bigger.
*/
int mpz_vector_cmp (mpz_vector_srcptr a, mpz_vector_srcptr b);
void mpz_vector_fprintf(FILE * file, mpz_vector_srcptr v);

/* Implementation of dot product and norm (skew and non-skew version) */
void mpz_vector_dot_product (mpz_t, mpz_vector_t, mpz_vector_t);
void mpz_vector_skew_dot_product (mpz_t, mpz_vector_t, mpz_vector_t, mpz_t);
void mpz_vector_norm (mpz_t, mpz_vector_t);
void mpz_vector_skew_norm (mpz_t, mpz_vector_t, mpz_t);

/* Convert from mpz_vector_t to mpz_poly */
void mpz_vector_get_mpz_poly (mpz_poly, mpz_vector_t);

/* Operations on vectors */
void mpz_vector_submul (mpz_vector_t r, mpz_t q, mpz_vector_t v);

/* Lagrange algo to reduce lattice of rank 2 */
void mpz_vector_Lagrange (mpz_vector_t, mpz_vector_t,
                          mpz_vector_t, mpz_vector_t, mpz_t);
void mpz_vector_reduce_with_max_skew (mpz_vector_t, mpz_vector_t, mpz_t,
                                      mpz_vector_t, mpz_vector_t, mpz_t, int);

#ifdef __cplusplus
}
#endif

#endif	/* MPZ_VECTOR_H_ */
