#ifndef GALOIS_UTILS_H_
#define GALOIS_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Galois automorphisms
   autom2.1: 1/x
   autom2.2: -x
   autom3.1: 1-1/x
   autom3.2: -1-1/x
   autom4.1: -(x+1)/(x-1)       [2 pops every two applications]
   autom6.1: -(2x+1)/(x-1)      [3 pops every two applications]

   the code is in galois_utils.c
*/

struct galois_automorphism_s {
    const char * name;
    int order;
    /* a and b get a common factor [factor] every [period] applications
     * of sigma. So in total we should have H^order =
     * factor^(order/period) where H is the homography below.  The code
     * currently wants the factor to be prime.
     */
    int factor;
    int period;
    int64_t H[4];       /* only because mpz_poly_homography wants int64_t's */
};
typedef struct galois_automorphism_s galois_automorphism[1];
typedef struct galois_automorphism_s * galois_automorphism_ptr;
typedef const struct galois_automorphism_s * galois_automorphism_srcptr;

galois_automorphism_srcptr galois_automorphism_get(const char * name);

/* returns how many times sigma->factor must by added to the
 * corresponding relation */
int galois_automorphism_apply_ab(galois_automorphism_srcptr sigma, uint64_t * a, int64_t * b);

void galois_automorphism_apply_root(galois_automorphism_srcptr sigma, mpz_ptr r1, mpz_srcptr r0, mpz_srcptr p);


void automorphism_init(int *ord, int mat[4], const char *galois_autom);
unsigned long automorphism_apply(residueul_t mat[4], unsigned long r, const modulusul_t mm, const unsigned long qq);

#ifdef __cplusplus
}
#endif


#endif /* GALOIS_UTILS_H_ */
