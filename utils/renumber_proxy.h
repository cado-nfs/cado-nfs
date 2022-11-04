#ifndef RENUMBER_PROXY_H_
#define RENUMBER_PROXY_H_

#ifndef __cplusplus
#include <stdbool.h>
#endif
#include <stddef.h>     // size_t
#include <stdint.h>     // uint64_t
#include "typedefs.h"
#include "mpz_poly.h"
#include "cado_poly.h"

/* This interface provides C-level access to the public interface of the
 * renumber structure.
 */

/*
 * To read a renumber table from a file, this goes as:

 renumber_proxy_t renumber;
 renumber_table_init(renumber, cpoly);
 renumber_table_read_from_file(renumber, renumberfilename);

 * There is currently no way to build a renumber table from the C-level
 * interface.
*/

struct renumber_proxy_s {
    void * x;
};
typedef struct renumber_proxy_s renumber_proxy_t[1];
typedef struct renumber_proxy_s * renumber_proxy_ptr;
typedef const struct renumber_proxy_s * renumber_proxy_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/* There is no "default constructor" for the renumber_proxy, because that
 * doesn't make sense in any of the uses that we have.
 */
extern void renumber_table_init(renumber_proxy_ptr, cado_poly_ptr);
extern void renumber_table_clear(renumber_proxy_ptr);
extern void renumber_table_set_lpb(renumber_proxy_ptr, const unsigned int *, size_t);
extern void renumber_table_read_from_file(renumber_proxy_ptr, const char * filename, int for_dl);

/* a few trivial accessors */
extern int renumber_table_get_format(renumber_proxy_srcptr);
extern unsigned int renumber_table_get_lpb(renumber_proxy_srcptr, int side);
extern unsigned int renumber_table_get_max_lpb(renumber_proxy_srcptr);
extern unsigned int renumber_table_get_min_lpb(renumber_proxy_srcptr);
extern uint64_t renumber_table_get_size(renumber_proxy_srcptr);
extern unsigned int renumber_table_get_nb_polys(renumber_proxy_srcptr);
extern mpz_poly_srcptr renumber_table_get_poly(renumber_proxy_srcptr, int side);
extern int renumber_table_get_poly_deg(renumber_proxy_srcptr, int side);
extern int renumber_table_get_rational_side(renumber_proxy_srcptr);
extern index_t renumber_table_get_max_index(renumber_proxy_srcptr);
extern index_t renumber_table_get_max_cached_index(renumber_proxy_srcptr);
extern index_t number_of_additional_columns(renumber_proxy_srcptr);
extern size_t renumber_table_get_sides_of_additional_columns(renumber_proxy_srcptr, int * sides, size_t * n);
extern index_t number_of_bad_ideals(renumber_proxy_srcptr);
extern size_t renumber_table_get_memory_size(renumber_proxy_srcptr);

/* most important outer-visible routines: lookups */

/* special lookups:
 *  a lookup of an additional column on side s returns {0, 0, s}
 *  a lookup of a bad ideal (given by (p,r)) returns the index of the
 *  _first_ ideal above this (p,r)
 * except for the non-injectivity of the mapping index->(p,r) for bad
 * ideals, the map is "almost" a bijection.
 */

/* return the number of bad ideals above x (and therefore zero if
 * x is not bad) ; likewise for index h.
 * If the ideal is bad, put in the reference [first] the
 * first index that corresponds to the bad ideals.
 */
extern int renumber_table_index_is_bad(renumber_proxy_srcptr, index_t * first_index, index_t);
extern int renumber_table_p_r_side_is_bad(renumber_proxy_srcptr, index_t *, p_r_values_t p, p_r_values_t r, int side);
extern bool renumber_table_index_is_additional_column(renumber_proxy_srcptr, index_t h);
extern index_t renumber_table_index_from_p_r (renumber_proxy_srcptr, p_r_values_t p, p_r_values_t r, int side);
extern bool renumber_table_p_r_from_index(renumber_proxy_srcptr, p_r_values_t *, p_r_values_t *, int *, index_t);
extern bool renumber_table_indices_from_p_a_b(renumber_proxy_srcptr R, index_t * first, int * exps, size_t * nexps, p_r_values_t p, p_r_values_t r, int side, int e, int64_t a, uint64_t b);

#ifdef __cplusplus
}
#endif

#endif /* RENUMBER_PROXY_H_ */
