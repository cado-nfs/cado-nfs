#ifndef RENUMBER_H_
#define RENUMBER_H_

#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "cado_poly.h"
#include "rootfinder.h"
#include "timing.h"
#include "math.h"
#include "typedefs.h"

#define MAX_LOG_CACHED 20
#define RENUMBER_MAXLINE 1024
#define RENUMBER_MAX_ABOVE_BADIDEALS 8
#define RENUMBER_SPECIAL_VALUE ((p_r_values_t) -1)
#define RENUMBER_DEFAULT_SIZE (1 << 22)
#define RENUMBER_ROOT_ON_RAT_SIDE ((p_r_values_t) -1)
/* Max number of char necessary to write all the roots in hexa for ONE prime p.
 * Must be greater than (max_nb_char_per_roots+1)*MAX_DEGREE*NB_POLYS_MAX.
 * max_nb_char_per_roots is 8 (a root is always smaller than a prime and a
 * prime is smaller than 2^64 (this bound is not tight at all, in practice a
 * prime is smaller 2^32~2^34 at most)), and the +1 is here to take into account
 * the '\n' at the end of each entry of the table
 */
#define RENUMBER_MAX_SIZE_PER_PRIME (1<<10)
#if RENUMBER_MAX_SIZE_PER_PRIME < (8+1)*MAX_DEGREE*NB_POLYS_MAX
  #error "RENUMBER_MAX_SIZE_PER_PRIME is too small."
#endif

struct bad_ideals_s
{
  int n; // number of p_r_values that correspond to more than one ideals.
  p_r_values_t * p;  // these p_r_values in two separate tables
  p_r_values_t * r;  //
  int * side;
  int * nb;            // the number of ideals for each
  p_r_values_t max_p;  // all badideals have p <= max_p
};

struct renumber_s
{
  FILE * file;                      // file containing the renumbering table
  struct bad_ideals_s bad_ideals; // the bad ideals
  p_r_values_t * table;             // renumbering table
  uint64_t size;                    // size of the renumbering table
  unsigned int nb_polys;            // Between 2 and NB_POLYS_MAX
  uint8_t nb_bits; // number of bits taken by an index in the file: 32 or 64
  int rat;         // At most 1 poly can have degree 1.
                   // If no poly has degree 1 => rat = -1
                   // If 1 poly has degree 1 => rat = <side of this poly>
                   // For now if nb_poly > 2, rat must be -1
  index_t *cached; // We cached the index for primes < 2^MAX_LOG_CACHED
  unsigned long *lpb;               // The large prime bounds
  unsigned long max_lpb;            // maximum value of lpb[0..nb_poly-1]
  p_r_values_t *biggest_prime_below_lpb;
  index_t *index_biggest_prime_below_lpb;
  p_r_values_t smallest_prime_not_cached;
  index_t index_smallest_prime_not_cached;
  uint64_t nonmonic; /* bit vector of which polynomials are non monic. */
  uint8_t naddcols; /* number of additional columns due to non monic poly */
};
typedef struct renumber_s renumber_t[1];
typedef struct renumber_s * renumber_ptr;
typedef const struct renumber_s * renumber_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void renumber_init_for_writing (renumber_ptr, unsigned int, int, int,
                                uint64_t, unsigned long *);
/* Last argument of renumber_write_open can be NULL. It will not print the
   polynomials on the file */
void renumber_write_open (renumber_ptr, const char *, const char *, cado_poly);
void renumber_write_p (renumber_ptr, unsigned long, unsigned long [][MAX_DEGREE],
                       int []);
size_t renumber_write_buffer_p (char *, renumber_ptr, unsigned long,
                                unsigned long [][MAX_DEGREE], int []);
void renumber_write_close (renumber_ptr, const char*);

void renumber_init_for_reading (renumber_ptr);
void renumber_read_table (renumber_ptr, const char *);

void renumber_clear (renumber_ptr);

int renumber_is_bad(int *, index_t*,renumber_srcptr, p_r_values_t, p_r_values_t, int);
int renumber_is_additional_column (renumber_srcptr, index_t);
index_t renumber_get_index_from_p_r (renumber_srcptr, p_r_values_t, p_r_values_t,int);
index_t renumber_get_random_index_from_p_side(renumber_srcptr renumber_info,
    p_r_values_t p, int side);
void renumber_get_p_r_from_index (renumber_srcptr, p_r_values_t *, p_r_values_t *,
                                                    int *, index_t, cado_poly);
int renumber_get_side_from_index (renumber_srcptr, index_t, cado_poly);
int renumber_badideal_get_p_r_below (renumber_srcptr, p_r_values_t *,
                                     p_r_values_t *, int *, index_t);

#ifdef __cplusplus
}
#endif

#endif
