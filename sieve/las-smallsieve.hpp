#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <stdarg.h>
#include "fb-types.h"

/* Structures for small sieves */

typedef struct {
    fbprime_t p;
    fbprime_t r;        // in [ 0, p [
    fbprime_t offset;   // in [ 0, p [
} ssp_t;

/* We currently *mandate* that this structure has the same size as ssp_t.
 * It would be possible to make it work with only a requirement on
 * identical alignment and smaller size. If extra fields are required, we
 * need to store them with the ssp_marker_t structure.
 */
typedef struct {
    fbprime_t g, q, U;
} ssp_bad_t;

#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD     (1u<<30)
#define SSP_END         (1u<<31)

typedef struct {
    int index;
    unsigned int event;
} ssp_marker_t;

typedef struct {
    ssp_marker_t * markers;
    ssp_t *ssp;
    int nb_ssp;
    unsigned char * logp;
} small_sieve_data_t;

/* Include this only now, as there's a cross dependency between the two
 * (our prototypes need sieve_info_t, which needs our datatypes...)
 */
#include "las-types.hpp"
#include "bucket.hpp"

extern void small_sieve_info(const char * what, int side, small_sieve_data_t * r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t * ssd);
extern void small_sieve_extract_interval(small_sieve_data_t * r, small_sieve_data_t * s, int bounds[2]);
extern void small_sieve_init(small_sieve_data_t *ssd, las_info& las, const std::vector<fb_general_entry> *fb,
                      sieve_info const & si, int side);
extern int64_t * small_sieve_copy_start(int64_t * base, int bounds[2]);
extern int64_t * small_sieve_start(small_sieve_data_t *ssd, unsigned int j0, sieve_info const & si);
extern void small_sieve_skip_stride(small_sieve_data_t *ssd, int64_t * ssdpos, unsigned int skip, sieve_info const & si);
extern void sieve_small_bucket_region(unsigned char *S, int N,
			       small_sieve_data_t * ssd, int64_t * ssdpos, sieve_info& si,
                               int side,
			       where_am_I& w MAYBE_UNUSED);

extern void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t *ssd, int64_t * ssdpos,
        sieve_info const & si, where_am_I& w MAYBE_UNUSED);


#endif	/* LAS_SMALLSIEVE_HPP_ */
