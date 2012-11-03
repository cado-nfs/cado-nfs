#ifndef __BUCKETS_H__
#define __BUCKETS_H__

#include <stdint.h>

#include "types.h"
#include "fb.h"


#ifdef BUCKET_RESIEVE
typedef struct {
    unsigned pos;
    unsigned hint;
} __replayable_update_struct;

typedef struct {
    __replayable_update_struct *b;
    unsigned n;
} __replayable_bucket_struct;

typedef       __replayable_bucket_struct  replayable_bucket_t[1];
typedef       __replayable_bucket_struct *replayable_bucket_ptr;
typedef const __replayable_bucket_struct *replayable_bucket_srcptr;
#else
typedef void *replayable_bucket_t;
typedef void *replayable_bucket_ptr;
typedef void *replayable_bucket_srcptr;
#endif


/* Array of buckets.
 *****************************************************************************/

// Forward declaration of bucket updates.
typedef struct __update_packed_struct update_packed_t;

typedef struct {
  unsigned          n;
  unsigned          max_size;
  unsigned          min_degp, max_degp;
#ifdef BUCKET_RESIEVE
  unsigned         *first_hint;
  unsigned         *current_hint;
#endif
  update_packed_t **start;
  update_packed_t **degp_end;
} __buckets_struct;

typedef       __buckets_struct  buckets_t[1];
typedef       __buckets_struct *buckets_ptr;
typedef const __buckets_struct *buckets_srcptr;

// Initialize structure and allocate buckets.
void buckets_init(buckets_ptr buckets, unsigned I, unsigned J,
                  double hits, unsigned min_degp, unsigned max_degp);

// Clean up memory.
void buckets_clear(buckets_ptr buckets);

// Return the size of a bucket region.
unsigned bucket_region_size();

// Print information about the buckets.
void print_bucket_info(buckets_t *buckets0, unsigned npol0,
                       buckets_t *buckets1, unsigned npol1);

// Fill the buckets with updates corresponding to divisibility by elements of
// the factor base.
void buckets_fill(buckets_ptr buckets, large_factor_base_srcptr FB,
        sublat_srcptr sublat, unsigned I, unsigned J, qlat_srcptr qlat);

// Apply all the updates from a given bucket to the sieve region S.
void bucket_apply(uint8_t *S, buckets_srcptr buckets, unsigned k);

#ifdef BUCKET_RESIEVE
void bucket_prepare_replay(replayable_bucket_ptr bb,
        buckets_srcptr buckets, uint8_t *S, unsigned k);
void bucket_apply_at_pos(fppol_ptr norm, ijpos_t pp, 
        replayable_bucket_srcptr buckets, large_factor_base_srcptr FB);
#endif


#endif  /* __BUCKETS_H__ */
