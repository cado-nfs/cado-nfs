#ifndef FILTER_IO_H_
#define FILTER_IO_H_

#define MAX_FILES 1000000
#include <unistd.h>
#include <time.h>

/* Size of relations buffer between parsing & processing.
 * CAREFUL! SIZE_BUF_REL must be greater (at least double) than (1<<(NNFR+1)).
 * Stores the sentences precomputed but not inserted. 
 * About 64K sentences for the optimal.
 */
// #define SIZE_BUF_REL MAX((1<<(NFR+NNFR+1+(NFR==0))),(1<<6))
#define SIZE_BUF_REL (1<<15)

/* we want the programs to specify in a completely explicit way whether
 * they want stuff in base 10 or 16 */
#define EARLYPARSE_NEED_LEGACY_AB_DECIMAL              1
// #define EARLYPARSE_NEED_LEGACY_AB_HEXA                 2
#define EARLYPARSE_NEED_LINE                    4
/* for reading ideals (e.g. output of las) */
#define EARLYPARSE_NEED_PRIMES                  8
/* for reading index (i.e. renumber ideal) */
#define EARLYPARSE_NEED_INDEX                  16
#define EARLYPARSE_NEED_SM                     32
#define EARLYPARSE_NEED_AB                     64


/* Initial size of primes_data array in earlyparsed_relation_s,
   If more than NB_PRIMES_OPT is needed (should be rare), *primes is
   allocated
*/
#define NB_PRIMES_OPT 31

#define EARLYPARSED_RELATION_MAX_AB 4


/* This field does not contain a proper relation, but only something
 * which has undergone quite limited parsing within a thread whose job is
 * to read data fast, and not to bother about the fine points of the
 * relation semantics. Because of this, there is no unique definition of
 * which are the valid fields below. This depends on how this field is
 * meant to be used downstream. Depending on the earlyparse_needed_data
 * bitmask argument fed to filter_relations, we may fill only some fields.
 * Which fields are filled is controlled by which of the
 * EARLYPARSE_NEED_* appear in the earlyparse_needed_data bitmask. The
 * callback thread function given to process_rels is then called for each
 * such "relation"
 */
struct earlyparsed_relation_s {
  int64_t ab[EARLYPARSED_RELATION_MAX_AB];
  prime_t *primes;      /*if nb_alloc <= NB_PRIME_OPT, primes == primes_data*/
  prime_t primes_data[NB_PRIMES_OPT];
  weight_t nb;           /* number of primes */
  weight_t nb_alloc;     /* allocated space for primes
                          * (if > NB_PRIMES_OPT: indirect addressing,
                          * otherwise primes == primes_data) */
  /* nb_above_min_index is counted only when ->primes is needed anyway,
   * so we defer it to the callback function instead.
   */
  // weight_t nb_above_min_index; /* nb of primes above min_index, must be <=nb */
  uint64_t num;          /* (absolute) relation number */
  char *line;           /* If not NULL, contains the relation with a '\n' at the end */
  mpz_t * sm;
  int sm_size;
  int sm_alloc;
};
typedef struct earlyparsed_relation_s earlyparsed_relation[1];
typedef struct earlyparsed_relation_s * earlyparsed_relation_ptr;
typedef const struct earlyparsed_relation_s * earlyparsed_relation_srcptr;

static const unsigned char ugly[256] = {
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   255, 255, 255, 255, 255, 255,
 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 
 };

#ifdef __cplusplus
extern "C" {
#endif

extern void realloc_buffer_primes (earlyparsed_relation_ptr buf);

static inline void realloc_buffer_primes_c (earlyparsed_relation_ptr buf)
{
  realloc_buffer_primes (buf);
}

extern int filter_rels_force_posix_threads;
/*
 * A pointer to such a structure must be provided to filter_rels, and
 * termination is indicated by f==NULL. The function specified in the
 * k-th member (starting from 1) in this array describes the operation
 * performed at level k in the process, level 0 being the implicit
 * production of relation from the input files.  Alongside with the
 * relation on which the thread is allowed to work, all threads at level
 * k receive the void* argument specified in the array member. n
 * specifies the number of worker threads to be used for each step.
 */

struct filter_rels_description {
    void * (*f)(void*, earlyparsed_relation_ptr);
    void * arg;
    int n;
};

typedef void *(*filter_rels_callback_t) (void *, earlyparsed_relation_ptr);

extern uint64_t filter_rels2(char ** input_files,
        struct filter_rels_description * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr);

static inline uint64_t filter_rels(char ** input_files,
        filter_rels_callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    /* Of course I would prefer designated intializers here.
     * Unfortunately this header is included by C++ code as well, which
     * makes them illegal (at least with gcc-4.6.3 ; gcc-4.8.1 groks
     * them). So I stick to dumb code.
     */
    struct filter_rels_description desc[2] = {
        { f, arg, 1, }, { NULL, NULL, 0, },
    };
    return filter_rels2(input_files, desc, earlyparse_needed_data, active, stats);
}


#ifdef __cplusplus
}
#endif

#endif /* FILTER_IO_H_ */

