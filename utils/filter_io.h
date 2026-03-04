#ifndef CADO_FILTER_IO_H
#define CADO_FILTER_IO_H

#include <stdint.h>     // for uint64_t, int64_t
#include <time.h>       // for NULL

#ifdef __cplusplus
#include <string>
#include <vector>
#endif

#include <gmp.h>        // for mpz_t
                        //
#include "bit_vector.h" // for bit_vector_srcptr
#include "timing.h"     // for timingstats_dict_ptr
#include "typedefs.h"   // for prime_t, weight_t

#define MAX_FILES 1000000

#define RELATION_MAX_BYTES 4096

/* Size of relations buffer between parsing & processing.
 * CAREFUL! SIZE_BUF_REL must be greater (at least double) than (1<<(NNFR+1)).
 * Stores the sentences precomputed but not inserted.
 * About 64K sentences for the optimal.
 */
// #define SIZE_BUF_REL MAX((1<<(NFR+NNFR+1+(NFR==0))),(1<<6))
#define SIZE_BUF_REL (1U << 15U)

/* we want the programs to specify in a completely explicit way whether
 * they want stuff in base 10 or 16 */
#define EARLYPARSE_NEED_AB_DECIMAL 1U
#define EARLYPARSE_NEED_AB_HEXA 2U
#define EARLYPARSE_NEED_LINE 4U
/* for reading ideals (e.g. output of las) */
#define EARLYPARSE_NEED_PRIMES 8U
/* for reading index (i.e. renumber ideal) */
#define EARLYPARSE_NEED_INDEX 16U
#define EARLYPARSE_NEED_SM 32U
#define EARLYPARSE_NEED_SORTED 64U
#define EARLYPARSE_NEED_INDEX_SORTED                                           \
    (EARLYPARSE_NEED_INDEX | EARLYPARSE_NEED_SORTED)

/* Initial size of primes_data array in earlyparsed_relation_s,
   If more than NB_PRIMES_OPT is needed (should be rare), *primes is
   allocated
*/
#define NB_PRIMES_OPT 31

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
    int64_t a;
    uint64_t b;
    int active_sides[2];
    prime_t * primes; /*if nb_alloc <= NB_PRIME_OPT, primes == primes_data*/
    prime_t primes_data[NB_PRIMES_OPT];
    weight_t nb;       /* number of primes */
    weight_t nb_alloc; /* allocated space for primes
                        * (if > NB_PRIMES_OPT: indirect addressing,
                        * otherwise primes == primes_data) */
    /* nb_above_min_index is counted only when ->primes is needed anyway,
     * so we defer it to the callback function instead.
     */
    // weight_t nb_above_min_index; /* nb of primes above min_index, must be
    // <=nb */
    uint64_t num; /* (absolute) relation number */
    char * line; /* If not NULL, contains the relation with a '\n' at the end */
    mpz_t * sm;
    int sm_size;
    int sm_alloc;
};
typedef struct earlyparsed_relation_s earlyparsed_relation[1];
typedef struct earlyparsed_relation_s * earlyparsed_relation_ptr;
typedef const struct earlyparsed_relation_s * earlyparsed_relation_srcptr;

struct earlyparsed_relation_mpz_s {
    mpz_t a;
    mpz_t b;
    int active_sides[2];
    prime_t * primes; /*if nb_alloc <= NB_PRIME_OPT, primes == primes_data*/
    prime_t primes_data[NB_PRIMES_OPT];
    weight_t nb;       /* number of primes */
    weight_t nb_alloc; /* allocated space for primes
                        * (if > NB_PRIMES_OPT: indirect addressing,
                        * otherwise primes == primes_data) */
    /* nb_above_min_index is counted only when ->primes is needed anyway,
     * so we defer it to the callback function instead.
     */
    // weight_t nb_above_min_index; /* nb of primes above min_index, must be
    // <=nb */
    uint64_t num; /* (absolute) relation number */
    char * line; /* If not NULL, contains the relation with a '\n' at the end */
    mpz_t * sm;
    int sm_size;
    int sm_alloc;
};
typedef struct earlyparsed_relation_mpz_s earlyparsed_relation_mpz[1];
typedef struct earlyparsed_relation_mpz_s * earlyparsed_relation_mpz_ptr;
typedef const struct earlyparsed_relation_mpz_s * earlyparsed_relation_mpz_srcptr;

static unsigned char const ugly[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   255, 255,
    255, 255, 255, 255, 255, 10,  11,  12,  13,  14,  15,  255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 10,  11,  12,  13,  14,  15,  255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255};

#ifdef __cplusplus
extern "C" {
#endif

void realloc_buffer_primes_c(earlyparsed_relation_ptr buf);

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
    void * (*f)(void *, earlyparsed_relation_ptr);
    void * arg;
    int n;
};

struct filter_rels_mpz_description {
    void * (*f)(void *, earlyparsed_relation_mpz_ptr);
    void * arg;
    int n;
};

extern uint64_t filter_rels2(char const ** input_files,
                             struct filter_rels_description * desc,
                             int earlyparse_needed_data,
                             bit_vector_srcptr active, timingstats_dict_ptr);

typedef void * (*filter_rels_callback_t)(void *, earlyparsed_relation_ptr);

extern uint64_t filter_rels(char const ** input_files, filter_rels_callback_t f,
                            void * arg, int earlyparse_needed_data,
                            bit_vector_srcptr active,
                            timingstats_dict_ptr stats);

typedef void * (*filter_rels_mpz_callback_t)(void *, earlyparsed_relation_mpz_ptr);
extern uint64_t filter_rels_mpz(
        char const ** input_files,
        filter_rels_mpz_callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

template<typename T>
concept filter_io_config = std::is_empty_v<T>
    && requires { typename T::rel_t; }
    && requires { typename T::rel_ptr; }
    && requires { typename T::rel_srcptr; }
    && requires { typename T::description_t; }
    && requires { typename T::callback_t; }
    && requires { T::ab_init; }
    && requires { T::ab_clear; };

struct filter_io_default_cfg
{
    using rel_t = earlyparsed_relation;
    using rel_ptr = earlyparsed_relation_ptr;
    using rel_srcptr = earlyparsed_relation_srcptr;
    using description_t = filter_rels_description;
    using callback_t = filter_rels_callback_t;
    static constexpr void (*ab_init)(uint64_t) = nullptr;
    static constexpr void (*ab_clear)(uint64_t) = nullptr;
};
static_assert(filter_io_config<filter_io_default_cfg>);

struct filter_io_large_ab_cfg
{
    using rel_t = earlyparsed_relation_mpz;
    using rel_ptr = earlyparsed_relation_mpz_ptr;
    using rel_srcptr = earlyparsed_relation_mpz_srcptr;
    using description_t = filter_rels_mpz_description;
    using callback_t = filter_rels_mpz_callback_t;
    static constexpr void (*ab_init)(mpz_ptr) = &mpz_init;
    static constexpr void (*ab_clear)(mpz_ptr) = &mpz_clear;
};
static_assert(filter_io_config<filter_io_large_ab_cfg>);

template<filter_io_config cfg>
void realloc_buffer_primes(typename cfg::rel_ptr buf);

template<filter_io_config cfg>
extern uint64_t filter_rels2(
        std::vector<std::string> const & input_files,
        typename cfg::description_t * desc,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats);

template<filter_io_config cfg>
uint64_t
filter_rels(
        std::vector<std::string> const & input_files,
        typename cfg::callback_t f,
        void * arg,
        int earlyparse_needed_data,
        bit_vector_srcptr active,
        timingstats_dict_ptr stats)
{
    typename cfg::description_t desc[2] = {
        { f, arg, 1, }, { nullptr, nullptr, 0, },
    };
    return filter_rels2<cfg>(input_files, desc, earlyparse_needed_data, active, stats);
}
#endif

#endif /* CADO_FILTER_IO_H */
