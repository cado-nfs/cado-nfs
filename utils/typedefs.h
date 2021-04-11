#ifndef	CADO_TYPEDEFS_H_
#define	CADO_TYPEDEFS_H_

// pragma no prototypes
#include <inttypes.h>   // PRIx32 etc
#include "ularith.h" /* NEEDED for LONG_BIT (32 or 64) */

/* data type to store the (p,r) values */
#ifndef SIZEOF_P_R_VALUES
#define SIZEOF_P_R_VALUES 4
#elif SIZEOF_P_R_VALUES != 4 && SIZEOF_P_R_VALUES != 8
#error "Defined constant SIZEOF_P_R_VALUES should be 4 or 8"
#endif

/* data type to store the renumber table */
#ifndef SIZEOF_INDEX
#define SIZEOF_INDEX 4
#elif SIZEOF_INDEX < 4 || 8 < SIZEOF_INDEX
#error "Defined constant SIZEOF_INDEX should be in [4..8]"
#endif

#if SIZEOF_INDEX > SIZEOF_P_R_VALUES
#error "SIZEOF_INDEX should be smaller or equal to SIZEOF_P_R_VALUES"
#endif

#if (SIZEOF_P_R_VALUES * 8) > LONG_BIT
#error "SIZEOF_P_R_VALUES cannot be greater than LONG_BIT / 8"
#endif

#if SIZEOF_P_R_VALUES == 4
typedef uint32_t p_r_values_t;
#define CADO_MPI_P_R_VALUES_T CADO_MPI_UINT32_T
#define PRpr PRIx32
#define SCNpr SCNx32
#else /* SIZEOF_P_R_VALUES == 8 */
typedef uint64_t p_r_values_t;
#define CADO_MPI_P_R_VALUES_T CADO_MPI_UINT64_T
#define PRpr PRIx64
#define SCNpr SCNx64
#endif

#if SIZEOF_INDEX == 4
typedef uint32_t index_t;
typedef int32_t index_signed_t;
#define PRid PRIx32
#define SCNid SCNx32
#else /* SIZEOF_INDEX == 8 */
typedef uint64_t index_t;
typedef int64_t index_signed_t;
#define PRid PRIx64
#define SCNid SCNx64
#endif 

/* The weight of ideals saturates at 255 */
/* For relations, we hope that there will never be more */
/* than 255 ideals per relation */
typedef uint8_t weight_t;
typedef int8_t exponent_t;
#define REL_MAX_SIZE 255

typedef struct {
  index_t h;
  p_r_values_t p;
  exponent_t e;
  uint8_t side;
} prime_t;

typedef struct {
  index_t id;
  int32_t e;
} ideal_merge_t;

typedef struct {
  uint64_t nrows;
  uint64_t ncols;
  double W; /* weight of the active part of the matrix */
} info_mat_t;

#endif	/* CADO_TYPEDEFS_H_ */
